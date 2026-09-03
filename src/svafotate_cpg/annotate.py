"""The ``annotate`` command: orchestrates loading, matching and writing."""

import logging
import multiprocessing
from argparse import Namespace
from pathlib import Path

import numpy as np
from cyvcf2 import VCF, Writer
from numpy.typing import NDArray

from .annotations import EXTRA_ANNOTATIONS, AnnotationPlan, required_columns
from .matching import SourceMatches, best_per_query, find_matches, group_by_query, passes_overlap_filter
from .query import Query, load_query
from .reference import MISSING, Reference, SourceData, column_type, load_reference
from .regions import RegionResults, compute_regions, load_targets, write_uniques_bed

logger = logging.getLogger(__name__)

# reciprocal overlap fraction used when -f is not provided
DEFAULT_MINF = 0.001
# SV size ceiling used when --lim is not provided
DEFAULT_SIZE_LIMIT = 250_000_000
WRITER_MODES = {'vcf': 'w', 'vcfgz': 'wz', 'bcf': 'wbu', 'bcfgz': 'wb'}


def resolve_cpus(requested: int) -> int:
    available = max(multiprocessing.cpu_count(), 1)
    if requested > available:
        logger.warning('Too many CPUs designated. Number of CPUs will be set to %d', available)
        return available
    return max(requested, 1)


def resolve_minfs(minf: list[float] | None, sources: list[str]) -> dict[str, float]:
    """One reciprocal overlap fraction per source: a single value applies to all, else one per -s source."""
    if not minf:
        logger.info('No reciprocal overlap fraction indicated; using default %s', DEFAULT_MINF)
        minf = [DEFAULT_MINF]
    if len(minf) == 1:
        return dict.fromkeys(sources, minf[0])
    if len(minf) != len(sources):
        raise SystemExit(
            'Please list either 1 overlap fraction for all sources or the same number of overlap fractions (-f)'
            ' as there are requested sources (-s)'
        )
    return dict(zip(sources, minf, strict=True))


class SourceAnnotations:
    """Per-source results in the shape the record writer needs, all keyed by query record index."""

    def __init__(self, source: SourceData, query: Query, minf: float, conv_ins: bool, plan: AnnotationPlan) -> None:
        self.source = source
        m: SourceMatches = find_matches(query, source, conv_ins)
        keep = passes_overlap_filter(m, query, minf, conv_ins)

        self.matches = group_by_query(m.q_idx[m.same_type & keep], m.r_idx[m.same_type & keep])
        self.best_match = best_per_query(m, m.same_type, source) if plan.best else {}
        if plan.mismatch:
            self.mismatches = group_by_query(m.q_idx[~m.same_type & keep], m.r_idx[~m.same_type & keep])
            self.best_mismatch = best_per_query(m, ~m.same_type, source)
        else:
            self.mismatches, self.best_mismatch = {}, {}

        if not self.matches:
            logger.info('There are no overlap matches for %s', source.name)
        if plan.mismatch and not self.mismatches:
            logger.info('There are no overlap mismatches for %s', source.name)

    def match_values(self, q: int, column: str) -> NDArray:
        rows = self.matches.get(q)
        if rows is None:
            return np.empty(0)
        return self.source.values(column, rows)

    def row_value(self, row: int, column: str):
        """A single reference cell as a VCF-typed value; missing becomes 0 (as upstream wrote)."""
        value = self.source.columns[column][row]
        if isinstance(value, float):
            if np.isnan(value):
                return 0
            return int(value) if column_type(column) == 'Integer' else float(value)
        return 0 if value in MISSING else value

    def full_string(self, q: int) -> str:
        assert self.source.raw_rows is not None
        return ','.join(
            '|'.join([str(self.source.sv_id[r]), str(self.source.svtype[r]), *self.source.raw_rows[r]])
            for r in self.matches[q].tolist()
        )


def write_record(v, q: int, plan: AnnotationPlan, per_source: list[SourceAnnotations], regions: RegionResults) -> None:
    """Fill every planned INFO field for one record."""
    # Max_* over all sources; 0 when nothing matched
    for column in plan.max_columns:
        values = np.concatenate([s.match_values(q, column) for s in per_source]) if per_source else np.empty(0)
        if len(values) == 0:
            v.INFO[f'Max_{column}'] = 0
        elif column_type(column) == 'Float':
            v.INFO[f'Max_{column}'] = float(np.max(values.astype(float)))
        else:
            v.INFO[f'Max_{column}'] = int(np.max(values.astype(float)))

    for s in per_source:
        name = s.source.name
        v.INFO[f'{name}_Count'] = len(s.matches.get(q, ()))

        if plan.best and q in s.best_match:
            row, ofp = s.best_match[q]
            v.INFO[f'Best_{name}_ID'] = str(s.source.sv_id[row])
            v.INFO[f'Best_{name}_OFP'] = ofp
            for column in plan.max_columns:
                v.INFO[f'Best_{name}_{column}'] = s.row_value(row, column)

        if plan.mismatch:
            if q in s.mismatches:
                rows = s.mismatches[q]
                v.INFO[f'{name}_Mismatches'] = ','.join(s.source.sv_id[rows].tolist())
                v.INFO[f'{name}_Mismatches_Count'] = len(rows)
                v.INFO[f'{name}_Mismatch_SVTYPEs'] = ','.join(sorted(set(s.source.svtype[rows].tolist())))
            if q in s.best_mismatch:
                row, ofp = s.best_mismatch[q]
                v.INFO[f'Best_{name}_Mismatch_ID'] = str(s.source.sv_id[row])
                v.INFO[f'Best_{name}_Mismatch_OFP'] = ofp
                v.INFO[f'Best_{name}_Mismatch_SVTYPE'] = str(s.source.svtype[row])
                for column in ('AF', 'Het', 'HomAlt'):
                    v.INFO[f'Best_{name}_Mismatch_{column}'] = s.row_value(row, column)

        if plan.full and q in s.matches:
            v.INFO[f'{name}_Matches'] = s.full_string(q)

        if plan.coverage:
            v.INFO['SV_Cov'] = float(regions.coverage_all[q])
            v.INFO[f'{name}_SV_Cov'] = float(regions.coverage_by_source[name][q])

    if plan.unique:
        v.INFO['SV_Uniq'] = len(regions.unique_pieces.get(q, ()))
    if plan.targets and q in regions.target_hits:
        v.INFO['Target_Overlaps'] = ','.join(sorted(regions.target_hits[q]))
    if plan.targets and plan.unique and q in regions.unique_target_hits:
        v.INFO['Unique_Targets'] = ','.join(sorted(set(regions.unique_target_hits[q])))


def run(args: Namespace) -> None:
    ncpus = resolve_cpus(args.cpu)
    extras: list[str] = args.ann or []
    conv_ins: bool = args.ins
    size_limit = args.lim if args.lim is not None else DEFAULT_SIZE_LIMIT

    logger.info('Annotating the following VCF: %s', args.vcf)
    keep_raw = 'full' in extras or 'all' in extras
    reference: Reference = load_reference(
        args.bed, args.sources, keep_raw=keep_raw, wanted_columns=None if keep_raw else required_columns(extras)
    )
    sources: list[str] = list(args.sources) if args.sources else reference.source_names
    minfs = resolve_minfs(args.minf, sources)

    logger.info('Annotating with the following sources and reciprocal overlap fractions:')
    for name in sources:
        logger.info('  %s\t%s', name, minfs[name])
    if extras:
        logger.info('Additional annotations requested:')
        for a in extras:
            logger.info('  %s', EXTRA_ANNOTATIONS[a])
    if args.cov is not None:
        logger.info('SV_Cov annotation requested; AF cutoff of: %s', args.cov)
    if args.uniq is not None:
        logger.info('Producing uniques.bed output and SV_Uniq annotation requested; AF cutoff of: %s', args.uniq)

    plan = AnnotationPlan(sources, extras, args.cov, args.uniq, args.target, args.bed)

    query = load_query(args.vcf, args, ncpus)
    targets = load_targets(args.target) if args.target else None

    per_source = []
    for name in sources:
        logger.info('Finding overlaps for: %s', name)
        per_source.append(SourceAnnotations(reference.sources[name], query, minfs[name], conv_ins, plan))

    regions = compute_regions(query, reference, sources, args.cov, args.uniq, size_limit, targets)
    if plan.unique:
        write_uniques_bed(Path(args.uniq_out), query, regions, with_targets=plan.targets)

    logger.info('Adding annotations to output VCF...')
    vcf = VCF(args.vcf, threads=ncpus)
    plan.add_headers(vcf)
    writer = Writer(args.out, vcf, WRITER_MODES[args.out_type])
    for q, v in enumerate(vcf):
        write_record(v, q, plan, per_source, regions)
        writer.write_record(v)
    writer.close()
    vcf.close()
    logger.info('Annotated VCF written to: %s', args.out)
    logger.info('DONE')
