"""Region-level annotations: observed coverage (-c), unique regions (-u) and target overlaps (-t)."""

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .intervals import Intervals, covered_fraction, overlap_pairs, subtract
from .matching import group_by_query
from .query import Query
from .reference import Reference, SourceData, open_text, strip_chr

logger = logging.getLogger(__name__)


@dataclass
class Targets:
    intervals: Intervals
    names: NDArray[np.object_]


def load_targets(path: str) -> Targets:
    """A CHROM START END ID BED; ``chr`` prefixes are stripped like everywhere else."""
    logger.info('Reading the following targets BED file: %s', path)
    chrom: list[str] = []
    start: list[int] = []
    end: list[int] = []
    names: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                raise SystemExit(f'Targets BED rows need CHROM START END ID; got: {line.rstrip()}')
            chrom.append(strip_chr(fields[0]))
            start.append(int(fields[1]))
            end.append(int(fields[2]))
            names.append(fields[3])
    return Targets(Intervals.from_columns(chrom, start, end), np.array(names, dtype=object))


@dataclass
class RegionResults:
    n: int
    coverage_all: NDArray[np.float64]
    coverage_by_source: dict[str, NDArray[np.float64]] = field(default_factory=dict)
    unique_pieces: dict[int, list[tuple[str, int, int]]] = field(default_factory=dict)
    target_hits: dict[int, list[str]] = field(default_factory=dict)
    unique_target_hits: dict[int, list[str]] = field(default_factory=dict)  # targets overlapping unique pieces
    piece_targets: dict[int, list[str]] = field(default_factory=dict)  # per unique piece (in output order)


def _eligible(source: SourceData, svtype: str, min_af: float, size_limit: int) -> NDArray[np.bool_]:
    """Reference rows of the given SVTYPE with AF above the cut-off and size within the limit."""
    with np.errstate(invalid='ignore'):
        af_ok = source.af > min_af
    return (source.svtype == svtype) & af_ok & (source.intervals.length <= size_limit)


def _concat(parts: list[Intervals]) -> Intervals:
    return Intervals(
        np.concatenate([p.chrom for p in parts]),
        np.concatenate([p.start for p in parts]),
        np.concatenate([p.end for p in parts]),
    )


def compute_regions(  # noqa: PLR0917
    query: Query,
    reference: Reference,
    sources: list[str],
    cov_af: float | None,
    uniq_af: float | None,
    size_limit: int,
    targets: Targets | None,
) -> RegionResults:
    results = RegionResults(n=len(query), coverage_all=np.zeros(len(query)))

    if targets is not None:
        q_idx, t_idx = overlap_pairs(query.intervals, targets.intervals)
        results.target_hits = {q: targets.names[rows].tolist() for q, rows in group_by_query(q_idx, t_idx).items()}

    if cov_af is None and uniq_af is None:
        return results

    if cov_af is not None:
        results.coverage_by_source = {s: np.zeros(len(query)) for s in sources}

    for svtype, q_rows in query.by_svtype().items():
        sub = query.subset(q_rows)

        if cov_af is not None:
            logger.info('Calculating coverage for: %s', svtype or '<no SVTYPE>')
            pooled: list[Intervals] = []
            for name in sources:
                source = reference.sources[name]
                rows = np.flatnonzero(_eligible(source, svtype, cov_af, size_limit))
                if len(rows) == 0:
                    logger.info('  %s not found in %s, skipping', svtype, name)
                    continue
                ref = source.intervals.take(rows)
                pooled.append(ref)
                results.coverage_by_source[name][q_rows] = covered_fraction(sub, ref)
            if pooled:
                results.coverage_all[q_rows] = covered_fraction(sub, _concat(pooled))

        if uniq_af is not None:
            logger.info('Calculating unique regions for: %s', svtype or '<no SVTYPE>')
            pooled = []
            for name in sources:
                source = reference.sources[name]
                rows = np.flatnonzero(_eligible(source, svtype, uniq_af, size_limit))
                if len(rows):
                    pooled.append(source.intervals.take(rows))
            if not pooled:
                logger.info('  %s not found in any sources, skipping', svtype)
                continue
            piece_q, piece_s, piece_e = subtract(sub, _concat(pooled))
            for a, s, e in zip(piece_q.tolist(), piece_s.tolist(), piece_e.tolist(), strict=True):
                q = int(q_rows[a])
                results.unique_pieces.setdefault(q, []).append((str(sub.chrom[a]), s, e))

    if uniq_af is not None and targets is not None and results.unique_pieces:
        owners = [q for q, pieces in results.unique_pieces.items() for _ in pieces]
        flat = [piece for pieces in results.unique_pieces.values() for piece in pieces]
        piece_iv = Intervals.from_columns([p[0] for p in flat], [p[1] for p in flat], [p[2] for p in flat])
        p_idx, t_idx = overlap_pairs(piece_iv, targets.intervals)
        results.piece_targets = {p: targets.names[rows].tolist() for p, rows in group_by_query(p_idx, t_idx).items()}
        for p, names in results.piece_targets.items():
            results.unique_target_hits.setdefault(owners[p], []).extend(names)

    return results


def write_uniques_bed(path: Path, query: Query, results: RegionResults, with_targets: bool) -> None:
    """The ``-u`` side file: one row per unique region, with the samples carrying the SV."""
    logger.info('Printing unique regions to %s', path)
    header = ['#CHROM', 'START', 'END', 'SVTYPE', 'SV_ID', 'Het_Samples', 'HomAlt_Samples']
    if with_targets:
        header.append('Targets')
    piece_no = 0
    with open(path, 'w', encoding='utf-8') as out:
        out.write('\t'.join(header) + '\n')
        for q, pieces in results.unique_pieces.items():
            for chrom, start, end in pieces:
                row = [
                    chrom,
                    str(start),
                    str(end),
                    str(query.svtype[q]),
                    str(query.sv_id[q]),
                    query.het_samples[q],
                    query.homalt_samples[q],
                ]
                if with_targets:
                    hits = results.piece_targets.get(piece_no)
                    row.append(','.join(hits) if hits else 'None')
                piece_no += 1
                out.write('\t'.join(row) + '\n')
