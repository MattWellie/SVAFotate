"""Reading the query VCF into arrays: one adjusted interval per record plus per-record genotype summaries."""

from __future__ import annotations

import logging
from argparse import Namespace
from dataclasses import dataclass

import numpy as np
from cyvcf2 import VCF
from numpy.typing import NDArray

from .intervals import Intervals
from .reference import strip_chr

logger = logging.getLogger(__name__)


@dataclass
class Query:
    intervals: Intervals
    sv_id: NDArray[np.object_]
    svlen: NDArray[np.int64]
    svtype: NDArray[np.object_]  # '' when the record has no SVTYPE
    het_samples: list[str]  # comma-joined sample names per record, 'None' when there are none
    homalt_samples: list[str]

    def __len__(self) -> int:
        return len(self.sv_id)

    def by_svtype(self) -> dict[str, NDArray[np.intp]]:
        """Record indices grouped by SVTYPE, in sorted SVTYPE order."""
        return {t: np.flatnonzero(self.svtype == t) for t in sorted(set(self.svtype.tolist()))}

    def subset(self, idx: NDArray[np.intp]) -> Intervals:
        return self.intervals.take(idx)


def _first_int(value) -> int | None:
    if value is None:
        return None
    if isinstance(value, tuple | list):
        value = value[0]
    return int(value)


def _ci(v, key: str, wanted: bool) -> tuple[int, int]:
    """A (lower, upper) confidence interval from INFO, or (0, 0) when absent or not requested."""
    if not wanted:
        return (0, 0)
    ci = v.INFO.get(key)
    if ci is None:
        logger.warning('%s: %s requested but absent from INFO; treating as (0, 0)', v.ID, key)
        return (0, 0)
    return (int(ci[0]), int(ci[1]))


def record_coordinates(v, args: Namespace) -> tuple[int, int]:
    """0-based half-open (start, end) for a record after every requested adjustment.

    This mirrors upstream exactly, including two quirks kept for output compatibility: END derived from
    SVLEN is ``start + |SVLEN| + 1`` (one base longer than the VCF convention implies), and confidence
    intervals are applied before ``--emb``/``--red``.
    """
    cipos = _ci(v, 'CIPOS', args.ci is not None)
    ciend = _ci(v, 'CIEND', args.ci is not None)
    cipos95 = _ci(v, 'CIPOS95', args.ci95 is not None)
    ciend95 = _ci(v, 'CIEND95', args.ci95 is not None)

    start = int(v.start)
    start += cipos[0] if args.ci == 'out' else cipos[1] if args.ci == 'in' else 0
    start += cipos95[0] if args.ci95 == 'out' else cipos95[1] if args.ci95 == 'in' else 0
    start -= args.emb or 0
    start += args.red or 0

    info_end = _first_int(v.INFO.get('END'))
    svlen = _first_int(v.INFO.get('SVLEN'))
    end = info_end if info_end is not None else int(v.POS)
    if info_end is None and svlen is not None:
        end = start + abs(svlen) + 1
    end += ciend[1] if args.ci == 'out' else ciend[0] if args.ci == 'in' else 0
    end += ciend95[1] if args.ci95 == 'out' else ciend95[0] if args.ci95 == 'in' else 0
    end += args.emb or 0
    end -= args.red or 0

    if (args.ci == 'in' or args.ci95 == 'in') and end < start:
        end = start + 1
    if args.red is not None and end < start:
        end = start + 1
    if start > end:
        logger.warning('%s: start > end (%d > %d); flipped', v.ID, start, end)
        start, end = end, start
    if start == end:
        logger.warning('%s: start == end (%d); end moved to %d', v.ID, start, end + 1)
        end += 1
    return start, end


def load_query(path: str, args: Namespace, threads: int) -> Query:
    """One pass over the VCF collecting adjusted coordinates, SVTYPE/SVLEN and het / hom-alt sample lists."""
    logger.info('Gathering SV coordinates from VCF file: %s', path)
    vcf = VCF(path, threads=threads, gts012=True)  # gts012: 0 hom-ref, 1 het, 2 hom-alt, 3 unknown
    samples = np.array(vcf.samples)

    chrom: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    sv_ids: list[str] = []
    svlens: list[int] = []
    svtypes: list[str] = []
    hets: list[str] = []
    homalts: list[str] = []

    for v in vcf:
        start, end = record_coordinates(v, args)
        svlen = _first_int(v.INFO.get('SVLEN'))
        chrom.append(strip_chr(v.CHROM))
        starts.append(start)
        ends.append(end)
        sv_ids.append(v.ID if v.ID is not None else f'{v.CHROM}:{v.POS}')
        svlens.append(abs(svlen) if svlen is not None else end - start)
        svtypes.append(v.INFO.get('SVTYPE') or '')
        gt_types = v.gt_types
        het = samples[gt_types == 1]
        homalt = samples[gt_types == 2]
        hets.append(','.join(sorted(het)) if len(het) else 'None')
        homalts.append(','.join(sorted(homalt)) if len(homalt) else 'None')
    vcf.close()

    return Query(
        intervals=Intervals.from_columns(chrom, starts, ends),
        sv_id=np.array(sv_ids, dtype=object),
        svlen=np.array(svlens, dtype=np.int64),
        svtype=np.array(svtypes, dtype=object),
        het_samples=hets,
        homalt_samples=homalts,
    )
