"""Overlap matching between query SVs and one reference source.

For every overlapping (query, reference) pair the reciprocal overlap fractions and their product, the
overlap fraction product (OFP), are computed once. Everything downstream (the ``-f`` filter, per-source
counts, best matches, mismatches) is a mask or a group-by over those arrays.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .intervals import overlap_pairs
from .query import Query
from .reference import SourceData


@dataclass
class SourceMatches:
    """All overlapping pairs between the query set and one source, in query order."""

    source: str
    q_idx: NDArray[np.intp]
    r_idx: NDArray[np.intp]
    fraction_q: NDArray[np.float64]  # overlap / query size
    fraction_r: NDArray[np.float64]  # overlap / reference size
    same_type: NDArray[np.bool_]  # query SVTYPE == reference SVTYPE

    @property
    def ofp(self) -> NDArray[np.float64]:
        return self.fraction_q * self.fraction_r

    def __len__(self) -> int:
        return len(self.q_idx)


def find_matches(query: Query, source: SourceData, conv_ins: bool) -> SourceMatches:
    """Overlap search plus overlap fractions.

    With ``--ins``, insertion queries (and the reference records they overlap) have their END extended
    to START + SVLEN before the fractions are computed. The overlap search itself always uses the
    coordinates as given, exactly as upstream did.
    """
    q_idx, r_idx = overlap_pairs(query.intervals, source.intervals)
    s1 = query.intervals.start[q_idx]
    e1 = query.intervals.end[q_idx]
    s2 = source.intervals.start[r_idx]
    e2 = source.intervals.end[r_idx]

    if conv_ins:
        ins = query.svtype[q_idx] == 'INS'
        grow = ins & ((e1 - s1) < query.svlen[q_idx])
        e1 = np.where(grow, s1 + query.svlen[q_idx], e1)
        grow_r = ins & ((e2 - s2) < source.svlen[r_idx])
        e2 = np.where(grow_r, s2 + source.svlen[r_idx], e2)

    overlap = (np.minimum(e1, e2) - np.maximum(s1, s2)).astype(np.float64)
    size_q = e1 - s1
    size_r = e2 - s2
    fraction_q = np.divide(overlap, size_q, out=np.zeros(len(overlap)), where=size_q != 0)
    fraction_r = np.divide(overlap, size_r, out=np.zeros(len(overlap)), where=size_r != 0)
    same_type = query.svtype[q_idx] == source.svtype[r_idx]
    return SourceMatches(source.name, q_idx, r_idx, fraction_q, fraction_r, same_type)


def passes_overlap_filter(m: SourceMatches, query: Query, minf: float, conv_ins: bool) -> NDArray[np.bool_]:
    """The ``-f`` reciprocal overlap filter. Insertions bypass it unless ``--ins`` widened them."""
    reciprocal = (m.fraction_q >= minf) & (m.fraction_r >= minf)
    ins_bypass = (query.svtype[m.q_idx] == 'INS') & (not conv_ins)
    return reciprocal | ins_bypass


def group_by_query(q_idx: NDArray[np.intp], r_idx: NDArray[np.intp]) -> dict[int, NDArray[np.intp]]:
    """Reference rows per query index, preserving pair order (pairs are already sorted by query)."""
    if len(q_idx) == 0:
        return {}
    starts = np.flatnonzero(np.r_[True, q_idx[1:] != q_idx[:-1]])
    ends = np.r_[starts[1:], len(q_idx)]
    return {int(q_idx[a]): r_idx[a:b] for a, b in zip(starts.tolist(), ends.tolist(), strict=True)}


def best_per_query(m: SourceMatches, mask: NDArray[np.bool_], source: SourceData) -> dict[int, tuple[int, float]]:
    """The best-matching reference row (and its OFP) for each query, among the pairs selected by ``mask``.

    Best means the highest OFP. Ties are broken by the highest AF, then the highest HomAlt count, then
    the first pair in reference order. Upstream intended the same tiebreak but compared the wrong columns.
    """
    sel = np.flatnonzero(mask)
    if len(sel) == 0:
        return {}
    ofp = m.ofp[sel]
    af = _numeric_or_zero(source, 'AF', m.r_idx[sel])
    homalt = _numeric_or_zero(source, 'HomAlt', m.r_idx[sel])
    # lexsort: last key is primary. Descending on ofp, af, homalt; ascending on position for stability.
    order = np.lexsort((np.arange(len(sel)), -homalt, -af, -ofp, m.q_idx[sel]))
    q_sorted = m.q_idx[sel][order]
    first = np.r_[True, q_sorted[1:] != q_sorted[:-1]]
    winners = order[first]
    return {int(m.q_idx[sel][w]): (int(m.r_idx[sel][w]), float(ofp[w])) for w in winners}


def _numeric_or_zero(source: SourceData, column: str, rows: NDArray[np.intp]) -> NDArray[np.float64]:
    col = source.columns.get(column)
    if col is None or col.dtype.kind != 'f':
        return np.zeros(len(rows))
    return np.nan_to_num(col[rows], nan=0.0)
