"""Interval operations on half-open ``[start, end)`` genomic intervals, implemented in numpy.

This module replaces the three pyranges operations SVAFotate relied on (``join``, ``coverage``,
``subtract``). Two intervals on the same chromosome overlap when ``a.start < b.end and b.start < a.end``,
which is the same strict test the NCLS backend of pyranges 0.x applied.

The join uses length-binned binary search: reference intervals are grouped by ``floor(log2(length))``,
so within a bin every interval is shorter than ``L = 2 ** (bin + 1)`` and any interval overlapping a query
must start in ``(q.start - L, q.end)``. Two ``searchsorted`` calls per bin give a tight candidate range
without an interval tree, and the whole thing is O((n + m) log m + candidates).
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class Intervals:
    """A column-oriented set of intervals. ``chrom`` is an object array of contig names."""

    chrom: NDArray[np.object_]
    start: NDArray[np.int64]
    end: NDArray[np.int64]

    def __post_init__(self) -> None:
        if not (len(self.chrom) == len(self.start) == len(self.end)):
            raise ValueError('chrom, start and end must have the same length')
        if len(self.start) and bool(np.any(self.end < self.start)):
            raise ValueError('every interval must satisfy start <= end')

    def __len__(self) -> int:
        return len(self.start)

    @classmethod
    def from_columns(cls, chrom: ArrayLike, start: ArrayLike, end: ArrayLike) -> 'Intervals':
        return cls(
            np.asarray(chrom, dtype=object),
            np.asarray(start, dtype=np.int64),
            np.asarray(end, dtype=np.int64),
        )

    def take(self, idx: NDArray[np.intp]) -> 'Intervals':
        return Intervals(self.chrom[idx], self.start[idx], self.end[idx])

    @property
    def length(self) -> NDArray[np.int64]:
        return self.end - self.start


def _chrom_groups(chrom: NDArray[np.object_]) -> dict[str, NDArray[np.intp]]:
    """Map each contig to the (ascending) indices of intervals on it."""
    if len(chrom) == 0:
        return {}
    names, inverse = np.unique(chrom.astype(str), return_inverse=True)
    order = np.argsort(inverse, kind='stable')
    bounds = np.searchsorted(inverse[order], np.arange(len(names) + 1))
    return {str(name): order[bounds[i] : bounds[i + 1]] for i, name in enumerate(names)}


def _expand_ranges(lo: NDArray[np.intp], hi: NDArray[np.intp]) -> tuple[NDArray[np.intp], NDArray[np.intp]]:
    """For per-query candidate ranges ``[lo, hi)`` return (query position, candidate position) pairs."""
    counts = np.maximum(hi - lo, 0)
    total = int(counts.sum())
    if total == 0:
        empty = np.empty(0, dtype=np.intp)
        return empty, empty
    q_pos = np.repeat(np.arange(len(lo), dtype=np.intp), counts)
    offsets = np.arange(total, dtype=np.intp) - np.repeat(np.cumsum(counts) - counts, counts)
    return q_pos, np.repeat(lo, counts) + offsets


def overlap_pairs(query: Intervals, ref: Intervals) -> tuple[NDArray[np.intp], NDArray[np.intp]]:
    """Every (query index, reference index) pair whose intervals overlap on the same contig.

    Pairs are returned sorted by query index, then reference start ascending, reference end descending
    (containing intervals before the intervals they contain) and reference index. This is the traversal
    order of the NCLS backend pyranges 0.x used, so list-valued annotations keep their historical order.
    """
    q_groups = _chrom_groups(query.chrom)
    r_groups = _chrom_groups(ref.chrom)
    q_out: list[NDArray[np.intp]] = []
    r_out: list[NDArray[np.intp]] = []

    for contig, qi in q_groups.items():
        ri = r_groups.get(contig)
        if ri is None:
            continue
        qs, qe = query.start[qi], query.end[qi]
        r_len = np.maximum(ref.end[ri] - ref.start[ri], 1)
        bins = np.floor(np.log2(r_len)).astype(np.int64)

        for b in np.unique(bins):
            rib = ri[bins == b]
            rib = rib[np.argsort(ref.start[rib], kind='stable')]
            rs, re_ = ref.start[rib], ref.end[rib]
            max_len = np.int64(2) ** (b + 1)
            lo = np.searchsorted(rs, qs - max_len, side='right')
            hi = np.searchsorted(rs, qe, side='left')
            q_pos, r_pos = _expand_ranges(lo, hi)
            keep = re_[r_pos] > qs[q_pos]
            q_out.append(qi[q_pos[keep]])
            r_out.append(rib[r_pos[keep]])

    if not q_out:
        empty = np.empty(0, dtype=np.intp)
        return empty, empty
    q_idx = np.concatenate(q_out)
    r_idx = np.concatenate(r_out)
    order = np.lexsort((r_idx, -ref.end[r_idx], ref.start[r_idx], q_idx))
    return q_idx[order], r_idx[order]


def overlap_bp(query: Intervals, ref: Intervals, q_idx: NDArray[np.intp], r_idx: NDArray[np.intp]) -> NDArray[np.int64]:
    """Number of overlapping bases for each (query, reference) pair."""
    return np.minimum(query.end[q_idx], ref.end[r_idx]) - np.maximum(query.start[q_idx], ref.start[r_idx])


def merge(ref: Intervals) -> Intervals:
    """Union of the intervals per contig: disjoint, sorted by contig then start. Touching intervals fuse."""
    chroms: list[NDArray[np.object_]] = []
    starts: list[NDArray[np.int64]] = []
    ends: list[NDArray[np.int64]] = []
    for contig, group in _chrom_groups(ref.chrom).items():
        ri = group[np.argsort(ref.start[group], kind='stable')]
        s, e = ref.start[ri], ref.end[ri]
        running_end = np.maximum.accumulate(e)
        new_cluster = np.ones(len(ri), dtype=bool)
        new_cluster[1:] = s[1:] > running_end[:-1]
        cluster_start = np.flatnonzero(new_cluster)
        cluster_end = np.append(cluster_start[1:], len(ri))
        starts.append(s[cluster_start])
        ends.append(running_end[cluster_end - 1])
        chroms.append(np.full(len(cluster_start), contig, dtype=object))
    if not chroms:
        return Intervals.from_columns([], [], [])
    return Intervals(np.concatenate(chroms), np.concatenate(starts), np.concatenate(ends))


def covered_fraction(query: Intervals, ref: Intervals) -> NDArray[np.float64]:
    """Fraction of each query interval covered by the union of the reference intervals (0.0 - 1.0)."""
    merged = merge(ref)
    q_idx, r_idx = overlap_pairs(query, merged)
    bp = overlap_bp(query, merged, q_idx, r_idx).astype(np.float64)
    covered = np.bincount(q_idx, weights=bp, minlength=len(query))
    length = query.length.astype(np.float64)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(length > 0, covered / length, 0.0)


def subtract(query: Intervals, ref: Intervals) -> tuple[NDArray[np.intp], NDArray[np.int64], NDArray[np.int64]]:
    """The parts of each query interval not covered by any reference interval.

    Returns ``(query index, start, end)`` for every non-empty residual piece, ordered by query index then
    start. A query with no overlaps is returned whole; a fully covered query contributes nothing.
    """
    merged = merge(ref)
    q_idx, r_idx = overlap_pairs(query, merged)  # sorted by q_idx then ref start; merged refs are disjoint
    n = len(query)

    if len(q_idx) == 0:
        keep = np.flatnonzero(query.length > 0)
        return keep.astype(np.intp), query.start[keep], query.end[keep]

    # first pair of each query, and whether a pair is the first / last for its query
    first_of_q = np.ones(len(q_idx), dtype=bool)
    first_of_q[1:] = q_idx[1:] != q_idx[:-1]
    last_of_q = np.ones(len(q_idx), dtype=bool)
    last_of_q[:-1] = q_idx[1:] != q_idx[:-1]

    # gap before each reference piece: from the previous ref end (or the query start) to this ref start
    prev_end = np.empty(len(q_idx), dtype=np.int64)
    prev_end[first_of_q] = query.start[q_idx[first_of_q]]
    prev_end[~first_of_q] = merged.end[r_idx[np.flatnonzero(~first_of_q) - 1]]
    gap_q = q_idx
    gap_s = prev_end
    gap_e = merged.start[r_idx]

    # tail after the last reference piece of each query
    tail_q = q_idx[last_of_q]
    tail_s = merged.end[r_idx[last_of_q]]
    tail_e = query.end[tail_q]

    # queries with no overlapping reference at all
    touched = np.zeros(n, dtype=bool)
    touched[q_idx] = True
    whole_q = np.flatnonzero(~touched).astype(np.intp)

    out_q = np.concatenate([gap_q, tail_q, whole_q])
    out_s = np.concatenate([gap_s, tail_s, query.start[whole_q]])
    out_e = np.concatenate([gap_e, tail_e, query.end[whole_q]])

    # clip to the query and drop empties
    out_s = np.maximum(out_s, query.start[out_q])
    out_e = np.minimum(out_e, query.end[out_q])
    keep = out_s < out_e
    out_q, out_s, out_e = out_q[keep], out_s[keep], out_e[keep]
    order = np.lexsort((out_s, out_q))
    return out_q[order], out_s[order], out_e[order]
