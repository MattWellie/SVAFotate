"""Randomised tests for the numpy interval engine against brute-force references (and pyranges 0.x if present)."""

import numpy as np
import pytest

from svafotate_cpg.intervals import Intervals, covered_fraction, merge, overlap_bp, overlap_pairs, subtract

CHROMS = ['1', '2', 'X']
SPAN = 300  # small coordinate space so coverage arrays are cheap and overlaps are dense


def random_intervals(rng: np.random.Generator, n: int, allow_empty: bool = True, huge: bool = False) -> Intervals:
    chrom = rng.choice(CHROMS, size=n)
    start = rng.integers(0, SPAN, size=n)
    length = rng.geometric(0.05, size=n) - (1 if allow_empty else 0)
    if huge:  # sprinkle a few intervals spanning most of the space, the pathological case for naive sweeps
        big = rng.random(n) < 0.05
        start[big] = rng.integers(0, 5, size=big.sum())
        length[big] = SPAN
    end = start + length
    return Intervals.from_columns(chrom, start, end)


def brute_pairs(q: Intervals, r: Intervals) -> set[tuple[int, int]]:
    return {
        (i, j)
        for i in range(len(q))
        for j in range(len(r))
        if q.chrom[i] == r.chrom[j] and q.start[i] < r.end[j] and r.start[j] < q.end[i]
    }


def brute_mask(r: Intervals, chrom: str) -> np.ndarray:
    mask = np.zeros(SPAN * 2, dtype=bool)
    for j in range(len(r)):
        if r.chrom[j] == chrom:
            mask[r.start[j] : r.end[j]] = True
    return mask


@pytest.mark.parametrize('seed', range(25))
def test_overlap_pairs_matches_brute_force(seed: int):
    rng = np.random.default_rng(seed)
    q = random_intervals(rng, rng.integers(0, 40), huge=True)
    r = random_intervals(rng, rng.integers(0, 120), huge=True)
    qi, ri = overlap_pairs(q, r)
    assert set(zip(qi.tolist(), ri.tolist(), strict=True)) == brute_pairs(q, r)
    assert len(set(zip(qi.tolist(), ri.tolist(), strict=True))) == len(qi), 'duplicate pairs'
    # ordering contract: by query, then ref start ascending, ref end descending, ref index
    keys = list(zip(qi.tolist(), r.start[ri].tolist(), (-r.end[ri]).tolist(), ri.tolist(), strict=True))
    assert keys == sorted(keys)
    assert np.all(overlap_bp(q, r, qi, ri) >= 0)


def test_zero_length_intervals_behave_like_points():
    q = Intervals.from_columns(['1', '1', '1'], [10, 10, 10], [10, 10, 10])
    r = Intervals.from_columns(['1', '1', '1'], [5, 10, 11], [15, 10, 20])
    qi, ri = overlap_pairs(q, r)
    # a point strictly inside [5,15) overlaps; a point equal to another point does not; [11,20) starts after
    assert sorted(set(ri.tolist())) == [0]
    assert len(qi) == 3


@pytest.mark.parametrize('seed', range(25))
def test_merge_is_a_disjoint_sorted_union(seed: int):
    rng = np.random.default_rng(seed)
    r = random_intervals(rng, rng.integers(0, 80))
    m = merge(r)
    for chrom in CHROMS:
        sel = np.flatnonzero(m.chrom == chrom)
        s, e = m.start[sel], m.end[sel]
        assert np.all(s[1:] > e[:-1]), 'merged intervals must be disjoint and non-touching'
        assert np.array_equal(brute_mask(m.take(sel), chrom), brute_mask(r, chrom))


@pytest.mark.parametrize('seed', range(25))
def test_covered_fraction_matches_brute_force(seed: int):
    rng = np.random.default_rng(seed)
    q = random_intervals(rng, rng.integers(0, 40))
    r = random_intervals(rng, rng.integers(0, 120), huge=True)
    got = covered_fraction(q, r)
    masks = {c: brute_mask(r, c) for c in CHROMS}
    for i in range(len(q)):
        length = q.end[i] - q.start[i]
        expected = masks[q.chrom[i]][q.start[i] : q.end[i]].sum() / length if length else 0.0
        assert got[i] == pytest.approx(expected)


@pytest.mark.parametrize('seed', range(25))
def test_subtract_matches_brute_force(seed: int):
    rng = np.random.default_rng(seed)
    q = random_intervals(rng, rng.integers(0, 40), allow_empty=False)
    r = random_intervals(rng, rng.integers(0, 120), huge=True)
    qi, s, e = subtract(q, r)
    assert np.all(s < e)
    keys = list(zip(qi.tolist(), s.tolist(), strict=True))
    assert keys == sorted(keys)
    masks = {c: brute_mask(r, c) for c in CHROMS}
    for i in range(len(q)):
        got = np.zeros(SPAN * 2, dtype=bool)
        for k in np.flatnonzero(qi == i):
            assert q.start[i] <= s[k] < e[k] <= q.end[i]
            got[s[k] : e[k]] = True
        expected = np.zeros(SPAN * 2, dtype=bool)
        expected[q.start[i] : q.end[i]] = ~masks[q.chrom[i]][q.start[i] : q.end[i]]
        assert np.array_equal(got, expected)


def test_empty_inputs():
    empty = Intervals.from_columns([], [], [])
    q = Intervals.from_columns(['1'], [0], [10])
    assert overlap_pairs(q, empty)[0].size == 0
    assert overlap_pairs(empty, q)[0].size == 0
    assert covered_fraction(q, empty).tolist() == [0.0]
    qi, s, e = subtract(q, empty)
    assert (qi.tolist(), s.tolist(), e.tolist()) == ([0], [0], [10])
    assert len(merge(empty)) == 0


def test_rejects_inverted_intervals():
    with pytest.raises(ValueError, match='start <= end'):
        Intervals.from_columns(['1'], [10], [5])
