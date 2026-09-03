"""Helpers for comparing an annotated VCF against a golden file.

Comparison is semantic rather than textual: records are matched by ID and their INFO fields compared
with a float tolerance, and INFO header lines are compared on ID/Number/Type. Descriptions are ignored
because they embed the reference BED path.
"""

import math
from pathlib import Path

from cyvcf2 import VCF

FLOAT_TOL = 1e-6

# header Type corrections this fork makes deliberately (upstream declared a float as String)
TYPE_CORRECTIONS = {'_Mismatch_OFP': ('String', 'Float')}

# comma-separated list annotations whose member order upstream inherited from NCLS traversal order; the
# members must match, the order is an implementation detail (this fork orders by reference start)
UNORDERED_LIST_SUFFIXES = ('_Matches', '_Mismatches')


def info_headers(vcf: VCF) -> dict[str, tuple[str, str]]:
    out = {}
    for h in vcf.header_iter():
        if h['HeaderType'] == 'INFO':
            out[h['ID']] = (h['Number'], h['Type'])
    return out


def records(vcf: VCF) -> dict[str, dict]:
    out = {}
    for v in vcf:
        assert v.ID not in out, f'duplicate ID {v.ID}'
        out[v.ID] = {'CHROM': v.CHROM, 'POS': v.POS, 'INFO': dict(v.INFO)}
    return out


def _close(a, b, key: str = '') -> bool:
    if key.endswith(UNORDERED_LIST_SUFFIXES) and isinstance(a, str) and isinstance(b, str):
        return sorted(a.split(',')) == sorted(b.split(','))
    if isinstance(a, float) or isinstance(b, float):
        try:
            return math.isclose(float(a), float(b), rel_tol=FLOAT_TOL, abs_tol=FLOAT_TOL)
        except (TypeError, ValueError):
            return False
    return a == b


def compare_vcfs(actual: Path, expected: Path) -> list[str]:
    """Return a list of human-readable differences; empty means identical."""
    diffs: list[str] = []
    va, ve = VCF(str(actual)), VCF(str(expected))

    ha, he = info_headers(va), info_headers(ve)
    for key in sorted(set(ha) | set(he)):
        if key not in ha:
            diffs.append(f'header INFO {key} missing from actual')
        elif key not in he:
            diffs.append(f'header INFO {key} unexpected in actual')
        elif ha[key] != he[key]:
            corrected = any(key.endswith(s) and (he[key][1], ha[key][1]) == c for s, c in TYPE_CORRECTIONS.items())
            if not corrected:
                diffs.append(f'header INFO {key}: actual {ha[key]} expected {he[key]}')

    ra, re_ = records(va), records(ve)
    if list(ra) != list(re_):
        diffs.append(f'record order/IDs differ: {list(ra)} vs {list(re_)}')
    for vid in re_:
        if vid not in ra:
            continue
        a, e = ra[vid], re_[vid]
        for field in ('CHROM', 'POS'):
            if a[field] != e[field]:
                diffs.append(f'{vid}: {field} {a[field]} != {e[field]}')
        for key in sorted(set(a['INFO']) | set(e['INFO'])):
            if key not in a['INFO']:
                diffs.append(f'{vid}: INFO {key} missing (expected {e["INFO"][key]!r})')
            elif key not in e['INFO']:
                diffs.append(f'{vid}: INFO {key} unexpected ({a["INFO"][key]!r})')
            elif not _close(a['INFO'][key], e['INFO'][key], key):
                diffs.append(f'{vid}: INFO {key} {a["INFO"][key]!r} != {e["INFO"][key]!r}')
    va.close()
    ve.close()
    return diffs


def compare_beds(actual: Path, expected: Path) -> list[str]:
    a = actual.read_text().splitlines() if actual.exists() else None
    e = expected.read_text().splitlines() if expected.exists() else None
    if a is None and e is None:
        return []
    if a is None:
        return [f'{actual.name} not produced']
    if e is None:
        return [f'{actual.name} produced but not expected']
    if a[0] != e[0]:
        return [f'header differs: {a[0]!r} vs {e[0]!r}']
    # upstream never populated HomAlt_Samples (it tested gt_types == 2, which is UNKNOWN in cyvcf2's default
    # encoding); this fork fixes that, so the column is excluded from the golden comparison
    homalt_col = a[0].split('\t').index('HomAlt_Samples')

    def drop(rows: list[str]) -> list[str]:
        return ['\t'.join(f for i, f in enumerate(r.split('\t')) if i != homalt_col) for r in rows]

    a, e = [a[0], *drop(a[1:])], [e[0], *drop(e[1:])]
    if sorted(a[1:]) != sorted(e[1:]):
        return [
            'rows differ:\n  actual:   '
            + '\n  actual:   '.join(a[1:])
            + '\n  expected: '
            + '\n  expected: '.join(e[1:])
        ]
    return []
