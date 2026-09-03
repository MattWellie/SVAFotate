"""Golden regression tests: the current code must reproduce the legacy oracle's output for every case."""

import sys
from pathlib import Path

import pytest
from golden_cases import CASES
from vcf_compare import compare_beds, compare_vcfs

from svafotate_cpg.main import main


@pytest.mark.parametrize('case', list(CASES))
def test_matches_golden(case: str, data_dir: Path, golden_dir: Path, tmp_path: Path, monkeypatch):
    vcf, args = CASES[case]
    expected_dir = golden_dir / case
    if not (expected_dir / 'out.vcf').exists():
        pytest.skip(f'no golden output for {case} (oracle failed?)')

    resolved = [str(data_dir / a) if (data_dir / a).exists() else a for a in args]
    out = tmp_path / 'out.vcf'
    monkeypatch.chdir(tmp_path)  # uniques.bed is written to the CWD
    monkeypatch.setattr(sys, 'argv', ['svafotate', 'annotate', '-v', str(data_dir / vcf), '-o', str(out), *resolved])
    main()

    diffs = compare_vcfs(out, expected_dir / 'out.vcf')
    diffs += compare_beds(tmp_path / 'uniques.bed', expected_dir / 'uniques.bed')
    assert not diffs, f'{len(diffs)} difference(s):\n' + '\n'.join(diffs)
