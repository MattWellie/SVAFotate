"""Regenerate the golden regression outputs from the legacy oracle.

The oracle is upstream SVAFotate (fakedrtom/SVAFotate @ 30b5004) installed with pandas<3 and
pyranges 0.1.4 in ``.venv-legacy`` (see CONTRIBUTING.md for the exact recipe). This script runs every
case in ``tests/golden_cases.py`` through it and stores ``out.vcf`` (and ``uniques.bed`` where produced)
under ``tests/golden/<case>/``.

Usage: ``python scripts/generate_golden.py [--oracle PATH] [case ...]``
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'tests'))
from golden_cases import CASES  # noqa: E402

DATA = ROOT / 'tests' / 'data'
GOLDEN = ROOT / 'tests' / 'golden'


def run_case(oracle: Path, name: str) -> None:
    vcf, args = CASES[name]
    out_dir = GOLDEN / name
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    resolved = [str(DATA / a) if (DATA / a).exists() else a for a in args]
    cmd = [str(oracle), 'annotate', '-v', str(DATA / vcf), '-o', 'out.vcf', *resolved]
    (out_dir / 'command.txt').write_text(' '.join(cmd) + '\n')
    proc = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True, check=False)
    (out_dir / 'oracle_stdout.txt').write_text(proc.stdout)
    if proc.returncode != 0:
        (out_dir / 'oracle_stderr.txt').write_text(proc.stderr)
        print(f'FAILED  {name}\n{proc.stderr[-2000:]}')
    else:
        print(f'ok      {name}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--oracle', default=str(ROOT / '.venv-legacy' / 'bin' / 'svafotate'))
    parser.add_argument('cases', nargs='*', default=list(CASES))
    ns = parser.parse_args()
    for name in ns.cases:
        run_case(Path(ns.oracle), name)


if __name__ == '__main__':
    main()
