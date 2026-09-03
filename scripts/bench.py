"""Benchmark `svafotate annotate` implementations on a realistic workload.

Builds (once) a synthetic query VCF by sampling records from a reference BED and jittering their
coordinates, then runs each given executable on it and reports wall time and peak memory.

Example:
    python scripts/bench.py --bed ~/code/talos/large_files/SVAFotate_reduced_gnomAD.bed.gz --n 10000 \\
        --exe .venv/bin/svafotate --exe .venv-legacy/bin/svafotate
"""

from __future__ import annotations

import argparse
import gzip
import random
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / 'bench'


def build_query(bed: Path, n: int, seed: int) -> Path:
    out = BENCH / f'query_{n}.vcf'
    if out.exists():
        return out
    BENCH.mkdir(exist_ok=True)
    rng = random.Random(seed)
    rows: list[tuple[str, int, int, str]] = []
    opener = gzip.open if str(bed).endswith('.gz') else open
    with opener(bed, 'rt') as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.split('\t', 6)
            if f[4] in ('DEL', 'DUP', 'INV', 'INS'):
                rows.append((f[0], int(f[1]), int(f[2]), f[4]))
    sample = rng.sample(rows, min(n, len(rows)))
    sample.sort(key=lambda r: (r[0], r[1]))
    with out.open('w') as fh:
        fh.write('##fileformat=VCFv4.2\n')
        for contig in sorted({r[0] for r in sample}):
            fh.write(f'##contig=<ID=chr{contig}>\n')
        fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n')
        fh.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n')
        fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n')
        for i, (chrom, start, end, svtype) in enumerate(sample):
            length = max(end - start, 1)
            jitter = int(length * rng.uniform(-0.1, 0.1))
            pos = max(start + 1 + jitter, 1)
            stop = max(end + int(length * rng.uniform(-0.1, 0.1)), pos + 1)
            gts = '\t'.join(rng.choice(['0/0', '0/1', '1/1']) for _ in range(3))
            fh.write(
                f'chr{chrom}\t{pos}\tbench_{i}\tN\t<{svtype}>\t.\tPASS\tSVTYPE={svtype};SVLEN={stop - pos + 1};END={stop}\tGT\t{gts}\n'
            )
    return out


def run(exe: str, query: Path, bed: Path, extra: list[str]) -> tuple[float, float, int]:
    out = BENCH / f'out_{re.sub(r"[^A-Za-z0-9]+", "_", exe)}.vcf'
    exe_path = Path(exe)
    cmd = [
        str(exe_path.resolve()) if exe_path.exists() else exe,
        'annotate',
        '-v',
        str(query),
        '-o',
        str(out),
        '-b',
        str(bed),
        *extra,
    ]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=BENCH, capture_output=True, text=True, check=False)
    wall = time.perf_counter() - t0
    # ru_maxrss is the max over all children so far, so run the heaviest candidate last if you want them all
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if proc.returncode != 0:
        print(proc.stderr[-2000:], file=sys.stderr)
    divisor = 1024**2 if sys.platform == 'darwin' else 1024  # bytes on macOS, KiB on Linux
    return wall, rss / divisor, proc.returncode


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bed', required=True, type=Path)
    p.add_argument('--n', type=int, default=10000)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--exe', action='append', required=True, help='svafotate executable (repeatable)')
    p.add_argument('extra', nargs='*', default=['-s', 'gnomAD', '-f', '0.5', '-a', 'best', '--ins'])
    ns = p.parse_args()
    query = build_query(ns.bed, ns.n, ns.seed)
    print(f'query: {query} ({ns.n} records)  reference: {ns.bed}  args: {" ".join(ns.extra)}')
    print(f'{"executable":60s} {"wall_s":>8s} {"peak_MB":>8s}  rc')
    for exe in ns.exe:
        wall, mb, rc = run(exe, query, ns.bed, ns.extra)
        print(f'{exe:60s} {wall:8.1f} {mb:8.0f}  {rc}')


if __name__ == '__main__':
    main()
