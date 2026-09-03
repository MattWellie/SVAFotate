# Contributing

## Development setup

```bash
uv sync --group dev          # creates .venv with the package (editable) and dev tools
uv run pytest                # golden regression suite + interval engine tests
uv run ruff check . && uv run ruff format --check .
uv run mypy src
pre-commit install
```

Dependencies must install from wheels on every supported Python (3.10-3.14). CI runs
`uv sync --no-build`, so adding a dependency that needs a compiler fails the build.

## Tests

`tests/test_golden.py` runs every case in `tests/golden_cases.py` through the CLI and compares the
result with the output of the **legacy oracle**: upstream SVAFotate at commit `30b5004`, running on
pandas 2.x and pyranges 0.1.4. Records are matched by ID and INFO values compared with a float
tolerance. Two deliberate deviations are allowed for (see `tests/vcf_compare.py`): the header type of
`Best_*_Mismatch_OFP`, and the member order of the `*_Matches` / `*_Mismatches` lists.

`tests/test_intervals.py` checks the numpy interval engine against brute-force implementations, and
against pyranges 0.x when it happens to be importable.

### Fixtures

| File | Origin |
|---|---|
| `tests/data/joint_sv.vcf.gz` | Talos `nextflow/inputs/joint_sv.vcf.bgz`: 8 adversarial SVs on chr1 |
| `tests/data/ref_gnomad_3rows.bed` | Talos `test/input/svafotate_slice.bed`: 3 gnomAD rows, full 177-column header |
| `tests/data/ref_multi_source.bed.gz` | All rows of `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz` within 50 kb of a fixture SV, all four sources |
| `tests/data/extended.vcf`, `tests/data/targets.bed` | Synthetic, written by `scripts/make_fixtures.py` |

### Regenerating the golden outputs

Only do this when a change is *meant* to alter output, and say so in the changelog.

```bash
uv venv --python 3.12 .venv-legacy
uv pip install --python .venv-legacy/bin/python 'pandas<3' 'numpy<2.3' pyranges==0.1.4 ncls 'cyvcf2>=0.31' \
    'git+https://github.com/fakedrtom/SVAFotate.git@30b5004a0f4d26959c6b9a82f165651585293626'
python scripts/generate_golden.py            # all cases, or name individual cases
```

`ncls` has no wheels, so this needs a C compiler. The oracle environment is gitignored.

To add a case, append to `CASES` in `tests/golden_cases.py`, run the generator for that case, and commit
the new `tests/golden/<case>/` directory. Upstream cannot process a record without `SVTYPE`, or a
`-ci`/`-ci95` run where any record lacks `CIPOS`/`CIEND`; keep such records out of oracle fixtures.

## Releasing

1. Update `CHANGELOG.md`, move the Unreleased section under the new version.
2. `uv run bump-my-version bump <major|minor|patch>` (or edit `pyproject.toml` and
   `src/svafotate_cpg/__init__.py` by hand; they must agree).
3. Tag `vX.Y.Z` and push the tag. The `Release` workflow builds and publishes to PyPI via trusted publishing.
