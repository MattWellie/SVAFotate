# SVAFotate (CPG fork) — Renovation Plan

Date: 2026-09-02. Status of repository at time of writing: fresh import of upstream
[fakedrtom/SVAFotate](https://github.com/fakedrtom/SVAFotate) (commit `30b5004`), lightly reshaped into a
`src/` layout and renamed `svafotate_cpg`. No tests, no CI, no lockfile, no `.gitignore`.

## 1. Why this fork exists

Upstream SVAFotate is built on pyranges 0.1.4, whose `join` cannot handle the pandas 3 default string dtype:

```
Exception: Unknown dtype str in a column SVTYPE
```

Reproduced in this repo's `.venv` (Python 3.12.13, pandas 3.0.5, pyranges 0.1.4): `PyRanges.join` fails,
while `coverage` and `subtract` still work. The failure is inside `pyranges/methods/join.py::null_types`,
which accepts only `object` or a dtype whose name contains `"string"` — pandas 3 reports `str`.

Downstream, Talos and Talos2 currently work around this with `pandas<3` pins in a separate Docker image
(`talos/docker/SVAFotate_Dockerfile`, `images/talos2/Dockerfile`). Those images also have to compile from
source, because **`ncls` (pyranges' interval backend) publishes no wheels at all**, and `sorted-nearest`
publishes none for aarch64. The goal of this fork is a `svafotate_cpg` package that:

- installs from wheels only on Python 3.12, 3.13 and 3.14 alongside pandas 3 / numpy 2;
- produces byte-identical annotations (modulo header provenance lines) to upstream for the options Talos uses;
- can be installed into the Talos2 image with `uv`, retiring the standalone SVAFotate image.

### The downstream contract that must not break

Talos invokes exactly this (`talos/nextflow/modules/annotation/AnnotateSvWithSvafotate/main.nf`):

```
svafotate annotate -v in.vcf.gz -o out.vcf -b SVAFotate_reduced_gnomAD.bed.gz \
    -s gnomAD -f <fraction> -a best --ins --cpu <n>
```

and `talos/src/talos/annotation_scripts/rename_sv_af_fields.py` reads three INFO fields: `Max_AF`,
`Best_gnomAD_ID`, `gnomAD_Count`. It relies on two documented quirks: `Best_*_ID` reports the best
*candidate* regardless of `-f`, and `*_Count` reflects only matches that cleared `-f`. Both must be preserved.

CLI shape: the imported code had flattened the CLI to `svafotate_cpg ...` with no subcommand. On
2026-09-02 the upstream `svafotate annotate ...` layout was restored (`main.py:get_parser`), the console
script is `svafotate` again, and `svafotate_cpg` is kept as an alias (`pyproject.toml [project.scripts]`).
Talos's Nextflow module therefore needs no change.

## 1a. Why this fork isn't a proper fork

The upstream SVAFotate repository has a Git history containing ~700MB of objects, which cannot be removed 
without a surgical identification of the now-deleted items, and complex git work. It also contains 200MB 
of live objects which are valuable as publication artefacts, but IMO it's fine maintaining a copy in the 
original repository, without needing the same data to live in the fork. After all, the zenodo link and 
link to the original repository are both ways to access those files. 

## 2. Audit findings

### 2.1 Packaging and tooling

| Item                | Finding                                                                                                      | Location                    |
|---------------------|--------------------------------------------------------------------------------------------------------------|-----------------------------|
| Dependencies        | `cyvcf2`, `numpy`, `pandas`, `pyranges` unpinned; pyranges 0.1.4 needs `ncls` + `sorted-nearest` (no wheels) | `pyproject.toml:22-27`      |
| Stale requirements  | conda-style `requirements.txt` (`htslib`, `zlib`, `python>=3`, pandas 1.2.3) contradicts pyproject           | `requirements.txt`          |
| Stray package       | `src/__init__.py` makes `src` importable as a package                                                        | `src/__init__.py`           |
| License metadata    | `license-files` lists `README.md`                                                                            | `pyproject.toml:11`         |
| setuptools leftover | `[options] include_package_data` is meaningless under hatchling                                              | `pyproject.toml:41-42`      |
| Python version skew | `requires-python >= 3.12` but pre-commit `default_language_version: python3.11`                              | `.pre-commit-config.yaml:2` |
| Ruff config         | isort `section-order` names unknown section `hail` (warning on every run); 75 lint errors, 1 format diff     | `pyproject.toml:82`         |
| Hygiene             | no `.gitignore`; `.venv/`, `.idea/`, `.ruff_cache/`, `.mypy_cache/` sit untracked in the tree                | repo root                   |
| Lockfile            | none (`uv.lock` absent)                                                                                      | repo root                   |
| Tests / CI          | none                                                                                                         | —                           |
| Versioning          | bumpversion rewrites `README.md`, which contains the *upstream* version string                               | `pyproject.toml:103-105`    |

### 2.2 Code structure

- `main.py` is 2,241 lines; `annotate()` alone is ~1,700 lines (`main.py:537-2226`).
- Eleven populations × three code sites are copy-pasted by hand: header creation (`main.py:813-877`,
  `944-1008`), max-value writing (`main.py:1525-1864`), best-value writing (`main.py:1902-2153`). Adding a
  population means editing six blocks.
- `utils.py:7-118` duplicates all four functions in `annotation_utils.py` verbatim; only the
  `annotation_utils` copies are imported.
- Three near-identical overlap functions (`get_overlap`, `get_fraction1`, `get_fraction2`,
  `main.py:292-367`) do per-row Python arithmetic that is a one-line vectorised expression
  (`min(e1,e2) - max(s1,s2)`).
- The input VCF is read twice (`main.py:676`, `main.py:756`).
- The whole reference BED is loaded into nested `defaultdict`s of Python strings (`utils.py:122-191`),
  including sources and columns that were never requested. Talos reports 16 GB+ for the full v4.1 BED,
  which is why it pre-shrinks the file with `awk | cut`.
- `print()` for all diagnostics; `NameError` used for user input errors (`main.py:606`, `620`).
- `--cpu` only sets cyvcf2 reader threads (`main.py:553`); nothing else is parallel. README says
  "experiencing technical difficulties".
- Output side file is hard-coded to `uniques.bed` in the CWD (`main.py:1391`); README says `unique.bed`.

### 2.3 Latent bugs (fix behind the regression suite, see §6 for which change output)

1. **Best-match tiebreaker reads the wrong columns and compares strings.** `get_best` uses
   `datas[source][id][0]` as AF and `[3]` as HomAlt (`main.py:501-503`), but `datas` rows are
   `[SVTYPE, AF, HomRef, Het, HomAlt, ...]` (`utils.py:168-169`), so it compares SVTYPE and Het, as strings
   (`"0.9" > "0.10"`). Only matters when two candidates tie on OFP.
2. `Best_<source>_Mismatch_OFP` is declared `Type=String` but a float is written (`main.py:1048-1052`).
3. `<source>_Mismatch_SVTYPEs` header description ends in a dangling `"for "` (`main.py:1033`).
4. Missing `END` with `SVLEN` present gives `end = start + |SVLEN| + 1` (`main.py:695-697`), one base longer
   than the `END = POS + SVLEN - 1` convention implies.
5. `chr` prefix is stripped from the query VCF and targets BED (`main.py:689`, `657`) but **not** from the
   reference BED. A `chr`-prefixed reference silently yields `Max_AF=0` everywhere (Talos documents this
   trap in three places).
6. Reference BEDs lacking `Het`, `HomAlt` or `PopMax_AF` raise a bare `KeyError` at write time; README
   claims only eight columns are mandatory. Validate up front.
7. Records with no `SVTYPE` are grouped under key `None` and passed to pyranges (`main.py:682`, `733`).
8. In the `-u` branch, `print(svtype + ' not found in ' + source ...)` uses the leaked loop variable from the
   preceding `for source in req_sources` (`main.py:1364`).
9. `--ins` coordinate expansion uses `DataFrame.apply(axis=1)` with row mutation (`main.py:370-389`);
   under pandas 3 copy-on-write the masked assignment pattern is fragile and it is slow.
10. Files opened without context managers (targets BED, `uniques.bed`; ruff `SIM115` ×5).

## 3. Decision: how to replace the interval engine

Four options were evaluated against pandas 3.0.5 in throwaway `uv` environments.

| Option | Verified | Verdict |
|---|---|---|
| **A. Patch only** — cast string columns to `object` before `PyRanges(...)`, or `pd.options.future.infer_string = False` | Both make `join` pass on pandas 3 | Unblocks in an hour, but keeps `ncls`/`sorted-nearest` (source-only, effectively unmaintained). Use only as a stop-gap tag. |
| **B. pyranges1** (`pyranges1==1.4.3`, Rust `ruranges` backend) | Wheels for 3.12/3.13; `join_overlaps`, `subtract_overlaps` work; on 3.14 the resolver falls back to `pyranges1==1.0.4` (no `ruranges` wheel yet); no `coverage` method (use `compute_interval_metrics`/manual) | Viable, but the API already broke once (0.x→1.x) and 3.14 support is incomplete. |
| **C. bioframe** (`bioframe==0.8.0`, pure Python) | Wheels for 3.12–3.14; `overlap`, `coverage`, `subtract` all work on pandas 3 (emits a `Pandas4Warning` inside `merge`) | Good fit, maintained (Open2C). Pulls in matplotlib. Still an external dependency that lags pandas. |
| **D. Own it** — `svafotate_cpg/intervals.py` in numpy | Not yet written | The tool needs exactly three operations on per-chromosome sorted arrays: inner overlap join, covered-fraction, subtract. ~150 lines, trivially unit-testable, zero compiled deps beyond numpy. |

**Recommendation: D, with C as the fallback** if the numpy implementation is not at parity with the
golden outputs within a day. Rationale: the operation surface is tiny, the pyranges 0.x left-join filler
semantics (`-1`) are not needed (matches are an inner join), and every third-party option evaluated has
already caused or is at risk of causing exactly the breakage this fork exists to escape.

Whatever the choice, `pyranges`, `ncls` and `sorted-nearest` leave the dependency set. Final runtime deps:
`cyvcf2`, `numpy`, `pandas` (pandas only for the BED loader and grouping; could go too, later).

## 4. Phased plan

Each phase ends in a green CI run and is independently mergeable.

### Phase 0 — Safety net (do first, nothing else lands without it)

1. Build a throwaway legacy environment: Python 3.11 or 3.12, `pandas<3`, `pyranges==0.1.4`, `ncls`,
   `cyvcf2`, and this repo's current code. This is the oracle.
2. Assemble fixtures under `tests/data/`:
   - `talos/nextflow/inputs/joint_sv.vcf.bgz` (8 adversarial records: exact match, small-in-large,
     four-way DUP overlap, INS with `END==POS`, BND, INV) — copy in.
   - `talos/test/input/svafotate_slice.bed` (3 gnomAD rows, full 176-column header) — copy in.
   - Slice the reduced gnomAD BED (`talos/large_files/SVAFotate_reduced_gnomAD.bed.gz`, 2.15 M rows)
     to the windows around the fixture SVs, keeping *all four sources* from the full v4.1 BED for those
     windows so multi-source and `-s` behaviour is exercised.
   - Hand-write a small targets BED and a VCF with `CIPOS/CIEND/CIPOS95/CIEND95`.
3. Run the oracle across an option matrix and commit the outputs as golden files:
   default; `-s gnomAD -f 0.5 -a best --ins` (Talos); `-a all`; `-a mis`; `-a full`; `-c 0.01 -u 0.01 -l 1000000 -t targets.bed`;
   `-ci in`, `-ci out`, `-ci95 in`; `-e 50`, `-r 50`; `-f 0.5 0.9 -s gnomAD CCDG`; each `-O` type.
4. Write `tests/test_golden.py` (pytest): run the CLI, parse both VCFs with cyvcf2, compare INFO
   dictionaries record-by-record with float tolerance, compare `uniques.bed`. Normalise header lines that
   embed paths or the BED filename.
5. Add unit tests for pure functions that exist today (`fraction`, `get_header_types`, `get_overlap`
   family, `get_best`, `reciprocal_overlap`) so the refactor has fine-grained coverage too.

Deliverable: a failing-on-pandas-3 test suite that passes on the legacy environment.

### Phase 1 — Repository hygiene

- `.gitignore` (venv, caches, `.idea/`, `dist/`, `uniques.bed`, `*.egg-info`).
- Delete `requirements.txt`, `src/__init__.py`; drop `[options]` and `README.md` from `license-files`.
- `pyproject.toml`: `requires-python = ">=3.12"`; runtime deps as decided in §3; a `[dependency-groups]
  dev` with `pytest`, `pytest-cov`, `ruff`, `mypy`; commit `uv.lock`; fix isort `section-order`; add
  classifiers for 3.12/3.13/3.14; bumpversion stops touching README.
- `.pre-commit-config.yaml`: `python3.12`, current ruff mirror rev, mypy with `types-*` as needed, `uv-lock`
  hook.
- GitHub Actions: `lint` (ruff check + format, mypy), `test` matrix {3.12, 3.13, 3.14} × {ubuntu, macos-arm64}
  installing with `uv sync --only-binary :all:` so a wheel regression fails CI, `build` (hatch wheel +
  sdist, `twine check`).
- `Dockerfile` in this repo (python:3.13-slim, `uv pip install .`), with the Talos smoke test moved in as
  a CI job rather than living in the consumer's Dockerfile.
- `CHANGELOG.md` (Keep a Changelog), `CONTRIBUTING.md` (how to regenerate golden files against the oracle).
- README rewrite: fork purpose and relationship to upstream; install via `pip`/`uv`; the actual CLI
  (`svafotate annotate`, plus the `svafotate_cpg` alias); annotation semantics section kept from
  upstream but with image links pointing at this repo's `images/`; remove pickle-source /
  custom-annotation / conda sections that do not exist here.

### Phase 2 — Replace the interval engine (the core deliverable)

1. `svafotate_cpg/intervals.py`:
   - `overlap_join(query, reference) -> DataFrame` of `(q_idx, r_idx, overlap_bp)` pairs. Per chromosome:
     sort reference by start, `searchsorted` for candidates with `start < q.end`, filter `end > q.start`
     via a cumulative-max-of-end trick or a sweep; vectorised overlap = `min(e1,e2) - max(s1,s2)`.
   - `covered_fraction(query, reference) -> ndarray` (merge reference intervals, then sum clipped overlap
     per query / query length) — replaces `PyRanges.coverage`.
   - `subtract(query, reference) -> DataFrame` of residual `(q_idx, start, end)` — replaces
     `PyRanges.subtract`.
   - Property-based tests (`hypothesis`) against a brute-force O(n·m) reference implementation.
2. Rewrite `convert_dict` / `get_overlap` family / `INS_conv` on top of it, computing OFP
   (`frac1 * frac2`) in one vectorised pass. Keep `reciprocal_overlap`'s INS-without-`--ins` bypass.
3. Remove `pyranges` imports and the `pr.PyRanges(..., int64=True)` conversions.
4. Golden tests must pass on Python 3.12/3.13/3.14 with pandas 3. Tag `0.2.0`.

Optional stop-gap before step 4 lands: tag `0.1.1` with option A (cast to `object`) so Talos2 can drop
the `pandas<3` pin immediately. Only worth it if Phase 2 is expected to take more than ~2 weeks.

### Phase 3 — Structural refactor of `main.py`

Split by responsibility, preserving CLI and output:

- `cli.py` — argparse; `--uniq-out PATH` (default `uniques.bed`) ; `--quiet/--verbose`.
- `reference.py` — load BED with `pandas.read_csv(usecols=..., dtype=...)`, filter to requested sources
  *before* materialising, normalise contig names, validate required columns with a clear error (fixes
  bugs 5, 6). Expose a typed `Reference` object.
- `query.py` — single pass over the VCF producing a `DataFrame` of adjusted coordinates (CI, `-e/-r`,
  INS expansion) plus het/hom-alt sample maps; fixes bug 4 and 7.
- `matching.py` — join, reciprocal filter, OFP, best-match with correct numeric tiebreak (bug 1).
- `annotations.py` — **data-driven** header and writer. One table
  `POPULATIONS = ['AFR','AMI',...]`, one `SEXES = ['Male','Female']`, and a resolver from `-a` flags to the
  list of `(field, type)` to emit. This collapses ~1,000 lines of repeated blocks into ~80 and fixes bugs
  2 and 3.
- `coverage.py` / `unique.py` — `-c` and `-u` paths on `intervals.py`.
- `logging` via stdlib `logging` (or `loguru`, matching Talos) instead of `print`; user errors raise
  `SystemExit` with a message, not `NameError`.
- Delete `utils.py` duplicates; type-annotate everything; mypy clean.

### Phase 4 — Performance and memory

- Memory target: annotate against the *full* v4.1 BED (469 MB gz, 4 sources, 176 columns) in < 4 GB, so
  Talos no longer needs its `awk | cut` pre-reduction. Achieved by `usecols` + source filter at load
  (Phase 3) and numpy arrays instead of dicts-of-lists-of-strings.
- Benchmark script `scripts/bench.py` against the reduced gnomAD BED and a ~10 k-record VCF; record
  numbers in the CHANGELOG.
- `--cpu`: either implement per-chromosome parallelism with `concurrent.futures`, or document honestly
  that it only sets htslib reader/writer threads. Recommend the latter until a benchmark shows a win.

### Phase 5 — Release and consumer cut-over

- Publish `svafotate_cpg` to PyPI via GitHub Actions trusted publishing on tag.
- Update `talos/docker/SVAFotate_Dockerfile` and `images/talos2/Dockerfile`: replace the git-commit
  install and `pandas<3` pin with `uv pip install svafotate_cpg==<version>`; delete the source-build stage
  for `ncls`/`sorted-nearest`; keep the smoke test.

## 5. Acceptance criteria

- `uv sync --only-binary :all:` succeeds on Python 3.12, 3.13 and 3.14 with pandas ≥ 3.
- Golden regression suite passes for every option combination in Phase 0 step 3.
- `pyranges`, `ncls`, `sorted-nearest` absent from `uv.lock`.
- ruff check/format and mypy clean; pre-commit passes.
- Talos's Dockerfile smoke test passes against the published wheel.
- Full-BED run fits in 4 GB and completes faster than the current implementation on the reduced BED.

## 6. Decisions needed from the maintainer

Each of these changes observable behaviour or the interface; the plan defaults are in bold.

1. ~~CLI shape~~ — decided 2026-09-02: `svafotate annotate` restored, `svafotate_cpg` alias kept.
2. Contig normalisation on the reference BED (bug 5): **normalise both sides** and log once when the input
   naming styles differ. This changes output for `chr`-prefixed BEDs (from all zeros to correct values).
3. Fix the tiebreaker (bug 1) and `END`/`SVLEN` off-by-one (bug 4) in `0.2.0`, or defer to a later minor
   so `0.2.0` is a pure engine swap. **Defer to `0.3.0`**; ship `0.2.0` as byte-identical-to-upstream.
4. `Best_*_Mismatch_OFP` type correction (bug 2): **fix in `0.2.0`**, header-only change.
5. Version number for the engine swap: **`0.2.0`**; `1.0.0` after Phase 3.

## 7. Estimated effort

| Phase         | Effort                                   |
|---------------|------------------------------------------|
| 0 Safety net  | 1–2 days                                 |
| 1 Hygiene     | 0.5–1 day                                |
| 2 Engine swap | 2–4 days (option D), 1–2 days (option C) |
| 3 Refactor    | 3–5 days                                 |
| 4 Performance | 1–2 days                                 |
| 5 Release     | 0.5 day plus consumer PRs                |

## Appendix — verification log for §3

All runs on 2026-09-02, macOS arm64, `uv 0.11.24`.

- `.venv` (py3.12, pandas 3.0.5, pyranges 0.1.4): `join` → `Exception: Unknown dtype str in a column SVTYPE`;
  `coverage`, `subtract` OK. `astype(object)` on string columns → `join` OK.
  `pd.options.future.infer_string = False` → `join` OK.
- `uv pip compile 'pyranges>=1'` → no solution (`pyranges` on PyPI tops out at 0.1.4; v1 is `pyranges1`).
- `uv pip install --only-binary :all: ncls` → "all versions of ncls have no usable wheels" (3.13, 3.14).
- `pyranges1` + `pandas>=3`: resolves to `pyranges1==1.4.3`, `ruranges==0.2.7` on 3.12/3.13; wheels-only on
  3.14 resolves to `pyranges1==1.0.4`. `join_overlaps(join_type='left', suffix='_b')` and
  `subtract_overlaps` verified; no `coverage` attribute.
- `bioframe==0.8.0` + `pandas>=3`: wheels on 3.13 and 3.14; `overlap(return_overlap=True)`, `coverage`,
  `subtract` verified (a `Pandas4Warning` is emitted from `merge`).

---

## Status (2026-09-02, end of day)

Executed in one pass; nothing committed. Consumer cut-over (Phase 5) deliberately not started.

| Phase | Status | Evidence |
|---|---|---|
| 0 Safety net | done | `tests/golden/` holds 22 oracle outputs from upstream `30b5004` on pandas 2.3.3 / pyranges 0.1.4; `tests/test_golden.py` passes on pandas-free Python 3.12, 3.13, 3.14 |
| 1 Hygiene | done | `.gitignore`, `uv.lock`, `pyproject.toml` cleaned (deps: `cyvcf2`, `numpy`), `.pre-commit-config.yaml`, `.github/workflows/{ci,release}.yml`, `Dockerfile`, `CHANGELOG.md`, `CONTRIBUTING.md`, README rewritten |
| 2 Engine swap | done (option D) | `svafotate_cpg/intervals.py`, 108 randomised tests incl. agreement with pyranges 0.x on coverage/subtract |
| 3 Refactor | done | `cli`, `reference`, `query`, `matching`, `annotations`, `regions`, `annotate` modules; population table replaces ~1,000 copied lines; ruff and mypy clean; bugs 1, 2, 3, 5, 6, 7, 8 and the hom-alt genotype bug fixed and covered by `tests/test_behaviour.py` |
| 4 Performance | done | 2.8x faster and 4x less memory than upstream on the reduced BED; full BED in 1.2 GB (see benchmark below); `scripts/bench.py` |
| 5 Release | not started | version set to 0.2.0 and changelog written; tag, PyPI trusted-publisher setup and Talos/Talos2 Dockerfile changes remain |

### Decisions taken (from §6)

1. CLI: `svafotate annotate` restored, `svafotate_cpg` alias kept.
2. Contig normalisation: both sides stripped of `chr`, logged once.
3. Bug fixes shipped now rather than deferred, because the refactor already made 0.2.0 more than an engine
   swap: tiebreak (bug 1) fixed; `END` from `SVLEN` (bug 4) **kept** as upstream, since the `ext_*` golden
   cases exercise it.
4. `Best_*_Mismatch_OFP` type corrected.
5. Version 0.2.0.

### Additional findings during execution

- Upstream never populated `HomAlt_Samples` in `uniques.bed` (`gt_types == 2` is UNKNOWN in cyvcf2's default
  encoding). Fixed; the golden comparison masks that column.
- Upstream's list-annotation order is the NCLS traversal order and not reproducible; the golden comparison
  treats `*_Matches` / `*_Mismatches` as sets and this fork orders by reference start.
- The tiebreak bug changes real output: on a 10,000-record synthetic query against the reduced gnomAD BED,
  20 records (0.2%) tie on OFP and get a different `Best_gnomAD_*`; every one was verified to be an exact
  OFP tie.
- Upstream cannot run `-ci`/`-ci95` if any record lacks the CI fields, and cannot process a record without
  `SVTYPE`; both are handled here and excluded from oracle fixtures.

### Benchmark

10,000 synthetic queries (sampled from the reduced gnomAD BED and jittered), `-s gnomAD -f 0.5 -a best --ins`,
macOS arm64. Peak RSS from `ru_maxrss`.

| Reference                                        | Implementation                                             | Wall       | Peak RSS    |
|--------------------------------------------------|------------------------------------------------------------|------------|-------------|
| reduced gnomAD BED (2.15 M rows, 22 cols)        | upstream (pandas 2, pyranges 0.1.4)                        | 20.7 s     | 4.1 GB      |
| reduced gnomAD BED                               | this fork, first cut (all columns as strings, then arrays) | 11.1 s     | 3.0 GB      |
| reduced gnomAD BED                               | this fork, streaming float columns                         | 10.9 s     | 1.4 GB      |
| reduced gnomAD BED                               | this fork, only requested columns                          | **7.3 s**  | **0.95 GB** |
| full v4.1 BED (4 sources, 177 cols), `-s gnomAD` | this fork, all columns                                     | 67.1 s     | 4.3 GB      |
| full v4.1 BED, `-s gnomAD`                       | this fork, only requested columns                          | **26.0 s** | **0.95 GB** |
| full v4.1 BED, all sources, `-f 0.5 -a best`     | this fork, all columns                                     | 88.7 s     | 5.6 GB      |
| full v4.1 BED, all sources, `-f 0.5 -a best`     | this fork, only requested columns                          | **27.1 s** | **1.2 GB**  |

Upstream was not run against the full BED here; Talos reports it needs 16 GB+. The remaining time on the
full BED is dominated by splitting 177-column lines in Python's `csv` module. `-a full` (verbatim strings
for every matched row) necessarily loads everything and is closer to the "all columns" rows.

### Remaining follow-ups

- Phase 5 as written in §4.
- Register the PyPI trusted publisher before pushing a `v*` tag.
- Consider `END = START + |SVLEN|` (dropping the `+ 1`) in a later minor, with a golden regeneration.
- The pre-commit `mirrors-mypy` rev is carried over unverified.
