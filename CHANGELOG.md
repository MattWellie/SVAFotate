# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.2.0] - unreleased

First release of the CPG fork of [fakedrtom/SVAFotate](https://github.com/fakedrtom/SVAFotate)
(forked at commit `30b5004`). The command line (`svafotate annotate ...`) and the annotations are
unchanged unless listed below; a golden regression suite checks 22 option combinations against upstream.

### Changed
- Package renamed `svafotate_cpg`; `svafotate_cpg` is also provided as a command alias for `svafotate`.
- Replaced pyranges 0.1.4 (and its source-only `ncls` / `sorted-nearest` backends) with an in-house
  numpy interval engine (`svafotate_cpg.intervals`). Runtime dependencies are now `cyvcf2` and `numpy`
  only, everything installs from wheels, and the tool runs on pandas-free Python 3.12-3.14.
- `main.py` (2,200 lines, one function) split into `cli`, `reference`, `query`, `matching`,
  `annotations`, `regions` and `annotate` modules. The eleven hand-copied population blocks are replaced
  by one table (`annotations.POPULATIONS`).
- The reference BED is loaded selectively: only rows of the requested sources, only the value columns the
  requested annotations read, as float64 arrays; verbatim strings only when `-a full` asks for them. On
  10,000 queries against the reduced gnomAD BED this takes 7 s and 0.95 GB (upstream: 21 s, 4.1 GB); the
  full four-source v4.1 BED with all sources takes 27 s and 1.2 GB and no longer needs pre-reduction.
- `chr` prefixes are stripped from reference BED contigs as well as query VCF contigs. Upstream only
  normalised the VCF, so a `chr`-prefixed BED silently produced `Max_AF=0` for every record.
- Diagnostics go to stderr through `logging` instead of `print`; `-q/--quiet` keeps only warnings.
- User errors (unknown source, malformed BED, missing columns) exit with a message instead of a
  `NameError` / `KeyError` traceback.
- The members of the `<source>_Matches` and `<source>_Mismatches` list annotations are ordered by
  reference start position. Upstream's order was whatever the NCLS index happened to yield.
- `Best_<source>_Mismatch_OFP` is declared `Type=Float` in the header (upstream: `String`).
- `--cpu` help text now says what it does (htslib read/write threads); nothing else was ever parallel.

### Fixed
- `HomAlt_Samples` in the `-u` BED was always `None`: upstream tested `gt_types == 2`, which is UNKNOWN in
  cyvcf2's default encoding. Hom-alt samples are now listed.
- Best-match tiebreak (equal OFP) compared the SVTYPE and Het columns as strings while meaning AF and
  HomAlt; it now compares AF, then HomAlt, numerically.
- A record without `SVTYPE` crashed the run (`TypeError` sorting `None`); it is now annotated with zeros
  and counted as a mismatch where it overlaps.
- A reference BED lacking a column needed by the requested annotations raised `KeyError` while writing
  the output; the columns are now validated before any work starts.
- `-ci`/`-ci95` on a record lacking `CIPOS`/`CIEND` crashed; it now warns and treats the interval as 0.
- The `-u` progress message named the wrong source (leaked loop variable).

### Added
- `--uniq-out PATH` to choose where `-u` writes its BED (default unchanged: `uniques.bed` in the CWD).
- Golden regression suite (`tests/test_golden.py`, `scripts/generate_golden.py`), randomised tests for
  the interval engine, behaviour tests for every fix above.
- CI on Python 3.12-3.14 (Linux and macOS) with wheel-only installs, Dockerfile smoke test, release
  workflow with PyPI trusted publishing.
- `Dockerfile`, `CONTRIBUTING.md`, `scripts/bench.py`.

### Removed
- `requirements.txt` (conda-era pins superseded by `pyproject.toml` / `uv.lock`).
- The `pickle-source` and `custom-annotation` subcommands (never functional upstream).

### Kept deliberately
- When a record has `SVLEN` but no `END`, the end is `START + |SVLEN| + 1`, one base longer than the
  VCF convention implies. Changing it would alter output for such records; it is documented instead.
