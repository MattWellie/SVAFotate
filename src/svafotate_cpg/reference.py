"""Loading the population-frequency reference BED.

The BED has a mandatory ``#``-prefixed header. Columns 1-7 are positional (CHROM, START, END, SVLEN,
SVTYPE, SOURCE, SV_ID); every column from the 8th onwards is a value column addressable by name
(``AF``, ``Het``, ``HomAlt``, ``PopMax_AF``, ``AFR_AF`` ...). ``NA`` marks a missing value.

Only rows from the requested sources are kept, and value columns are stored as float64 arrays (``NaN``
for ``NA``) rather than Python strings, which is what keeps the full four-source BED within a few GB.
Columns that turn out to be non-numeric (e.g. ``In_Pop``) stay as object arrays. The verbatim strings
are retained only when the ``full`` annotation asks for them.
"""

import csv
import gzip
import io
import logging
import sys
from array import array
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from .intervals import Intervals

logger = logging.getLogger(__name__)

POSITIONAL = ('CHROM', 'START', 'END', 'SVLEN', 'SVTYPE', 'SOURCE', 'SV_ID')
MISSING = frozenset({'NA', 'None', ''})


def open_text(path: str) -> io.TextIOBase:
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8')  # type: ignore[return-value]
    return open(path, encoding='utf-8')


def strip_chr(contig: str) -> str:
    return contig[3:] if contig.startswith('chr') else contig


def column_type(name: str) -> str:
    """VCF INFO type for a value column: anything named ``AF`` or ``*_AF`` is a Float, the rest Integers."""
    return 'Float' if name == 'AF' or '_AF' in name else 'Integer'


def to_int(value: str) -> int:
    try:
        return int(value)
    except ValueError:
        return int(float(value))


@dataclass
class SourceData:
    """Everything known about one data source (e.g. gnomAD)."""

    name: str
    intervals: Intervals
    svtype: NDArray[np.object_]
    sv_id: NDArray[np.object_]
    svlen: NDArray[np.int64]
    af: NDArray[np.float64]  # the AF column, used for the --cov/--uniq frequency cut-offs
    columns: dict[str, NDArray] = field(default_factory=dict)  # value column name -> per-row array
    raw_rows: list[list[str]] | None = None  # verbatim value strings, only when --ann full is requested

    def __len__(self) -> int:
        return len(self.sv_id)

    def values(self, column: str, rows: NDArray[np.intp]) -> NDArray:
        """Values of one column for the given rows, with missing (NaN / 'NA') entries removed."""
        col = self.columns[column][rows]
        if col.dtype.kind == 'f':
            return col[~np.isnan(col)]
        return col[[v not in MISSING for v in col]]


@dataclass
class Reference:
    path: str
    value_columns: list[str]  # names of columns 8..N in file order
    sources: dict[str, SourceData]

    @property
    def source_names(self) -> list[str]:
        return sorted(self.sources)


class ColumnBuffer:
    """Accumulates one value column while streaming the BED.

    Values are parsed to float64 as they arrive and kept in a compact ``array('d')`` (8 bytes each). If a
    non-missing value turns out not to be a number, the column switches to string mode for good; the
    values seen so far are re-materialised as strings (``NA`` for missing).
    """

    __slots__ = ('numbers', 'strings')

    def __init__(self) -> None:
        self.numbers: array | None = array('d')
        self.strings: list[str] | None = None

    def append(self, value: str) -> None:
        if self.strings is not None:
            self.strings.append(value)
            return
        assert self.numbers is not None
        if value in MISSING:
            self.numbers.append(float('nan'))
            return
        try:
            self.numbers.append(float(value))
        except ValueError:
            self.strings = ['NA' if np.isnan(x) else repr(x) for x in self.numbers]
            self.strings.append(value)
            self.numbers = None

    def finish(self) -> NDArray:
        if self.strings is not None:
            return np.array(self.strings, dtype=object)
        assert self.numbers is not None
        return np.frombuffer(self.numbers, dtype=np.float64).copy()


def _clean_af(value: str) -> str:
    # gnomAD multi-allelic CNVs list AF as lowAF:highAF; use the high end, as upstream did for --cov/--uniq
    return value.split(':')[1] if ':' in value else value


def load_reference(
    path: str,
    wanted_sources: list[str] | None,
    keep_raw: bool,
    wanted_columns: list[str] | None = None,
) -> Reference:
    """Read the BED, keeping only rows of ``wanted_sources`` (all when None) and only ``wanted_columns``
    (all when None) as numeric arrays. ``AF`` is always kept. Missing wanted columns are reported at once."""
    logger.info('Reading the following data source file: %s', path)
    header: list[str] = []
    keep_names: list[str] = []
    keep_idx: list[int] = []
    buckets: dict[str, dict[str, list]] = {}
    seen_sources: set[str] = set()
    chr_prefixed = False

    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for fields in reader:
            if not fields:
                continue
            if fields[0].startswith('#'):
                header = [fields[0].lstrip('#'), *fields[1:]]
                if tuple(header[:7]) != POSITIONAL:
                    raise SystemExit(
                        f'The reference BED must start with the columns {" ".join(POSITIONAL)};'
                        f' found {" ".join(header[:7])}'
                    )
                if 'AF' not in header[7:]:
                    raise SystemExit('The reference BED has no AF column')
                keep_names = header[7:] if wanted_columns is None else ['AF', *[c for c in wanted_columns if c != 'AF']]
                missing = [c for c in keep_names if c not in header[7:]]
                if missing:
                    raise SystemExit(
                        f'The reference BED {path} lacks the column(s) required for the requested annotations: '
                        + ', '.join(missing)
                        + f'. Available value columns: {", ".join(header[7:])}'
                    )
                keep_idx = [header.index(c) for c in keep_names]
                continue
            if not header:
                raise SystemExit('A header is required in the bed source file')

            source = fields[5]
            seen_sources.add(source)
            if wanted_sources is not None and source not in wanted_sources:
                continue
            if len(fields) != len(header):
                raise SystemExit(f'Row for {fields[6]} has {len(fields)} columns; the header has {len(header)}')

            contig = fields[0]
            if contig.startswith('chr'):
                chr_prefixed = True
            b = buckets.get(source)
            if b is None:
                b = buckets[source] = {k: [] for k in ('chrom', 'start', 'end', 'svlen', 'svtype', 'sv_id', 'raw')}
                b['columns'] = [ColumnBuffer() for _ in keep_idx]
            b['chrom'].append(sys.intern(strip_chr(contig)))
            b['start'].append(int(fields[1]))
            b['end'].append(int(fields[2]))
            b['svlen'].append(to_int(fields[3]) if fields[3] not in MISSING else 0)
            b['svtype'].append(sys.intern(fields[4]))
            b['sv_id'].append(fields[6])
            fields[7] = _clean_af(fields[7])
            for buffer, j in zip(b['columns'], keep_idx, strict=True):
                buffer.append(fields[j])
            if keep_raw:
                b['raw'].append(fields[7:])

    if not header:
        raise SystemExit('A header is required in the bed source file')
    if chr_prefixed:
        logger.info('Reference BED contigs carry a "chr" prefix; it has been stripped to match the query VCF')

    value_columns = header[7:]
    if wanted_sources is not None:
        unknown = [s for s in wanted_sources if s not in seen_sources]
        if unknown:
            raise SystemExit(
                f'{", ".join(unknown)} is an unexpected source; acceptable sources include (case-sensitive): '
                + ','.join(sorted(seen_sources))
            )

    sources: dict[str, SourceData] = {}
    for name, b in buckets.items():
        columns: dict[str, NDArray] = {
            col_name: buffer.finish() for col_name, buffer in zip(keep_names, b['columns'], strict=True)
        }
        columns['SVTYPE'] = np.array(b['svtype'], dtype=object)
        af = columns['AF']
        if af.dtype.kind != 'f':
            raise SystemExit(f'The AF column of source {name} contains non-numeric values')
        sources[name] = SourceData(
            name=name,
            intervals=Intervals.from_columns(b['chrom'], b['start'], b['end']),
            svtype=columns['SVTYPE'],
            sv_id=np.array(b['sv_id'], dtype=object),
            svlen=np.array(b['svlen'], dtype=np.int64),
            af=af,
            columns=columns,
            raw_rows=b['raw'] if keep_raw else None,
        )
        logger.info('  %s: %d records', name, len(b['sv_id']))

    return Reference(path=path, value_columns=value_columns, sources=sources)
