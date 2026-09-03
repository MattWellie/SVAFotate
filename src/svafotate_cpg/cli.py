"""Command line entry point: ``svafotate annotate ...``."""

from __future__ import annotations

import logging
import sys
from argparse import ArgumentParser, ArgumentTypeError

from . import __version__
from .annotate import DEFAULT_MINF, run
from .annotations import EXTRA_ANNOTATIONS


def fraction(value: str) -> float:
    """argparse type for values which must sit between 0.0 and 1.0."""
    try:
        number = float(value)
    except ValueError:
        raise ArgumentTypeError(f'{value} is not a number') from None
    if number < 0 or number > 1:
        raise ArgumentTypeError(f'{value} must be between 0.0 and 1.0')
    return number


def get_parser() -> tuple[ArgumentParser, ArgumentParser]:
    """The top-level parser and the ``annotate`` sub-parser (returned so validation can use its ``error``)."""
    parser = ArgumentParser(prog='svafotate', description='SVAFotate: Structural Variant Allele Frequency annotator')
    parser.add_argument(
        '-v', '--version', help='Installed version', action='version', version=f'%(prog)s {__version__}'
    )
    parser.add_argument('-q', '--quiet', action='store_true', help='Only report warnings and errors')

    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', required=True)
    annotate_parser = subparsers.add_parser(
        'annotate',
        help='Annotate SV VCF File',
        description='Annotate your SV VCF files with Population Allele Frequency Info',
    )

    req = annotate_parser.add_argument_group('Required Arguments')
    opt = annotate_parser.add_argument_group('Optional Arguments')

    req.add_argument(
        '-v', '--vcf', metavar='INPUT VCF', required=True, help='Path and/or name of the VCF file to annotate.'
    )
    req.add_argument(
        '-o', '--out', metavar='OUTPUT VCF', required=True, help='Path and/or name of the output VCF file.'
    )
    req.add_argument(
        '-b', '--bed', metavar='SOURCE BED', required=True, help='Path and/or name of the combined sources AF bed file.'
    )

    opt.add_argument(
        '-f',
        '--minf',
        metavar='MINIMUM OVERLAP FRACTION',
        nargs='*',
        type=fraction,
        help='A space separated list of minimum reciprocal overlap fractions required between SVs for each source'
        ' listed with the `-s` option. If `-s` is not used, only the first minf will be used and will be applied'
        f' to all sources. minf values must be between 0.0 and 1.0 (Default = {DEFAULT_MINF}).',
    )
    opt.add_argument(
        '-s',
        '--sources',
        metavar='SOURCES TO ANNOTATE',
        nargs='*',
        help="Space separated list of data sources to use for annotation. If '-s' is not used, all sources available in the"
        " source bed file will be used (Example: ' -s CCDG gnomAD ' ).",
    )
    opt.add_argument(
        '-a',
        '--ann',
        metavar='EXTRA ANNOTATIONS',
        nargs='*',
        default=[],
        choices=list(EXTRA_ANNOTATIONS),
        help='By default, only the Max_AF, Max_Het and Max_HomAlt counts, and Max_PopMax_AF are annotated in the output VCF file.'
        ' `-a` can be used to add additional annotations, with each annotation separated by a space'
        f" (Example ' -a mf best pops ' ). Choices = [{', '.join(EXTRA_ANNOTATIONS)}]",
    )
    opt.add_argument(
        '-c',
        '--cov',
        metavar='OBSERVED SV COVERAGE',
        type=fraction,
        help='Add an annotation reflecting how much of the queried SV genomic space has been previously observed with the same SVTYPE.'
        ' Uses the data sources listed with -s as the previously observed SVs.'
        ' Please provide minimum AF to exclude all SVs from data sources with a total AF below that value (must be between 0 and 1.0).',
    )
    opt.add_argument(
        '-u',
        '--uniq',
        metavar='UNIQUE SV REGIONS',
        type=fraction,
        help='Generate a file of unique SV regions (see --uniq-out).'
        ' These regions reflect genomic space within the queried SV region that have not been previously observed with the same SVTYPE.'
        ' This will also add an annotation regarding the number of unique regions within a given SV.'
        ' Please provide minimum AF to exclude all SVs from data sources with a total AF below that value (must be between 0 and 1.0).',
    )
    opt.add_argument(
        '--uniq-out',
        metavar='UNIQUE BED',
        default='uniques.bed',
        help='Where to write the unique regions BED (Default = uniques.bed in the working directory).',
    )
    opt.add_argument(
        '-l',
        '--lim',
        metavar='SV SIZE LIMIT',
        type=int,
        help='Only include previously observed SVs from data sources with a size less than or equal to this value (only available when using --cov or --uniq).',
    )
    opt.add_argument(
        '-t',
        '--target',
        metavar='TARGETS BED FILE',
        help='Path to target regions BED file. Expected format is a tab delimited file listing CHROM START END ID'
        ' where ID is a genomic region identifier that will be listed as an annotation if an overlap exists between a given SV and the target regions.',
    )

    ci_group = opt.add_mutually_exclusive_group()
    ci_group.add_argument(
        '-ci',
        '--ci',
        metavar='USE CI BOUNDARIES',
        choices=['in', 'out'],
        help='Expects CIPOS and CIEND to be included in the INFO field of the input VCF (--vcf).'
        " If argument is selected, use 'inner' or 'outer' confidence intervals (CIPOS, CIEND) for SV boundaries. Choices = [in, out]",
    )
    ci_group.add_argument(
        '-ci95',
        '--ci95',
        metavar='USE CI BOUNDARIES',
        choices=['in', 'out'],
        help='Expects CIPOS95 and CIEND95 to be included in the INFO field of the input VCF (--vcf).'
        " If argument is selected, use 'inner' or 'outer' confidence intervals (CIPOS95, CIEND95) for SV boundaries. Choices = [in, out]",
    )

    size_group = opt.add_mutually_exclusive_group()
    size_group.add_argument(
        '-e',
        '--emb',
        metavar='EMBIGGEN THE SV SIZE',
        type=int,
        help='Increase the size of the SV coordinates in the input VCF (--vcf) by a single integer; Subtract that value from the start and add it to the end of each set of coordinates.',
    )
    size_group.add_argument(
        '-r',
        '--red',
        metavar='REDUCE THE SV SIZE',
        type=int,
        help='Reduce the size of the SV coordinates in the input VCF (--vcf) by a single integer; Add that value to the start and subtract it from the end of each set of coordinates.',
    )

    opt.add_argument(
        '--cpu',
        metavar='CPU Count',
        type=int,
        default=1,
        help='The number of threads htslib may use to read and write the VCF (Default = 1).',
    )
    opt.add_argument(
        '-O',
        '--out_type',
        metavar='OUTPUT FILE TYPE',
        choices=['vcf', 'vcfgz', 'bcf', 'bcfgz'],
        default='vcf',
        help='Specify output type as VCF (vcf), compressed VCF (vcfgz), BCF (bcf), or compressed BCF (bcfgz).',
    )
    opt.add_argument(
        '--ins',
        action='store_true',
        help='Use SVLEN to adjust insertions coordinates, allowing for reciprocal overlap matching with insertions',
    )
    return parser, annotate_parser


def main(argv: list[str] | None = None) -> int:
    parser, annotate_parser = get_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO, format='%(message)s', stream=sys.stderr, force=True
    )

    if args.command == 'annotate':
        # --lim only means anything alongside the coverage/unique calculations
        if args.lim is not None and args.cov is None and args.uniq is None:
            annotate_parser.error(
                'The argument --lim must be used alongside --cov or --uniq; please include either --cov or --uniq'
            )
        run(args)
    return 0
