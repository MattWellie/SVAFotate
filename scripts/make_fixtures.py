"""Generate the hand-written regression fixtures under tests/data.

The 8-record ``joint_sv.vcf.gz`` and the reference BED slices are copied from real data (see
CONTRIBUTING.md); everything written here is synthetic and exists to exercise code paths those files
do not reach: confidence-interval boundaries, a missing END, START == END, a second chromosome, and a
targets BED with ``chr``-prefixed contigs.

Run from the repository root: ``python scripts/make_fixtures.py``
"""

from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / 'tests' / 'data'

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Breakend">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description="95% confidence interval around POS">
##INFO=<ID=CIEND95,Number=2,Type=Integer,Description="95% confidence interval around END">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband\tfather\tmother
"""

CI = 'CIPOS=-500,500;CIEND=-500,500;CIPOS95=-100,100;CIEND95=-100,100'

# (chrom, pos, id, alt, info, genotypes). Every record carries CIs because upstream subscripts
# CIPOS/CIEND unconditionally once -ci/-ci95 is given.
RECORDS = [
    # exact gnomAD DEL match, with CIs so -ci in/out shrink/grow it
    ('chr1', 509968, 'ci_common_del', '<DEL>', f'SVTYPE=DEL;SVLEN=180032;END=690000;{CI}', ('0/1', '0/1', '0/0')),
    # DUP with CIs, overlapping four gnomAD DUPs of very different AF
    ('chr1', 21489967, 'ci_dup', '<DUP>', f'SVTYPE=DUP;SVLEN=4889;END=21494856;{CI}', ('0/1', '0/0', '1/1')),
    # no END at all, only SVLEN: exercises the END-from-SVLEN fallback
    ('chr1', 17372196, 'no_end_del', '<DEL>', f'SVTYPE=DEL;SVLEN=29503;{CI}', ('1/1', '0/0', '0/0')),
    # START == END after coordinate conversion (INS with END == POS and no SVLEN)
    ('chr1', 66340, 'ins_no_svlen', '<INS>', f'SVTYPE=INS;END=66340;{CI}', ('0/1', '0/0', '0/1')),
    # NOTE: a record with no SVTYPE is deliberately absent - upstream crashes on it (TypeError sorting
    # None against str), so there is no oracle behaviour to match. See tests/test_query.py instead.
    # a second chromosome with nothing in the reference slice: must annotate as zero, not crash
    ('chr2', 100000, 'chr2_del', '<DEL>', f'SVTYPE=DEL;SVLEN=5000;END=105000;{CI}', ('0/1', '0/1', '0/1')),
    # deletion nested inside the large gnomAD DEL, all three samples hom-alt -> only HomAlt_Samples populated
    ('chr1', 600001, 'nested_del_homalt', '<DEL>', f'SVTYPE=DEL;SVLEN=200;END=600200;{CI}', ('1/1', '1/1', '1/1')),
]


def write_extended_vcf() -> None:
    lines = [HEADER]
    for chrom, pos, vid, alt, info, gts in RECORDS:
        lines.append('\t'.join([chrom, str(pos), vid, 'N', alt, '.', 'PASS', info, 'GT', *gts]) + '\n')
    (DATA / 'extended.vcf').write_text(''.join(lines))


def write_targets_bed() -> None:
    # chr-prefixed on purpose: the tool strips the prefix from targets, so these must still hit
    rows = [
        ('chr1', 17372196, 17401699, 'PADI6'),
        ('chr1', 38874069, 38881587, 'GJA9'),
        ('chr1', 685679, 686673, 'OR4F16'),
        ('chr1', 600050, 600100, 'tiny_in_nested_del'),
        ('chr1', 1180000, 1190000, 'inside_inv'),
        ('chr2', 100500, 100600, 'chr2_target'),
    ]
    text = '#CHROM\tSTART\tEND\tID\n' + ''.join('\t'.join(map(str, r)) + '\n' for r in rows)
    (DATA / 'targets.bed').write_text(text)


if __name__ == '__main__':
    write_extended_vcf()
    write_targets_bed()
    print(f'wrote {DATA / "extended.vcf"} and {DATA / "targets.bed"}')
