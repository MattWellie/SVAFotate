"""Behaviour this fork defines beyond the upstream oracle (bug fixes and new options)."""

import gzip
from pathlib import Path

import pytest
from cyvcf2 import VCF

from svafotate_cpg.cli import main

VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband\tfather\tmother
"""

BED_HEADER = '#CHROM\tSTART\tEND\tSVLEN\tSVTYPE\tSOURCE\tSV_ID\tAF\tHomRef\tHet\tHomAlt\tPopMax_AF\n'


def write_vcf(path: Path, records: list[str]) -> Path:
    path.write_text(VCF_HEADER + ''.join(r + '\n' for r in records))
    return path


def write_bed(path: Path, rows: list[str], gz: bool = False) -> Path:
    text = BED_HEADER + ''.join(r + '\n' for r in rows)
    if gz:
        with gzip.open(path, 'wt') as f:
            f.write(text)
    else:
        path.write_text(text)
    return path


def run(vcf: Path, out: Path, bed: Path, *extra: str) -> dict[str, dict]:
    main(['-q', 'annotate', '-v', str(vcf), '-o', str(out), '-b', str(bed), *extra])
    return {v.ID: dict(v.INFO) for v in VCF(str(out))}


def test_record_without_svtype_is_annotated_with_zeros(tmp_path: Path):
    """Upstream crashed sorting a None SVTYPE; here the record simply matches nothing."""
    vcf = write_vcf(
        tmp_path / 'q.vcf', ['chr1\t1000\tno_type\tN\t<DEL>\t.\tPASS\tSVLEN=500;END=1500\tGT\t0/1\t0/0\t0/0']
    )
    bed = write_bed(tmp_path / 'r.bed', ['1\t999\t1500\t501\tDEL\tgnomAD\tref1\t0.5\t10\t10\t10\t0.6'])
    info = run(
        vcf,
        tmp_path / 'out.vcf',
        bed,
        '-a',
        'best',
        'mis',
        '-c',
        '0.0',
        '-u',
        '0.0',
        '--uniq-out',
        str(tmp_path / 'u.bed'),
    )
    assert info['no_type']['Max_AF'] == 0
    assert info['no_type']['gnomAD_Count'] == 0
    # it does count as a mismatch against the DEL it overlaps
    assert info['no_type']['gnomAD_Mismatches'] == 'ref1'


def test_homalt_samples_are_reported_in_uniques_bed(tmp_path: Path):
    """Upstream tested gt_types == 2 (UNKNOWN in cyvcf2's default encoding), so HomAlt_Samples was always None."""
    vcf = write_vcf(
        tmp_path / 'q.vcf', ['chr1\t1000\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=1000;END=2000\tGT\t1/1\t0/1\t0/0']
    )
    bed = write_bed(tmp_path / 'r.bed', ['1\t1499\t2000\t501\tDEL\tgnomAD\tref1\t0.5\t10\t10\t10\t0.6'])
    out_bed = tmp_path / 'my_uniques.bed'
    run(vcf, tmp_path / 'out.vcf', bed, '-u', '0.0', '--uniq-out', str(out_bed))
    rows = out_bed.read_text().splitlines()
    assert rows[0].split('\t') == ['#CHROM', 'START', 'END', 'SVTYPE', 'SV_ID', 'Het_Samples', 'HomAlt_Samples']
    assert rows[1].split('\t') == ['1', '999', '1499', 'DEL', 'del', 'father', 'proband']


def test_missing_reference_column_is_a_clear_error(tmp_path: Path):
    vcf = write_vcf(
        tmp_path / 'q.vcf', ['chr1\t1000\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=1000;END=2000\tGT\t0/1\t0/0\t0/0']
    )
    bed = tmp_path / 'r.bed'
    bed.write_text('#CHROM\tSTART\tEND\tSVLEN\tSVTYPE\tSOURCE\tSV_ID\tAF\n1\t999\t2000\t1001\tDEL\tgnomAD\tref1\t0.5\n')
    with pytest.raises(SystemExit, match='Het, HomAlt, PopMax_AF'):
        run(vcf, tmp_path / 'out.vcf', bed)


def test_chr_prefixed_reference_bed_matches(tmp_path: Path):
    """Upstream stripped 'chr' from the VCF only, so a chr-prefixed BED silently produced all-zero AFs."""
    vcf = write_vcf(
        tmp_path / 'q.vcf', ['chr1\t1000\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=1000;END=2000\tGT\t0/1\t0/0\t0/0']
    )
    bed = write_bed(tmp_path / 'r.bed.gz', ['chr1\t999\t2000\t1001\tDEL\tgnomAD\tref1\t0.25\t10\t10\t10\t0.6'], gz=True)
    info = run(vcf, tmp_path / 'out.vcf', bed, '-f', '0.9', '-a', 'best')
    assert info['del']['Max_AF'] == pytest.approx(0.25)
    assert info['del']['Best_gnomAD_ID'] == 'ref1'
    assert info['del']['gnomAD_Count'] == 1


def test_best_match_tiebreak_is_numeric_af_then_homalt(tmp_path: Path):
    """Two identical-coordinate references tie on OFP; the higher AF wins (compared as numbers, not strings)."""
    vcf = write_vcf(
        tmp_path / 'q.vcf', ['chr1\t1000\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=1000;END=2000\tGT\t0/1\t0/0\t0/0']
    )
    bed = write_bed(
        tmp_path / 'r.bed',
        [
            '1\t999\t2000\t1001\tDEL\tgnomAD\tlow_af_first\t0.10\t10\t10\t1\t0.1',
            '1\t999\t2000\t1001\tDEL\tgnomAD\thigh_af\t0.9\t10\t10\t1\t0.9',
            '1\t999\t2000\t1001\tDEL\tgnomAD\tsame_af_more_homalt\t0.9\t10\t10\t5\t0.9',
        ],
    )
    info = run(vcf, tmp_path / 'out.vcf', bed, '-a', 'best')
    assert info['del']['Best_gnomAD_ID'] == 'same_af_more_homalt'
    assert info['del']['Best_gnomAD_OFP'] == pytest.approx(1.0)
    assert info['del']['gnomAD_Count'] == 3


@pytest.mark.parametrize('out_type', ['vcfgz', 'bcf', 'bcfgz'])
def test_compressed_output_types_are_readable(tmp_path: Path, out_type: str):
    vcf = write_vcf(
        tmp_path / 'q.vcf', ['chr1\t1000\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=1000;END=2000\tGT\t0/1\t0/0\t0/0']
    )
    bed = write_bed(tmp_path / 'r.bed', ['1\t999\t2000\t1001\tDEL\tgnomAD\tref1\t0.25\t10\t10\t10\t0.6'])
    suffix = {'vcfgz': 'vcf.gz', 'bcf': 'bcf', 'bcfgz': 'bcf'}[out_type]
    info = run(vcf, tmp_path / f'out.{suffix}', bed, '-O', out_type)
    assert info['del']['Max_AF'] == pytest.approx(0.25)


def test_unknown_source_is_rejected(tmp_path: Path):
    vcf = write_vcf(
        tmp_path / 'q.vcf', ['chr1\t1000\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=1000;END=2000\tGT\t0/1\t0/0\t0/0']
    )
    bed = write_bed(tmp_path / 'r.bed', ['1\t999\t2000\t1001\tDEL\tgnomAD\tref1\t0.25\t10\t10\t10\t0.6'])
    with pytest.raises(SystemExit, match='CCDG is an unexpected source'):
        run(vcf, tmp_path / 'out.vcf', bed, '-s', 'CCDG')


def test_lim_requires_cov_or_uniq(capsys):
    with pytest.raises(SystemExit):
        main(['annotate', '-v', 'x.vcf', '-o', 'y.vcf', '-b', 'z.bed', '-l', '100'])
    assert '--lim must be used alongside --cov or --uniq' in capsys.readouterr().err
