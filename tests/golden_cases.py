"""The regression case matrix shared by scripts/generate_golden.py (oracle) and tests/test_golden.py.

Each case is the argument list handed to ``svafotate annotate`` after ``-v <vcf> -o <out>``. Paths in
``ARGS`` are relative to ``tests/data`` and are resolved by the runner. Every case writes ``out.vcf`` to
its own directory; cases using ``-u`` also produce ``uniques.bed`` in that directory.
"""

MULTI = 'ref_multi_source.bed.gz'
GNOMAD3 = 'ref_gnomad_3rows.bed'
JOINT = 'joint_sv.vcf.gz'
EXT = 'extended.vcf'
TARGETS = 'targets.bed'

# name -> (query vcf, [annotate args])
CASES: dict[str, tuple[str, list[str]]] = {
    # the Talos production invocation, against the reduced-style single-source BED and the multi-source one
    'talos_gnomad_3rows': (JOINT, ['-b', GNOMAD3, '-s', 'gnomAD', '-f', '0.5', '-a', 'best', '--ins', '--cpu', '1']),
    'talos_multi': (JOINT, ['-b', MULTI, '-s', 'gnomAD', '-f', '0.5', '-a', 'best', '--ins']),
    # defaults: every source, minf 0.001, no extras
    'defaults_multi': (JOINT, ['-b', MULTI]),
    # extras
    'all_extras': (JOINT, ['-b', MULTI, '-a', 'all', '-f', '0.1']),
    'mis_full': (JOINT, ['-b', MULTI, '-a', 'mis', 'full']),
    'pops_mf_best': (JOINT, ['-b', MULTI, '-a', 'pops', 'mf', 'best', '-f', '0.5']),
    'single_pops': (JOINT, ['-b', MULTI, '-a', 'AFR', 'EAS', 'NFE']),
    # per-source overlap fractions and source subsetting
    'per_source_minf': (JOINT, ['-b', MULTI, '-s', 'gnomAD', 'CCDG', '-f', '0.5', '0.9', '-a', 'best']),
    'sources_topmed_thousg': (JOINT, ['-b', MULTI, '-s', 'TOPMed', 'ThousG']),
    # coverage / unique / size limit / targets
    'cov_uniq_targets': (JOINT, ['-b', MULTI, '-c', '0.01', '-u', '0.01', '-l', '1000000', '-t', TARGETS]),
    'cov_only': (JOINT, ['-b', MULTI, '-c', '0.0']),
    'uniq_only_lim': (JOINT, ['-b', MULTI, '-u', '0.05', '-l', '50000']),
    'targets_only': (JOINT, ['-b', MULTI, '-t', TARGETS]),
    # coordinate adjustments on the extended VCF
    'ext_defaults': (EXT, ['-b', MULTI, '-a', 'best']),
    'ext_ci_in': (EXT, ['-b', MULTI, '-ci', 'in', '-a', 'best']),
    'ext_ci_out': (EXT, ['-b', MULTI, '-ci', 'out', '-a', 'best']),
    'ext_ci95_in': (EXT, ['-b', MULTI, '-ci95', 'in', '-a', 'best']),
    'ext_ci95_out': (EXT, ['-b', MULTI, '-ci95', 'out', '-a', 'best']),
    'ext_emb': (EXT, ['-b', MULTI, '-e', '500', '-f', '0.5', '-a', 'best']),
    'ext_red': (EXT, ['-b', MULTI, '-r', '500', '-f', '0.5', '-a', 'best']),
    'ext_ins': (EXT, ['-b', MULTI, '--ins', '-f', '0.5', '-a', 'best', 'mis']),
    'ext_uniq_targets': (EXT, ['-b', MULTI, '-u', '0.0', '-t', TARGETS, '-c', '0.0']),
}
