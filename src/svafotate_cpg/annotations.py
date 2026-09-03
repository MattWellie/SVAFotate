"""Which INFO fields an ``annotate`` run adds, and how each is filled.

The optional annotations requested with ``-a`` expand to *groups* (``Male``, ``AFR``, ``AFR_Female`` ...)
and every group contributes an ``AF``, ``Het`` and ``HomAlt`` column. Everything else is derived from that
one table, so adding a population is a one-line change.
"""

from dataclasses import dataclass, field

from cyvcf2 import VCF

from .reference import column_type

POPULATIONS = ('AFR', 'AMI', 'AMR', 'ASJ', 'EAS', 'EUR', 'FIN', 'MID', 'NFE', 'OTH', 'SAS')
SEXES = ('Male', 'Female')
BASE_COLUMNS = ('AF', 'Het', 'HomAlt', 'PopMax_AF')
GROUP_METRICS = ('AF', 'Het', 'HomAlt')

# the extra annotations that can be requested with -a, and what each one adds
EXTRA_ANNOTATIONS = {
    'all': 'All additional frequencies and counts and best matches',
    'mf': 'Male/Female frequencies and counts',
    'best': 'Best match ID, frequencies and counts',
    'pops': 'All populations frequencies and counts',
    **{pop: f'{pop} population frequencies and counts' for pop in POPULATIONS},
    'full': 'Full annotations for all SV matches',
    'mis': 'Mismatch SV counts and Best mismatch SV ID, frequencies and counts',
}


def requested_groups(extras: list[str]) -> list[str]:
    """Expand ``-a`` choices into the ordered list of column groups beyond the base columns."""
    everything = 'all' in extras
    mf = everything or 'mf' in extras
    groups: list[str] = list(SEXES) if mf else []
    for pop in POPULATIONS:
        if everything or 'pops' in extras or pop in extras:
            groups.append(pop)
            if mf:
                groups.extend(f'{pop}_{sex}' for sex in SEXES)
    return groups


def required_columns(extras: list[str]) -> list[str]:
    """Reference value columns a run with these ``-a`` choices reads (``AF`` is always needed for -c/-u)."""
    needed = list(BASE_COLUMNS) + [f'{g}_{m}' for g in requested_groups(extras) for m in GROUP_METRICS]
    if 'mis' in extras or 'all' in extras:
        needed += ['AF', 'Het', 'HomAlt']
    return list(dict.fromkeys(needed))


@dataclass
class AnnotationPlan:
    sources: list[str]
    extras: list[str]
    cov_af: float | None  # -c
    uniq_af: float | None  # -u
    targets_path: str | None  # -t
    bed_path: str
    max_columns: list[str] = field(init=False)

    def __post_init__(self) -> None:
        self.max_columns = list(BASE_COLUMNS) + [
            f'{g}_{m}' for g in requested_groups(self.extras) for m in GROUP_METRICS
        ]

    @property
    def best(self) -> bool:
        return 'best' in self.extras or 'all' in self.extras

    @property
    def mismatch(self) -> bool:
        return 'mis' in self.extras or 'all' in self.extras

    @property
    def full(self) -> bool:
        return 'full' in self.extras or 'all' in self.extras

    @property
    def coverage(self) -> bool:
        return self.cov_af is not None

    @property
    def unique(self) -> bool:
        return self.uniq_af is not None

    @property
    def targets(self) -> bool:
        return self.targets_path is not None

    def required_reference_columns(self) -> list[str]:
        return required_columns(self.extras)

    def add_headers(self, vcf: VCF) -> None:
        """Declare every INFO field this run may write, in the order upstream declared them."""
        source_list = ', '.join(self.sources)

        def add(field_id: str, description: str, field_type: str, number: str = '1') -> None:
            vcf.add_info_to_header({'ID': field_id, 'Description': description, 'Type': field_type, 'Number': number})

        for column in self.max_columns:
            add(
                f'Max_{column}',
                f'The maximum {column} from all matching SVs across all specified data sources ({source_list})',
                column_type(column),
            )

        for source in self.sources:
            add(f'{source}_Count', f'The number of matching SVs with {source}', 'Integer')

            if self.best:
                add(f'Best_{source}_ID', f'The {source} ID of the best matching SV', 'String')
                add(f'Best_{source}_OFP', f'The {source} OFP of the best matching SV', 'Float')
                for column in self.max_columns:
                    add(f'Best_{source}_{column}', f'The best {column} match for {source}', column_type(column))

            if self.mismatch:
                add(
                    f'{source}_Mismatches',
                    f'Comma-separated list of the {source} IDs of overlapping SVs with different SVTYPEs',
                    'String',
                    '.',
                )
                add(
                    f'{source}_Mismatches_Count',
                    f'The number of {source} overlapping SVs with different SVTYPEs',
                    'Integer',
                )
                add(
                    f'{source}_Mismatch_SVTYPEs',
                    f'Comma-separated list of the other overlapping SVTYPEs for {source}',
                    'String',
                    '.',
                )
                add(
                    f'Best_{source}_Mismatch_ID',
                    f'The {source} ID of the best overlapping SV with different SVTYPE',
                    'String',
                )
                add(
                    f'Best_{source}_Mismatch_OFP',
                    f'The {source} OFP of the best overlapping SV with different SVTYPE',
                    'Float',
                )
                add(
                    f'Best_{source}_Mismatch_SVTYPE',
                    f'The {source} SVTYPE of the best overlapping SV with different SVTYPE',
                    'String',
                )
                add(
                    f'Best_{source}_Mismatch_AF',
                    f'The {source} AF for the best overlapping SV with different SVTYPE',
                    'Float',
                )
                add(
                    f'Best_{source}_Mismatch_Het',
                    f'The {source} Het count for the best overlapping SV with different SVTYPE',
                    'Integer',
                )
                add(
                    f'Best_{source}_Mismatch_HomAlt',
                    f'The {source} HomAlt count for the best overlapping SV with different SVTYPE',
                    'Integer',
                )

            if self.full:
                add(
                    f'{source}_Matches',
                    f'Comma-separated list of each SV match with {source}, where all annotations are listed as SV_ID|SVTYPE|<value columns in BED order> as found in {self.bed_path}',
                    'String',
                    '.',
                )

            if self.coverage:
                add(
                    'SV_Cov',
                    f'The amount of the SV covered by matching SVs from {source_list} that have an AF greater than {self.cov_af}',
                    'Float',
                )
                add(
                    f'{source}_SV_Cov',
                    f'The amount of the SV covered by matching SVs from {source} that have an AF greater than {self.cov_af}',
                    'Float',
                )

            if self.unique:
                add(
                    'SV_Uniq',
                    f'The number of unique regions within the SV that are not found in {source_list} that have an AF greater than {self.uniq_af}',
                    'Integer',
                )

            if self.targets:
                add(
                    'Target_Overlaps',
                    f'The target overlaps found for the SV corresponding to targets listed in {self.targets_path}',
                    'String',
                    '.',
                )

            if self.targets and self.unique:
                add(
                    'Unique_Targets',
                    f'The target overlaps found for the unique region within the SV corresponding to targets listed in {self.targets_path}',
                    'String',
                    '.',
                )
