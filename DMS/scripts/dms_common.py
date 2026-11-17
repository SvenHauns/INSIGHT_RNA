"""Shared bulk DMS settings, mutation counters and CLI validation."""

DEFAULT_MIN_BASE_QUALITY = 20
CONSTRAINT_THRESHOLDS = ('0.04', '0.06', '0.1', '0.25')
COVERAGE_SUFFIXES = ('99_20', '95_20', '90_20', 'all')
CANONICAL_CHROMOSOMES = frozenset(f'chr{i}' for i in range(1, 23)) | {'chrX', 'chrY'}


def mutation_counters():
    """Return independent match and substitution counters in transcript orientation."""
    stable = dict.fromkeys('ACGT', 0)
    changed = {f'{a}->{b}': 0 for a in 'ACGT' for b in 'ACGT' if a != b}
    return stable, changed


def add_counts(aggregate, counts):
    for key, value in counts.items():
        aggregate[key] += value


def add_counting_arguments(parser):
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--checkpoint-dir', default='.dms-checkpoints')
    parser.add_argument('--no-resume', action='store_true', help='Recompute targets without reading/writing checkpoints')
    parser.add_argument('--min-base-quality', type=int, default=DEFAULT_MIN_BASE_QUALITY)


def validate_counting_arguments(parser, args):
    if args.workers < 1:
        parser.error('--workers must be positive')
    if not 0 <= args.min_base_quality <= 93:
        parser.error('--min-base-quality must be between 0 and 93')
