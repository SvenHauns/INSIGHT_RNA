"""Extract bulk DMS profiles using independent transcript workers."""
import argparse

from base_counts import count_interval
from checkpoints import checkpointed_results
from dms_common import (CONSTRAINT_THRESHOLDS, DEFAULT_MIN_BASE_QUALITY,
                        add_counting_arguments, add_counts, mutation_counters,
                        validate_counting_arguments)
from parallel_targets import initialize_bam, get_bam
from profile_io import load_profile_tasks, write_profiles


def mutation_fractions(coverage, mutations):
    """Return per-position mutation fractions (zero for unobserved positions)."""
    return [m / n if n else 0.0 for n, m in zip(coverage, mutations)]


def format_values(values):
    """Preserve trailing spaces expected by the downstream profile readers."""
    return ''.join(str(value) + ' ' for value in values) + '\n'


def build_profile(bam, chrom, starts, ends, strand, name, sequences,
                  min_base_quality=DEFAULT_MIN_BASE_QUALITY):
    """Count one transcript and return its profile, constraints and counters."""
    if len(starts) != len(ends) or not (len(starts) == len(sequences)):
        raise ValueError('Mismatched exon/sequence counts')
    intervals = sorted(zip(starts, ends, sequences), key=lambda exon: exon[0],
                       reverse=(strand == '-'))
    stable, changed = mutation_counters()
    rates, raw, depths, sequence_parts = [], [], [], []
    for start, end, sequence in intervals:
        if len(sequence) != end - start:
            raise ValueError('Sequence length disagrees with BED interval')
        counts, mutations, exon_stable, exon_changed = count_interval(
            bam, chrom, start, end, strand, min_base_quality)
        add_counts(stable, exon_stable)
        add_counts(changed, exon_changed)
        if strand == '-':
            counts.reverse()
            mutations.reverse()
        rates.extend(mutation_fractions(counts, mutations))
        raw.extend(mutations)
        depths.extend(counts)
        # The FASTA extractor already reverse-complements each minus-strand exon.
        sequence_parts.append(sequence)
    sequence = ''.join(sequence_parts)
    header = f'>{name}\n{sequence}\n'
    profile = header + format_values(rates) + format_values(raw) + format_values(depths)
    constraints = {}
    for threshold in CONSTRAINT_THRESHOLDS:
        calls = ['x' if rate > float(threshold) and base in 'AC' else '.'
                 for rate, base in zip(rates, sequence)]
        constraints[threshold] = header + format_values(calls)
    return profile, constraints, stable, changed


def evaluate_profile(task):
    chrom, starts, ends, strand, name, sequences = task[:6]
    min_base_quality = task[6] if len(task) > 6 else DEFAULT_MIN_BASE_QUALITY
    try:
        if len(starts) != 1 or len(ends) != 1:
            raise ValueError('Expected one transcript per profile task')
        return name, build_profile(get_bam(), chrom, starts[0], ends[0], strand,
                                   name, sequences, min_base_quality)
    except Exception as exc:
        raise RuntimeError(f'Profile target {name} ({chrom}): {exc}') from exc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coverage_list', required=True)
    parser.add_argument('--extracted_sequences', required=True)
    parser.add_argument('--bed_file', '--fixed_bed_file', dest='bed_file', required=True,
                        help='BED6 with 0-based, half-open coordinates')
    parser.add_argument('--bam_file', required=True)
    # Retain the existing CLI argument; plotting is no longer part of bulk profiling.
    parser.add_argument('--save_path', help='Legacy plotting path (accepted for compatibility)')
    parser.add_argument('--output_file', required=True)
    parser.add_argument('--stable_dict', required=True)
    parser.add_argument('--changed_dict', required=True)
    add_counting_arguments(parser)
    args = parser.parse_args()
    validate_counting_arguments(parser, args)
    tasks = load_profile_tasks(args.coverage_list, args.bed_file, args.extracted_sequences,
                               args.min_base_quality)
    results = checkpointed_results(
        evaluate_profile, tasks, args.workers, initialize_bam, (args.bam_file,),
        stage='bulk-profiles', inputs=[args.bam_file, args.bam_file + '.bai',
                                      args.coverage_list, args.bed_file, args.extracted_sequences],
        checkpoint_dir=args.checkpoint_dir, resume=not args.no_resume)
    write_profiles(results, args.output_file, args.stable_dict, args.changed_dict,
                   args.bam_file, len(tasks))


if __name__ == '__main__':
    main()
