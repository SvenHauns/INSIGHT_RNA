"""Evaluate candidate coverage using independent transcript workers."""
import argparse
from collections import defaultdict
from pathlib import Path
import pickle

from annotation_regions import extract_refseq_utr
from base_counts import count_interval
from checkpoints import checkpointed_results
from dms_common import (CANONICAL_CHROMOSOMES, COVERAGE_SUFFIXES,
                        DEFAULT_MIN_BASE_QUALITY, add_counting_arguments,
                        add_counts, mutation_counters, validate_counting_arguments)
from parallel_targets import initialize_bam, get_bam, atomic_output


def metrics(stable, changed):
    def rate(bases):
        mutations = sum(value for key, value in changed.items() if key[0] in bases)
        total = mutations + sum(stable[base] for base in bases)
        return mutations / total if total else 0.0
    return tuple(rate(base) for base in 'ACGT') + (rate('AC'),)


def coverage_passes(depths):
    """Require <1%, <5% or <10% of bases below depth 20, respectively."""
    fraction_below = sum(depth < 20 for depth in depths) / len(depths)
    return tuple(fraction_below < limit for limit in (0.01, 0.05, 0.10))


def evaluate_candidate(task):
    """Return the historical four candidate dictionaries for one transcript."""
    name, chrom, strand, starts, ends, gene = task[:6]
    min_base_quality = task[6] if len(task) > 6 else DEFAULT_MIN_BASE_QUALITY
    collection = [defaultdict(list) for _ in COVERAGE_SUFFIXES]
    stable, changed = mutation_counters()
    depths = []
    try:
        for start, end in zip(starts, ends):
            if start == end:
                continue
            counts, _, exon_stable, exon_changed = count_interval(
                get_bam(), chrom, start, end, strand, min_base_quality)
            depths.extend(min(depth, 20) for depth in counts)
            add_counts(stable, exon_stable)
            add_counts(changed, exon_changed)
        if depths:
            record = [[chrom, [starts, ends]], name, gene, name, strand,
                      *metrics(stable, changed)]
            for selected, passed in zip(collection, (*coverage_passes(depths), True)):
                if passed:
                    selected[name].append(record)
        return collection
    except Exception as exc:
        raise RuntimeError(f"Candidate target {name} ({chrom}): {exc}") from exc


def candidate_tasks(annotation, run_type, min_base_quality):
    regions, genes, strands, chromosomes = extract_refseq_utr(annotation, run_type)
    for name, exons in regions.items():
        if name and exons and chromosomes[name] in CANONICAL_CHROMOSOMES:
            yield (name, chromosomes[name], strands[name], [s for s, _ in exons],
                   [e for _, e in exons], genes[name], min_base_quality)


def write_candidates(results, target_folder, bamfile):
    collection = [defaultdict(list) for _ in COVERAGE_SUFFIXES]
    for number, result in enumerate(results, 1):
        for aggregate, local in zip(collection, result):
            for name, records in local.items():
                aggregate[name].extend(records)
        if number % 100 == 0:
            print(f"{bamfile}: evaluated {number} candidate targets", flush=True)
    prefix = Path(bamfile).name.split('.bam')[0]
    for suffix, records in zip(COVERAGE_SUFFIXES, collection):
        with atomic_output(Path(target_folder) / f'{prefix}_{suffix}.pkl', 'wb') as stream:
            pickle.dump(records, stream)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff_path', required=True)
    parser.add_argument('--bamfile', required=True)
    parser.add_argument('--target_folder', required=True)
    parser.add_argument('--run_type', required=True, choices=('full', '5UTR', '3UTR'))
    add_counting_arguments(parser)
    args = parser.parse_args()
    validate_counting_arguments(parser, args)
    tasks = candidate_tasks(args.gff_path, args.run_type, args.min_base_quality)
    results = checkpointed_results(
        evaluate_candidate, tasks, args.workers, initialize_bam, (args.bamfile,),
        stage='bulk-candidates', inputs=[args.bamfile, args.bamfile + '.bai', args.gff_path],
        checkpoint_dir=args.checkpoint_dir, resume=not args.no_resume)
    write_candidates(results, args.target_folder, args.bamfile)


if __name__ == '__main__':
    main()
