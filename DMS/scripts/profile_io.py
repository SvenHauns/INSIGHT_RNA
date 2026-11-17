"""Validate aligned bulk inputs and publish profiles, counters and constraints."""
from pathlib import Path
import os
import pickle
import tempfile

from dms_common import CONSTRAINT_THRESHOLDS, add_counts, mutation_counters
from parallel_targets import atomic_output


def read_sequences(path):
    # Keep Biopython at the FASTA IO boundary so calculations import independently.
    from Bio import SeqIO
    return [str(record.seq).upper() for record in SeqIO.parse(path, 'fasta')]


def read_coverage_targets(path):
    """Read the historical TSV layout: ID, metric, chrom, ranges, blank, gene."""
    targets = []
    previous = None
    version = 1
    seen = set()
    with open(path) as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 5:
                raise ValueError(f'Coverage line {number}: missing target fields')
            original, _, chrom = fields[:3]
            version = version + 1 if original == previous else 1
            name = original if version == 1 else f'{original}_v{version}'
            previous = original
            if name in seen or Path(name).name != name or name in {'', '.', '..'}:
                raise ValueError(f'Duplicate or unsafe profile target ID: {name}')
            seen.add(name)
            try:
                intervals = [tuple(map(int, item.split('-'))) for item in fields[3:-1] if item]
                if not intervals or any(len(exon) != 2 or not 0 <= exon[0] < exon[1]
                                        for exon in intervals):
                    raise ValueError('invalid interval')
            except ValueError as exc:
                raise ValueError(f'Coverage line {number}: expected BED0 start-end ranges') from exc
            targets.append((name, chrom, intervals))
    return targets


def load_profile_tasks(coverage_list, bed_file, fasta_file, min_base_quality):
    targets = read_coverage_targets(coverage_list)
    sequences = read_sequences(fasta_file)
    with open(bed_file) as stream:
        bed = [line.rstrip('\r\n').split('\t') for line in stream if line.strip()]
    count = sum(len(intervals) for _, _, intervals in targets)
    if len(sequences) != count or len(bed) != count:
        raise ValueError('Coverage list, BED and FASTA exon counts disagree')

    tasks = []
    offset = 0
    for name, chrom, intervals in targets:
        exon_rows = bed[offset:offset + len(intervals)]
        if any(len(row) != 6 for row in exon_rows):
            raise ValueError(f'Target {name}: expected BED6 rows')
        strand = exon_rows[0][5]
        if strand not in {'+', '-'}:
            raise ValueError(f'Invalid strand for {name}: {strand}')
        exon_sequences = sequences[offset:offset + len(intervals)]
        for (start, end), row, sequence in zip(intervals, exon_rows, exon_sequences):
            if (row[0], int(row[1]), int(row[2]), row[5]) != (chrom, start, end, strand):
                raise ValueError('Coverage list and BED exon coordinates/strands disagree')
            if len(sequence) != end - start:
                raise ValueError(f'Target {name}: sequence length disagrees with BED interval')
        starts, ends = zip(*intervals)
        tasks.append((chrom, [list(starts)], [list(ends)], strand, name,
                      exon_sequences, min_base_quality))
        offset += len(intervals)
    return tasks


def publish_constraints(staged, output):
    """Replace current constraints and remove obsolete files for this output only."""
    for folder in staged.iterdir():
        destination = output.parent / folder.name
        destination.mkdir(exist_ok=True)
        current_names = {item.name for item in folder.iterdir()}
        for previous in destination.iterdir():
            if (previous.name.startswith(output.stem + '_')
                    and previous.name.endswith('_' + folder.name[4:] + '.txt')
                    and previous.name not in current_names):
                previous.unlink()
        for item in folder.iterdir():
            os.replace(item, destination / item.name)


def write_profiles(results, output_file, stable_file, changed_file, bam_label, total):
    """Stage aggregates until all target workers finish; preserve five-line layout."""
    stable, changed = mutation_counters()
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.dms-targets-', dir=output.parent) as temporary:
        staged = Path(temporary)
        for threshold in CONSTRAINT_THRESHOLDS:
            (staged / ('dms_' + threshold)).mkdir()
        with atomic_output(output) as stream:
            for number, (name, result) in enumerate(results, 1):
                profile, constraints, local_stable, local_changed = result
                stream.write(profile)
                add_counts(stable, local_stable)
                add_counts(changed, local_changed)
                for threshold, text in constraints.items():
                    filename = f'{output.stem}_{name}_{threshold}.txt'
                    (staged / ('dms_' + threshold) / filename).write_text(text)
                if number % 100 == 0:
                    print(f'{bam_label}: evaluated {number}/{total} profile targets', flush=True)
            publish_constraints(staged, output)
            for filename, counters in ((stable_file, stable), (changed_file, changed)):
                with atomic_output(filename, 'wb') as counters_file:
                    pickle.dump(counters, counters_file)
