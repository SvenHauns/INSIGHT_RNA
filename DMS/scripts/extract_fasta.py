"""Extract strand-oriented BED6 sequences from the indexed preprocessing FASTA."""
import argparse

from parallel_targets import atomic_output
from seismic_clustering import reverse_complement


def extract_sequences(reference, bed, output):
    with open(bed) as intervals, atomic_output(output) as sequences:
        for number, line in enumerate(intervals, 1):
            if not line.strip():
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 6:
                raise ValueError(f'{bed}:{number}: expected BED6')
            chrom, start, end, _, _, strand = fields
            start, end = int(start), int(end)
            if start < 0 or end <= start or strand not in {'+', '-'}:
                raise ValueError(f'{bed}:{number}: invalid BED0 interval/strand')
            sequence = reference.fetch(chrom, start, end).upper()
            if len(sequence) != end - start:
                raise ValueError(f'{bed}:{number}: interval extends outside FASTA')
            if strand == '-':
                sequence = reverse_complement(sequence)
            sequences.write(f'>{chrom}:{start}-{end}({strand})\n{sequence}\n')


def main():
    import pysam
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--beds', nargs='+', required=True)
    parser.add_argument('--outputs', nargs='+', required=True)
    args = parser.parse_args()
    if len(args.beds) != len(args.outputs):
        parser.error('Each BED requires one output')
    with pysam.FastaFile(args.fasta) as reference:
        for bed, output in zip(args.beds, args.outputs):
            extract_sequences(reference, bed, output)


if __name__ == '__main__':
    main()
