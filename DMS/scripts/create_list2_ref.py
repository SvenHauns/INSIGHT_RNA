"""Export coverage lists without discarding distinct transcript targets."""
import argparse
import pickle
from target_records import target_records


def main_func(input_string, output_string, formated_output):
    with open(input_string, 'rb') as stream:
        records = pickle.load(stream)
    with open(output_string, 'w') as raw, open(formated_output, 'w') as formatted:
        for record in target_records(records):
            chrom, (starts, ends) = record[0]
            ranges = [f'{s}-{e}' for s, e in zip(starts, ends)]
            identity, gene, rate = str(record[1]), str(record[2]), float(record[9])
            raw.write('\t'.join([identity, f'{rate:.4f}', chrom, *ranges, '', gene]) + '\n')
            formatted.write(f'{identity:<20}{gene:<20}{rate:<20.4f}{chrom:<20}{" ".join(ranges):<20}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-f', '--input', required=True, help='Candidate pickle')
    parser.add_argument('-o', '--output', required=True, help='Coverage list')
    parser.add_argument('-x', '--formated_output', required=True, help='Formatted coverage list')
    args = parser.parse_args()
    main_func(args.input, args.output, args.formated_output)


if __name__ == '__main__':
    main()
