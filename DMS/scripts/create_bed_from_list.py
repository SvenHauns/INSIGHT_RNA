"""Write BED6 in exactly the same target order as the coverage list."""
import argparse
import pickle
from pathlib import Path
from target_records import target_records


def main_func(input_string, input_folder):
    with open(input_string, 'rb') as stream:
        records = pickle.load(stream)
    output = Path(input_folder) / (Path(input_string).stem + '.bed')
    with output.open('w') as stream:
        for record in target_records(records):
            chrom, (starts, ends) = record[0]
            for start, end in zip(starts, ends):
                stream.write(f'{chrom}\t{start}\t{end}\t{record[1]}\t0\t{record[4]}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-f', '--folder', required=True, help='Output BED folder')
    parser.add_argument('-i', '--input', required=True, help='Candidate pickle')
    args = parser.parse_args()
    main_func(args.input, args.folder)


if __name__ == '__main__':
    main()
