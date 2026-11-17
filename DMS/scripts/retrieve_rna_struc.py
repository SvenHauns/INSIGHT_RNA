"""Fold bulk profiles with and without A/C mutation-fraction constraints."""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile

from parallel_targets import atomic_output


def read_profiles(path):
    with open(path) as stream:
        while header := stream.readline():
            name = header.strip().removeprefix('>')
            sequence = stream.readline().strip().upper().replace('T', 'U')
            rates = [float(value) for value in stream.readline().split()]
            raw = stream.readline().split()
            coverage = stream.readline().split()
            if (not header.startswith('>') or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name)
                    or not sequence or not set(sequence) <= set('ACGUN')
                    or not len(sequence) == len(rates) == len(raw) == len(coverage)
                    or any(not 0 <= value <= 1 for value in rates)):
                raise ValueError(f'Invalid five-line profile: {name}')
            yield name, sequence, rates


def fold_profiles(profiles, output_dir, summary, threshold=0.06):
    if not 0 <= threshold <= 1:
        raise ValueError('Constraint threshold must be between 0 and 1')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    names = set()
    with atomic_output(summary) as table:
        table.write('target\tmode\tstructure\tenergy_kcal_mol\n')
        # Isolate any RNAfold auxiliary files from parallel jobs.
        with tempfile.TemporaryDirectory(prefix='dms-fold-') as temporary:
            for name, sequence, rates in read_profiles(profiles):
                if name in names:
                    raise ValueError(f'Duplicate profile ID: {name}')
                names.add(name)
                constraint = ''.join('x' if rate > threshold and base in 'AC' else '.'
                                     for rate, base in zip(rates, sequence))
                for mode in ('sequence', 'constrained'):
                    data = f'>{name}\n{sequence}\n'
                    command = ['RNAfold', '--noPS']
                    if mode == 'constrained':
                        data += constraint + '\n'
                        command.append('-C')
                    result = subprocess.run(command, input=data, text=True, capture_output=True,
                                            cwd=temporary, check=True)
                    match = re.search(r'^([.()]+)\s+\(\s*(-?\d+(?:\.\d+)?)\s*\)',
                                      result.stdout, re.MULTILINE)
                    if not match or len(match[1]) != len(sequence):
                        raise RuntimeError(f'Invalid RNAfold output for {name} ({mode})')
                    with atomic_output(output_dir / f'{name}.{mode}.dbn') as structure:
                        structure.write(result.stdout)
                    table.write(f'{name}\t{mode}\t{match[1]}\t{match[2]}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profiles', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--summary', required=True)
    parser.add_argument('--threshold', type=float, default=0.06)
    args = parser.parse_args()
    fold_profiles(args.profiles, args.output_dir, args.summary, args.threshold)


if __name__ == '__main__':
    main()
