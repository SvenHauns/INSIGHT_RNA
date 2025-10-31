import pysam
mismatch_count = 0
import pandas as pd
import collections
import pickle
import numpy as np
import argparse


def main_func(input_file, output_file):
    with open(input_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            parts = line.strip().split("\t")
            start = int(parts[1]) - 1
            fout.write("\t".join([parts[0], str(start), parts[2], parts[3], parts[4], parts[5]]) + "\n")



if __name__ == '__main__':
    cmdline_parser = argparse.ArgumentParser('coverage check')

    cmdline_parser.add_argument('-f', '--input_file',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-o', '--output_file',
                                default="",
                                help='output file',
                                required = True,
                                type=str)
    args, unknowns = cmdline_parser.parse_known_args()
    main_func(args.input_file, args.output_file)


    
    
    
    
