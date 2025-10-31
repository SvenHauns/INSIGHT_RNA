import argparse
import pickle
import numpy as np
import pickle
from pathlib import Path

def main_func(input_string, input_folder):

    with open(input_string, "rb") as f:
        loaded_dict = pickle.load(f)

    out_path = Path(input_folder) / (Path(input_string).stem + ".bed")

    last_sub = None
    with open(out_path, "w") as fout:
        for key in loaded_dict: 
            for sub in loaded_dict[key]:

                chrom = str(sub[0][0])
                starts = sub[0][1][0]
                ends   = sub[0][1][1]
                strand = sub[4]

                total_len = sum(abs(int(s) - int(e)) for s, e in zip(starts, ends))

                if total_len > 30 and sub[0] != last_sub:
                    last_sub = sub[0]
                    for s, e in zip(starts, ends):
                        fout.write("\t".join([chrom, str(s), str(e), "NaN", "0", strand]) + "\n")
                     
if __name__ == '__main__':

    cmdline_parser = argparse.ArgumentParser('coverage check')

    cmdline_parser.add_argument('-f', '--folder',
                                default="",
                                help='input pkl folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-i', '--input',
                                default="",
                                help='input pkl file',
                                required = True,
                                type=str)
                                
    args, unknowns = cmdline_parser.parse_known_args()   
    main_func(args.input, args.folder)
    