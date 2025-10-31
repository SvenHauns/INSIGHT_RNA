import pysam
import pandas as pd
import pickle
import numpy as np
import argparse

def main_func(input_string, output_string, formated_output):

    with open(input_string, "rb") as f:
        loaded_dict = pickle.load(f)

    last_st, last_en = None, None
    with open(output_string, "w") as fout:
        for key in loaded_dict: 
            for sub in loaded_dict[key]:
                starts = sub[0][1][0]
                ends   = sub[0][1][1]

                if sum(abs(int(s) - int(e)) for s, e in zip(starts, ends)) <= 30:
                    continue

                if starts != last_st or ends != last_en:
                    header = [sub[1], f"{float(sub[8]):.4f}", str(sub[0][0])]
                    ranges = [f"{s}-{e}" for s, e in zip(starts, ends)]

                    line_parts = header + ranges + ["", str(sub[2])]
                    fout.write("\t".join(line_parts) + "\n")

                    last_st, last_en = starts, ends

    last_range_string = None
    with open(formated_output, "w") as f_out, open(output_string, "r") as f_in:
        for line in f_in:
            s = line.split()
            if not s:
                continue
            gene_name = s[-1]                 
            core = s[:-1]                     
            range_string = " ".join(core[3:]) 
            if last_range_string != range_string:
                f_out.write(f"{core[0]:<20}{gene_name:<20}{core[1]:<20}{core[2]:<20}{range_string:<20}\n")
            last_range_string = range_string
          
    
if __name__ == '__main__':


    cmdline_parser = argparse.ArgumentParser('coverage check')

    cmdline_parser.add_argument('-f', '--input',
                                default="",
                                help='input pkl file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-o', '--output',
                                default="",
                                help='output',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-x', '--formated_output',
                                default="",
                                help='formated output',
                                required = True,
                                type=str)

    args, unknowns = cmdline_parser.parse_known_args()
    main_func(args.input, args.output, args.formated_output)


    
    
    
