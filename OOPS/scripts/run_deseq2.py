
mismatch_count = 0
import pandas as pd
import collections
import pickle
import numpy as np
import argparse
import subprocess as sp
import glob
import os
import pickle as pkl

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data



def run_deseq2(output_folder, input_folder, input_meta_file, run_checks = False):

    OUTPUT_PATH = output_folder
    deseq_list = []

    files = glob.glob(input_folder + "*_data.csv")
    meta_file_ = input_meta_file
    metadata_file = pd.read_csv(meta_file_, index_col=0)

    print(files)

    for file_ in files:
     
        #if "IGFL1" not in file_: continue
        
        try:
            print(file_)
            name = file_.split("/")[-1].split("_data.csv")[0]

    
            counts_df = pd.read_csv(file_, index_col=0)

            samples_to_keep = ~metadata_file.condition.isna()

            counts_df = counts_df.loc[samples_to_keep]

            metadata = metadata_file.loc[samples_to_keep]

            if run_checks == True:
                keep_sum = counts_df[counts_df.mean(axis=1) > 5].index
                if len(list(keep_sum))==0: continue
                common_elements1 = list(set(["NC1", "NC2", "NC3"]) & set(list(keep_sum)))
                common_elements2 = list(set(["CL1", "CL2", "CL3"]) & set(list(keep_sum)))
    
                if len(common_elements1) == 0: continue
                if len(common_elements2) == 0: continue

                counts_df = counts_df.loc[keep_sum]
                metadata = metadata.loc[keep_sum]


                genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 5]
                counts_df = counts_df[genes_to_keep]

            counts_df = counts_df.astype(int)
            print(counts_df)
            print("meta data file")
            print(metadata)


            inference = DefaultInference(n_cpus=8)
            dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design_factors="condition",
            refit_cooks=True,
            inference=inference,
            )


            dds.deseq2()


            
            stat_res = DeseqStats(dds, inference=inference, contrast = ["condition", "CL", "NC"], quiet = True)
            stat_res.summary()
            stat_res.results_df.to_csv(OUTPUT_PATH + name +  "_stat_results.csv")
            deseq_list.append(file_ + "\n")

        except:
            print("Exception")
    return deseq_list

if __name__ == '__main__':
    cmdline_parser = argparse.ArgumentParser('extract bed file windows')

    cmdline_parser.add_argument('-m', '--input_meta_file',
                                default="",
                                help='input meta file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-i', '--input_folder',
                                default="",
                                help='input folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-o', '--output_folder',
                                default="",
                                help='output folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-d', '--deseq_summary',
                                default="",
                                help='deseq summary file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-c', '--checks',
                                default=True,
                                help='deseq summary file',
                                required = False,
                                type=str)

    args, unknowns = cmdline_parser.parse_known_args()
    deseq_list = run_deseq2(args.output_folder, args.input_folder, args.input_meta_file, args.checks)
    
    print("WRITE")
    print(deseq_list)
    print(args.deseq_summary)
    file_ = open(args.deseq_summary, "w")
    file_.write("".join(deseq_list))
    file_.close()
    print("END")
