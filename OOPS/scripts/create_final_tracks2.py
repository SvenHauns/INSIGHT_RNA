
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

def create_track_file(deseq_folder, track_file, track_table):

    header = 'track name=window description="window" color=0,255,0,\n#chrom chromStart chromEnd\n'
    track_file = open(track_file, "w")
    track_file.write(header)

    table_list = []
    miss_counter = 0

    for file_ in glob.glob(deseq_folder + "/*"):

        file_read = pd.read_csv(file_)
        columns = file_read.columns
        region = file_read["Unnamed: 0"].iloc[-1]
        if region.split("_")[-1] != "True": raise NotImplementedError

        result = file_read.iloc[-1]
        if result["padj"] < 0.05 and result["baseMean"]>50:
            table_list.append(result)
        else: 
            miss_counter = miss_counter + 1

        track_file.write(result["Unnamed: 0"].split("_")[0] + "\t" + result["Unnamed: 0"].split("_")[1] + "\t" + result["Unnamed: 0"].split("_")[2] + "\n")



    print("miss_counter")
    print(miss_counter)

    df = pd.DataFrame(table_list, columns = columns)
    df.to_csv(track_table)


    return

if __name__ == '__main__':
    cmdline_parser = argparse.ArgumentParser('extract bed file windows')

    cmdline_parser.add_argument('-a', '--deseq_folder_tracks',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-b', '--deseq_folder_tracks_p1',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-c', '--deseq_folder_tracks_p2',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-d', '--deseq_folder_positive',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    

    cmdline_parser.add_argument('-e', '--track_file',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-f', '--track_file_p1',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-g', '--track_file_p2',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-i', '--track_file_positive',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    
    cmdline_parser.add_argument('-j', '--track_table',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-k', '--track_table_p1',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-l', '--track_table_p2',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-m', '--track_table_positive',
                                default="",
                                help='output html',
                                required = True,
                                type=str)


    args, unknowns = cmdline_parser.parse_known_args()
    

    
    create_track_file(args.deseq_folder_tracks, args.track_file, args.track_table)
    create_track_file(args.deseq_folder_tracks_p1, args.track_file_p1, args.track_table_p1)
    create_track_file(args.deseq_folder_tracks_p2, args.track_file_p2, args.track_table_p2)
    create_track_file(args.deseq_folder_positive, args.track_file_positive, args.track_table_positive)


