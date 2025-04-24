
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

def run_file(file_, file_out, table):

    print("run_file")
    info_list = []
    offset = 0

    for line_num, line_ in enumerate(file_):
        print(line_num)
        if line_[0] == "t" or line_[0] == "#":
            file_out.write(line_)
            offset = offset + 1
            continue

        id_ = line_.split("\t")[-1][:-1]
        print(args.deseq_output + id_ + "_stat_results.csv")
        deseq_data = pd.read_csv(args.deseq_output + id_ + "_stat_results.csv")

        chr_ = line_.split("\t")[0]
        start_ = line_.split("\t")[1]
        end_ =  line_.split("\t")[2]

        print(deseq_data.columns)

        regions = deseq_data['Unnamed: 0']
        start_found = False

        for enum, region in enumerate(regions):
            print("compare")
            print(region-1)
            print(start_)
            if int(region) - 1 == int(start_):
                info_list.append(deseq_data.loc[enum])
                start_found = True
            
        if start_found == False: raise NotImplementedError

        file_out.write(chr_)
        file_out.write("\t")
        file_out.write(start_)
        file_out.write("\t")
        file_out.write(end_)
        file_out.write("\n")

    print(file_[0:2])
    print(len(info_list))
    print(line_num)
    print(offset)
    print(line_num - 2 - offset)
    assert len(info_list)-1 == line_num - offset

    df = pd.DataFrame(info_list, columns = ["chr", "window_start", "window_end", "transcript_id", "gene", "baseMean", "padj", "lf"])
    df.to_csv(table)


    return


if __name__ == '__main__':
    cmdline_parser = argparse.ArgumentParser('extract bed file windows')

    cmdline_parser.add_argument('-a', '--deseq_output',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-b', '--track_output_p2',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-c', '--track_output_p',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-d', '--positive_tracks',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-e', '--track_output',
                                default="",
                                help='output html',
                                required = True,
                                type=str)

    cmdline_parser.add_argument('-f', '--track_output_p2_table',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-g', '--track_output_p_table',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-t', '--positive_tracks_table',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-i', '--track_output_table',
                                default="",
                                help='output html',
                                required = True,
                                type=str)

    cmdline_parser.add_argument('-o', '--output_folder',
                                default="",
                                help='output html',
                                required = True,
                                type=str)

    args, unknowns = cmdline_parser.parse_known_args()
    
    print("RUNNING")
    print("RUNNING")
    print("RUNNING")




    #file_ = open(args.track_output_p2).readlines()
    #file_out = open(args.output_folder + args.track_output_p2.split("/")[-1], "w")
    #run_file(file_, file_out, args.track_output_p2_table) 

    #### did i move the offset?

    #file_ = open(args.track_output_p).readlines()
    #file_out = open(args.output_folder + args.track_output_p.split("/")[-1], "w")
    #run_file(file_, file_out, args.track_output_p_table) 


    #file_ = open(args.positive_tracks).readlines()
    #file_out = open(args.output_folder + args.positive_tracks.split("/")[-1], "w")
    #run_file(file_, file_out, args.positive_tracks_table) 


    file_ = open(args.track_output).readlines()
    file_out = open(args.output_folder + args.track_output.split("/")[-1], "w")
    run_file(file_, file_out, args.track_output_table) 
