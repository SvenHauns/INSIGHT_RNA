
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

def get_exons(file_):

    file_ = open(file_).readlines()
    
    strand_dict =  {}
    exon_dict = {}
    
    for line_ in file_:
        id_ = line_.split("\t")[-1][1:-3]
        id_ = line_.split("\t")[3]

        
        if id_ not in strand_dict.keys():
            exon_dict[id_] = []
            strand_dict[id_] = line_.split("\t")[-1][:-1]

        exon_dict[id_].append(line_.split("\t")[1] + "-" + line_.split("\t")[2])


    return strand_dict, exon_dict

def write_bed_file(possible_regions, file):

    for region in possible_regions:
        file.write(f"{region[0]}\t{region[1]}\t{region[2]}\t{region[3]}_{region[4]}\t0\t{region[-1]}\n")


    return

def create_second_bed_file(input_file, exon_dict, strand_dict):

    possible_regions = []


    for enum, line_ in enumerate(input_file):
        

        if line_[0] == "t" or line_[0] == "#": continue
        chr_ = line_.split("\t")[0]
        start = line_.split("\t")[1]
        end = line_.split("\t")[2]

        print(line_)

        if int(start) > int(end):
            end = start
            start = int(line_.split("\t")[2])


        region_length = np.abs(int(end)-int(start))
        if region_length == 0: continue
        id_ = line_.split("\t")[-1][:-1]
        tid_ = id_.split("|")[0]
        valid_regions = exon_dict[tid_]
        strand = strand_dict[tid_]

        ### get possible regions
        print(valid_regions)
        for region in valid_regions:
            
            exon_start = int(region.split("-")[0])
            exon_end = int(region.split("-")[1])

            if exon_start > exon_end:
                exon_start = int(region.split("-")[1])
                exon_end = int(region.split("-")[0])


            start_region = exon_start

            while start_region < int(exon_end): 

                print("REGIONS")
                print(start_region)
                print(int(region.split("-")[1]))
                print(region_length)

                possible_region = [chr_, start_region, min(exon_end , start_region + region_length), enum, id_, strand]
                if possible_region[1] > int(end) or possible_region[2] > int(start):
                    possible_regions.append(possible_region)
                start_region =  start_region + min(20, region_length)

        possible_regions.append([chr_, start,end, enum, id_ + "_hit", strand])


    return possible_regions

if __name__ == '__main__':
    cmdline_parser = argparse.ArgumentParser('extract bed file windows')

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
    cmdline_parser.add_argument('-f', '--exons',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-g', '--track_output_cov',
                                default="",
                                help='track_output_cov',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-i', '--track_output_p_cov',
                                default="",
                                help='track_output_p_cov',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-j', '--track_output_p2_cov',
                                default="",
                                help='track_output_2_cov',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-k', '--positive_tracks_cov',
                                default="",
                                help='positive_tracks_cov',
                                required = True,
                                type=str)





    args, unknowns = cmdline_parser.parse_known_args()

    strand_dict , exon_dict = get_exons(args.exons)


    file_ = open(args.track_output_p2).readlines()

    print(file_)
    print(args.track_output_p2)
    regions = create_second_bed_file(file_, exon_dict, strand_dict)
    file_out = open(args.track_output_p2_cov, "w")
    write_bed_file(regions, file_out)


    file_ = open(args.positive_tracks).readlines()
    regions = create_second_bed_file(file_, exon_dict, strand_dict)
    file_out = open(args.positive_tracks_cov, "w")
    write_bed_file(regions, file_out)

    file_ = open(args.track_output_p).readlines()
    regions = create_second_bed_file(file_, exon_dict, strand_dict)
    file_out = open(args.track_output_p_cov, "w")
    write_bed_file(regions, file_out)



    file_ = open(args.track_output).readlines()
    regions = create_second_bed_file(file_, exon_dict, strand_dict)
    file_out = open(args.track_output_cov, "w")
    write_bed_file(regions, file_out)
