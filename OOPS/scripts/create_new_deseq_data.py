
import pandas as pd
import argparse
import glob
import os
from create_deseq2_data import initialize_lists, to_csv

def clean_dicts(dcits, columns):
    column_key = list(columns.keys())
    for key in column_key:
        if key.split("_")[-1] == "hit":
            del columns[key]
            add_key = "_".join(key.split("_")[:-1])
            columns[add_key] = []
            
    for subdict in dcits:
        subdict_key_list = list(subdict.keys())
        for subdict_key in subdict_key_list:
            if subdict_key.split("_")[-1] == "hit":
                del subdict[subdict_key]
                add_key = "_".join(subdict_key.split("_")[:-1])
                subdict[add_key] = []


    return dcits, columns

def create_counts_region(coverage_file, dict, columns, set_columns):

    hisat1 = open(coverage_file).readlines()
    print(dict.keys())
    for enum, line in enumerate(hisat1):
        chr_ = line.split("\t")[0]
        start = str(line.split("\t")[1])
        end = int(line.split("\t")[2])
        tid = str(line.split("\t")[3].split("|")[0] + "|" + line.split("\t")[3].split("|")[1])
        print(tid)
        hit_found = False
        if tid.split("_")[-1] == "hit":
            print("found hit")
            hit_found = True
            line_number = tid.split("_")[0]
            tid = str(line_number) + "_" + "_".join(tid.split("_")[1:-1])
            print(tid)


        coverage = float(line.split("\t")[-1])
        dict[tid].append(coverage)
        if set_columns: columns[tid].append(chr_ + "_" + str(start) +"_" + str(end) + "_" + str(hit_found))


    return dict, columns

if __name__ == '__main__':

    cmdline_parser = argparse.ArgumentParser('coverage check')


    cmdline_parser.add_argument('-c', '--crosslinked_files',
                                default="",
                                help='crosslinked files',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-n', '--non_crosslinked_files',
                                default="",
                                help='not crosslinked files',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-m', '--meta_data_file',
                                default="",
                                help='meta data file output',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-d', '--data_folder',
                                default="",
                                help='data folder',
                                required = True,
                                type=str)
                                
    args, unknowns = cmdline_parser.parse_known_args()
    
    
    print("input files")
    print(args.crosslinked_files)
    print(args.non_crosslinked_files)
    
    hisat_list = args.crosslinked_files.split(",")
    nc_list = args.non_crosslinked_files.split(",")
    
    print(hisat_list)
    print(nc_list)
    all_lists = hisat_list.copy()
    all_lists.extend(nc_list)
    rownames = []
    condition = []
    dcits, columns = initialize_lists(all_lists)
    dcits, columns = clean_dicts(dcits, columns)

    
    counter = 0
    set_columns = False
    
    for coverage_file in hisat_list:
        
        if counter == 0: set_columns = True
        else: set_columns = False
        
        hisat1_counts, columns = create_counts_region(coverage_file, dcits[counter], columns, set_columns = set_columns)
        dcits[counter] = hisat1_counts

        counter = counter + 1
    
        rownames.append("CL" + str(counter))
        condition.append("CL")
    
    
    nc_counter = 0
    for coverage_file in nc_list:
    
        if counter == 0: set_columns = True
        else: set_columns = False
        
        hisat1_counts, columns = create_counts_region(coverage_file, dcits[counter], columns, set_columns = set_columns)
        dcits[counter] = hisat1_counts

        counter = counter + 1
        nc_counter = nc_counter + 1
        
        rownames.append("NC" + str(nc_counter))
        condition.append("NC")


    print(rownames)
    print(dcits)
    print(columns)
    to_csv(columns, dcits, rownames, args.data_folder)
    
    
    df = pd.DataFrame(condition, rownames, ["condition"])
    df.to_csv(args.meta_data_file)
    
