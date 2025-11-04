import pandas as pd
import argparse
import glob
import os
    
def initialize_lists(num_lists):

    hisat1 = open(num_lists[0])
    dicts_ = [{} for _ in range(len(num_lists))]
    columns = {}
    
    for enum, line in enumerate(hisat1):

        tid = str(line.split("\t")[3].split("|")[0] + "|" + line.split("\t")[3].split("|")[1])
        for dict_ in dicts_:
            dict_[tid] =  []
            columns[tid] =  []
            
    hisat1.close()
        
    return dicts_, columns
    
def create_counts(coverage_file, hisat1_counts, columns, set_columns=False, window_size=40):
    junction_list_left = []
    junction_list_right = []
    junction_list_middle = []
    end_last = False

    hisat1 = open(coverage_file).readlines()
    last_tid = None
    skip = False
    last_start = None

    for enum, line in enumerate(hisat1):
        fields = line.strip().split("\t")
        chr_ = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        count = float(fields[6])
        tid = f"{fields[3].split('|')[0]}|{fields[3].split('|')[1]}"
        tid1 = fields[3].split("|")[0]
        junc = fields[3].split("|")[-1]

        current_length = end - start

        if junc == "intra":
            if junction_list_left or junction_list_right:
                if not junction_list_left:
                    for half in junction_list_right:
                        junc_tid = half[3]
                        if set_columns:

                            if half[1] + half[0] < half[1]: raise NotImplementedError
                            if "." in str(half[1] + half[0]): raise NotImplementedError
                            columns[junc_tid].append(str(half[1]) + "-" + str(half[1] + half[0]))
                        hisat1_counts[junc_tid].append(half[2])
                elif not junction_list_right:
                    for half in junction_list_left:
                        junc_tid = half[3]
                        if set_columns:

                            if half[1] + half[0] < half[1]: raise NotImplementedError
                            if "." in str(half[1] + half[0]): raise NotImplementedError
                            columns[junc_tid].append(str(half[1]) + "-" + str(half[1] + half[0]))
                        hisat1_counts[junc_tid].append(half[2])
                else:
                    for half in junction_list_right:
                        junc_tid = half[3]

                        middle_window = [w for w in junction_list_middle if w[-2] == half[-1] and w[-3] == half[-2] + 1 and half[3] == w[3]] if junction_list_middle else []
                        left_window = [w for w in junction_list_left if w[-1] == half[-1] and w[-2] == half[-2] + 1 and half[3] == w[3]] if not middle_window else \
                                      [w for w in junction_list_left if w[-1] == middle_window[0][-1] and w[-2] == middle_window[0][-3] + 1 and w[3] == middle_window[0][3]]

                        if not left_window:
                            if set_columns:

                                if half[1] + half[0] < half[1]: raise NotImplementedError
                                if "." in str(half[1] + half[0]): raise NotImplementedError
                                columns[junc_tid].append((str(half[1]) + "-" + str(half[1] + half[0])))
                            hisat1_counts[junc_tid].append(half[0])
                            continue

                        combined_coverage = half[2] * half[0]/window_size + left_window[0][2] * left_window[0][0]/window_size
                        if middle_window:
                            combined_coverage += middle_window[0][2] * middle_window[0][0]/window_size

                        if set_columns:

                            start_pos = half[1]
                            end_pos = left_window[0][1] + left_window[0][0]

                            if int(end_pos) < int(start_pos): raise NotImplementedError
                            if "." in str(end_pos): raise NotImplementedError
                            columns[half[3]].append(str(start_pos)  + "-"  + str(end_pos))

                        hisat1_counts[half[3]].append(combined_coverage)

                        if middle_window:
                            assert half[0] + left_window[0][0] + middle_window[0][0] == window_size
                        else:
                            assert half[0] + left_window[0][0] == window_size

                        junction_list_left.remove(left_window[0])
                        if middle_window:
                            junction_list_middle.remove(middle_window[0])

            last_tid = tid
            if set_columns:

                if int(end) < int(start): raise NotImplementedError
                if "." in str(end): raise NotImplementedError
                columns[tid].append(str(start) + "-" + str(end))
            hisat1_counts[tid].append(count)
            junction_list_left = []
            junction_list_right = []
            junction_list_middle = []

        else:
            if last_start is None:
                last_start = start

            if "left" in junc:
                junc_id = int(junc.split("_")[-1])
                exon_id = int(junc.split("_")[2])
                junction_list_left.append([current_length, start, count, tid, exon_id, junc_id])
            elif "right" in junc:
                junc_id = int(junc.split("_")[-1])
                exon_id = int(junc.split("_")[2])
                junction_list_right.append([current_length, start, count, tid, exon_id, junc_id])
            elif "middle" in junc:
                junc_id_right = int(junc.split("_")[-1])
                junc_id_left = int(junc.split("_")[-2])
                exon_id = int(junc.split("_")[2])
                junction_list_middle.append([current_length, start, count, tid, exon_id, junc_id_left, junc_id_right])

            last_tid = tid
            last_start = start

    if junction_list_right:
        for half in junction_list_right:
            junc_tid = half[3]
            if set_columns:
                if half[1] + half[0] < half[1]: raise NotImplementedError
                
                if "." in str(half[1] + half[0]): raise NotImplementedError
                columns[junc_tid].append(str(half[1]) + "-" + str(half[1] + half[0]))
            hisat1_counts[junc_tid].append(half[2])


    return hisat1_counts, columns
     
    

    
def to_csv(columns, dcits, rownames, data_folder):

    final_stack = []
    
    for key in list(columns.keys()):
    

        for dict_ in dcits:
            final_stack.append(dict_[key])


        data = final_stack
        column_names = columns[key]
        name = key
        
        import numpy as np
        # Create the DataFrame


        df = pd.DataFrame(data, rownames, column_names)
        df.to_csv(data_folder + name + "_data.csv")
        final_stack = []
    
    return
    
    
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
    

    
    hisat_list = args.crosslinked_files.split(",")
    nc_list = args.non_crosslinked_files.split(",")

    all_lists = hisat_list.copy()
    all_lists.extend(nc_list)
    rownames = []
    condition = []
    dcits, columns = initialize_lists(all_lists)

    
    counter = 0
    set_columns = False
    
    for coverage_file in hisat_list:
        
        if counter == 0: set_columns = True
        else: set_columns = False
        
        hisat1_counts, columns = create_counts(coverage_file, dcits[counter], columns, set_columns = set_columns, window_size = 40)
        dcits[counter] = hisat1_counts

        counter = counter + 1
    
        rownames.append("CL" + str(counter))
        condition.append("CL")
    
    
    nc_counter = 0
    for coverage_file in nc_list:
    
        if counter == 0: set_columns = True
        else: set_columns = False
        
        hisat1_counts, columns = create_counts(coverage_file, dcits[counter], columns, set_columns = set_columns, window_size = 40)
        dcits[counter] = hisat1_counts

        counter = counter + 1
        nc_counter = nc_counter + 1
        
        rownames.append("NC" + str(nc_counter))
        condition.append("NC")

    
    
    to_csv(columns, dcits, rownames, args.data_folder)
    
    
    df = pd.DataFrame(condition, rownames, ["condition"])
    df.to_csv(args.meta_data_file)
    
