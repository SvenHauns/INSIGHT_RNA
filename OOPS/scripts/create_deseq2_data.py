
import pandas as pd
import argparse
import glob
import os
    
    
    
"""
 XDG_CACHE_HOME=/mnt/tmp/ CONDA_PKGS_DIRS=/mnt/tmp/conda_pkgs snakemake --cores 1 --use-conda  --config samples=NC-1-1.hisat_mine,NC-1-2.hisat_mine,NC-1-3.hisat_mine,NC-1-4.hisat_mine,400-1-1.hisat_mine,400-1-2.hisat_mine,400-1-3.hisat_mine,400-1-4.hisat_mine
    
    XDG_CACHE_HOME=/mnt/tmp/ CONDA_PKGS_DIRS=/mnt/tmp/conda_pkgs snakemake --cores 1 --use-conda  --config samples=OOPS-IL-1,OOPS-IL-2,OOPS-IL-3,OOPS-NOIL-1,OOPS-NOIL-2,OOPS-NOIL-3,RNA-IL-1,RNA-IL-2,RNA-IL-3,RNA-NOIL-1,RNA-NOIL-2,RNA-NOIL-3
    

XDG_CACHE_HOME=/mnt/tmp/ CONDA_PKGS_DIRS=/mnt/tmp/conda_pkgs snakemake --cores 1 --use-conda  --config samples=OOPS-NOIL-1,OOPS-NOIL-2,OOPS-NOIL-3,RNA-NOIL-1,RNA-NOIL-2,RNA-NOIL-3 nc=OOPS-NOIL-1,OOPS-NOIL-2,OOPS-NOIL-3 cl=RNA-NOIL-1,RNA-NOIL-2,RNA-NOIL-3


"""
    
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
def create_counts(coverage_file, hisat1_counts, columns, set_columns = False, window_size = 40):

    junction_list_left = []
    junction_list_right = []
    junction_list_middle = []
    end_last = False
    
    hisat1 = open(coverage_file).readlines()
    last_tid = None
    skip = False
    last_start = None
    
    print("create coutns")
    print(len(hisat1))
    
    for enum, line in enumerate(hisat1):
    
        #print(enum)
        chr_ = line.split("\t")[0]
        start = int(line.split("\t")[1])
        end = int(line.split("\t")[2])
        count = float(line.split("\t")[6][:-1])
        #print(line.split("\t")[3])
        tid = str(line.split("\t")[3].split("|")[0] + "|" + line.split("\t")[3].split("|")[1])
        tid1 = str(line.split("\t")[3].split("|")[0])

        #if  "XM_047439107.1|SPHK2" == tid: raise NotImplementedError
        #if tid == "NM_001305275.2|AGRN": continue
        #if tid == "XM_024446454.2|ARHGEF16": continue
        #if tid == "XM_024446455.2|ARHGEF16": continue
        #if  "XM_047439107.1|SPHK2" != tid:continue
        #if enum > 10000: continue
        # NM_001401365.1|TTC7B
        #if  "XM_047439107.1|SPHK2" == tid:end_last = True
        #print(line)
        #if end_last == True and tid != "XM_047439107.1|SPHK2": raise NotImplementedError
        #print("junction_list_right")
        #print(junction_list_right)
        #print("junction_list_left")
        #print(junction_list_left)
        junc = line.split("\t")[3].split("|")[-1]
        
        if junc == "intra":
        
            if junction_list_left != [] or junction_list_right != []:
                
                
                if junction_list_left == []:
                
                    for enum, half in enumerate(junction_list_right):
                         junc_tid = half[3]
                         if set_columns: columns[junc_tid].append(half[1])
                         #print("87")
                         #print(columns[junc_tid])
                         hisat1_counts[junc_tid].append(half[2])

                elif junction_list_right == []:
                    for enum, half in enumerate(junction_list_left):
                         junc_tid = half[3]
                         if set_columns: columns[junc_tid].append(half[1])
                         #print("91")
                         #print(columns[junc_tid])
                         hisat1_counts[junc_tid].append(half[2])
                         
                         
                         

                else:
                    for enum, half in enumerate(junction_list_right):
                         
                         ## fetch fitting window 
                         #print("all")
                         #print(junction_list_right)
                         #print(junction_list_middle)
                         #print(junction_list_left)
                         junc_tid = half[3]
                         if junction_list_left == []:
                             if set_columns: columns[junc_tid].append(half[1])
                             #print("108")
                             #print(columns[junc_tid])
                             hisat1_counts[junc_tid].append(half[2])
                             continue
                         
                         #print("right")
                         #print(half)
                         
                         middle_window = []
                         if junction_list_middle != []:
                              middle_window = [window for window in junction_list_middle if window[-2] == half[-1] and window[-3] == half[-2]+1 and half[3] == window[3]]
                         
                         #print(middle_window)
                         if middle_window == []: left_window = [window for window in junction_list_left if window[-1] == half[-1] and window[-2] == half[-2]+1 and half[3] == window[3]]
                         else: left_window = [window for window in junction_list_left if window[-1] == middle_window[0][-1] and window[-2] == middle_window[0][-3]+1 and middle_window[0][3] == window[3]]
                         
                         
                         #print("left")
                         #print(left_window)
                         
                         if left_window == []: 
                             junc_tid = half[3]
                             if set_columns: columns[junc_tid].append(half[1])
                             #print("129")
                             #print(columns[junc_tid])
                             hisat1_counts[junc_tid].append(half[0])
                             continue ##### catch error here
                         
                         
                         if middle_window == []: combined_coverage = half[2] * half[0]/window_size + left_window[0][2] * left_window[0][0]/window_size
                         else: combined_coverage = half[2] * half[0]/window_size + left_window[0][2] * left_window[0][0]/window_size + middle_window[0][2] *middle_window[0][0]/window_size
                     
                         if set_columns: columns[half[3]].append(half[1])
                         #print("136")
                         #print(half)
                         #print(columns[half[3]])
                         hisat1_counts[half[3]].append(combined_coverage)

                         skip = False
                         last_start = None
                         
                         
                         if middle_window == []: assert half[0] + left_window[0][0] == window_size
                         else: 
                             #print(half)
                             #print(middle_window)
                             #print(left_window)
                             #print(half[0])
                             #print(left_window[0][0])
                             #print(middle_window[0][0])
                             assert half[0] + left_window[0][0] + middle_window[0][0] == window_size
                         
                         junction_list_left.remove(left_window[0])
                         if middle_window != []: junction_list_middle.remove(middle_window[0])

            last_tid = tid
                
            if set_columns: columns[tid].append(int(line.split("\t")[1]))
            #print(columns[tid])
            #print("160")
            
            hisat1_counts[tid].append(float(line.split("\t")[6][:-1]))
            junction_list_left = []
            junction_list_right = []
            junction_list_middle = []
        else:
            
             #if (last_tid == tid or last_tid == None):
             
            if last_start == None: last_start = int(line.split("\t")[1])
            current_length = end - start
            if junc.split("_")[1] == "left":
                junc_id =  int(junc.split("_")[-1])
                exon_id =  int(junc.split("_")[2])
                junction_list_left.append([current_length, int(line.split("\t")[1]), float(line.split("\t")[6][:-1]), tid, exon_id, junc_id])
            elif junc.split("_")[1] == "right":
               junc_id =  int(junc.split("_")[-1])
               exon_id =  int(junc.split("_")[2])
               junction_list_right.append([current_length, int(line.split("\t")[1]), float(line.split("\t")[6][:-1]),tid,  exon_id,  junc_id])
               #print("junction_list_right appended")
               #print(junction_list_right)
            elif junc.split("_")[1] == "middle":
                junc_id_right =  int(junc.split("_")[-1])
                junc_id_left =  int(junc.split("_")[-2])
                exon_id =  int(junc.split("_")[2])
                    
                junction_list_middle.append([current_length, int(line.split("\t")[1]), float(line.split("\t")[6][:-1]), tid, exon_id,  junc_id_left, junc_id_right])
                    
            last_tid = tid
            last_start = int(line.split("\t")[1])
             
            """
            else:

                for enum, half in enumerate(junction_list_left):
                    junc_id = half[3]
                    if set_columns: columns[junc_id].append(half[1])
                    print("193")
                    print(last_tid)
                    print(tid)
                    print(columns[junc_id])
                    hisat1_counts[junc_id].append(half[2])
                    skip = False
                    last_start = None
                    raise NotImplementedError
                         
                junction_list_left = []
                junction_list_right = []
                junction_list_middle = []
                last_tid_save = last_tid
                last_tid = tid
             """
        #if  "XM_047439107.1|SPHK2" == tid:print(columns[last_tid_save])
        #if  "XM_047439107.1|SPHK2" == tid:print(columns[tid])
        #if  "XM_047439107.1|SPHK2" == tid:print(last_tid_save)
        #if  "XM_047439107.1|SPHK2" == tid:raise NotImplementedError
        # NM_000979.4|RPL18
        
        
    if junction_list_right != []:
                
        for enum, half in enumerate(junction_list_right):
            junc_tid = half[3]
            if set_columns: columns[junc_tid].append(half[1])
            #print("87")
            #print(columns[junc_tid])
            hisat1_counts[junc_tid].append(half[2])
        
    return hisat1_counts, columns
    
    
    
def create_counts_old(coverage_file, hisat1_counts, columns, set_columns = False, window_size = 40):

    junction_list = []
    hisat1 = open(coverage_file)
    last_tid = None
    skip = False
    last_start = None
    
    for enum, line in enumerate(hisat1):
    
    
        chr_ = line.split("\t")[0] 
        start = int(line.split("\t")[1])
        end = int(line.split("\t")[2])
        count = int(line.split("\t")[4])

        tid = str(line.split("\t")[3].split("|")[0] + "|" + line.split("\t")[3].split("|")[1])
        
        
        junc = line.split("\t")[3].split("|")[-1]
        
        if junc == "intra":
        
            if junction_list != []:
                #print(len(junction_list)/2)
                second_half = junction_list[int(len(junction_list)/2):]
                if last_tid == tid: 
                    for enum, half in enumerate(junction_list[:int(len(junction_list)/2)]):
                
                         combined_coverage = half[2] * half[0]/window_size + second_half[enum][2] * second_half[enum][0]/window_size
                     
                         if set_columns: columns[tid].append(half[1])
                         hisat1_counts[tid].append(combined_coverage)
                         print(tid)
                         print(junction_list)
                         skip = False
                         last_start = None
                         assert half[0] + second_half[enum][0] == window_size
                else:
                
                    for enum, half in enumerate(junction_list):
                    
                         if set_columns: columns[tid].append(half[1])
                         hisat1_counts[tid].append(half[2])
                         skip = False
                         last_start = None
                
            last_tid = tid
                
            if set_columns: columns[tid].append(int(line.split("\t")[1]))
            hisat1_counts[tid].append(int(line.split("\t")[4]))
            print("junction info")
            print(tid)
            print(int(line.split("\t")[1]))
            
            junction_list = []
            
        else:


            if (last_tid == tid or last_tid == None) and skip == False or len(junction_list)%2 != 0:
                current_length = end - start
                print("append")
                print(current_length)
                print([current_length, int(line.split("\t")[1]), int(line.split("\t")[4])])
                if current_length <= 2: 
                    print("continue")
                    continue
                if current_length == 38 or current_length == 39 or current_length == 40:
                    if set_columns: columns[tid].append(half[1])
                    hisat1_counts[tid].append(half[2])
                else:
                    print("append")
                    if last_start == None: last_start = int(line.split("\t")[1])
                    if int(line.split("\t")[1]) > last_start + 40:
                        skip = True
                    
                    junction_list.append([current_length, int(line.split("\t")[1]), int(line.split("\t")[4])])
                    last_tid = tid
                    last_start = int(line.split("\t")[1])

            elif skip == True and len(junction_list)%2 != 0:
                for enum, half in enumerate(junction_list[:int(len(junction_list)/2)]):
                    combined_coverage = half[2] * half[0]/window_size + second_half[enum][2] * second_half[enum][0]/window_size
                    if set_columns: columns[tid].append(half[1])
                    hisat1_counts[tid].append(combined_coverage)
                    skip = False
                    last_start = None
                    assert half[0] + second_half[enum][0] == window_size
            
                junction_list = []
                last_start = None
                
                if current_length <= 2: 
                    print("continue")
                    continue
                if current_length == 38 or current_length == 39 or current_length == 40:
                    if set_columns: columns[tid].append(half[1])
                    hisat1_counts[tid].append(half[2])
                else:
                    print("append")
                    
                    if last_start == None: last_start = int(line.split("\t")[1])
                    if int(line.split("\t")[1]) > last_start + 40:
                        skip = True
                    
                    junction_list.append([current_length, int(line.split("\t")[1]), int(line.split("\t")[4])])
                    last_tid = tid
                    last_start = int(line.split("\t")[1])
                    
                    
            else:

                for enum, half in enumerate(junction_list):
                    if set_columns: columns[tid].append(half[1])
                    hisat1_counts[tid].append(half[2])
                    skip = False
                    last_start = None
                         
                junction_list = []
                last_tid = tid
            
    hisat1.close()
    
    if junction_list != []:
        for enum, half in enumerate(junction_list):
            if set_columns: columns[tid].append(half[1])
            hisat1_counts[tid].append(half[2])
            skip = False
            last_start = None
            
        junction_list = []
            
    return
    
def to_csv(columns, dcits, data_folder):

    final_stack = []
    
    for key in list(columns.keys()):
    

        for dict_ in dcits:
            final_stack.append(dict_[key])


        data = final_stack
        column_names = columns[key]
        name = key
        

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
    
    print("output")
    #print(dcits)
    print(all_lists)
    print(len(dcits))
    print(hisat_list)
    print(nc_list)
    print(len(all_lists))
    print(len(hisat_list))
    print(len(nc_list))

    
    
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

    
    
    to_csv(columns, dcits, args.data_folder)
    
    
    df = pd.DataFrame(condition, rownames, ["condition"])
    df.to_csv(args.meta_data_file)
    

    """
    check error:

    NM_001305275.2|AGRN|intra
    ['chr1', 1050771, 1050811, 0, 'NM_001305275.2|AGRN|intra']
    NM_001305275.2|AGRN|intra
    ['chr1', 1050791, 1050831, 0, 'NM_001305275.2|AGRN|intra']
    NM_001305275.2|AGRN|junction_right_0
    ['chr1', 1050811, 1050837, 0, 'NM_001305275.2|AGRN|junction_right_0']
    NM_001305275.2|AGRN|junction_right_1
    ['chr1', 1050831, 1050837, 0, 'NM_001305275.2|AGRN|junction_right_1']
    NM_001305275.2|AGRN|junction_right_0
    ['chr1', 1051031, 1051043, 0, 'NM_001305275.2|AGRN|junction_right_0']
    NM_001305275.2|AGRN|junction_left_0
    ['chr1', 1051252, 1051266, 0, 'NM_001305275.2|AGRN|junction_left_0']
    NM_001305275.2|AGRN|junction_left_1
    ['chr1', 1051252, 1051286, 0, 'NM_001305275.2|AGRN|junction_left_1']
    NM_001305275.2|AGRN|junction_left_2
    ['chr1', 1051266, 1051294, 0, 'NM_001305275.2|AGRN|junction_left_2']
    NM_001305275.2|AGRN|intra
    ['chr1', 1051274, 1051314, 0, 'NM_001305275.2|AGRN|intra']

NM_001305275.2|AGRN
XM_024446454.2|ARHGEF16

['chr1', 3481092, 3481113, 0, 'XM_024446454.2|ARHGEF16|junction_right_0']
XM_024446454.2|ARHGEF16|junction_right_1
['chr1', 3481112, 3481113, 0, 'XM_024446454.2|ARHGEF16|junction_right_1']
XM_024446455.2|ARHGEF16|junction_right_0
['chr1', 3467330, 3467350, 0, 'XM_024446455.2|ARHGEF16|junction_right_0']
XM_024446455.2|ARHGEF16|junction_right_1
['chr1', 3467350, 3467350, 0, 'XM_024446455.2|ARHGEF16|junction_right_1']
XM_024446455.2|ARHGEF16|junction_left_0
['chr1', 3468879, 3468899, 0, 'XM_024446455.2|ARHGEF16|junction_left_0']
XM_024446455.2|ARHGEF16|junction_left_1
['chr1', 3468879, 3468919, 0, 'XM_024446455.2|ARHGEF16|junction_left_1']
XM_024446455.2|ARHGEF16|junction_right_0
['chr1', 3468899, 3468936, 0, 'XM_024446455.2|ARHGEF16|junction_right_0']
XM_024446455.2|ARHGEF16|junction_right_1
['chr1', 3468919, 3468936, 0, 'XM_024446455.2|ARHGEF16|junction_right_1']
XM_024446455.2|ARHGEF16|junction_left_0
['chr1', 3469432, 3469435, 0, 'XM_024446455.2|ARHGEF16|junction_left_0']
XM_024446455.2|ARHGEF16|junction_left_1
['chr1', 3469432, 3469455, 0, 'XM_024446455.2|ARHGEF16|junction_left_1']
115248

NM_001025242.2|IRAK1
    """
