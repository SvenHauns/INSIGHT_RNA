import subprocess as sp
import glob
import os
import argparse
import numpy as np
import shlex

def read_dps_file(file, seq_len):

    file_ = open(file).readlines()
    


    start_seq = "%start of base pair probability data\n"
    end_seq = "showpage\n"
    
    start = 0
    end = 0
    for enum, line in enumerate(file_):
    
        if line == start_seq:
            start = enum +1
        elif line == end_seq:
            end = enum -1
            
    probabilities = file_[start:end+1]
    

    
    prob_dict = {}

    for line in probabilities:
        if line.split(" ")[-1] == "lbox\n": continue
        if int(line.split(" ")[0])-1 not in prob_dict.keys():
            prob_dict[int(line.split(" ")[0])-1] = []
        prob_dict[int(line.split(" ")[0])-1].append([int(line.split(" ")[1])-1, float(line.split(" ")[2])])
        
        if int(line.split(" ")[1])-1 not in prob_dict.keys():
            prob_dict[int(line.split(" ")[1])-1] = []
        prob_dict[int(line.split(" ")[1])-1].append([int(line.split(" ")[0])-1, float(line.split(" ")[2])])
            
    ### create matrix
    
    matrix = np.zeros((seq_len, seq_len))
    
    for key in prob_dict:
        for line in prob_dict[key]:
            matrix[key-1][line[0]-1] = line[1]
            
    


    return matrix
    
    
def elementwise_euclidean_distance(matrix1, matrix2):

    
    return np.mean(np.sqrt((matrix1 - matrix2) ** 2))
    
if __name__ == "__main__":


    cmdline_parser = argparse.ArgumentParser('retrieve RNA seq varna structure')

    cmdline_parser.add_argument('-f', '--target_folder',
                                default="",
                                help='target folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-y', '--target_folder_3utr',
                                default="",
                                help='target folder 3utr',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-z', '--target_folder_5utr',
                                default="",
                                help='target folder 5utr',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-o', '--output_summary',
                                default="",
                                help='output summary',
                                required = True,
                                type=str)
                                
                                
    args, unknowns = cmdline_parser.parse_known_args()
    print("starting run")
    
    
    target_folder = glob.glob(args.target_folder_5utr + "/*")
    varna_files_created = []
    
    mean1 = []
    mean2 = []
    
    evaluation = []
    
    for folder in target_folder:

        if os.path.isfile(folder) == True: continue
        
        plain_file = folder + "/plain_dp.ps"
        dms_file = folder + "/dms_dp.ps"
        oops_file = folder + "/oops_single_dp.ps"
        eclip_eval_file =  folder + "/eclip_dp.ps"
        
        
        plain_file_log = folder + "/plain_RNAfold.log"
        dms_file_log = folder + "/dms_RNAfold.log"
        eclip_file_log = folder + "/eclip_RNAfold.log"
        
        plain_struc = open(plain_file_log).readlines()[2].split(" ")[0]
        dms_struc = open(dms_file_log).readlines()[2].split(" ")[0]
        eclip_struc = open(eclip_file_log).readlines()[2].split(" ")[0]
        
        
        

        if os.path.isfile(plain_file) == False: continue
        if os.path.isfile(dms_file) == False: continue
        if os.path.isfile(oops_file) == False: continue
        if os.path.isfile(eclip_eval_file) == False: continue
        
        eclip_file = open(folder + "/" + "eclip_regions_5utr.txt").readlines()
        info_file = open(folder + "/" + "complete_info_5utr.txt").readlines()
        seq = info_file[1][:-1]
        dms_values =  [1 if float(i) > 0.06 else 0 for i in info_file[1][:-1].split(" ")]


        
        #### get relevant regions
        relevant_regions = []
        
        for line_ in eclip_file:
            relevant_regions.append([line_.split("\t")[5], line_.split("\t")[6]])
            
        ### create mask            
        mask = [0 for _ in range(0, len(seq))]
        

        
        for relevant_region in relevant_regions:
            mask[int(relevant_region[0])-1: int(relevant_region[1])-1] = [1 for _ in range(int(relevant_region[0])-1, int(relevant_region[1])-1)]

        ### get matrices
        
        plain_matrix = read_dps_file(plain_file, len(seq))[mask,:]
        dms_matrix = read_dps_file(dms_file, len(seq))[mask,:]
        oops_matrix = read_dps_file(oops_file, len(seq))[mask,:]
        
        eclip_full_matrix = read_dps_file(eclip_eval_file, len(seq))
        plain_full_matrix = read_dps_file(plain_file, len(seq))
        dms_full_matrix = read_dps_file(dms_file, len(seq))
        
        distance_plain_dms = elementwise_euclidean_distance(plain_full_matrix, dms_full_matrix) 
        distance_eclip_dms = elementwise_euclidean_distance(eclip_full_matrix, dms_full_matrix) 
        
        evaluation.append(str(folder) + "\t" + str(distance_plain_dms) + "\t" + str(distance_eclip_dms))
        
        if distance_plain_dms != 0:
            mean1.append(distance_plain_dms)
            mean2.append(distance_eclip_dms)
            
        if distance_plain_dms != 0:
            if distance_plain_dms > distance_eclip_dms:
                print(str(folder) + "\t" + str(distance_plain_dms) + "\t" + str(distance_eclip_dms))
                print(distance_plain_dms - distance_eclip_dms)
            else:
                print("wrong way")
                print(str(folder) + "\t" + str(distance_plain_dms) + "\t" + str(distance_eclip_dms))
                
                
        ### evaluate pure counts
        
        plain_sum = sum([1 for enum,e in enumerate(dms_values) if e == 1 and plain_struc[enum] == "."])
        eclip_sum = sum([1 for enum, e in enumerate(dms_values) if e == 1 and eclip_struc[enum] == "."])
        
        if eclip_sum > plain_sum:
            print("sum distance")
            print(str(folder) + "\t" + str(plain_sum) + "\t" + str(eclip_sum))
        
        
        ### calculate distances
        distance1 = elementwise_euclidean_distance(plain_matrix, oops_matrix) 
        distance2 = elementwise_euclidean_distance(dms_matrix, oops_matrix)
        
        evaluation.append(str(folder) + "\t" + str(distance1) + "\t" + str(distance2))
        
        

        
#print(evaluation)
#print(np.mean([float(e.split("\t")[1]) for e in evaluation]))
#print(np.mean([float(e.split("\t")[2]) for e in evaluation]))
        
print("means")
print(np.mean(mean1))
print(np.mean(mean2))
        

file_ = open(args.output_summary, "w")
file_.write("\n".join(evaluation))
file_.close()
