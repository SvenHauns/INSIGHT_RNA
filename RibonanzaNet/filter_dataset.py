import pandas as pd
import torch
import matplotlib.pyplot as plt
import numpy as np
import torch
import random
import time

def calculate_dms(seq, dms):
    
    dms_val = [float(d) for d in dms.split(" ") if d != ""]
    AC_values = np.quantile([val for enum, val in enumerate(dms_val) if (seq[enum] == "A" or seq[enum] == "C") and (val > 0)], 0.9)
    AC_normd = [val/AC_values if seq[enum] == "A" or seq[enum] == "C" else 0 for enum, val in enumerate(dms_val)]
    AC_normd = [min(u, 1) for u in AC_normd]
           
    return AC_normd

def check_coverage(cov, min_cov = 20):

    cov = [int(c) for c in cov[:-1].split(" ") if c != ""]
    check = sum([1 for c in cov if c < min_cov])
    
    if check/len(cov)  > 0.1: return False
    else: return True
    

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="cutting dataset")
    
    # Define arguments
    parser.add_argument("--not_normalized_dataset", type=str, required=True, help="path to dataset")
    parser.add_argument("--normalized_output", type=str, required=True, help="path to dataset")
    parser.add_argument("--min_cov", type=int, required=True, help="minimum coverage in sequences")
    # Parse arguments
    args = parser.parse_args()

    file_ = open(args.not_normalized_dataset)
    file_out = open(args.normalized_output, "w")
    
    counter = 0
    check_cov = False
    sample_counter = 0
    len_list = []

    for line_ in file_:

        line = line.strip()
    
        # Start new record
        if line.startswith(">"):
            id_ = line
            counter = 0
            check_cov = False
            continue
    
        # Parse fields in order
        if counter == 0:
            seq = line
        elif counter == 1:
            dms = line
        elif counter == 2:
            counts = line
        elif counter == 3:
            cov = line
            check_cov = True
    
        counter += 1

        if check_cov: 
            check = check_coverage(cov, min_cov =args.min_cov)
            if check: sample_counter = sample_counter + 1
            if check: len_list.append(len(seq))
            if check:
                ac = calculate_dms(seq[:-1], dms[:-1])
                ac_string = ""
            
                for a in ac:
                    ac_string = ac_string + str(a) + " "
            
                file_out.write(id_)
                file_out.write(seq)
                file_out.write(str(ac_string))
                file_out.write("\n")
                #file_out.write(cov)
