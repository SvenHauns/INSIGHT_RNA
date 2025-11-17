import pandas as pd
import torch
import matplotlib.pyplot as plt
import numpy as np
import torch
import random
import time
import glob

    
def create_sliding_windows(sequence, window_size, step_length):
    """
    Create subsequences using a sliding window.

    Args:
        sequence (list or str): The input sequence.
        window_size (int): The size of each window.
        step_length (int): The step length for the sliding window.

    Returns:
        list: List of subsequences.
    """
    return [sequence[i:i + window_size] for i in range(0, len(sequence) - 50*3, step_length)]
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="cutting dataset")
    
    # Define arguments
    parser.add_argument("--full_dms_dataset", type=str, required=True, help="path to dataset")
    parser.add_argument("--cut_dataset", type=str, required=True, help="path to dataset")
    parser.add_argument("--window_size", type=int, required=True, help="window size of sequences")
    parser.add_argument("--step_size", type=int, required=True, help="step size of sequences")
    parser.add_argument("--coverage_limit", type=int, required=True, help="setting the limit of the coverage")
    
    # Parse arguments
    args = parser.parse_args()
    
    file_ = open(args.full_dms_dataset).readlines()
    file_out = open(args.cut_dataset, "w")

    check_cov = False
    len_list = []
    window_size = args.window_size
    step_length = args.step_size

    for enum, line_ in enumerate(file_):
        if line_[0] == ">":
    
            id_ = line_[:-1]
            seq = file_[enum+1][:-1]
            dms = file_[enum+2][1:-2]
            cov = file_[enum+3][:-1]
        
            dms = [float(d) for d in dms.split(", ") if d != ""]
            cov = [int(d) for d in cov.split(" ") if d != ""]
        
            seqs = create_sliding_windows(seq, window_size, step_length)
            dms_values = create_sliding_windows(dms, window_size, step_length)
            coverages = create_sliding_windows(cov, window_size, step_length)
        
            for enum2, s in enumerate(seqs):
                dms_string = ""
                for v in dms_values[enum2]:
                    dms_string = dms_string + str(v) + " "
                
                cont_val = False    
                cov_string = ""
                for v in coverages[enum2]:
                    cov_string = cov_string + str(v) + " "
                    if int(v) < args.coverage_limit: cont_val = True

                if cont_val: continue
                file_out.write(id_ + "_" +str(enum2) + "\n")
                file_out.write(s + "\n")
                file_out.write(dms_string + "\n")
