import pysam
import pandas as pd
import pickle
import numpy as np
import pandas as pd
import argparse



file_ = open("./5UTR_all/treated_5UTR.fa").readlines()
#file_ = open("./5UTR_all/untreated_full_5UTR.fa").readlines()

filtered_samples = []

for enum, line_ in enumerate(file_):

    if line_[0] == ">":
        seq = file_[enum+1]
        cov = [int(c) for c in file_[enum+4].split(" ") if c != " " and c != "\n" ]


        dms_val =  [int(c) for c in file_[enum+3].split(" ") if c != " " and c != "\n" ]
        
        
        cov_check = [True if c >=20 else False for c in cov]
        if 0 in cov: continue

        if False in cov_check: continue
        


        AC_values = np.quantile([ val for enum, val in enumerate(dms_val) if seq[enum] == "A" or seq[enum] == "C"], 0.9)
        
        A_values = np.quantile([ val for enum, val in enumerate(dms_val) if seq[enum] == "A"], 0.9)
        C_values = np.quantile([ val for enum, val in enumerate(dms_val) if seq[enum] == "C"], 0.9)
        
        
        AC_normd = [val/AC_values if seq[enum] == "A" or seq[enum] == "C" else 0 for enum, val in enumerate(dms_val)]
        unique_normd = []

        
        
        for enum, s in enumerate(seq[:-1]):
        
            if s == "A":
                unique_normd.append(dms_val[enum]/max(A_values,1))
            elif s == "C":
                unique_normd.append(dms_val[enum]/max(C_values,1))
            else:
                unique_normd.append(0.0)
        

        assert len(seq)-1 == len(unique_normd)
        AC_normd = [min(u,1) for u in AC_normd] 
        unique_normd = [min(u,1) for u in unique_normd] 
        unique_normd = " ".join([str(c) for c in unique_normd])

        


        filtered_samples.append([line_, seq, dms_val, cov, AC_normd, unique_normd])



file_ = open("./5UTR_all/5UTR_treated_normd.fa", "w")


for sample in filtered_samples:
    file_.write(sample[0])
    file_.write(sample[1])
    file_.write(sample[-1])
    file_.write("\n")
    
    
