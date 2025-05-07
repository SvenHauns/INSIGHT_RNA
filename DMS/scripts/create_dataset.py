import pysam
import pandas as pd
import pickle
import numpy as np
import pandas as pd
import argparse
import random

def main_func(input_file, output_file, output_file_500_train, output_file_500_test, output_file_800_train, output_file_800_test, cutoff_cov, cov_perc):

    file_ = open(input_file).readlines()

    filtered_samples = []

    for enum, line_ in enumerate(file_):
        
        if line_[0] == ">":
            seq = file_[enum+1]
            print("counter")
            print(len(file_))
            print(enum)
            cov = [int(c) for c in file_[enum+4].split(" ") if c != " " and c != "\n" ]


            dms_val =  [int(c) for c in file_[enum+3].split(" ") if c != " " and c != "\n" ]
        
        
            cov_check = [True if c >=cutoff_cov else False for c in cov]

            #if False in cov_check: continue
            if sum(cov_check)/len(cov_check) < cov_perc: continue
        


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



    file_ = open(output_file, "w")


    for sample in filtered_samples:
        file_.write(sample[0])
        file_.write(sample[1])
        file_.write(sample[-1])
        file_.write("\n")
    

    ###### write in cutoffs ######
    step = 200
    dataset_500 = []
    dataset_800 = []
    seq_500 = []
    seq_800 = []

    for sample in filtered_samples:
        for x in range(0, len(sample[1]), step):
            #if  sample[1][:-1][x:x+500] not in seq_500:
                dms_string = " ".join(sample[-1][:-1].split(" ")[x:x+500])
                #print(dms_string)
                #if sum([float(f) for f in dms_string]) == 0: continue
                dataset_500.append([sample[0][:-1] + "_" + str(x) + "\n", sample[1][:-1][x:x+500] + "\n", dms_string])
                seq_500.append(sample[1][:-1][x:x+500])
            #if  sample[1][:-1][x:x+800] not in seq_800:
                dms_string = " ".join(sample[-1][:-1].split(" ")[x:x+700])
                dataset_800.append([sample[0][:-1] + "_" + str(x) + "\n", sample[1][:-1][x:x+700] + "\n", dms_string])
                seq_800.append(sample[1][:-1][x:x+700])



    random.shuffle(dataset_500)
    random.shuffle(dataset_800)

    print("unique 500")

    uniq_keys, idx = np.unique(seq_500, return_index=True)
    dataset_500 = np.asarray(dataset_500, dtype=object)[idx]

    print("unique 800")
    uniq_keys, idx = np.unique(seq_800, return_index=True)
    dataset_800 = np.asarray(dataset_800, dtype=object)[idx]


    file_ = open(output_file_500_train, "w")
    file_test = open(output_file_500_test, "w")


    for sample in dataset_500[:int(len(dataset_500)*0.9)]:
        file_.write(sample[0])
        file_.write(sample[1])
        file_.write(sample[-1])
        file_.write("\n")


    for sample in dataset_500[int(len(dataset_500)*0.9):]:
        file_test.write(sample[0])
        file_test.write(sample[1])
        file_test.write(sample[-1])
        file_test.write("\n")


    file_ = open(output_file_800_train, "w")
    file_test = open(output_file_800_test, "w")


    for sample in dataset_800[:int(len(dataset_800)*0.9)]:
        file_.write(sample[0])
        file_.write(sample[1])
        file_.write(sample[-1])
        file_.write("\n")


    for sample in dataset_800[int(len(dataset_800)*0.9):]:
        file_test.write(sample[0])
        file_test.write(sample[1])
        file_test.write(sample[-1])
        file_test.write("\n")

    
if __name__ == '__main__':


    cmdline_parser = argparse.ArgumentParser('coverage check')

    cmdline_parser.add_argument('-f', '--input_file',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-o', '--output_file',
                                default="",
                                help='output',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-x', '--output_file_500_train',
                                default="",
                                help='output_file_500_train',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-y', '--output_file_800_train',
                                default="",
                                help='output_file_800_train',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-m', '--output_file_500_test',
                                default="",
                                help='output_file_500_test',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-n', '--output_file_800_test',
                                default="",
                                help='output_file_800_test',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-c', '--cutoff_cov',
                                default="",
                                help='cutoff_cov',
                                required = True,
                                type=int)
    cmdline_parser.add_argument('-g', '--cov_perc',
                                default="",
                                help='cov_perc',
                                required = True,
                                type=float)

    args, unknowns = cmdline_parser.parse_known_args()
    

    main_func(args.input_file, args.output_file, args.output_file_500_train, args.output_file_500_test, args.output_file_800_train, 
              args.output_file_800_test, args.cutoff_cov, args.cov_perc)


    
    