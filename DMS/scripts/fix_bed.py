import pysam
mismatch_count = 0
import pandas as pd
import collections
import pickle
import numpy as np
import argparse


    



def main_func(input_file, output_file):

    # 1l beta
    
    
    file_ = open(input_file)
    file_2  = open(output_file, "w")
    
    for line_ in file_:
    
        
        file_2.write(line_.split("\t")[0])
        #file_2.write(line_.split("\t")[0].split("chr")[1])
        file_2.write("\t")

        
        if  line_.split("\t")[5][0] == "-":
        
            

            start = int(line_.split("\t")[1]) - 1
            file_2.write(str(start))
            file_2.write("\t")
            file_2.write(line_.split("\t")[2])
            
            
        elif line_.split("\t")[5][0] == "+":

            start = int(line_.split("\t")[1]) - 1
            file_2.write(str(start))
            file_2.write("\t")
            file_2.write(line_.split("\t")[2])
        
         
        file_2.write("\t")
        file_2.write(line_.split("\t")[3])
        file_2.write("\t")
        file_2.write(line_.split("\t")[4])
        file_2.write("\t")
        file_2.write(line_.split("\t")[5])
        
    file_.close()
    file_2.close()




if __name__ == '__main__':
    cmdline_parser = argparse.ArgumentParser('coverage check')

    cmdline_parser.add_argument('-f', '--input_file',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-o', '--output_file',
                                default="",
                                help='output file',
                                required = True,
                                type=str)
    args, unknowns = cmdline_parser.parse_known_args()
    main_func(args.input_file, args.output_file)


    
    
    
    
