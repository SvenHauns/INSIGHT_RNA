import pysam
import inspect
mismatch_count = 0
import argparse


import pandas as pd
import collections
import pickle
import numpy as np




def main_func(input_string, input_folder):



    with open(input_string, 'rb') as f:
        loaded_dict = pickle.load(f)





    file_ = open(input_folder + "/" + input_string.split("/")[-1].split(".pkl")[0] + ".bed", "w")   
    last_sub = None
        
           
    
    trancript_id_list = []
    
    for key_enum, key in enumerate(list(loaded_dict.keys())):

        key_list = loaded_dict[key]
        
        for sublist in key_list:
            trancript_id_list.append(sublist[2])
        
        




    for key_enum, key in enumerate(list(loaded_dict.keys())):

        key_list = loaded_dict[key]
        
        for sublist in key_list:
        
            t_id = sublist[2]
            t_id_index = trancript_id_list.index(t_id)
            

            length = [abs(int(start) - int(sublist[0][1][1][enum])) for enum, start in enumerate(sublist[0][1][0])]
            



            if sum(length) > 30 and sublist[0] != last_sub:

                last_sub = sublist[0]

                for enum, start in enumerate(sublist[0][1][0]):
                    file_.write(str(sublist[0][0]))
                    file_.write("\t")
                    file_.write(str(start))                
                    file_.write("\t")
                    file_.write(str(sublist[0][1][1][enum]))
                    file_.write("\t")
                    file_.write("NaN")
                    file_.write("\t")
                    file_.write("0")
                    file_.write("\t")
                    file_.write(sublist[4])
                    file_.write("\n")
                    
                    
    file_.close()

        

             
if __name__ == '__main__':

    cmdline_parser = argparse.ArgumentParser('coverage check')

    cmdline_parser.add_argument('-f', '--folder',
                                default="",
                                help='input pkl folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-i', '--input',
                                default="",
                                help='input pkl file',
                                required = True,
                                type=str)
                                
    args, unknowns = cmdline_parser.parse_known_args()   
    main_func(args.input, args.folder)


    
    
    
