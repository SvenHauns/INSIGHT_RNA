import pysam
import pandas as pd
import pickle
import numpy as np
import pandas as pd
import argparse



def main_func(input_string, output_string, formated_output):


    

    with open(input_string, 'rb') as f:
        loaded_dict = pickle.load(f)



    last_st = ""
    last_en = ""

    file_ = open(output_string, "w")   
    
    
    
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

            
            if sum(length) > 30:
           # if True:

                if sublist[0][1][0] != last_st or sublist[0][1][1] != last_en:
                    


                    file_.write(sublist[1])
                    file_.write("\t")
                    file_.write(str("{:.4f}".format(float(sublist[8]))))
                    file_.write("\t")                
                    file_.write(str(sublist[0][0]))
                    file_.write("\t")                
                

                

                
                    for enum, start in enumerate(sublist[0][1][0]):
                        file_.write(str(start))                
                        file_.write("-")
                        file_.write(str(sublist[0][1][1][enum]))
                        file_.write("\t")


                    file_.write("\t")
                    file_.write(str(t_id))
                    file_.write("\n")



                last_st = sublist[0][1][0]
                last_en = sublist[0][1][1]
    
    file_.close()
    
    max_len = 0
    with open(output_string, 'r') as f:
        for line in f:
            s = line.split()[:-1]
            if len(s) > max_len:
                max_len = len(s)
                

    max_len = 13
    max_len = 8
    last_list = []
    last_range_string = ""
    
    file_ = open(formated_output, "w")   
    with open(output_string, 'r') as f:
        for line in f:
            s = line.split()
            

            gene_name = s[-1]

            s = s[:-1]
            range_string = ""
            
            for x in range(3, len(s)):
                range_string = range_string + " " +s[x]
                

            
            
            s[3] =range_string
            

            if last_range_string != range_string:

            

                file_.write(f'{s[0]:<20}{gene_name:<20}{s[1]:<20}{s[2]:<20}{s[3]:<20}')
                file_.write("\n")
            
            last_range_string = range_string
            
            

    
if __name__ == '__main__':


    cmdline_parser = argparse.ArgumentParser('coverage check')

    cmdline_parser.add_argument('-f', '--input',
                                default="",
                                help='input pkl file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-o', '--output',
                                default="",
                                help='output',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-x', '--formated_output',
                                default="",
                                help='formated output',
                                required = True,
                                type=str)
                                ######## require input #####################################



    args, unknowns = cmdline_parser.parse_known_args()
    

    main_func(args.input, args.output, args.formated_output)


    
    
    
