
import pandas as pd
import argparse
import glob
import os
import numpy as np
    
def write_to_output(data, name, header):
    file_ = open(name, "w")

    for num in range(len(data[0])):
        print(data)
        print(data[0])
        print(len(data[0]))
        print(num)

        for a in data:
            print(a)
            print(num)
            print(a[num])
        file_.write(header[num] + ",")
        file_.write(",".join([str(a[num]) for a in data]))
        if [str(a[num]) for a in data][-1][:-1] != "\n": file_.write("\n")


    file_.close()


    return
    
if __name__ == '__main__':

    cmdline_parser = argparse.ArgumentParser('coverage check')


    cmdline_parser.add_argument('-c', '--track_file',
                                default="",
                                help='track_file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-n', '--deseq_data',
                                default="",
                                help='deseq data folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-e', '--deseq_data_corrected',
                                default="",
                                help='deseq data folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-o', '--meta_output',
                                default="",
                                help='meta_output',
                                required = True,
                                type=str)
                                
    args, unknowns = cmdline_parser.parse_known_args()
    
    fixed_coverage = open(args.track_file).readlines()
    coverage_dict = {}

    for enum, line_ in enumerate(fixed_coverage):
        print(line_)
        if enum <= 1: continue
        if line_.split("\t")[-1][:-1] not in coverage_dict.keys():
            coverage_dict[line_.split("\t")[-1][:-1]] = []

        coverage_dict[line_.split("\t")[-1][:-1]].append([line_.split("\t")[0],line_.split("\t")[1],line_.split("\t")[2]])


    for key in coverage_dict.keys():

        data_load = open(args.deseq_data + key + "_data.csv").readlines()
        data_load_header = data_load[0].split(",")
        file_header = [line_.split(",")[0] for line_ in data_load]
        file_header2 = []
        fixed_data = []
        header_dict = {}
        print(file_header)
        print(data_load_header)
        for num, data in enumerate(data_load_header):
            if data == "":continue
            found = False
            for region in coverage_dict[key]:
                if region[1] not in header_dict.keys(): header_dict[region[1]] = []
                print("region")
                print(header_dict)
                print(region)
                print(data)
                if region[1] == '8870850': raise NotImplementedError
                if region[1] == '8870450': raise NotImplementedError
                if region[1] == '8870849': raise NotImplementedError
                print(int(data) >= int(region[1]))
                print(int(data) <= int(region[2]))
                if int(data) >= int(region[1]) and int(data) <= int(region[2]): 
                    found = True
                    header_dict[region[1]].append([val.split(",")[num] for val in data_load])

            if found == False: 
                fixed_data.append([val.split(",")[num] for val in data_load])
                file_header2.append(data_load_header[num])


        for header in header_dict.keys():
            print("header merge")
            print(header_dict[header])
            if header_dict[header] == []: continue
            print(header_dict[header][0][0])
            to_append = [header_dict[header][0][0]]
            header_dict[header] =  np.array([h[1:] for h in header_dict[header]]).astype(np.float64)
            if len(header_dict[header]) > 1:
                print(to_append)
                print(header_dict[header])
                print(np.shape(header_dict[header]))
                print(np.mean(header_dict[header], axis = 0 ))


            if len(header_dict[header]) > 1: to_append.extend(list(np.mean(header_dict[header], axis =0 )))
            else:
                print("APPEND")
                print(list(header_dict[header]))
                to_append.extend(list(header_dict[header])[0])
                print(to_append)
            print(to_append)
            if len(to_append) != 7: raise NotImplementedError
            print("header append")
            fixed_data.append(to_append)
            file_header2.append(str(int(header)+1))
        print(header_dict)
        print(file_header2)
        print([f[0] for f in fixed_data])

        fixed_data_corrected = []
        for enum, data in enumerate(fixed_data):
            to_append = [file_header2[enum]]
            to_append.extend(data[1:])
            fixed_data_corrected.append(to_append)
        
        fixed_data = fixed_data_corrected

        fixed_data = sorted(fixed_data, key = lambda l:int(l[0]))
        #file_header2 = sorted(file_header2, key=lambda l: int(l))
        print(fixed_data)

        write_to_output(fixed_data, args.deseq_data_corrected + key + "_data.csv", file_header)


    print(file_header)
    rownames = file_header
    condition = [d[:-1] for d in file_header]


    df = pd.DataFrame(condition, rownames, ["condition"])
    df.to_csv(args.meta_output)



    """
    can the values be correct? in deseq_data
    0.2 10 epochs
    """
    #### create bed tool files
    #### sample random windows
    #### create deseq files
