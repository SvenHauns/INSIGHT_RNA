import subprocess as sp
import glob
import os
import argparse
import shlex


def run_varna(sequence, structure, output, highlight_region, dms_color_map):


	
	if highlight_region != []:

	    cmd = 'java -Djava.awt.headless=true -Djava.io.tmpdir=/mnt/tmp/ -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "{sequnence}" -structureDBN "{structure}" -exportFormat png -o {output} -resolution 15.0 -auxBPs "0.6:type=B,anchor=7,size=5,color=#00FF00" -colorMapStyle "0:#FFFFFF;1:#FF0000"  -colorMapMax 1.0 -colorMap "{dms_color_map}" -highlightRegion {highlight_region}'
	    cmd = cmd.format(sequnence=sequence, structure=structure, highlight_region = highlight_region, dms_color_map = dms_color_map, output = output)
	else:
	
	    cmd = 'java -Djava.awt.headless=true -Djava.io.tmpdir=/mnt/tmp/ -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "{sequnence}" -structureDBN "{structure}" -exportFormat png -o {output} -resolution 15.0 -auxBPs "0.6:type=B,anchor=7,size=5,color=#00FF00" -colorMapStyle "0:#FFFFFF;1:#FF0000"  -colorMapMax 1.0 -colorMap "{dms_color_map}"'
	    cmd = cmd.format(sequnence=sequence, structure=structure, dms_color_map = dms_color_map, output = output)
	
	log_file = "varna.log"
	print("COMMAND")
	print(cmd)
	with open(log_file, 'w') as lf:
		sp.call(shlex.split(cmd), stdout=lf)
        
    
	return
    
    

    
    
if __name__ == "__main__":


    cmdline_parser = argparse.ArgumentParser('retrieve RNA seq varna structure')

    cmdline_parser.add_argument('-f', '--target_folder',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-y', '--output_summary',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-z', '--version',
                                default="",
                                help='output html',
                                required = True,
                                type=int)
                                
    args, unknowns = cmdline_parser.parse_known_args()
    print("starting run")
    target_folder = glob.glob(args.target_folder + "/*")
    varna_files_created = []
    
    print(target_folder)
    
    append = ""
    if args.version == 2:
         append = "_5utr"
    elif args.version == 3:
         append = "_3utr"
    
    print(append)
    
    for folder in target_folder:
        folder_files = glob.glob(folder + "/*log*")
        print(folder_files)
    
        for log_file_name in folder_files:
            print(log_file_name)
            log_file = open(log_file_name).readlines()
            print("LOG_FILE")
            print("LOG_FILE")
            print("LOG_FILE")
            print("LOG_FILE")
            print("LOG_FILE")
            print(log_file)
            if log_file == []: continue
            id_ = log_file[0][:-1]
            seqes = log_file[1][:-1]
            print("seq")
            print(seqes)
            print(len(seqes))
            struc = log_file[2][:-1].split(" ")[0]
            print(struc)
            print(len(struc))
            output_file = log_file_name.split(".log")[0] + "_varna.png"
            print(folder + "/complete_info" + str(append) + ".txt")
            info_file = open(folder  + "/complete_info" + str(append) + ".txt").readlines()
            print("dms")
            print(info_file[0][:-1].split(" "))
            print(len(info_file[0][:-1].split(" ")))

            dms_color_map = ",".join(info_file[0][:-1].split(" "))
            highlight_region = ""
            
            
        
        
            signal_list = []     
            sequence_list = []     
            print("info_file")
            print(folder  + "/complete_info" + str(append) + ".txt")
            print(info_file)
            print(info_file[1][:-1].split(" "))
            print(len(info_file[1][:-1].split(" ")))
            for enum, signal in enumerate(info_file[1][:-2].split(" ")):
                if int(signal) == 0:
                    signal_list.append(enum)
                elif int(signal) == 1:
                    if signal_list != []:
                        sequence_list.append(signal_list)
                        signal_list = []
        
        
            high_lights = []
            for seq in sequence_list:
                start = seq[0]
                end = seq[-1]
                high_lights.extend([str(start+1) + "-" + str(end+1)])

            if high_lights != []: high_lights_str = ";".join(high_lights)
            else : high_lights_str = high_lights
            print(high_lights_str)
            varna_files_created.append(output_file)
            print("INPUT SEQ")
            print("INPUT SEQ")
            print("INPUT SEQ")
            print("INPUT SEQ")
            print("INPUT SEQ")
            print(seqes)
            run_varna(seqes, struc, output_file, high_lights_str, dms_color_map)

file_ = open(args.output_summary, "w")
file_.write("\n".join(varna_files_created))
file_.close()
