import subprocess as sp
import glob
import os
import argparse
import numpy as np
import shlex
import pickle

def run_varna(sequence, structure, output, highlight_region, dms_color_map, oops_seq_region, annotation_command):

	if highlight_region != []:

	    cmd = 'java -Djava.awt.headless=true -Djava.io.tmpdir=/mnt/tmp/ -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "{sequnence}" -structureDBN "{structure}" -resolution 35.0  -exportFormat png  -o {output} -auxBPs "0.6:type=B,anchor=7,size=5,color=#00FF00" -colorMapStyle "0:#FFFFFF;1:#FF0000"  -colorMapMin 0.0 -colorMapMax 0.2 -colorMap "{dms_color_map}" -highlightRegion {highlight_region}'
	    
	    if dms_color_map == "":
	    
	        cmd = 'java -Djava.awt.headless=true -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd \
-sequenceDBN "{sequnence}" \
-structureDBN "{structure}" \
-exportFormat png \
-o {output} \
-resolution 35.0 \
-highlightRegion "{highlight_region}" \
-applyBasesStyle2on "{highlight_oops}" -basesStyle2 "fill=#000080,outline=#000000,label=#000000" \
-annotations "{annotation_command}" -flat true' 

	        cmd = cmd.format(sequnence=sequence, structure=structure, highlight_region = highlight_region, dms_color_map = dms_color_map, output = output, highlight_oops = oops_seq_region, annotation_command = annotation_command)
	        
	        
	    else: 
	    
	        if dms_color_map != "":
	    
	    
	            cmd = 'java -Djava.awt.headless=true -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd \
-sequenceDBN "{sequnence}" \
-structureDBN "{structure}" \
-exportFormat png \
-resolution 35.0 \
-o {output} \
-highlightRegion "{highlight_region}" \
-applyBasesStyle2on "{highlight_oops}" -basesStyle2 "fill=#000080,outline=#000000,label=#000000" \
-applyBasesStyle1on "{dms_color_map}" -basesStyle1 "fill=#FF1000,outline=#000000,label=#000000" \
-annotations "{annotation_command}" -flat true'

	            cmd = cmd.format(sequnence=sequence, structure=structure, highlight_region = highlight_region, dms_color_map = dms_color_map, output = output, highlight_oops = oops_seq_region, annotation_command = annotation_command)
	    
	        else:
	            cmd = 'java -Djava.awt.headless=true -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd \
-sequenceDBN "{sequnence}" \
-structureDBN "{structure}" \
-exportFormat png \
-resolution 35.0 \
-o {output} \
-highlightRegion "{highlight_region}" \
-applyBasesStyle2on "{highlight_oops}" -basesStyle2 "fill=#000080,outline=#000000,label=#000000" \
 -basesStyle1 "fill=#FF1000,outline=#000000,label=#000000" \
-annotations "{annotation_command}" -flat true'

	            cmd = cmd.format(sequnence=sequence, structure=structure, highlight_region = highlight_region, output = output, highlight_oops = oops_seq_region, annotation_command = annotation_command)
	            
	            
	else:
	    
	    if dms_color_map != "":

	        cmd = 'java -Djava.awt.headless=true -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd \
-sequenceDBN "{sequnence}" \
-structureDBN "{structure}" \
-exportFormat png \
-resolution 35.0 \
-o {output} \
-applyBasesStyle1on "{dms_color_map}" -basesStyle1 "fill=#FF1000,outline=#000000,label=#000000" -flat true'
	        cmd = cmd.format(sequnence=sequence, structure=structure, dms_color_map = dms_color_map, output = output)
	    else:

	        cmd = 'java -Djava.awt.headless=true -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd \
-sequenceDBN "{sequnence}" \
-structureDBN "{structure}" \
-exportFormat png \
-resolution 35.0 \
-o {output} \
 -basesStyle1 "fill=#FF1000,outline=#000000,label=#000000" -flat true'
	        cmd = cmd.format(sequnence=sequence, structure=structure, output = output)

	log_file = "varna.log"

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
    target_folder = glob.glob(args.target_folder + "/*")
    varna_files_created = []
    
    
    append = ""
    if args.version == 2:
         append = "_3utr"
    elif args.version == 3:
         append = "_5utr"


    
    for folder in target_folder:
    
        folder_files = glob.glob(folder + "/*log*")

    
        for log_file_name in folder_files:

            log_file = open(log_file_name).readlines()

            if log_file == []: continue
            id_ = log_file[0][:-1]
            seqes = log_file[1][:-1]
            seq = "".join([s if s != "U" else "T" for s in seqes])

            struc = log_file[2][:-1].split(" ")[0]

            output_file = log_file_name.split(".log")[0] + "_varna.png"
            info_file = open(folder  + "/complete_info" + str(append) + ".txt").readlines()
            coverage = [int(f) for f in info_file[-1].split(" ")]
            enhance_signal = np.quantile(coverage, 0.6)
            dms_color_map = [str(f) if float(f) >= 0.25  else str(0) for f in info_file[1][:-1].split(" ")]

            dms_color_pred = [float(f) for f in info_file[-3].split(" ")]
            pred_cutoff = np.quantile(dms_color_pred, 0.95)

            
            if len(log_file_name.split("/")[-1].split("_")) > 3:
                if log_file_name.split("/")[-1].split("_")[2] == "dms":
                    dms_color_map = [str(1) if dms_color_pred[enum] >= pred_cutoff and coverage[enum] <= enhance_signal else str(f) for enum, f in enumerate(dms_color_map)]
                    
            elif "ml" in log_file_name.split("/")[-1]:
                if log_file_name.split("/")[-1].split("_")[1] != "dms": 
                    dms_color_map = [str(1) if dms_color_pred[enum] >= pred_cutoff and coverage[enum] <= enhance_signal else str(0) for enum, f in enumerate(dms_color_map)]
                else: 
                    dms_color_map = [str(1) if dms_color_pred[enum] >= pred_cutoff and coverage[enum] <= enhance_signal else str(f) for enum, f in enumerate(dms_color_map)]
            dms_color_map =  ",".join(dms_color_map)
                

            
            highlight_region = ""           
            signal_list = []     
            sequence_list = []     
            sequence_list_oops = []     
            

                    
            for enum, signal in enumerate(info_file[-2][:-3].split(" ")):
                if signal == "0" or signal == "X" or signal == "U":
                    sequence_list_oops.append(str(enum+1))
                    

            for enum, signal in enumerate(info_file[2][:-3].split(" ")):
                if signal == "0" or signal == "X" or signal == "U":
                    signal_list.append(enum)
                elif signal == "1":
                    if signal_list != []:
                        sequence_list.append(signal_list)
                        signal_list = []
                        
            if signal_list != []: sequence_list.append(signal_list)
            
            
            dms_high_light = []
            for enum, seq in enumerate(dms_color_map.split(",")):
                if float(seq) >0: dms_high_light.append(str(enum+1))


            if dms_high_light != []: dms_high_light_str = ",".join(dms_high_light)
            else: dms_high_light_str = ""
            
            
            sequence_list_oops = [s for s in sequence_list_oops if s not in dms_high_light]

            oops_high_lights_str = ",".join(sequence_list_oops)
            

            high_lights = []

            for seq in sequence_list:
                start = seq[0]
                end = seq[-1]
                high_lights.extend([str(start+1) + "-" + str(end+1)])

            if high_lights != []: high_lights_str = ";".join(high_lights)
            else : high_lights_str = high_lights


            if args.version == 3: eclip_info_file = open(folder  + "/" + "eclip_regions_5utr.txt").readlines()
            else: eclip_info_file = open(folder  + "/" + "eclip_regions_3utr.txt").readlines()

            varna_files_created.append(output_file)
            

            annotation_command = ""
            anchor_list = []
            command_lits = []
            last_anchor = - np.inf
            name = None
            
            for line in eclip_info_file:
                new_name = line.split("\t")[3]
                anchor = int(line.split("\t")[5])


                if anchor > last_anchor + 5:
                    if name != None: annotation_command += f"{name}:type=B,anchor={last_anchor},angle=0,size=12,color=#000000;"
                    name = new_name
                else:
                    name = name + "/" + new_name


                last_anchor = anchor
            
            annotation_command += f"{name}:type=B,anchor={last_anchor},angle=0,size=12,color=#000000;"

            run_varna(seqes, struc, output_file, high_lights_str, dms_high_light_str, oops_high_lights_str, annotation_command)

file_ = open(args.output_summary, "w")
file_.write("\n".join(varna_files_created))
file_.close()
