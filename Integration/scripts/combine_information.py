import subprocess as sp
import glob
import os
import pandas as pd
import argparse
import numpy as np


def get_DMS_values(file_name, clamp_value = 0.2):
    
    file_ = open(file_name)
    file_lines = file_.readlines()
    
    value_list = file_lines[1].split(" ")[:-1]
    value_list = [str(min(float(v),clamp_value)) + "\n" for v in value_list]

    return value_list


def color_ps_file(file_name, dms_vaues_path, function_file):

    file2_ = open(function_file)
    lines = file2_.readlines()
    file2_.close()
    replace_int = None
    

    for enum,line in enumerate(lines):

        if str(line) == "REPLACE_ME_WITH_VALUES\n":
            replace_int = enum
            
            
    values = get_DMS_values(dms_vaues_path)
    print(values)
    
    lines[replace_int:replace_int+1] = values
    
    
    file_ = open(file_name)
    lines_read = file_.readlines()
    at_line = None
    draw_line = None
    
    for enum, line in enumerate(lines_read):
    
        if str(line) == "% switch off outline pairs or bases by removing these lines\n":
            at_line = enum
            print(at_line)
        if str(line) == "drawoutline\n":
            draw_line = enum
            
    file_.close()
    draw_list = ["/invert false def\n", "drawreliability\n", "0.1 0.1 colorbar\n"]
    lines_read[at_line:at_line] = lines
    lines_read[draw_line + len(lines):draw_line + len(lines)] = draw_list
    
    file_ = open(file_name, "w")
    for line in lines_read:
        file_.write(line)
        
    file_.close()
    

    return


def RNAfold(prodigal_cmd, fasta_file, shape_file, target_file, id_):
    
	fasta_file_preffix = fasta_file.rsplit('.', 1)[0]
	output_pdf = fasta_file_preffix + '_proteins.fa'
	prodigal_cmd += ' {input_fasta} --filename-full -p'
	if shape_file!= None: prodigal_cmd += " --shape {shape_file}" 
        
	    
	prodigal_cmd = prodigal_cmd.format(prodigal=prodigal_cmd, input_fasta=fasta_file, shape_file = shape_file)
	
	print(prodigal_cmd)
	
	log_file = target_file + '_RNAfold.log'
	
	with open(log_file, 'w') as lf:
		sp.call(prodigal_cmd.split(), stdout=lf)
        
	os.rename(id_ + "_dp.ps", target_file + "_dp.ps")
	os.rename(id_ + "_ss.ps", target_file + "_ss.ps")
    
	return



def get_concat_v(im1, im2, im3):
    dst = Image.new('RGB', (im1.width, im1.height + im2.height + im3.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    dst.paste(im3, (0, im1.height + im2.height))
    return dst



def expand_white_space(im1):
    dst = Image.new('RGB', (im1.width, im1.height + 80), color = "white")
    dst.paste(im1, (0, 80))
    
    return dst
    
    
def get_exons(file_):

    file_ = open(file_).readlines()
    
    dict_ =  {}
    exon_dict = {}
    
    for line_ in file_:

        id_ = line_.split("\t")[3]
        
        if id_ not in dict_.keys():
            exon_dict[id_] = []
            dict_[id_] = line_.split("\t")[0]

        exon_dict[id_].append(line_.split("\t")[1] + "-" + line_.split("\t")[2])


    return dict_, exon_dict
    
    
def create_shape_file(id_, seq, dms, coverage, oops_signal, shape_file1, shape_file2, shape_file3, input_fasta):

    input_fasta = open(input_fasta, "w")
    input_fasta.write(id_)
    input_fasta.write(seq)
    input_fasta.close()
    
    shape_file1 = open(shape_file1, "w")
    for enum, nuc in enumerate(seq):
        shape_file1.write(str(enum + 1))
        shape_file1.write(" ")
        shape_file1.write(nuc)
        shape_file1.write(" ")
        if coverage[enum] > 20: shape_file1.write(str(dms[enum]))
        shape_file1.write("\n")
        
    shape_file1.close()

    shape_file2 = open(shape_file2, "w")
    for enum, nuc in enumerate(seq):
        shape_file2.write(str(enum + 1))
        shape_file2.write(" ")
        shape_file2.write(nuc)
        shape_file2.write(" ")
        if coverage[enum] > 20: 
            if oops_signal[enum] == 1:
                shape_file2.write(str(dms[enum]))
            else:
                shape_file2.write(str("0.0"))

                
        shape_file2.write("\n")
        
    shape_file2.close()
    
    
    shape_file3 = open(shape_file3, "w")
    for enum, nuc in enumerate(seq):
        shape_file3.write(str(enum + 1))
        shape_file3.write(" ")
        shape_file3.write(nuc)
        shape_file3.write(" ")
        if coverage[enum] > 20: 
            if oops_signal[enum] == 1:
                shape_file3.write(str(dms[enum]))
            else:
                shape_file3.write(str("1.0"))

                
        shape_file3.write("\n")
        
    shape_file3.close()
    
    return
    
    
if __name__ == "__main__":


    cmdline_parser = argparse.ArgumentParser('retrieve RNA seq structure')

    cmdline_parser.add_argument('-f', '--dms_analysis_file',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-t', '--target_file',
                                default="",
                                help='target file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-g', '--target_folder',
                                default="",
                                help='target folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-s', '--oops_seq_folder',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-e', '--exon_file',
                                default="",
                                help='exon_file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-a', '--summary',
                                default="",
                                help='output file summary',
                                required = True,
                                type=str)

                                
    args, unknowns = cmdline_parser.parse_known_args()
    print(args.target_file)
    print(open(args.target_file).readlines())
    targets = [t.split("\t")[-1][:-1] for t in open(args.target_file).readlines()]
    
    dms_file = open(args.dms_analysis_file).readlines()
    dms_targets_inds = [[d[1:-1], enum] for enum, d in enumerate(dms_file) if d[0] == ">"]
    dms_targets = [d[0] for d in dms_targets_inds]
    
    oops_seq_files = glob.glob(args.oops_seq_folder + "*")
    oops_seq_folder_ids = [o.split("|")[0].split("/")[-1] for o in oops_seq_files]
    dict_, exon_dict = get_exons(args.exon_file)
    
    print(targets)
    created = []
    
    for target in targets:
        ind = dms_targets.index(target)
        id_ = dms_file[dms_targets_inds[ind][1]]
        seq = dms_file[dms_targets_inds[ind][1]+1][:-1]
        dms = np.array([float(c) for c in dms_file[dms_targets_inds[ind][1]+2][:-1].split(" ") if c != ""])
        cov = [int(c) for c in dms_file[dms_targets_inds[ind][1]+4][:-1].split(" ") if c != ""]
        
        dms = (dms - dms.min()) / (dms.max() - dms.min())

            
            
        oops_search = target
        
        if oops_search not in oops_seq_folder_ids: continue
        oops_data = oops_seq_files[oops_seq_folder_ids.index(oops_search)]
        
        oops_data = pd.read_csv(oops_data)
        regions = oops_data['Unnamed: 0']
        mean = oops_data['baseMean']
        lf = oops_data['log2FoldChange']
        padj_val = oops_data['padj']
        
        ################################################################################
        ############################# create selection mask ###########################
        ################################################################################
        
        selected_region = [r for enum, r in enumerate(regions) if mean[enum] >= 50 and padj_val[enum] > 0.05] # select-non significant regions for analysis
            
        exons = sorted(exon_dict[id_[1:-1]], key=lambda l:int(l.split("-")[0]))
            
        lengths = [int(e.split("-")[1]) - int(e.split("-")[0]) for e in exons]

            
        signal_fade_out = [1 for _ in range(len(seq))]
        for region_selected in selected_region:
            for exon_num, exon_range in enumerate(exons):
               if (int(region_selected) >= int(exon_range.split("-")[0])) and (int(region_selected) < int(exon_range.split("-")[1])) :
                   offset = int(region_selected) - int(exon_range.split("-")[0])
                   if exon_num>0:
                       offset = offset + lengths[exon_num-1]
                   signal_fade_out[offset:offset+40] = np.zeros(len(signal_fade_out[offset:offset+40]))
                   
                   
        shape_file1 = args.target_folder + id_[1:-1] + "1.shape"
        shape_file2 = args.target_folder + id_[1:-1] + "2.shape"
        shape_file3 = args.target_folder + id_[1:-1] + "3.shape"
        input_fasta = args.target_folder + id_[1:-1] + ".fasta"
        info_file = args.target_folder + id_[1:-1] + "/" + "complete_info.txt"
        
        create_shape_file(id_, seq, dms, cov, signal_fade_out, shape_file1, shape_file2, shape_file3, input_fasta)

        
        target_folder = args.target_folder + id_[1:-1] + "/"
        if os.path.isdir(target_folder) == False: os.mkdir(target_folder)
        info_file = open(info_file, "w")
        info_file.write(str(" ".join([str(d) for d in dms])))
        info_file.write("\n")
        info_file.write(str(" ".join([str(int(s)) for s in signal_fade_out])))
        info_file.close()
        
        RNAfold("RNAfold", input_fasta, shape_file=None, target_file = target_folder + "plain", id_ = id_[1:-1])
        RNAfold("RNAfold", input_fasta, shape_file1, target_file = target_folder + "dms", id_ = id_[1:-1])
        RNAfold("RNAfold", input_fasta, shape_file2, target_file = target_folder + "oops_double", id_ = id_[1:-1])
        RNAfold("RNAfold", input_fasta, shape_file3, target_file = target_folder + "oops_single", id_ = id_[1:-1])
        
        os.remove(shape_file1)
        os.remove(shape_file2)
        os.remove(shape_file3)
        os.remove(input_fasta)
        
        created.append(id_)


    summary = open(args.summary, "w")
    summary.write("".join(id_))
    summary.close()
    
