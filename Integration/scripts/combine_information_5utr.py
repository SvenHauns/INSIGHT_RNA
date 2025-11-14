import subprocess as sp
import glob
import os
import pandas as pd
import argparse
import numpy as np
import sys
from Network import *
import yaml
import pickle

class Config:
    def __init__(self, **entries):
        self.__dict__.update(entries)
        self.entries=entries

    def print(self):
        print(self.entries)

def load_config_from_yaml(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return Config(**config)

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

    
    lines[replace_int:replace_int+1] = values
    
    
    file_ = open(file_name)
    lines_read = file_.readlines()
    at_line = None
    draw_line = None
    
    for enum, line in enumerate(lines_read):
    
        if str(line) == "% switch off outline pairs or bases by removing these lines\n":
            at_line = enum
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
	if shape_file== None: prodigal_cmd += " -C" 
	if shape_file!= None: prodigal_cmd += " --shape {shape_file}" 
        
	    
	prodigal_cmd = prodigal_cmd.format(prodigal=prodigal_cmd, input_fasta=fasta_file, shape_file = shape_file)
	
	
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


def create_shape_file_constraint2(id_, seq, dms, dms_predicted, coverage, oops_signal, shape_file1, shape_file2, shape_file3, shape_file4, shape_file5, input_fasta, cutoff = 0.06):

    input_fasta = open(input_fasta, "w")
    input_fasta.write(id_)
    input_fasta.write(seq+ "\n")
    input_fasta.close()
    
    shape_file1 = open(shape_file1, "w")
    shape_file1.write(id_)
    shape_file1.write(seq+ "\n")
    constraint = "".join(["x" if dms[enum] >= cutoff and coverage[enum] > 20 else "." for enum, _ in enumerate(dms)])
    shape_file1.write(constraint)
    shape_file1.close()
    
                
    shape_file2 = open(shape_file2, "w")
    shape_file2.write(id_)
    shape_file2.write(seq+ "\n")
    
    
    constraint = ""
    for enum, _ in enumerate(dms):
        if coverage[enum] > 20: 
            if oops_signal[enum] == 1:
                if dms[enum] >= cutoff: constraint = constraint + "x"
                else: constraint = constraint + "."
            else:
                constraint = constraint + "x"
        else:
            constraint = constraint + "."

    shape_file2.write(constraint)
    shape_file2.close()
                    
                
    shape_file3 = open(shape_file3, "w")
    shape_file3.write(id_)
    shape_file3.write(seq+ "\n")
    
    constraint = ""
    for enum, _ in enumerate(dms):
        if coverage[enum] > 20: 
            if oops_signal[enum] == 1:
                if dms[enum] >= cutoff: constraint = constraint + "x"
                else: constraint = constraint + "."
            else:
                constraint = constraint + "|"
        else:
            constraint = constraint + "."
            

    shape_file3.write(constraint)
    shape_file3.close()
    
    shape_file4 = open(shape_file4, "w")
    shape_file4.write(id_)
    shape_file4.write(seq+ "\n")

    
    constraint = ""
    pred_cutoff = np.quantile(dms_predicted, 0.9)
    for enum, _ in enumerate(dms):
        if coverage[enum] > 20: 
            if oops_signal[enum] == 1:
                if dms[enum] >= cutoff and dms_predicted[enum] >= pred_cutoff: constraint = constraint + "x"
                else: constraint = constraint + "."
            else:
                constraint = constraint + "x"
        else:
            constraint = constraint + "."

    shape_file4.write(constraint)
    shape_file4.close()
    
    
    
    shape_file5 = open(shape_file5, "w")
    shape_file5.write(id_)
    shape_file5.write(seq+ "\n")

    
    constraint = ""

    pred_cutoff = np.quantile(dms_predicted, 0.9)
    for enum, _ in enumerate(dms):
        if coverage[enum] > 20: 
            if oops_signal[enum] == 1:
                if dms[enum] >= cutoff and dms_predicted[enum] >= pred_cutoff: constraint = constraint + "x"
                else: constraint = constraint + "."
            else:
                constraint = constraint + "|"
        else:
            constraint = constraint + "."

    shape_file5.write(constraint)
    shape_file5.close()

    
    return



def create_shape_file_constraint(id_, seq, dms, dms_predicted, coverage, oops_signal, shape_file1, shape_file2, shape_file3, shape_file4, shape_file5,shape_file6, shape_file7, shape_file8, input_fasta, enhance_signal, cutoff = 0.25):

    input_fasta = open(input_fasta, "w")
    input_fasta.write(id_)
    input_fasta.write(seq+ "\n")
    input_fasta.close()
    
    shape_file1 = open(shape_file1, "w")
    shape_file1.write(id_)
    shape_file1.write(seq+ "\n")
    constraint = "".join(["x" if dms[enum] >= cutoff else "." for enum, _ in enumerate(dms)])
    shape_file1.write(constraint)
    shape_file1.close()
    
    # oops-single oops-double oops-single-dms oops_double_dms
    """
    0 in oops_signal: single stranded
    X in oops_signal: double stranded
    U in oops_signal: either way
    """
    
    shape_file2 = open(shape_file2, "w")
    shape_file2.write(id_)
    shape_file2.write(seq+ "\n")
    
    constraint = ""
    for enum, _ in enumerate(dms):
        #if coverage[enum] > 20: 
            if oops_signal[enum] == "1":
                if dms[enum] >= cutoff: constraint = constraint + "x"
                else: constraint = constraint + "."
            elif oops_signal[enum] == "0":
                constraint = constraint + "x"
            elif oops_signal[enum] == "X":
                constraint = constraint + "|"
            elif oops_signal[enum] == "U":
                constraint = constraint + "x"
        #else:
        #    constraint = constraint + "."

    shape_file2.write(constraint)
    shape_file2.close()
                    
                
    shape_file3 = open(shape_file3, "w")
    shape_file3.write(id_)
    shape_file3.write(seq+ "\n")
    
    constraint = ""
    for enum, _ in enumerate(dms):
        #if coverage[enum] > 20: 
            if oops_signal[enum] == "1":
                if dms[enum] >= cutoff: constraint = constraint + "x"
                else: constraint = constraint + "."
            elif oops_signal[enum] == "0":
                constraint = constraint + "x"
            elif oops_signal[enum] == "X":
                constraint = constraint + "|"
            elif oops_signal[enum] == "U":
                constraint = constraint + "|"
        #else:
        #    constraint = constraint + "."
            

    shape_file3.write(constraint)
    shape_file3.close()
    
    shape_file4 = open(shape_file4, "w")
    shape_file4.write(id_)
    shape_file4.write(seq+ "\n")

    
    constraint = ""
    pred_cutoff = np.quantile(dms_predicted, 0.95)
    for enum, _ in enumerate(dms):
        #if coverage[enum] > 20: 
            if oops_signal[enum] == "1":
                if coverage[enum] > enhance_signal:
                    #if dms[enum] >= cutoff or dms_predicted[enum] >= pred_cutoff: constraint = constraint + "x"
                    if dms[enum] >= cutoff: constraint = constraint + "x"
                    else: constraint = constraint + "."

                else:
                    if dms[enum] >= cutoff or dms_predicted[enum] >= pred_cutoff: constraint = constraint + "x"
                    else: constraint = constraint + "."
            elif oops_signal[enum] == "0":
                constraint = constraint + "x"
            elif oops_signal[enum] == "X":
                constraint = constraint + "|"
            elif oops_signal[enum] == "U":
                constraint = constraint + "x"
        #else:
        #    constraint = constraint + "."

    shape_file4.write(constraint)
    shape_file4.close()
    
    
    
    shape_file5 = open(shape_file5, "w")
    shape_file5.write(id_)
    shape_file5.write(seq+ "\n")

    
    constraint = ""

    pred_cutoff = np.quantile(dms_predicted, 0.95)
    for enum, _ in enumerate(dms):
        #if coverage[enum] > 20: 
            if oops_signal[enum] == "1":
                if coverage[enum] > enhance_signal:
                     #if dms[enum] >= cutoff or dms_predicted[enum] >= pred_cutoff: constraint = constraint + "x"
                     if dms[enum] >= cutoff: constraint = constraint + "x"
                     else: constraint = constraint + "."
                else:
                     if dms[enum] >= cutoff or dms_predicted[enum] >= pred_cutoff: constraint = constraint + "x"
                     else: constraint = constraint + "."
            elif oops_signal[enum] == "0":
                constraint = constraint + "x"
            elif oops_signal[enum] == "X":
                constraint = constraint + "|"
            elif oops_signal[enum] == "U":
                constraint = constraint + "|"
        #else:
        #   constraint = constraint + "."

    shape_file5.write(constraint)
    shape_file5.close()
    
    """
    eclip only file
    """
    shape_file6 = open(shape_file6, "w")
    shape_file6.write(id_)
    shape_file6.write(seq+ "\n")
    pred_cutoff = np.quantile(dms_predicted, 0.95)
    constraint = ""
    for enum, _ in enumerate(dms):
            if oops_signal[enum] == "0":
                constraint = constraint + "x"
            elif oops_signal[enum] == "X":
                constraint = constraint + "x"
            elif oops_signal[enum] == "U":
                constraint = constraint + "x"
            else:
                constraint = constraint + "."

    shape_file6.write(constraint)
    shape_file6.close()
    
    
    """
    ML only file
    """
    pred_cutoff = np.quantile(dms_predicted, 0.95)
    shape_file7 = open(shape_file7, "w")
    shape_file7.write(id_)
    shape_file7.write(seq+ "\n")
    constraint = "".join(["x" if dms_predicted[enum] >= pred_cutoff and coverage[enum] < enhance_signal else "." for enum, _ in enumerate(dms_predicted)])
    shape_file7.write(constraint)
    shape_file7.close()
    
    """
    ML+DMS only file
    """

    pred_cutoff = np.quantile(dms_predicted, 0.95)
    shape_file8 = open(shape_file8, "w")
    shape_file8.write(id_)
    shape_file8.write(seq+ "\n")
    
    constraint = ""
    for enum, _ in enumerate(dms):
        if coverage[enum] > enhance_signal:
            #if dms[enum] >= cutoff or dms_predicted[enum] >= pred_cutoff: constraint = constraint + "x"
            if dms[enum] >= cutoff: constraint = constraint + "x"
            else: constraint = constraint + "."
        else:
            #if dms_predicted[enum] >= pred_cutoff: constraint = constraint + "x"
            if dms[enum] >= cutoff or dms_predicted[enum] >= pred_cutoff: constraint = constraint + "x"
            else: constraint = constraint + "."
    #constraint = "".join(["x" if dms[enum] >= cutoff or dms_predicted[enum] >= pred_cutoff else "." for enum, _ in enumerate(dms)])
    shape_file8.write(constraint)
    shape_file8.close()
    
    
    return



def create_shape_file(id_, seq, dms, dms_predicted, coverage, oops_signal, shape_file1, shape_file2, shape_file3, shape_file4, shape_file5, input_fasta):

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

    shape_file4 = open(shape_file4, "w")
    for enum, nuc in enumerate(seq):
        shape_file4.write(str(enum + 1))
        shape_file4.write(" ")
        shape_file4.write(nuc)
        shape_file4.write(" ")
        if coverage[enum] > 20: 
            if oops_signal[enum] == 1:
                shape_file4.write(str(min(dms[enum], dms_predicted[enum])))
            else:
                shape_file4.write(str("0.0"))

                
        shape_file4.write("\n")
        
    shape_file4.close()


    shape_file5 = open(shape_file5, "w")
    for enum, nuc in enumerate(seq):
        shape_file5.write(str(enum + 1))
        shape_file5.write(" ")
        shape_file5.write(nuc)
        shape_file5.write(" ")
        if coverage[enum] > 20: 
            if oops_signal[enum] == 1:
                shape_file5.write(str(min(dms[enum], dms_predicted[enum])))
            else:
                shape_file5.write(str("1.0"))

                
        shape_file5.write("\n")
        
    shape_file5.close()
    
    
    return


def etract_refseq_utr(gff_path, run_type = "full"):


    hg_file = open(gff_path)
    exon_dict = {}
    start_dict = {}
    stop_dict = {}
    transcript_lists = []

    strand_dict = {}
    gene_name_dict = {}
    counter = 0
    chr_dict = {}
    haCDSs_dict = {}
    last_tid = ""
    for enum, line in enumerate(hg_file):
        line_data = line.split()
        
        if enum > 4:
            if line[0] == "#":continue
            if "transcript_id" not in line: print(line)
            transcript_id = line.split("transcript_id")[1].split(";")[0][2:-1]

            
            if transcript_id != last_tid:
                transcript_lists.append(transcript_id)
                exon_dict[transcript_id] = []
                start_dict[transcript_id] = []
                stop_dict[transcript_id] = []

                gene_name_dict[transcript_id] = ""
                counter = counter + 1
                last_tid = transcript_id
                
            gen_type_ = line.split("\t")[2]

            
            if gen_type_ == "exon":
                gene_name_dict[transcript_id] = line.split("gene_id")[1].split(";")[0][2:-1]
                start = line.split("\t")[3]
                end = line.split("\t")[4]
                exon_dict[transcript_id].append([int(start), int(end)])
                strand_dict[transcript_id] = line.split("\t")[6]
                chr_dict[transcript_id] = line.split("\t")[0]
            if  gen_type_ == "start_codon":
                start = line.split("\t")[3]
                end = line.split("\t")[4]
                start_dict[transcript_id]=[int(start), int(end)]
            elif  gen_type_ == "stop_codon":
                start = line.split("\t")[3]
                end = line.split("\t")[4]
                stop_dict[transcript_id]=[int(start), int(end)]

            if gen_type_ == "CDS":
                haCDSs_dict[transcript_id] = True
            
    

    
    
    
            
    hg_file.close()

    #utr_regions = extract_utr_region(transcript_lists, exon_dict, start_dict, strand_dict)
    utr5_regions = extract_utr_region(transcript_lists, exon_dict, start_dict, strand_dict)
    utr3_regions = extract_3utr_region(transcript_lists, exon_dict, stop_dict, strand_dict)

    for key in utr5_regions.keys():
        size = 0
        for val in utr5_regions[key]:
            size = size + np.abs(val[0] - val[1])
        utr5_regions[key] = size


    for key in utr3_regions.keys():
        size = 0
        for val in utr3_regions[key]:
            size = size + np.abs(val[0] - val[1])
        utr3_regions[key] = size

            
    return utr5_regions, utr3_regions, gene_name_dict

def extract_full_exon(transcript_lists, exon_dict, start_dict, strand_dict):
    utrs = {}
    for key in transcript_lists:
    
        exon_dict_list = []
        
        if exon_dict[key] == []:
            continue
        if start_dict[key] == []:
            continue
            
        for exon in exon_dict[key]:
                exon_dict_list.append(exon)
        utrs[key] = exon_dict_list
    
    return utrs


    
def extract_3utr_region(transcript_lists, exon_dict, start_dict, strand_dict):

    utrs = {}
  
    for key in transcript_lists:
    
        if exon_dict[key] == []:
            continue
        if start_dict[key] == []:
            continue
    
        strand = strand_dict[key]
        exon_dict_list = []

        if strand == "+":

            for exon in exon_dict[key]:
                if exon[1] > start_dict[key][1]:
            
                    exon_dict_list.append(exon)
                
                    if exon[0] <= start_dict[key][0]:

                        exon[0] = start_dict[key][1] + 1 
                        
        elif strand == "-":
            for exon in exon_dict[key]:
            
                if exon[0] < start_dict[key][0]:
                
                    if exon[1] >= start_dict[key][1]:
                        exon[1] = start_dict[key][0] - 1 
                   
                    exon_dict_list.append(exon)
                
        utrs[key] = exon_dict_list
        exon_dict_list = []

    return utrs
    
    
    
def extract_utr_region(transcript_lists, exon_dict, start_dict, strand_dict):

    utrs = {}
  
    for key in transcript_lists:
    
    
        if exon_dict[key] == []:
            continue
        if start_dict[key] == []:
            continue

    
        strand = strand_dict[key]
        exon_dict_list = []
        if strand == "+":
        
            offset = 0
            extra_nuc = 0

            for exon in exon_dict[key]:
                
                if exon[0] < start_dict[key][0] +offset:
            
              
                    if exon[1] >= start_dict[key][0] + offset:

                        exon[1] = start_dict[key][0] - 1 + offset

                    exon_dict_list.append(exon)
        
        
        elif strand == "-":
        
            offset = 0
            extra_nuc = 0
            
            for exon in exon_dict[key]:
            
                if exon[1] > start_dict[key][1] - offset:
                
                    if exon[0] <= start_dict[key][1] -offset:
                        exon[0] = start_dict[key][1] + 2 - offset
                   
                    if exon[0] <= exon[1]:exon_dict_list.append(exon)
                    
                    
        utrs[key] = exon_dict_list
        exon_dict_list = []
    

    return utrs


def load_dataset(datafile, length_limit = 850, length_min = 100):

    dataset = []
    file_ = open(datafile).readlines()

    for enum, line_ in enumerate(file_):

        if line_[0] == ">":

            if len(file_[enum+1][:-1]) > length_limit: continue
            if len(file_[enum+1][:-1]) < length_min: continue
            dataset.append([line_[:-1],file_[enum+1][:-1], file_[enum+2][:-1]])

    return dataset


def load_eClip_data(paths_to_data):

    eClip_dict = {"chr1":[], "chr2":[], "chr3":[], "chr4":[], "chr5":[], "chr6":[], "chr7":[], "chr8":[], "chr9":[], "chr10":[], 
    "chr11":[], "chr12":[], "chr13":[], "chr14":[], "chr15":[], "chr16":[], "chr17":[], "chr18":[], "chr19":[], "chr20":[], "chr21":[], 
    "chr22":[], "chrX":[], "chrY":[], "chrM": []}


    for path_to_data in glob.glob(paths_to_data + "/*"):
        folders = glob.glob(path_to_data + "/*")
    
        for folder in folders:
            bedfiles = glob.glob(folder + "/*.bed")
            for bedfile in bedfiles:
                bedfile = open(glob.glob(bedfile+"/*.bed")[0]).readlines()
            
                for line in bedfile:
                    chr_ = line.split("\t")[0]
                    start = line.split("\t")[1]
                    end = line.split("\t")[2]
                    target = line.split("\t")[3]
                    if chr_ not in eClip_dict.keys():continue
                
                    eClip_dict[chr_].append([int(start), int(end), target, folder.split("/")[-1]])
                
    return eClip_dict
    
def load_eClip_data2(paths_to_data):

    eClip_dict = {"chr1":[], "chr2":[], "chr3":[], "chr4":[], "chr5":[], "chr6":[], "chr7":[], "chr8":[], "chr9":[], "chr10":[], 
    "chr11":[], "chr12":[], "chr13":[], "chr14":[], "chr15":[], "chr16":[], "chr17":[], "chr18":[], "chr19":[], "chr20":[], "chr21":[], 
    "chr22":[], "chrX":[], "chrY":[], "chrM": []}


    for path_to_data in glob.glob(paths_to_data + "/*"):
        folders = glob.glob(path_to_data + "/*.bed")
    
        for folder in folders:
            #bedfiles = glob.glob(folder + "/*.bed")
            #for bedfile in bedfiles:
            print(folder)
            print(glob.glob(folder+"/*.bed"))
            bedfile = glob.glob(folder+"/*.bed")
            
            if bedfile == []: continue
            bedfile = open(bedfile[0]).readlines()
            
            
            for line in bedfile:
                chr_ = line.split("\t")[0]
                start = line.split("\t")[1]
                end = line.split("\t")[2]
                target = line.split("\t")[3]
                if chr_ not in eClip_dict.keys():continue
                
                eClip_dict[chr_].append([int(start), int(end), target, folder.split("/")[-1]])
                
    return eClip_dict
    
    
def load_motif_data():


    return
    
class RNA_Dataset2(torch.utils.data.Dataset):
    def __init__(self,data):
        self.data=data
        self.tokens={nt:i for i,nt in enumerate('ACGT')}
        self.counter = 0

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sequence=[self.tokens[nt] for nt in self.data[idx][1]]
        sequence=np.array(sequence)
        sequence=torch.tensor(sequence)
        target = self.data[idx][2].split(" ")
        target = [float(t) for t in target  if t != ""]
        target=torch.tensor(target)

        mask = [1 if i=="A" or i =="C" else 0 for i in self.data[idx][1]]
        mask=torch.tensor(mask)

        if torch.isnan(target).sum() >0:
            target = torch.zeros(target.size())
            self.counter = self.counter + 1

        assert torch.isnan(target).sum() == 0, "NaN values found in target data!"

        return {'sequence':sequence, 'id':self.data[idx][0], 'target': target, 'mask':mask}


class RNA_Dataset(torch.utils.data.Dataset):
    def __init__(self,data):
        self.data=data
        self.tokens={nt:i for i,nt in enumerate('ACGU')}

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sequence=[self.tokens[nt] for nt in self.data.loc[idx,'sequence']]
        sequence=np.array(sequence)
        sequence=torch.tensor(sequence)

        return {'sequence':sequence}
    
def get_model_prediction(model, sequence, step_size = 200):
    tokens={nt:i for i,nt in enumerate('ACGT')}
    sequence=[tokens[nt] for nt in sequence]


    model.eval()
    prediction = []
    with torch.no_grad():
        for i in range(0, len(sequence), step_size):
            subseq=torch.tensor(sequence[i:i+step_size]).cpu().unsqueeze(0)
            res = model(subseq,torch.ones_like(subseq))
            prediction.extend(res[..., 1:2].squeeze(dim=0).squeeze(dim=-1))


    return prediction
    
def find_region(eClip_list, start, end):

    detected_region = []
    for eclip in eClip_list:
        if eclip[0] >= start:
            if eclip[1] <= end:
                detected_region.append(eclip)

    return detected_region


if __name__ == "__main__":


    cmdline_parser = argparse.ArgumentParser('retrieve RNA seq structure')

    cmdline_parser.add_argument('-f', '--dms_analysis_file',
                                default="",
                                help='input file',
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
    cmdline_parser.add_argument('-b', '--gff_path',
                                default="",
                                help='gff_path',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-u', '--target_folder_5utr',
                                default="",
                                help='target folder 5utr',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-i', '--target_folder_3utr',
                                default="",
                                help='target folder 3utr',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-j', '--model_path',
                                default="",
                                help='model_path',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-c', '--eclip_path',
                                default="",
                                help='path to eclip data',
                                required = True,
                                type=str)


    args, unknowns = cmdline_parser.parse_known_args()

    
    eClip_dict = load_eClip_data(args.eclip_path)
    
    
    
    strandedness = {"AKAP1": "single", 
    "AQR": "unknown",
    "FTO": "single",
    "FXR2": "unknown",
    "G3BP1": "single",
    "KHSRP": "single",
     "NIP7": "unknown",
     "PABPN1": "single",
     "PRPF4": "double",
     "SAFB": "unknown",
     "SDAD1": "unknown",
     "SSB": "single",
     "TIAL1": "unknown",
     "UPF1": "single",
     "YBX3": "single",
     "ZNF800": "unknown",
     "ABCF1": "unknown",
     "CPEB4": "unknown",
     "MATR3": "single",
     "PCBP1": "single"}


    utr5, utr3, gene_names = etract_refseq_utr(args.gff_path, run_type = "full")

    
    dms_file = open(args.dms_analysis_file).readlines()
    dms_targets_inds = [[d[1:-1], enum] for enum, d in enumerate(dms_file) if d[0] == ">"]
    dms_targets = [d[0] for d in dms_targets_inds]

    oops_seq_files = glob.glob(args.oops_seq_folder + "*")

    
    oops_seq_folder_ids = [o.split("|")[0].split("/")[-1] for o in oops_seq_files]
    dict_, exon_dict = get_exons(args.exon_file)
    

    created = []

    model=RibonanzaNet(load_config_from_yaml("./Integration/scripts/pairwise.yaml")).cpu()
    model.load_state_dict(torch.load(args.model_path,map_location='cpu'))

    cont_ = True
    
    test_count = 0
    
    for target_z in dms_targets_inds:
        target = target_z[0]

        gene_name = gene_names[target]
        #if "/mnt/StrucIntegration/Integration/targets/output_5utr/" + target in glob.glob("mnt/StrucIntegration/Integration/targets/output_5utr/*"):continue
        #if target.split("_")[0] == "XM" or target.split("_")[0] == "XP" or target.split("_")[0] == "XN": continue


        #### load normalization files
        
        a_baseline_file = "Integration/a_baseline"
        c_baseline_file = "Integration/c_baseline"
        with open(a_baseline_file,"rb") as f: 
            baseline = pickle.load(f)
        
        with open(c_baseline_file,"rb") as f: 
            c_baseline = pickle.load(f)
        
        baseline.update(c_baseline)

        ind = dms_targets.index(target)
        id_ = dms_file[dms_targets_inds[ind][1]]
        seq = dms_file[dms_targets_inds[ind][1]+1][:-1]
        dms = np.array([float(c) for c in dms_file[dms_targets_inds[ind][1]+2][:-1].split(" ") if c != ""])
        cov = [int(c) for c in dms_file[dms_targets_inds[ind][1]+4][:-1].split(" ") if c != ""]
        dms_counts = [int(c) for c in dms_file[dms_targets_inds[ind][1]+3][:-1].split(" ") if c != ""]
        dms_raw_signal = [int(c)/cov[enum] if cov[enum] >20 else 0 for enum, c in enumerate(dms_counts)]
        
        dms = (dms - dms.min()) / (dms.max() - dms.min())

        inds = [enum for enum, c in enumerate(dms_raw_signal) if c > 0.06 ]
        

        dms_raw_signal = [c/baseline[seq[max(0,enum-1):min(enum+2, len(seq))]] if seq[enum] == "A" or seq[enum] == "C" else 0 for enum, c in enumerate(dms_raw_signal)]


        if len(seq) > 1500: continue

            
            
        oops_search = target

        utr5_len = utr5[target]
        #utr3_len = utr3[target]
        


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
        

        #selected_region = [int(r.split("-")[0]) for enum, r in enumerate(regions) if mean[enum] >= 50 and padj_val[enum] < 0.05 and lf[enum] > 3] # select-non significant regions for analysis
        selected_region = [int(r.split("-")[0]) for enum, r in enumerate(regions) if mean[enum] >= 10 and padj_val[enum] < 0.1] # select-non significant regions for analysis
        #selected_region = [r for enum, r in enumerate(regions) if mean[enum] >= 50 and padj_val[enum] > 0.05] # select-non significant regions for analysis
        
        if len(selected_region) == 0: continue
        exons = sorted(exon_dict[id_[1:-1]], key=lambda l:int(l.split("-")[0]))
        lengths = [int(e.split("-")[1]) - int(e.split("-")[0]) for e in exons]
        chr_ = dict_[id_[1:-1]]
        chr_eClip = eClip_dict[chr_]

        
        found_regions = []
        signal_fade_out = ["1" for _ in range(len(seq))]
        signal_fade_out_relevant_oops = ["1" for _ in range(len(seq))]
        selected_eclips = []
        oops_seq_region = []
        
        for region_selected in selected_region:
        

            regions = find_region(chr_eClip, region_selected, region_selected+40)
            
            if len(regions) == 0: continue

            for region in regions: 


                ss_strand = "single"
                

                

                for exon_num, exon_range in enumerate(exons):
            
                    if (int(region[0]) >= int(exon_range.split("-")[0])) and (int(region[0]) < int(exon_range.split("-")[1])):

                        offset = int(region[0]) - int(exon_range.split("-")[0])
                        offset2 = int(region_selected) - int(exon_range.split("-")[0])

                   
                        if exon_num>0:

                   
                            offset = offset + sum(lengths[:exon_num])
                            offset2 = offset2 + sum(lengths[:exon_num])
                            #offset = offset + lengths[exon_num-1]

                        dms_raw_signal_in_region = dms_raw_signal[offset:offset+(region[1]-region[0])]
                        
                        sig_dms_signals = sum([1 for d in dms_raw_signal_in_region if d > 0.25])
                        
                        if sig_dms_signals != 0 or ss_strand == "single":
                            signal_fade_out[offset:offset+(region[1]-region[0])] = ["0" for _ in range(0,len(signal_fade_out[offset:offset+(region[1]-region[0])]))]
                        elif ss_strand == "double":
                            signal_fade_out[offset:offset+(region[1]-region[0])] = ["X" for _ in range(0,len(signal_fade_out[offset:offset+(region[1]-region[0])]))]
                        elif ss_strand == "unknown":
                            signal_fade_out[offset:offset+(region[1]-region[0])] = ["U" for _ in range(0,len(signal_fade_out[offset:offset+(region[1]-region[0])]))]
                        
                        signal_fade_out_relevant_oops[offset2:offset2+40] = ["0" for _ in range(0,len(signal_fade_out[offset2:offset2+40]))]
                        
                        
                        start = offset
                        end = offset+(region[1]-region[0])
                        to_append = region.copy()
                        to_append.extend([ss_strand, start, end])
                        found_regions.append(to_append)
                        oops_seq_region.append([ss_strand, offset2, offset2+40])

                        
        if found_regions == []: continue
        dms_prediction = get_model_prediction(model, seq)#
        


        """           
            
            
        
            for exon_num, exon_range in enumerate(exons):
            
            
            
               if (int(region_selected) >= int(exon_range.split("-")[0])) and (int(region_selected) < int(exon_range.split("-")[1])):
                   print("found")
                   print(region_selected)
                   print(exon_range)
                   offset = int(region_selected) - int(exon_range.split("-")[0])
                   print(offset)
                   print("exon_num")
                   print(exon_num)
                   
                   if exon_num>0:
                       print(lengths)
                       print(sum(lengths[:exon_num]))
                   
                       offset = offset + sum(lengths[:exon_num])
                       #offset = offset + lengths[exon_num-1]
                       print(len(signal_fade_out))
                       print("diff")
                       print(lengths[exon_num-1])
                       print(sum(lengths[:exon_num]))
                       print("offset res")
                       
                       print(offset)
                       
                   signal_fade_out[offset:offset+40] = np.zeros(len(signal_fade_out[offset:offset+40]))
        """           



        ################################################################################
        ############################# create 5 UTR data ###########################
        ################################################################################
        
        
        

        #seq_5utr = seq[:utr5_len+1]
        #dms_5utr = dms_raw_signal[:utr5_len+1]
        #cov_5utr = cov[:utr5_len+1]
        #dms_pred_5utr = dms_prediction[:utr5_len+1]
        #signal_fade_out_5utr = signal_fade_out[:utr5_len+1]
        #signal_fade_out_relevant_oops_5utr = signal_fade_out_relevant_oops[:utr5_len+1]
        
        
        """
        fixed 5 utr
        
        """
        
        seq_5utr = seq
        dms_5utr = dms_raw_signal[:len(seq)]
        cov_5utr = cov[:len(seq)]
        dms_pred_5utr = dms_prediction
        signal_fade_out_5utr = signal_fade_out[:len(seq)]
        signal_fade_out_relevant_oops_5utr = signal_fade_out_relevant_oops[:len(seq)]
        
        enhance_signal = np.quantile(cov_5utr, 0.6)

        found_signal = False
        if "0" in signal_fade_out_5utr:
            found_signal = True
        elif "X" in signal_fade_out_5utr:
            found_signal = True
        elif "U" in signal_fade_out_5utr:
            found_signal = True
            
        if found_signal == False: continue
        
        if os.path.isdir(args.target_folder_5utr) == False: os.mkdir(args.target_folder_5utr)
        if os.path.isdir(args.target_folder_5utr + id_[1:-1] + "/") == False: os.mkdir(args.target_folder_5utr + id_[1:-1] + "/")
        
        eclip_info_file = open(args.target_folder_5utr + id_[1:-1] + "/" + "eclip_regions_5utr.txt", "w")
        for eClip_num, e_reg in enumerate(found_regions):
            
            #if e_reg[-2] > utr5_len: continue
            if e_reg[-2] > len(seq): continue
            eclip_info_file.write(str(e_reg[0]))
            eclip_info_file.write(str("\t"))
            eclip_info_file.write(str(e_reg[1]))
            eclip_info_file.write(str("\t"))
            eclip_info_file.write(str(e_reg[2]))
            eclip_info_file.write(str("\t"))
            eclip_info_file.write(str(e_reg[3]))
            eclip_info_file.write(str("\t"))
            eclip_info_file.write(str(e_reg[4]))
            eclip_info_file.write(str("\t"))
            eclip_info_file.write(str(e_reg[5]))
            eclip_info_file.write(str("\t"))
            eclip_info_file.write(str(e_reg[6]))
            eclip_info_file.write(str("\t"))
            eclip_info_file.write(str(oops_seq_region[eClip_num][1]))
            eclip_info_file.write(str("\t"))
            eclip_info_file.write(str(oops_seq_region[eClip_num][2]))
            eclip_info_file.write(str("\n"))
            
            
            
        eclip_info_file.close()


        shape_file1 = args.target_folder_5utr + id_[1:-1] + "1.shape.fa"
        shape_file2 = args.target_folder_5utr + id_[1:-1] + "2.shape.fa"
        shape_file3 = args.target_folder_5utr + id_[1:-1] + "3.shape.fa"
        shape_file4 = args.target_folder_5utr + id_[1:-1] + "4.shape.fa"
        shape_file5 = args.target_folder_5utr + id_[1:-1] + "5.shape.fa"
        shape_file6 = args.target_folder_5utr + id_[1:-1] + "6.shape.fa"
        shape_file7 = args.target_folder_5utr + id_[1:-1] + "7.shape.fa"
        shape_file8 = args.target_folder_5utr + id_[1:-1] + "8.shape.fa"
        
        input_fasta = args.target_folder_5utr + id_[1:-1] + ".fasta"
        info_file = args.target_folder_5utr + id_[1:-1] + "/" + "complete_info_5utr.txt"
        
        create_shape_file_constraint(id_, seq_5utr, dms_5utr, dms_pred_5utr, cov_5utr, signal_fade_out_5utr, shape_file1,
                           shape_file2, shape_file3, shape_file4, shape_file5,shape_file6,shape_file7, shape_file8, input_fasta, enhance_signal)

        
        target_folder = args.target_folder_5utr + id_[1:-1] + "/"
        if os.path.isdir(target_folder) == False: os.mkdir(target_folder)
        info_file = open(info_file, "w")
        info_file.write(gene_name + " " + target + "\n")
        info_file.write(str(" ".join([str(d) for d in dms_5utr])))
        info_file.write("\n")
        info_file.write(str(" ".join([str(s) for s in signal_fade_out_5utr])))
        info_file.write("\n")
        info_file.write(str(" ".join([str(d.item()) for d in dms_pred_5utr])))
        info_file.write("\n")
        info_file.write(str(" ".join([str(d) for d in signal_fade_out_relevant_oops_5utr])))
        info_file.write("\n")
        info_file.write(str(" ".join([str(d) for d in cov_5utr])))
        info_file.close()
        """
        ML only and ML + DMS !
        
        RNAfold("RNAfold", input_fasta, shape_file=None, target_file = target_folder + "plain", id_ = id_[1:-1])
        RNAfold("RNAfold", input_fasta, shape_file1, target_file = target_folder + "dms", id_ = id_[1:-1])
        RNAfold("RNAfold", input_fasta, shape_file2, target_file = target_folder + "oops_double", id_ = id_[1:-1])
        RNAfold("RNAfold", input_fasta, shape_file3, target_file = target_folder + "oops_single", id_ = id_[1:-1])
        RNAfold("RNAfold", input_fasta, shape_file4, target_file = target_folder + "oops_double_dms", id_ = id_[1:-1])
        RNAfold("RNAfold", input_fasta, shape_file5, target_file = target_folder + "oops_single_dms", id_ = id_[1:-1])
        """
        ###### rnafold eclip only
        RNAfold("RNAfold", shape_file7, shape_file=None, target_file = target_folder + "ml", id_ = id_[1:-1])
        RNAfold("RNAfold", shape_file8, shape_file=None, target_file = target_folder + "ml_dms", id_ = id_[1:-1])
        
        
        RNAfold("RNAfold", shape_file6, shape_file=None, target_file = target_folder + "eclip", id_ = id_[1:-1])
        RNAfold("RNAfold", input_fasta, shape_file=None, target_file = target_folder + "plain", id_ = id_[1:-1])
        RNAfold("RNAfold", shape_file1, shape_file=None, target_file = target_folder + "dms", id_ = id_[1:-1])
        RNAfold("RNAfold", shape_file2, shape_file=None, target_file = target_folder + "oops_single", id_ = id_[1:-1])
        RNAfold("RNAfold", shape_file3, shape_file=None, target_file = target_folder + "oops_double", id_ = id_[1:-1])
        RNAfold("RNAfold", shape_file4, shape_file=None, target_file = target_folder + "oops_single_dms", id_ = id_[1:-1])
        RNAfold("RNAfold", shape_file5, shape_file=None, target_file = target_folder + "oops_double_dms", id_ = id_[1:-1])
        
        #raise NotImplementedError
        os.remove(shape_file1)
        os.remove(shape_file2)
        os.remove(shape_file3)
        os.remove(shape_file4)
        os.remove(shape_file5)
        os.remove(shape_file6)
        os.remove(shape_file7)
        os.remove(shape_file8)
        os.remove(input_fasta)
        
        test_count = test_count + 1



    summary = open(args.summary, "w")
    summary.write("".join(id_))
    summary.close()
    
