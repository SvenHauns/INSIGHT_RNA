import pysam
import inspect
mismatch_count = 0
import time
import psutil
import multiprocessing
import time
from multiprocessing import Pool
import collections
import pickle
import glob
import random
import numpy as np
import matplotlib.pyplot as plt
import pickle
from Bio import SeqIO
import argparse


def autolabel(rects, ax, sum_):

    for rect in rects:
        height = rect.get_height()

        ax.text(rect.get_x() + rect.get_width()/2., 1.02*height,
            str(height)[0:4] + " total: " + str(height*sum_),
            ha='center', va='bottom', rotation=90, fontsize= 4)

def create_plot(ax, data, title, biotype, exon_number):

    r1=ax.plot([x for x in range(len(data))], data, color ='red'
            )


    ax.set_ylim(0, 10)
    #ytick = ax.get_yticks()
    #ax.set_yticks(ytick, fontsize=4)
    ax.set_xlabel("Base", fontsize = 4)
    ax.set_ylabel("Coverage", fontsize = 4)
    ax.set_title(title + ""+ str(exon_number) + "_" + biotype , fontsize= 4)

    ax.tick_params(axis = "both", which = "major", labelsize = 4)
    ax.tick_params(axis = "both", which = "minor", labelsize = 4)
    

    return



def save_bed(chr_, start, end, filename):

    file_ = open(filename, "a+")
    
    file_.write("chr" + (str(chr_)))
    file_.write("\t")
    file_.write(str(start))
    file_.write("\t")
    file_.write(str(end))
    file_.write("\n")
    
    file_.close()

    return

def metrics(stat_stable_dict, stat_changed_dict):


                        
    change_a = (stat_changed_dict["A->C"] + stat_changed_dict["A->T"] + stat_changed_dict["A->G"]) / max(stat_stable_dict["A"],1)
    change_c = (stat_changed_dict["C->A"] + stat_changed_dict["C->T"] + stat_changed_dict["C->G"]) / max(stat_stable_dict["C"],1)
    change_g = (stat_changed_dict["G->C"] + stat_changed_dict["G->T"] + stat_changed_dict["G->C"]) / max(stat_stable_dict["G"],1)    
    change_t = (stat_changed_dict["T->C"] + stat_changed_dict["T->A"] + stat_changed_dict["T->G"]) / max(stat_stable_dict["T"],1)     
     
    change_dms = (stat_changed_dict["A->C"] + stat_changed_dict["A->T"] + stat_changed_dict["A->G"] + stat_changed_dict["C->A"] + stat_changed_dict["C->T"] + stat_changed_dict["C->G"]) / (max(stat_stable_dict["A"] + stat_stable_dict["C"],1))      
                       
    return change_a, change_c, change_g,change_t, change_dms
    
    
def analysis_normd2(bamfile, starts, ends, strand, chr_, name_, values, values2, coverage, coverage2, gene_name, transcript_id, transcript_name, collection, flag = "complete"):


    dict_for_genes_99_20 = collection[0]
    dict_for_genes_95_20 = collection[1]
    dict_for_genes_90_20 = collection[2]
    dict_all = collection[3]

    saved_gene = [0,0,0,0,0]
    saved_gene2 = [0,0,0,0,0]
    
    #bamfile = pysam.AlignmentFile("./data/DMS_1lb.bam")
    bamfile = pysam.AlignmentFile(bamfile)


    
    b_val = values[0]
    s_val = values[1]
    t_val = values[2]
    f_val = values[3]
    fi_val = values[4]
    

    
    b_val2 = values2[0]
    s_val2 = values2[1]
    t_val2 = values2[2]
    f_val2 = values2[3]
    fi_val2 = values2[4]
    
    to_save_list = []
    last_transcript_name = transcript_name[0]
        
    

    for enum_start, exon_starts in enumerate(starts):
    
        exon_ends = ends[enum_start]


        
        best_coverage, second_coverage, third_coverage, fourth_coverage, fifth_coverage = [], [], [], [], []
        
        dms_data = []
        dms_data_peaked = []

        
        full_length_exons = len(exon_starts)
        full_data_save = []
        full_data_forward_reverse = []
        full_data_reverse_save = []
        
        
        stat_stable_dict = {"A": 0, "C": 0, "T": 0, "G": 0}
        stat_changed_dict = {"A->C": 0, "A->T": 0, "A->G": 0, 
                       "C->A": 0, "C->T": 0, "C->G": 0, 
                       "G->A": 0, "G->T": 0, "G->C": 0, 
                        "T->A": 0, "T->C": 0, "T->G": 0}
        

        for ex_num, exon_start in enumerate(exon_starts):

        
            occupied_list = collections.defaultdict(list)
            mutated_list = collections.defaultdict(list)
            occupied_list_forward = collections.defaultdict(list)
            occupied_list_reverse = collections.defaultdict(list)            
            mutated_list_forward = collections.defaultdict(list)
            mutated_list_reverse = collections.defaultdict(list)

            for i in range(exon_start,exon_ends[ex_num]):
                occupied_list[i] = []
                mutated_list[i] = []
                occupied_list_forward[i] = []
                occupied_list_reverse[i] = []
                mutated_list_forward[i] = []
                mutated_list_reverse[i] = []
                
                
            if exon_start == exon_ends[ex_num]:
              #  print(exon_start)
              #  print(exon_ends[ex_num])
                full_length_exons = full_length_exons - 1
                continue

            reads = bamfile.fetch(chr_, start = exon_start, stop = exon_ends[ex_num])
                    
            for read in reads:
    
                read_pairs = read.get_aligned_pairs(matches_only = True, with_seq = True)
                flag = False
                


                for r in read_pairs:
                

                    if r[1] < exon_ends[ex_num] and r[1] >= exon_start:
                    
                    
                    
                         if occupied_list[r[1]] == []: occupied_list[r[1]] = 0
                         occupied_list[r[1]] = occupied_list[r[1]] + 1
                    
                    
                         if read.is_reverse == True:
                             if occupied_list_reverse[r[1]] == []: occupied_list_reverse[r[1]] = 0
                             occupied_list_reverse[r[1]] = occupied_list_reverse[r[1]] + 1
                         else:
                             if occupied_list_forward[r[1]] == []: occupied_list_forward[r[1]] = 0
                             occupied_list_forward[r[1]] = occupied_list_forward[r[1]] + 1
                         
                         

                         
                         if r[2].islower() == True:

                             if read.is_reverse == False:
                         
                                 if mutated_list[r[1]] == []: mutated_list[r[1]] = 0
                                 mutated_list[r[1]] = mutated_list[r[1]] + 1
                                 
                                 if mutated_list_forward[r[1]] == []: mutated_list_forward[r[1]] = 0
                                 mutated_list_forward[r[1]] = mutated_list_forward[r[1]] + 1
                                 
                             if read.is_reverse == True:
                         
                                 if mutated_list[r[1]] == []: mutated_list[r[1]] = 0
                                 mutated_list[r[1]] = mutated_list[r[1]] + 1
                                 
                                 if mutated_list_reverse[r[1]] == []: mutated_list_reverse[r[1]] = 0
                                 mutated_list_reverse[r[1]] = mutated_list_reverse[r[1]] + 1
                                 
                         stat_stable_dict, stat_changed_dict = count_mutation(stat_stable_dict, stat_changed_dict, r, read, exon_ends[ex_num], exon_start)


            full_data_forward_mutated = [0 if mutated_list_forward[x] == [] else mutated_list_forward[x] for x in range(exon_start, exon_ends[ex_num])]
            full_data_reverse_mutated = [0 if mutated_list_reverse[x] == [] else mutated_list_reverse[x] for x in range(exon_start, exon_ends[ex_num])]
            
            
            full_data_foward2 = [1 if occupied_list_forward[x] == [] else occupied_list_forward[x] for x in range(exon_start, exon_ends[ex_num])]
            full_data_reverse2 = [1 if occupied_list_reverse[x] == [] else occupied_list_reverse[x] for x in range(exon_start, exon_ends[ex_num])]
            
            
            full_data_forward_mutated = [m/full_data_foward2[enum] for enum, m in enumerate(full_data_forward_mutated)]
            

            
            full_data_reverse_mutated = [m/full_data_reverse2[enum] for enum, m in enumerate(full_data_reverse_mutated)]

            
            for enum, data in enumerate(full_data_forward_mutated):
                full_data_forward_mutated[enum] = full_data_forward_mutated[enum] + full_data_reverse_mutated[enum]
            
            full_data_mutated = full_data_forward_mutated
            

            
            full_data = [0 if occupied_list[x] == [] else occupied_list[x] for x in range(exon_start, exon_ends[ex_num])]
            full_data_forward = [0 if occupied_list_forward[x] == [] else occupied_list_forward[x] for x in range(exon_start, exon_ends[ex_num])]
            full_data_reverse = [0 if occupied_list_reverse[x] == [] else occupied_list_reverse[x] for x in range(exon_start, exon_ends[ex_num])]
            
            full_data_mutated = [m  if full_data_forward[enum] > 0 else "NA" for enum, m in enumerate(full_data_mutated)]
            
            full_data = [min(20, x) for x in full_data]
            
            if len(full_data) == 0: 

                raise NotImplementedError
    


            # a99_full_coverage_2, a95_full_coverage_2, a90_full_coverage_2
            a99_full_coverage_20, a95_full_coverage_20, a90_full_coverage_20 = coverage_metric_check(full_data)
            full_data_save.extend(full_data)
            dms_data.append(full_data_mutated)
            


        if full_data_save == []:
            continue

        #a95_full_coverage_5, a90_full_coverage_5, a95_full_coverage_2, a90_full_coverage_2 = coverage_metric_check(full_data_save)
        a99_full_coverage_20, a95_full_coverage_20, a90_full_coverage_20 = coverage_metric_check(full_data_save)
        
        change_a, change_c, change_g,change_t, change_dms = metrics(stat_stable_dict, stat_changed_dict)
        


        
        if a99_full_coverage_20: dict_for_genes_99_20[transcript_id[enum_start]].append([[chr_, [exon_starts, ends[enum_start]]],gene_name[0], transcript_name, transcript_id[enum_start][0], strand, change_a, change_c, change_g, change_t, change_dms])
        if a95_full_coverage_20: dict_for_genes_95_20[transcript_id[enum_start]].append([[chr_, [exon_starts, ends[enum_start]]],gene_name[0], transcript_name, transcript_id[enum_start][0], strand, change_a, change_c, change_g, change_t, change_dms]) 
        if a90_full_coverage_20: dict_for_genes_90_20[transcript_id[enum_start]].append([[chr_, [exon_starts, ends[enum_start]]],gene_name[0], transcript_name, transcript_id[enum_start][0], strand, change_a, change_c, change_g, change_t, change_dms])
        dict_all[transcript_id[enum_start]].append([[chr_, [exon_starts, ends[enum_start]]],gene_name[0], transcript_name, transcript_id[enum_start][0],strand, change_a, change_c, change_g, change_t, change_dms])
        
    values1 = [b_val, s_val, t_val, f_val, fi_val]
    values2 = [b_val2, s_val2, t_val2, f_val2, fi_val2]
        
    collection = [dict_for_genes_99_20, dict_for_genes_95_20, dict_for_genes_90_20, dict_all]
    

        
    return values1,values2, saved_gene, saved_gene2, collection
    
    
    
    
def save_dms(chr_, exon_start,exon_end, dms_data, file_name):

    file_ =  open(file_name, "a+")
    
    file_.write("> chr" + chr_ + " " + str(exon_start) + " " + str(exon_end))
    file_.write("\n")
    
    for val in dms_data:
    
        file_.write(str(val))
        file_.write(" ")
                
    file_.write("\n")
    file_.close()
    

    return
    
    
    
    
    
def save_conjunction(file_name, in_list):

    file_ =  open(file_name, "a+")
    
    for val in in_list:
    
        file_.write(str(val))
        file_.write(" ")
                
    file_.close()
    


    return
    
    
    
    
def count_mutation(stat_stable_dict, stat_changed_dict, r, read, exon_end, exon_start):




                
    if r[1] < exon_end and r[1] >= exon_start:
                
        forward_string = read.get_forward_sequence()
             
        
                     
        if read.is_reverse: forward_string = forward_string[::-1]                
        if r[2].islower() == False:
             stat_stable_dict[forward_string[r[0]]] += 1
        else:

                        
            if read.is_reverse == True:
                        
                if r[2] == "t": complement = "A"
                if r[2] == "g": complement = "C"
                if r[2] == "c": complement = "G"
                if r[2] == "a": complement = "T"
                            
                            
            else:
                complement = r[2].upper()

                        
                        
            if complement != "N" and forward_string[r[0]] != "N":
                        
                ke =  complement + "->" + forward_string[r[0]]
                stat_changed_dict[ke] += 1


    return stat_stable_dict, stat_changed_dict
    
    
    
    
def create_gene_sublist(starts, ends, exon_number):


    full_gene_exons_start = []
    full_gene_exons_end = []
    
    for enum, start in enumerate(starts):
        exons = exon_number[enum]
        
        exon_list_start_full = []
        exon_list_end_full = []
        
        exon_list_start = []
        exon_list_end = []
        
        last_ex = 1
        
        for ex_num, ex in enumerate(exons):
            if ex > last_ex:
                exon_list_start.append(starts[enum][ex_num])
                exon_list_end.append(ends[enum][ex_num])
                
                
                last_ex = ex
                
            else:
                if exon_list_start != []: exon_list_start_full.append(exon_list_start)
                if exon_list_end != []: exon_list_end_full.append(exon_list_end)
                
                exon_list_start = []
                exon_list_end = []
    
                exon_list_start.append(starts[enum][ex_num])
                exon_list_end.append(ends[enum][ex_num])
                
                last_ex = ex
                
                
        if exon_list_start != []:
            exon_list_start_full.append(exon_list_start)
            exon_list_end_full.append(exon_list_end)
                
        full_gene_exons_start.append(exon_list_start_full)
        full_gene_exons_end.append(exon_list_end_full)
                

                
    return full_gene_exons_start, full_gene_exons_end
    
    
    
def coverage_metric_check(exon_list):
   
    a99_full_coverage_20 = True
    a95_full_coverage_20 = True
    a90_full_coverage_20 = True
    
    smaller_20 = 0
    
    
    for val in exon_list:

        if val < 20:
            smaller_20 = smaller_20 + 1


    if smaller_20/len(exon_list) >= 0.01:
         a99_full_coverage_20 = False
    if smaller_20/len(exon_list) >= 0.05:
         a95_full_coverage_20 = False
    if smaller_20/len(exon_list) >= 0.10:
         a90_full_coverage_20 = False
    
        
    return a99_full_coverage_20, a95_full_coverage_20, a90_full_coverage_20
    
    
    
def coverage_metric(avg_list):


    full_coverage = [[],[],[]]
    two_coverage = [[],[],[]]
    five_coverage = [[],[],[]]
    ten_coverage = [[],[],[]]


    for l in avg_list:
        

        i = 1
        for val in l:
            if val < 10:
                i = 0
        full_coverage[0].append(i)
        
        i = 1
        for val in l:
            if val < 5:
                i = 0
        full_coverage[1].append(i)
        
        i = 1
        for val in l:
            if val < 2:
                i = 0
        full_coverage[2].append(i)
        
            
        smaller_10 = 0
        smaller_5 = 0
        smaller_2 = 0
            
        for val in l:
            if val < 10:
                smaller_10 = smaller_10 + 1
            if val < 5:
                smaller_5 = smaller_5 + 1
            if val < 2:
                smaller_2 = smaller_2 + 1
                
                
        if smaller_10/len(l) >= 0.02:
             two_coverage[0].append(0)
        else:
             two_coverage[0].append(1)
        
        if smaller_5/len(l) >= 0.02:
             two_coverage[1].append(0)
        else:
             two_coverage[1].append(1)
             
        if smaller_2/len(l) >= 0.02:
             two_coverage[2].append(0)
        else:
             two_coverage[2].append(1)
            
        if smaller_10/len(l) >= 0.05:
             five_coverage[0].append(0)
        else:
             five_coverage[0].append(1)
        
        if smaller_5/len(l) >= 0.05:
             five_coverage[1].append(0)
        else:
             five_coverage[1].append(1)
             
        if smaller_2/len(l) >= 0.05:
             five_coverage[2].append(0)
        else:
             five_coverage[2].append(1)
             
        if smaller_10/len(l) >= 0.10:
             ten_coverage[0].append(0)
        else:
             ten_coverage[0].append(1)
        
        if smaller_5/len(l) >= 0.10:
             ten_coverage[1].append(0)
        else:
             ten_coverage[1].append(1)
             
        if smaller_2/len(l) >= 0.10:
             ten_coverage[2].append(0)
        else:
             ten_coverage[2].append(1)
             


    return full_coverage, two_coverage, five_coverage, ten_coverage
    
    
    
    
    
def uniques(ind_list, chr_list, name_list):

    lastval = None
    
    new_ind_list = []
    new_chr_list = []
    new_name_list = []
    
    for enum, val in enumerate(ind_list):
        if val != lastval:
            lastval = val
            new_ind_list.append(val)
            new_chr_list.append(chr_list[enum])
            new_name_list.append(name_list[enum])
    


    return new_ind_list, new_chr_list, new_name_list



    
def etract_refseq_utr_genome_gff(annotation="/media/sven/Intenso/95_20_stuff/data/genome.gff"):


    hg_file = open(annotation)
    exon_dict = {}
    start_dict = {}
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
            if "transcript_id" not in line: continue
            transcript_id = line.split("transcript_id")[1].split(";")[0][1:-1]
            
            
            if True:
            
                if transcript_id != last_tid:
                    transcript_lists.append(transcript_id)
                    exon_dict[transcript_id] = []
                    start_dict[transcript_id] = []
                    gene_name_dict[transcript_id] = ""
                    counter = counter + 1
                    last_tid = transcript_id
                
                gen_type_ = line.split("\t")[2]

            
                if gen_type_ == "exon":

                    gene_name_dict[transcript_id] = line.split("gene=")[1].split(";")[0]


                    start = line.split("\t")[3]
                    end = line.split("\t")[4]
                    exon_dict[transcript_id].append([int(start), int(end)])
                    strand_dict[transcript_id] = line.split("\t")[6]
                    chr_dict[transcript_id] = line.split("\t")[0]
                if gen_type_ == "start_codon":
                    start = line.split("\t")[3]
                    end = line.split("\t")[4]
                    start_dict[transcript_id]=[int(start), int(end)]
                

                if gen_type_ == "CDS":
                    haCDSs_dict[transcript_id] = True
            
    

    
    
    
            
    hg_file.close()

    #utr_regions = extract_utr_region(transcript_lists, exon_dict, strand_dict)
    utr_regions = extract_utr_region(transcript_lists, exon_dict, start_dict, strand_dict)


    #print(strand_dict)
    
    return_chr_type = {"NC_000001":"chr1","NC_000002":"chr2","NC_000003":"chr3","NC_000004":"chr4","NC_000005":"chr5","NC_000006":"chr6","NC_000007":"chr7","NC_000008":"chr8","NC_000009":"chr9"
    ,"NC_000010":"chr10","NC_000011":"chr11","NC_000012":"chr12","NC_000013":"chr13","NC_000014":"chr14","NC_000015":"chr15","NC_000016":"chr16","NC_000017":"chr17","NC_000018":"chr18","NC_000019":"chr19",
    "NC_000020":"chr20","NC_000021":"chr21","NC_000022":"chr22","NC_000023":"chrX","NC_000024":"chrY","NC_012920":"chrMT" }
    
    for key in list(chr_dict.keys()):

    
        if chr_dict[key].split(".")[0] in list(return_chr_type.keys()): chr_dict[key] = return_chr_type[chr_dict[key].split(".")[0]]
        
        if chr_dict[key] == "chrMT":
            print("chrMT")
        
   # print(chr_dict)
   

    

    
    
    return utr_regions, gene_name_dict, strand_dict, chr_dict
    
    
def etract_refseq_utr(gff_path):


    hg_file = open(gff_path)
    exon_dict = {}
    start_dict = {}
    transcript_lists = []
    strand_dict = {}
    gene_name_dict = {}
    counter = 0
    chr_dict = {}
    haCDSs_dict = {}
    last_tid = ""
    for enum, line in enumerate(hg_file):
        line_data = line.split()
        
        if counter%10000==0: print("counter: " + str(counter))
            
        
        if enum > 4:
            if line[0] == "#":continue
            if "transcript_id" not in line: print(line)
            transcript_id = line.split("transcript_id")[1].split(";")[0][2:-1]

            
            if transcript_id != last_tid:
                transcript_lists.append(transcript_id)
                exon_dict[transcript_id] = []
                start_dict[transcript_id] = []
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
            if gen_type_ == "start_codon":
                start = line.split("\t")[3]
                end = line.split("\t")[4]
                start_dict[transcript_id]=[int(start), int(end)]
                

            if gen_type_ == "CDS":
                haCDSs_dict[transcript_id] = True
            
    

    
    
    
            
    hg_file.close()

    #utr_regions = extract_utr_region(transcript_lists, exon_dict, start_dict, strand_dict)
    utr_regions = extract_full_exon(transcript_lists, exon_dict, start_dict, strand_dict)


    #print(strand_dict)
    
    return_chr_type = {"NC_000001":"chr1","NC_000002":"chr2","NC_000003":"chr3","NC_000004":"chr4","NC_000005":"chr5","NC_000006":"chr6","NC_000007":"chr7","NC_000008":"chr8","NC_000009":"chr9"
    ,"NC_000010":"chr10","NC_000011":"chr11","NC_000012":"chr12","NC_000013":"chr13","NC_000014":"chr14","NC_000015":"chr15","NC_000016":"chr16","NC_000017":"chr17","NC_000018":"chr18","NC_000019":"chr19",
    "NC_000020":"chr20","NC_000021":"chr21","NC_000022":"chr22","NC_000023":"chrX","NC_000024":"chrY","NC_012920":"chrMT" }
    
    for key in list(chr_dict.keys()):

    
        if chr_dict[key].split(".")[0] in list(return_chr_type.keys()): chr_dict[key] = return_chr_type[chr_dict[key].split(".")[0]]
        
        if chr_dict[key] == "chrMT":
            print("chrMT")
        
   # print(chr_dict)
   

    

    
    
    return utr_regions, gene_name_dict, strand_dict, chr_dict

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

    
def extract_utr_region(transcript_lists, exon_dict, start_dict, strand_dict):

    utrs = {}
    
    #print(transcript_lists["NM_001354840.3"])

  
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
    
    print("return")
    print(utrs)
    print("end")
    return utrs


    
    """
    
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
                
                
                    if exon[1] >= start_dict[key][0]:
                    
                    
                        last_nuc_number = extra_nuc

                        if exon[0] < start_dict[key][0]: 

                             extra_nuc =  extra_nuc + np.abs(exon[1] - start_dict[key][0])
                        else: 

                             extra_nuc =  extra_nuc + np.abs(exon[1] - exon[0])
                        


                        
                        
                        if extra_nuc >= 200:
                        
                            offset = (exon[0] + (200-last_nuc_number)) - start_dict[key][0]
                            
                            break;
                            
                        
                         
            for exon in exon_dict[key]:
                
                if exon[0] < start_dict[key][0] +offset:
            
              
                    if exon[1] >= start_dict[key][0] + offset:

                        exon[1] = start_dict[key][0] - 1 + offset

                    exon_dict_list.append(exon)
        
        
        
        
        
        
        elif strand == "-":
        
            offset = 0
            extra_nuc = 0
            
            for exon in exon_dict[key]:
                
                
                    if exon[0] <= start_dict[key][1]:
                    
                    
                        last_nuc_number = extra_nuc
                        
                        
                        if exon[1] > start_dict[key][1]: 

                             extra_nuc =  extra_nuc + np.abs(start_dict[key][1] - exon[0])
                        else: 

                             extra_nuc =  extra_nuc + np.abs(exon[1] - exon[0])
                        
                        if extra_nuc >= 200:
                        
                            offset = start_dict[key][0] - (exon[1] - (200-last_nuc_number))
                            
                            break;
        
        
            for exon in exon_dict[key]:
            
            

                if exon[1] > start_dict[key][1] - offset:
                
                    if exon[0] <= start_dict[key][1] -offset:
                        exon[0] = start_dict[key][1] + 2 - offset
                   
        
                    if exon[0] <= exon[1]:exon_dict_list.append(exon)
        
        
        utrs[key] = exon_dict_list
        exon_dict_list = []



    return utrs
    
    """

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="loads potential DMS candidates")
    
    # Define arguments
    parser.add_argument("--gff_path", type=str, required=True, help="Path to gff file")
    parser.add_argument("--bamfile", type=str, required=True, help="Path to sample bam file")
    parser.add_argument("--target_folder", type=str, required=True, help="Path to output folder")
    # Parse arguments
    args = parser.parse_args()



    # 1l beta
    utr_regions, gene_name_dict, strand_dict, chr_dict = etract_refseq_utr(args.gff_path)

    return_chr_type = {"NC_000001":"chr1","NC_000002":"chr2","NC_000003":"chr3","NC_000004":"chr4","NC_000005":"chr5","NC_000006":"chr6","NC_000007":"chr7","NC_000008":"chr8","NC_000009":"chr9"
    ,"NC_000010":"chr10","NC_000011":"chr11","NC_000012":"chr12","NC_000013":"chr13","NC_000014":"chr14","NC_000015":"chr15","NC_000016":"chr16","NC_000017":"chr17","NC_000018":"chr18","NC_000019":"chr19",
    "NC_000020":"chr20","NC_000021":"chr21","NC_000022":"chr22","NC_000023":"chrX","NC_000024":"chrY","NC_012920":"chrMT" }

    #return_chr_type = {"NC_000001":"1","NC_000002":"2","NC_000003":"3","NC_000004":"4","NC_000005":"5","NC_000006":"6","NC_000007":"7","NC_000008":"8","NC_000009":"9"
    #,"NC_000010":"10","NC_000011":"11","NC_000012":"12","NC_000013":"13","NC_000014":"14","NC_000015":"15","NC_000016":"16","NC_000017":"17","NC_000018":"18","NC_000019":"19",
    #"NC_000020":"20","NC_000021":"21","NC_000022":"22","NC_000023":"X","NC_000024":"Y","NC_012920":"MT" }
    
    
    dict_for_genes_99_20 = collections.defaultdict(list)
    dict_for_genes_95_20 = collections.defaultdict(list)
    dict_for_genes_90_5 = collections.defaultdict(list)
    dict_all = collections.defaultdict(list)
    # utr_regions, gene_name_dict

    collection = [dict_for_genes_99_20, dict_for_genes_95_20, dict_for_genes_90_5, dict_all]
    flags = ["reverse"]
    

    
    for flag in flags:
    
        iterator = 10
        counter = 0


        values1 = [0,0,0,0,0]
        coverage1 = [[],[],[],[], []]    
        saved_gene_list = []
        saved_gene2_list = []
        values2 = [0,0,0,0,0]
        coverage2 = [[],[],[],[], []]    
        

    
        avg_protein_list = []

       
        for en, key in enumerate(list(utr_regions.keys())):
               
                chr_ =chr_dict[key]
                
                if chr_ == "chrMT": continue

                if en % 1000 == 0: 
                    print(en)

                strand = strand_dict[key]
                
                
                starts_sublists = [x[0] for x in utr_regions[key]]
                end_sublists = [x[1] for x in utr_regions[key]]

                
                if chr_ not in list(return_chr_type.values()): continue
                    
                values1,values2, saved_gene, saved_gene2, collection = analysis_normd2(args.bamfile, [starts_sublists], [end_sublists], strand, chr_, [key], values1, values2, coverage1, coverage2, [key], [key], gene_name_dict[key], collection, flag)
                

                    
                    
                saved_gene_list.append(saved_gene)
                saved_gene2_list.append(saved_gene2)
                    
                    
        dict_for_genes_99_20 = collection[0]
        dict_for_genes_95_20 = collection[1]
        dict_for_genes_90_20 = collection[2]
        dict_all = collection[3]
        
        output_name = "/" + args.bamfile.split("/")[-1].split(".bam")[0] + "_"


        with open(args.target_folder + output_name + "99_20.pkl", 'wb') as f:
            pickle.dump(dict_for_genes_99_20, f)
            
        with open(args.target_folder + output_name + "95_20.pkl", 'wb') as f:
            pickle.dump(dict_for_genes_95_20, f)
            
        with open(args.target_folder + output_name + "90_20.pkl", 'wb') as f:
            pickle.dump(dict_for_genes_90_20, f)
            
            
        with open(args.target_folder + output_name + "all.pkl", 'wb') as f:
            pickle.dump(dict_all, f)
   
    
    
