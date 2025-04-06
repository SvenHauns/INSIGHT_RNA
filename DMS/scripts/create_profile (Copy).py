import pysam
import inspect
mismatch_count = 0
import time
import psutil
import multiprocessing
import time
from multiprocessing import Pool
import pandas as pd
import collections
import pickle
import glob
import random
import numpy as np
import matplotlib.pyplot as plt
import pickle
from Bio import SeqIO




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



    
    

    
    

def save_bed(chr_, start, end, strand, filename):

    name_tag = chr_ + "_" + str(start) + "_" + str(end) 
    file_ = open(filename, "a+")
    
    file_.write("chr" + (str(chr_)))
    file_.write("\t")
    file_.write(str(start))
    file_.write("\t")
    file_.write(str(end))
    file_.write("\t")
    file_.write(name_tag)
    file_.write("\t")
    file_.write("0")
    file_.write("\t")
    file_.write(str(strand))
    file_.write("\n")
    
    file_.close()

    return

# videmp@informatik.uni-freiburg.de




def rna_seq_mutation(chr_, exon_start, exon_end):

    bamfile = pysam.AlignmentFile("/media/sven/Intenso/95_20_stuff/data/RNA-NOIL-3.bam")
    
    
    occupied_list = collections.defaultdict(list)
    mutated_list = collections.defaultdict(list)
    occupied_list_forward = collections.defaultdict(list)
    occupied_list_reverse = collections.defaultdict(list)            
    mutated_list_forward = collections.defaultdict(list)
    mutated_list_reverse = collections.defaultdict(list)
            


    for i in range(exon_start,exon_end):
        occupied_list[i] = []
        mutated_list[i] = []
        occupied_list_forward[i] = []
        occupied_list_reverse[i] = []
        mutated_list_forward[i] = []
        mutated_list_reverse[i] = []
                
                

            
            
            
    reads = bamfile.fetch(chr_, start = exon_start, stop = exon_end)

            

    for read in reads:
    
            read_pairs = read.get_aligned_pairs(matches_only = True, with_seq = True)
            

                
            for r in read_pairs:
                
                       
                if r[1] < exon_end and r[1] >= exon_start:
                    
                     if occupied_list[r[1]] == []: occupied_list[r[1]] = 0
                     occupied_list[r[1]] = occupied_list[r[1]] + 1
                    
                    
                     if read.is_reverse == True:
                         if occupied_list_reverse[r[1]] == []: occupied_list_reverse[r[1]] = 0
                         occupied_list_reverse[r[1]] = occupied_list_reverse[r[1]] + 1
                     else:
                         if occupied_list_forward[r[1]] == []: occupied_list_forward[r[1]] = 0
                         occupied_list_forward[r[1]] = occupied_list_forward[r[1]] + 1
                         

                         
                     if r[2].islower() == True:

                         if read.is_reverse == False and r[2].islower() == True:
                         
                             if mutated_list[r[1]] == []: mutated_list[r[1]] = 0
                             mutated_list[r[1]] = mutated_list[r[1]] + 1
                                 
                             if mutated_list_forward[r[1]] == []: mutated_list_forward[r[1]] = 0
                             mutated_list_forward[r[1]] = mutated_list_forward[r[1]] + 1
                                
                         if read.is_reverse == True and r[2].islower() == True and (r[2] == "g" or r[2] == "t"  or r[2] == "a" or r[2] == "c"):
                         
                             if mutated_list[r[1]] == []: mutated_list[r[1]] = 0
                             mutated_list[r[1]] = mutated_list[r[1]] + 1
                                 
                             if mutated_list_reverse[r[1]] == []: mutated_list_reverse[r[1]] = 0
                             mutated_list_reverse[r[1]] = mutated_list_reverse[r[1]] + 1
                                 
                                 
                                 

                             
    full_data_forward_mutated = [0 if mutated_list_forward[x] == [] else mutated_list_forward[x] for x in range(exon_start, exon_end)]
    full_data_reverse_mutated = [0 if mutated_list_reverse[x] == [] else mutated_list_reverse[x] for x in range(exon_start, exon_end)]
    occupied_list = [0 if occupied_list[x] == [] else occupied_list[x] for x in range(exon_start, exon_end)]

            
            
    full_mutations = [d + full_data_reverse_mutated[enum] for enum, d in enumerate(full_data_forward_mutated)]
        
        
    ############## LIMIT TO OCCUPATION #################
        
        
    full_mutations = [x if occupied_list[en] >= 5 else 0 for en, x in enumerate(full_mutations)]
        

        


            
    return occupied_list, full_mutations



def norm_in_window(full_data_occupied, full_data_both_mutated_not_normd):



    binsize = int(len(full_data_occupied)/5) + 1
    normd_window = [[full_data_occupied[r*binsize:(r+1)*binsize]] for r in range(0, 5)]
            

    normd_window = [m for m in normd_window if m[0] != [] ]
    normd_window = [max(max(m[0]),1) for m in normd_window]
            
    normd_window_forward = []
    for enum, window in enumerate(normd_window):
        t_values = [t/window for t in full_data_both_mutated_not_normd[enum*binsize:(enum+1)*binsize]]               
        normd_window_forward.extend(t_values)


 
    return normd_window_forward



def analysis_normd2(starts, ends, chr_, strand, name_, values, values2, coverage, coverage2, flag, sequences, bamfile_path, is_rna = False, bam_num = 0):

    saved_gene = [0,0,0,0,0]
    saved_gene2 = [0,0,0,0,0]
    


    
    bamfile = pysam.AlignmentFile(bamfile_path)

    best_coverage_comb = coverage[0]
    second_coverage_comb = coverage[1]
    third_coverage_comb = coverage[2]
    fourth_coverage_comb = coverage[3]
    fifth_coverage_comb = coverage[4]
    
    b_val = values[0]
    s_val = values[1]
    t_val = values[2]
    f_val = values[3]
    fi_val = values[4]
    

    best_coverage_comb2 = coverage2[0]
    second_coverage_comb2 = coverage2[1]
    third_coverage_comb2 = coverage2[2]
    fourth_coverage_comb2 = coverage2[3]
    fifth_coverage_comb2 = coverage2[4]
    
    b_val2 = values2[0]
    s_val2 = values2[1]
    t_val2 = values2[2]
    f_val2 = values2[3]
    fi_val2 = values2[4]
    






    full_data_both_mutated_normd = []
    full_data_both_mutated_not_normd2_list = []
    full_data_occupied = []
    full_data_occupied_list = []
    constraint_1_list = []
    constraint_2_list = []
    constraint_3_list = []
    constraint_4_list = []
    
    
    

    
            
    for enum, exon_starts in enumerate(starts):
    

        normd_gene_data = []
    
        exon_ends = ends[enum]
        
        best_coverage, second_coverage, third_coverage, fourth_coverage, fifth_coverage = [], [], [], [], []
        
        dms_data = []
        dms_data_peaked = []

        
        full_length_exons = len(exon_starts)
        #print(full_length_exons)
        full_data_forward_save = []
        full_data_forward_reverse = []
        full_data_reverse_save = []
        
        
        
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
                full_length_exons = full_length_exons - 1
                continue


            reads = bamfile.fetch(chr_, start = exon_start, stop = exon_ends[ex_num])

            

            for read in reads:
    
                read_pairs = read.get_aligned_pairs(matches_only = True, with_seq = True)

                
                for r in read_pairs:
                
                       
                    #if r[1] < exon_ends[ex_num] + 1 and r[1] >= exon_start:
                    
                    
                    #if r[1] < exon_ends[ex_num] and r[1] >= exon_start -1:
                    
                    if r[1] < exon_ends[ex_num] and r[1] >= exon_start:
                    
                         l = [r for r in read_pairs if r[1] in range(exon_start, exon_ends[ex_num])]


                         
                         
                    
                         if occupied_list[r[1]] == []: occupied_list[r[1]] = 0
                         occupied_list[r[1]] = occupied_list[r[1]] + 1
                    
                    
                         if read.is_reverse == True:
                             if occupied_list_reverse[r[1]] == []: occupied_list_reverse[r[1]] = 0
                             occupied_list_reverse[r[1]] = occupied_list_reverse[r[1]] + 1
                         else:
                             if occupied_list_forward[r[1]] == []: occupied_list_forward[r[1]] = 0
                             occupied_list_forward[r[1]] = occupied_list_forward[r[1]] + 1
                         

                         
                         if r[2].islower() == True:

                             if read.is_reverse == False and r[2].islower() == True:
                         
                                 if mutated_list[r[1]] == []: mutated_list[r[1]] = 0
                                 mutated_list[r[1]] = mutated_list[r[1]] + 1
                                 
                                 if mutated_list_forward[r[1]] == []: mutated_list_forward[r[1]] = 0
                                 mutated_list_forward[r[1]] = mutated_list_forward[r[1]] + 1
                                 
                             if read.is_reverse == True and r[2].islower() == True:
                         
                                 if mutated_list[r[1]] == []: mutated_list[r[1]] = 0
                                 mutated_list[r[1]] = mutated_list[r[1]] + 1
                                 
                                 if mutated_list_reverse[r[1]] == []: mutated_list_reverse[r[1]] = 0
                                 mutated_list_reverse[r[1]] = mutated_list_reverse[r[1]] + 1
                                 
                                 
                                 


            full_data_forward_mutated = [0 if mutated_list_forward[x] == [] else mutated_list_forward[x] for x in range(exon_start, exon_ends[ex_num])]
            full_data_reverse_mutated = [0 if mutated_list_reverse[x] == [] else mutated_list_reverse[x] for x in range(exon_start, exon_ends[ex_num])]

            
            

            #full_data_reverse_mutated = [0 if full_data_reverse_mutated[x] > 100 else full_data_reverse_mutated[x] for x in range(0, len(full_data_reverse_mutated))]
            


            full_data_foward2 = [0 if occupied_list_forward[x] == [] else occupied_list_forward[x] for x in range(exon_start, exon_ends[ex_num])]
            full_data_reverse2 = [0 if occupied_list_reverse[x] == [] else occupied_list_reverse[x] for x in range(exon_start, exon_ends[ex_num])]
            

            
            full_data_forward_mutated_not_normd = [m for enum, m in enumerate(full_data_forward_mutated)]
            full_data_forward_mutated = [m/max(full_data_foward2[enum],1) for enum, m in enumerate(full_data_forward_mutated)]
            

            
            full_data_reverse_mutated_not_normd = [m for enum, m in enumerate(full_data_reverse_mutated)]
            full_data_reverse_mutated = [m/max(full_data_reverse2[enum],1) for enum, m in enumerate(full_data_reverse_mutated)]
            
            


            
            for enum, data in enumerate(full_data_forward_mutated):
            #    if full_data_reverse_mutated[enum] != 0 and full_data_forward_mutated[enum]  != 0: raise NotImplementedError
                full_data_forward_mutated[enum] = full_data_forward_mutated[enum] + full_data_reverse_mutated[enum]


            full_data_mutated = full_data_forward_mutated
            
            

            


            
            full_data = [0 if occupied_list[x] == [] else occupied_list[x] for x in range(exon_start, exon_ends[ex_num])]
            full_data_forward = [0 if occupied_list_forward[x] == [] else occupied_list_forward[x] for x in range(exon_start, exon_ends[ex_num])]
            full_data_reverse = [0 if occupied_list_reverse[x] == [] else occupied_list_reverse[x] for x in range(exon_start, exon_ends[ex_num])]
            
            full_data_mutated = [m  if full_data_forward[enum] > 0 or full_data_reverse[enum] >0 else "NA" for enum, m in enumerate(full_data_mutated)]
            
            
            
            
            
            full_data_both_mutated_not_normd_ = [full_data_reverse_mutated_not_normd[enum] + m for enum, m in enumerate(full_data_forward_mutated_not_normd)]
            
            
            full_data_occupied_ = [full_data_foward2[enum] + m for enum, m in enumerate(full_data_reverse2)]
            
            
            
            full_data_both_mutated_not_normd_ = [d if full_data_occupied_[enum] > 5 else 0 for enum, d in enumerate(full_data_both_mutated_not_normd_)]
            full_data_both_mutated_not_normd2_ = [d for enum, d in enumerate(full_data_both_mutated_not_normd_)]
            
            
            full_data_normd = norm_in_window(full_data_occupied_, full_data_both_mutated_not_normd_)


 
            #is_rna = True
            """
            if is_rna == False:
            
                rna_mutations_occupied, rna_mutations = rna_seq_mutation(chr_,exon_starts[ex_num], exon_ends[ex_num])
            

                #### have to norm first before you substract ####
            
            
                rna_mutations_normd =  norm_in_window(rna_mutations_occupied, rna_mutations)
                        
                full_data_normd = [max(x - rna_mutations_normd[enum],0) for enum, x in enumerate(full_data_normd)]
            
            
            """
            #full_data_occupied_ = [x if x <=5 else 5 for x in full_data_occupied_]
            full_data_occupied_list.append(full_data_occupied_)
            full_data_both_mutated_normd.append(full_data_normd)
            full_data_both_mutated_not_normd2_list.append(full_data_both_mutated_not_normd2_)

            
            constraint_1 = calculate_constraint(full_data_normd, sequences, 0.04, ex_num, strand)
            constraint_2 = calculate_constraint(full_data_normd, sequences, 0.06, ex_num, strand)
            constraint_3 = calculate_constraint(full_data_normd, sequences, 0.25, ex_num, strand)
            constraint_4 = calculate_constraint(full_data_normd, sequences, 0.1, ex_num, strand)
            


            
            
            constraint_1_list.extend(constraint_1)
            constraint_2_list.extend(constraint_2)
            constraint_3_list.extend(constraint_3)
            constraint_4_list.extend(constraint_4)
            
            
            
        ######################## also try mean !


        normd_window_forward = full_data_both_mutated_normd
        window_forward = full_data_both_mutated_not_normd2_list
        
        
        if strand == "-":
            

            for norm2 in window_forward:
                norm2.reverse()
        
        
            for norm in normd_window_forward:
                norm.reverse()
        
            #normd_window_forward.reverse()
            for occ in full_data_occupied_list:
                occ.reverse()
                

            
        normd_window_forward2,full_data_occupied_list2, window_forward2 = [], [], []
        for m in normd_window_forward: normd_window_forward2.extend(m)
        for m in full_data_occupied_list:full_data_occupied_list2.extend(m)
        for m in window_forward:window_forward2.extend(m)
        
        
        normd_window_forward = normd_window_forward2
        full_data_occupied_list = full_data_occupied_list2
        window_forward = window_forward2

             
    


    bamfile.close()
    
    
    
    ############################# write dms to fasta file ###########################################
    
    output_files = ["untreated_3UTR.fa", "untreated2.fa", "rna_seq2.fa"]
    #output_files = ["all_not_dedup_small.fa"]
    output_files = ["all_not_dedup_small_g.fa"]
    
    output_file = output_files[bam_num]
    
    coverage_save = ""
    
    
    for cov_l in full_data_occupied_list:

        coverage_save = coverage_save + str(cov_l) + str(" ")
        


    
    
    
    x_list = ""
    for seq in sequences:
        x_list = x_list + seq    
    
    
    full_data_dms = ""
    for dms in normd_window_forward:
        full_data_dms = full_data_dms + str(dms) + " "
        
        
    raw_data_dms = ""
    for dms in window_forward:
        raw_data_dms = raw_data_dms + str(dms) + " "
    
    file_ = open("./coverage_full/"  + output_file, "a")
    file_.write(">")
    file_.write(name_)
    file_.write("\n")
    file_.write(str(x_list))
    file_.write("\n")
    file_.write(str(full_data_dms))
    file_.write("\n")
    file_.write(str(raw_data_dms))
    file_.write("\n")
    file_.write(str(coverage_save))
    file_.write("\n")
    file_.close()
    


        
        
    constraint_list1 = ""
    for con in constraint_1_list:
        constraint_list1 = constraint_list1 + con + " "
        
    constraint_list2 = ""
    for con in constraint_2_list:
        constraint_list2 = constraint_list2 + con + " "
        
    constraint_list3 = ""
    for con in constraint_3_list:
        constraint_list3 = constraint_list3 + con + " "
        
    constraint_list4 = ""
    for con in constraint_4_list:
        constraint_list4 = constraint_list4 + con + " "
    ############################# write to fasta file ###############################################

    output_files = ["constraints_untreated_3UTR", "constraints_untreated_2", "constraints_rna_seq_2"]
    
    output_file = output_files[bam_num]
    
    
    """


    file_ = open("./95_res/"  + output_file + "_0.04.fa", "a")
    file_.write(">")
    file_.write(name_)
    file_.write("\n")
    file_.write(str(x_list))
    file_.write("\n")
    file_.write(str(constraint_list1))
    file_.write("\n")
    file_.close()
    
    file_ = open("./95_res/"   + output_file + "_0.06.fa", "a")
    file_.write(">")
    file_.write(name_)
    file_.write("\n")
    file_.write(str(x_list))
    file_.write("\n")
    file_.write(str(constraint_list2))
    file_.write("\n")
    file_.close()
    
    
    file_ = open("./95_res/"   + output_file + "_0.1.fa", "a")
    file_.write(">")
    file_.write(name_)
    file_.write("\n")
    file_.write(str(x_list))
    file_.write("\n")
    file_.write(str(constraint_list4))
    file_.write("\n")
    file_.close()
    
    
    file_ = open("./95_res/"   + output_file + "_0.25.fa", "a")
    file_.write(">")
    file_.write(name_)
    file_.write("\n")
    file_.write(str(x_list))
    file_.write("\n")
    file_.write(str(constraint_list3))
    file_.write("\n")
    file_.close()
    """
    
    print("DONE: " + str(name_))
            
    return normd_window_forward, full_data_occupied_list, sequences, name_
    
    # PFN2-201
# PFN2-201
    
def create_fig(sequences, normd_windows, coverage, name_, strand):
    

    
    x_list = ""
    for seq in sequences:
        x_list = x_list + seq
           

    
    color = []
    x_list_ = []
    for enum, x_ in enumerate(x_list):
            
        if x_ == "C":
            color.append("blue")
        elif x_ == "G":
            color.append("yellow")
        elif x_ == "T":
            color.append("grey")
        elif x_ == "A":
            color.append("red")
        x_list_.append(x_)
                
    x_list = x_list_
    
    


    
    
    
    tiltes = ["1lbeta-treated","untreated","RNA seq"]
    
    fig, axes = plt.subplots(1, 3, sharey=False, figsize = (30, 5))

    
    for enum, window in enumerate(normd_windows):
            
        axes[enum].bar(range(0, len(x_list)), window, color = color)
        axes[enum].set_xticks(range(0, len(x_list)))
        axes[enum].set_xticklabels(x_list, fontsize = 6)


        axes[enum].set_title(tiltes[enum])           

        axes[enum].set_ylabel('DMS signal')           
        axes[enum].set_ylim(bottom=0, top=0.35)        
        
        ax_twin = axes[enum].twinx()
        ax_twin.plot(coverage[enum], color='green', linewidth=5)
        ax_twin.grid(None)
        ax_twin.set_ylabel('coverage')
        ax_twin.set_ylim(0, 5)
        ax_twin.tick_params(axis = "y", labelsize = 6)
        axes[enum].tick_params(axis = "y", labelsize = 6)
        
        labels = ["A", "C", "G", "T", "coverage"]
        colors = ["red", "blue", "yellow", "grey", "green"]
            
        handles = [plt.Rectangle((0,0),1,1, color=colors[en]) for en,label in enumerate(labels)]
        axes[enum].legend(handles, labels)
        
        
    fig.suptitle(name_)
    

        
    plt.savefig("/media/sven/Elements/RNA_STAR/deduplicated_files/figs_result/combined/99_5/" + str(name_) +"_" +str(strand) + ".svg")
    plt.close()
        

    return
    
    
def create_multiple(sequences, normd_windows, coverage, name_, strand):

    
    x_list = ""
    for seq in sequences:
        x_list = x_list + seq
           
    color = []
    x_list_ = []
    for enum, x_ in enumerate(x_list):
            
        if x_ == "C":
            color.append("blue")
        elif x_ == "G":
            color.append("yellow")
        elif x_ == "T":
            color.append("grey")
        elif x_ == "A":
            color.append("red")
        x_list_.append(x_)
                
    x_list = x_list_
    
    


    
    
    
    tiltes = ["1lbeta-treated","untreated","RNA seq"]
    
    fig, axes = plt.subplots(20, 3, sharey=False, figsize = (30, 5))
    
    
    
    
    for en, normd_windows in enumerate(normd_windows_list):
    
    
        iterate_axis(axes[en], normd_windows, tiltes, color)
        
        
    
    for enum, window in enumerate(normd_windows):
            
        axes[enum].bar(range(0, len(x_list)), window, color = color)
        axes[enum].set_xticks(range(0, len(x_list)))
        axes[enum].set_xticklabels(x_list, fontsize = 6)


        axes[enum].set_title(tiltes[enum])           

        axes[enum].set_ylabel('DMS signal')           
        axes[enum].set_ylim(bottom=0, top=0.35)        
        
        ax_twin = axes[enum].twinx()
        ax_twin.plot(coverage[enum], color='green', linewidth=5)
        ax_twin.grid(None)
        ax_twin.set_ylabel('coverage')
        ax_twin.set_ylim(0, 5)
        ax_twin.tick_params(axis = "y", labelsize = 6)
        axes[enum].tick_params(axis = "y", labelsize = 6)
        
        labels = ["A", "C", "G", "T", "coverage"]
        colors = ["red", "blue", "yellow", "grey", "green"]
            
        handles = [plt.Rectangle((0,0),1,1, color=colors[en]) for en,label in enumerate(labels)]
        axes[enum].legend(handles, labels)
        
        
    fig.suptitle(name_)



    return
    
    
def iterate_axis(axes, normd_windows, tiltes, color):


    for enum, window in enumerate(normd_windows):
            
        axes[enum].bar(range(0, len(x_list)), window, color = color)
        axes[enum].set_xticks(range(0, len(x_list)))
        axes[enum].set_xticklabels(x_list, fontsize = 6)


        axes[enum].set_title(tiltes[enum])           

        axes[enum].set_ylabel('DMS signal')           
        axes[enum].set_ylim(bottom=0, top=0.35)        
        
        ax_twin = axes[enum].twinx()
        ax_twin.plot(coverage[enum], color='green', linewidth=5)
        ax_twin.grid(None)
        ax_twin.set_ylabel('coverage')
        ax_twin.set_ylim(0, 5)
        ax_twin.tick_params(axis = "y", labelsize = 6)
        axes[enum].tick_params(axis = "y", labelsize = 6)
        
        labels = ["A", "C", "G", "T", "coverage"]
        colors = ["red", "blue", "yellow", "grey", "green"]
            
        handles = [plt.Rectangle((0,0),1,1, color=colors[en]) for en,label in enumerate(labels)]
        axes[enum].legend(handles, labels)
        

    return
    
    
    
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
                
        if  exon_list_start != []:exon_list_start_full.append(exon_list_start)
        if  exon_list_end != []: exon_list_end_full.append(exon_list_end)
                
        full_gene_exons_start.append(exon_list_start_full)
        full_gene_exons_end.append(exon_list_end_full)
                

                
    return full_gene_exons_start, full_gene_exons_end
    
    
    
def coverage_metric_check(exon_list):

    full_coverage_10 = True
    a98_full_coverage_10 = True
    a95_full_coverage_10 = True
    a90_full_coverage_10 = True
    
    full_coverage_5 = True
    a98_full_coverage_5 = True
    a95_full_coverage_5 = True
    a90_full_coverage_5 = True
    
    full_coverage_2 = True
    a98_full_coverage_2 = True
    a95_full_coverage_2 = True
    a90_full_coverage_2 = True
    
    a80_full_coverage_2 = True
    
    smaller_10 = 0
    smaller_5 = 0
    smaller_2 = 0
    
    
    for val in exon_list:

        if val < 10:
            full_coverage_10 = False
            smaller_10 = smaller_10 + 1
        if val < 5:
            full_coverage_5 = False
            smaller_5 = smaller_5 + 1
        if val < 2:
            full_coverage_2 = False
            smaller_2 = smaller_2 + 1



    if smaller_10/len(exon_list) >= 0.02:
         a98_full_coverage_10 = False
    if smaller_10/len(exon_list) >= 0.05:
         a95_full_coverage_10 = False
    if smaller_10/len(exon_list) >= 0.10:
         a90_full_coverage_10 = False
    
    
    if smaller_5/len(exon_list) >= 0.02:
         a98_full_coverage_5 = False
    if smaller_5/len(exon_list) >= 0.05:
         a95_full_coverage_5 = False
    if smaller_5/len(exon_list) >= 0.10:
         a90_full_coverage_5 = False
    
    
    if smaller_2/len(exon_list) >= 0.02:
         a98_full_coverage_2 = False
    if smaller_2/len(exon_list) >= 0.05:
         a95_full_coverage_2 = False
    if smaller_2/len(exon_list) >= 0.10:
         a90_full_coverage_2 = False
    if smaller_2/len(exon_list) >= 0.99:
         a80_full_coverage_2 = False
    
    return a95_full_coverage_5, a90_full_coverage_5, a95_full_coverage_2, a90_full_coverage_2, a80_full_coverage_2
    
    
    
    
    
    
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



def save_gene_fasta(chr_list, gene_start, gene_end):



    with open("./dedup/Galaxy293-[Extract_Genomic_DNA_on_data_291].fasta") as handle:
        for record in SeqIO.parse(handle, "fasta"):

            s1 = str(record.id).split(":")
            chr_ = s1[0].split("chr")[1]
            start_ = s1[1].split("-")[0]
            end_ = s1[1].split("-")[1].split("(")[0]
            
            if chr_ in chr_list:

                for enum, s in enumerate(gene_start):
                
                    if start_ in s:
                
                        if end_ in gene_end[enum]:
                        
                            seq_ = record.seq
            
            

    return
    
    
    
def calculate_constraint(dms_values_in, sequences, threshold, ex_num, strand):
     


    dms_values =dms_values_in.copy()
    sequence = sequences[ex_num]
    

    
    
    
    if strand == "-":
    


        dms_values.reverse()



    constraint = ""
    


    
    for en, value in enumerate(dms_values):
    
        if value > threshold and (sequence[en].upper() == "A" or sequence[en].upper() == "C"):
            constraint = constraint + "x" 
        else:
            constraint = constraint + "."
            

    
    
    
    return constraint




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="extract final DMS profiles")
    
    # Define arguments
    parser.add_argument("--coverage_list", type=str, required=True, help="Path to gff file")
    # Parse arguments
    args = parser.parse_args()


### read in genes, chromosom and 5'UTR regions


    #files_ = open("./coverage_gene_lists_275_x20_5UTR/treated_refseq/treated_genes_99_5.txt")
    files_ = open("./coverage_full/all.txt")
    transcript_id_list = []
    chr_list = []
    starts = []
    ends = []
    
    for line in files_:

    
    
        mod_line = list(filter(None, line.split("\t")))[:-1]
        transcript_id_list.append(mod_line[0])
        chr_list.append(mod_line[2])
        start_ends = mod_line[3:]

        
        start_val = []
        end_val = []
        
        for start_end in start_ends:

            start_val.append(int(start_end.split("-")[0]))
            end_val.append(int(start_end.split("-")[1]))
        
        starts.append([start_val])
        ends.append([end_val])
        

    
    files_.close()


    


    
    version_counter = 2
    prev = "None"
    for en, x in enumerate(transcript_id_list):
        if x == prev:
            transcript_id_list[en] = x + "_v" + str(version_counter)
            version_counter = version_counter + 1

        else:
            version_counter = 2 
        prev = x
                


    ind_ = 1
    

    ############################ read in proteins #########################################################
    
    sequences = []


    #for record in SeqIO.parse("./extracted_seq/Galaxy507-[Extract_Genomic_DNA_on_data_497].fasta", "fasta"):
    for record in SeqIO.parse("./coverage_full/Galaxy125-[Extract_Genomic_DNA_on_data_124].fasta", "fasta"):
    #for record in SeqIO.parse("./coverage_3utr/Galaxy738-[Extract_Genomic_DNA_on_data_736].fasta", "fasta"):
    
         sequences.append(str(record.seq))
         

    seq_sub_list = []
    iterator = 0
    for enum, sub in enumerate(starts):
         
         subs = []

         
         lenght_counter = 0
         
         for start_enum, subs in enumerate(sub[0]):
         
             dist = ends[enum][0][start_enum] - subs
             #if dist != 0: lenght_counter = lenght_counter + 1
             lenght_counter = lenght_counter + 1
         

         subs = sequences[iterator:iterator+lenght_counter]
         
         for en, seq in enumerate(subs):
             subs[en] = seq.upper()
         
         
         

         iterator = iterator + lenght_counter
         seq_sub_list.append(subs)
         

    strands = []
    
    
    #file_ = open("./coverage_3utr/dict_for_genes_90_5_fixed_3utr.bed")
    file_ = open("./coverage_full/all_fixed.bed")
    for line in file_:
    
        strands.append(line.split("\t")[-1].split("\n")[0])
    

    file_.close()
    
    strand_sub = []
    iterator = 0

    for sub in starts:
         
         subs = []
         length_ = len(sub[0])
         subs = strands[iterator:iterator+length_]
         iterator = iterator + length_
         strand_sub.append(subs)
         

    flags = ["combined"]
    

    chr_list = chr_list[ind_-1:]
    starts = starts[ind_-1:]
    ends = ends[ind_-1:]
    strand_sub = strand_sub[ind_-1:]
    transcript_id_list = transcript_id_list[ind_-1:]
    seq_sub_list = seq_sub_list[ind_-1:]
    
    
    
    
    run_from = False
    start_en = None

    
    for flag in flags:
    
        iterator = 10
        counter = 0
        values1 = [0,0,0,0,0]
        coverage1 = [[],[],[],[], []]    
        saved_gene_list = []
        saved_gene2_list = []
        values2 = [0,0,0,0,0]
        coverage2 = [[],[],[],[], []]    
        # 28796101
        # 28796174

        


        #bam_file_list = ["./data/DMS_1lb.bam", "./data/DMS.bam", "./data/RNA-NOIL-3.bam"]
        bam_file_list = ["./data/DMS_1lb.bam", "./data/DMS.bam", "./data/RNA-NOIL-3.bam"]
        
        #bam_file_list = ["/media/sven/Intenso/95_20_stuff/data/DMS.bam"]
        #bam_file_list = ["/media/sven/Intenso/95_20_stuff/data/DMS_1lb.bam"]
        bam_file_list = ["/mnt/Zubradt_extract/Galaxy91-[Samtools_view_on_data_89__filtered_alignments].bam"]
        #bam_file_list = ["/media/sven/Intenso/Zubradt_data/Galaxy91-[Samtools_view_on_data_89__filtered_alignments].bam"]
        is_rna_list = [True, True, True]
        #Galaxy91-[Samtools_view_on_data_89__filtered_alignments].bam
        
        cont = True     
        for en, chr_ in enumerate(chr_list):

            normd_windows_for_bam = []
            coverage_bam_list = []
            
            #if transcript_id_list[en] == "NM_001330460.1": 
            
            #	run_from = True
            #	start_en = en
            	# XM_047422382.1
            #if transcript_id_list[en] == "NM_001375406.1": cont = False
            print(en)
            if transcript_id_list[en] == "NM_001278244.1": cont = False
            if transcript_id_list[en] == "XM_017017591.2": cont = True
            if cont == True: continue
        		
            
            #if run_from == False: continue
            #if en > start_en + run_to or en < start_en + run_from_num:  continue
            
            
            for enum, exon_starts in enumerate(starts[en]): 
    
    
                for ex_num, exon_start in enumerate(exon_starts):
    
                    if strand_sub[en][0] == "+":
            
                        starts[en][enum][ex_num] = starts[en][enum][ex_num] -1


                    elif strand_sub[en][0] == "-":
            
                        #ends[en][enum][ex_num] = ends[en][enum][ex_num] + 1
                        starts[en][enum][ex_num] = starts[en][enum][ex_num] -1
                    else:
        
                        raise NotImplementedError
   
            
            for bam_num, bam_file in enumerate(bam_file_list):

                normd_window, full_data_occupied_list, sequence, name_= analysis_normd2(starts[en], ends[en], chr_, strand_sub[en][0] ,transcript_id_list[en],values1, values2,coverage1, coverage2, flag, seq_sub_list[en], bam_file, is_rna = is_rna_list[bam_num], bam_num = bam_num)
                
                
                

