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
    

def analysis_normd2(
    bamfile, starts, ends, strand, chr_, name_, values, values2,
    coverage, coverage2, gene_name, transcript_id, transcript_name,
    collection, flag="complete"
):
    dict_for_genes_99_20, dict_for_genes_95_20, dict_for_genes_90_20, dict_all = collection

    saved_gene  = [0, 0, 0, 0, 0]
    saved_gene2 = [0, 0, 0, 0, 0]

    bam = pysam.AlignmentFile(bamfile)

    b_val,  s_val,  t_val,  f_val,  fi_val  = values
    b_val2, s_val2, t_val2, f_val2, fi_val2 = values2
    values1 = [b_val,  s_val,  t_val,  f_val,  fi_val]
    values2 = [b_val2, s_val2, t_val2, f_val2, fi_val2]

    for enum_start, exon_starts in enumerate(starts):
        exon_ends = ends[enum_start]

        full_data_save = []
        stat_stable_dict  = {"A": 0, "C": 0, "T": 0, "G": 0}
        stat_changed_dict = {
            "A->C": 0, "A->T": 0, "A->G": 0,
            "C->A": 0, "C->T": 0, "C->G": 0,
            "G->A": 0, "G->T": 0, "G->C": 0,
            "T->A": 0, "T->C": 0, "T->G": 0
        }

        for ex_num, exon_start in enumerate(exon_starts):
            exon_end = exon_ends[ex_num]
            if exon_start == exon_end:
                continue

            L = exon_end - exon_start
            occ_total   = np.zeros(L, dtype=np.int32)
            occ_forward = np.zeros(L, dtype=np.int32)
            occ_reverse = np.zeros(L, dtype=np.int32)
            mut_total   = np.zeros(L, dtype=np.int32)
            mut_forward = np.zeros(L, dtype=np.int32)
            mut_reverse = np.zeros(L, dtype=np.int32)

            for read in bam.fetch(chr_, start=exon_start, stop=exon_end):
                pairs = read.get_aligned_pairs(matches_only=True, with_seq=True)
                if not pairs:
                    continue

                gpos = np.fromiter((p[1] for p in pairs), dtype=np.int64, count=len(pairs))
                in_span = (gpos >= exon_start) & (gpos < exon_end)
                if not in_span.any():
                    continue

                idx = (gpos[in_span] - exon_start).astype(np.int64)
                bases = [p[2] for p in pairs]  # list of chars
                is_mut = np.fromiter((b.islower() for b in bases), dtype=bool, count=len(bases))[in_span]

                np.add.at(occ_total, idx, 1)
                if read.is_reverse:
                    np.add.at(occ_reverse, idx, 1)
                else:
                    np.add.at(occ_forward, idx, 1)

                if is_mut.any():
                    idxm = idx[is_mut]
                    np.add.at(mut_total, idxm, 1)
                    if read.is_reverse:
                        np.add.at(mut_reverse, idxm, 1)
                    else:
                        np.add.at(mut_forward, idxm, 1)

                # Per-pair stable/changed stats (depends on strand/ref; keep loop)
                for r in pairs:
                    stat_stable_dict, stat_changed_dict = count_mutation(
                        stat_stable_dict, stat_changed_dict, r, read, exon_end, exon_start
                    )

            cov_fwd_safe = np.maximum(occ_forward, 1)  # avoid div-by-zero
            cov_rev_safe = np.maximum(occ_reverse, 1)
            frac_sum = (mut_forward / cov_fwd_safe) + (mut_reverse / cov_rev_safe)

            cov_total_capped = np.minimum(occ_total, 20)
            if cov_total_capped.size == 0:
                raise NotImplementedError
            full_data_save.extend(cov_total_capped.tolist())

        if not full_data_save:
            continue

        a99_full_coverage_20, a95_full_coverage_20, a90_full_coverage_20 = coverage_metric_check(full_data_save)
        change_a, change_c, change_g, change_t, change_dms = metrics(stat_stable_dict, stat_changed_dict)

        record = [
            [chr_, [exon_starts, ends[enum_start]]],
            gene_name[0],
            transcript_name,
            transcript_id[enum_start][0],
            strand,
            change_a, change_c, change_g, change_t, change_dms
        ]

        if a99_full_coverage_20:
            dict_for_genes_99_20[transcript_id[enum_start]].append(record)
        if a95_full_coverage_20:
            dict_for_genes_95_20[transcript_id[enum_start]].append(record)
        if a90_full_coverage_20:
            dict_for_genes_90_20[transcript_id[enum_start]].append(record)
        dict_all[transcript_id[enum_start]].append(record)

    collection = [dict_for_genes_99_20, dict_for_genes_95_20, dict_for_genes_90_20, dict_all]
    return values1, values2, saved_gene, saved_gene2, collection


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
    pos, gpos, base = r 

    if not (exon_start <= gpos < exon_end):
        return stat_stable_dict, stat_changed_dict

    seq = read.get_forward_sequence()
    if read.is_reverse:
        seq = seq[::-1]

    if not (0 <= pos < len(seq)):
        return stat_stable_dict, stat_changed_dict

    if not base.islower():
        b = seq[pos]
        stat_stable_dict[b] = stat_stable_dict.get(b, 0) + 1
        return stat_stable_dict, stat_changed_dict

    if read.is_reverse:
        comp_map = {"t": "A", "g": "C", "c": "G", "a": "T"}
        complement = comp_map.get(base, base.upper())
    else:
        complement = base.upper()

    obs = seq[pos]
    if complement != "N" and obs != "N":
        key = f"{complement}->{obs}"
        stat_changed_dict[key] = stat_changed_dict.get(key, 0) + 1

    return stat_stable_dict, stat_changed_dict
    
    
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
    thresholds = [10, 5, 2]
    fracs = [0.02, 0.05, 0.10]

    full_coverage = [[] for _ in thresholds]
    two_coverage  = [[] for _ in thresholds]
    five_coverage = [[] for _ in thresholds]
    ten_coverage  = [[] for _ in thresholds]

    for l in avg_list:
        n = len(l)
        for i, T in enumerate(thresholds):
            full_coverage[i].append(1 if all(v >= T for v in l) else 0)

        below_counts = [sum(v < T for v in l) for T in thresholds]
        below_fracs = [c / n for c in below_counts] 
        
        for i, frac in enumerate(below_fracs):
            two_coverage[i].append(0 if frac >= fracs[0] else 1)
            five_coverage[i].append(0 if frac >= fracs[1] else 1)
            ten_coverage[i].append(0 if frac >= fracs[2] else 1)

    return full_coverage, two_coverage, five_coverage, ten_coverage   
    
    
def extract_refseq_utr(gff_path, run_type="full"):

    exon_dict, start_dict, strand_dict, gene_name_dict, chr_dict, has_cds = (
        {}, {}, {}, {}, {}, {}
    )
    transcripts = []
    last_tid = None

    with open(gff_path) as f:
        for i, line in enumerate(f):
            if i <= 4 or line.startswith("#"):
                continue
            if "transcript_id" not in line:
                print(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # skip malformed lines

            chrom, _, feature, start, end, _, strand, _, attrs = fields

            # Extract transcript and gene IDs
            try:
                transcript_id = attrs.split("transcript_id")[1].split(";")[0].strip().strip('"')
            except IndexError:
                continue

            if transcript_id != last_tid:
                transcripts.append(transcript_id)
                exon_dict[transcript_id] = []
                start_dict[transcript_id] = []
                gene_name_dict[transcript_id] = ""
                last_tid = transcript_id

            # Extract gene name if present
            if "gene_id" in attrs:
                gene_name_dict[transcript_id] = attrs.split("gene_id")[1].split(";")[0].strip().strip('"')

            start, end = int(start), int(end)

            # Record exon and strand
            if feature == "exon":
                exon_dict[transcript_id].append([start, end])
                strand_dict[transcript_id] = strand
                chr_dict[transcript_id] = chrom

            # Record start/stop codons based on run_type
            if run_type != "3UTR" and feature == "start_codon":
                start_dict[transcript_id] = [start, end]
            elif run_type == "3UTR" and feature == "stop_codon":
                start_dict[transcript_id] = [start, end]

            # Mark transcripts with coding sequence
            if feature == "CDS":
                has_cds[transcript_id] = True

    # Choose UTR extraction function
    if run_type == "full":
        utr_regions = extract_full_exon(transcripts, exon_dict, start_dict, strand_dict)
    elif run_type == "5UTR":
        utr_regions = extract_utr_region(transcripts, exon_dict, start_dict, strand_dict)
    elif run_type == "3UTR":
        utr_regions = extract_3utr_region(transcripts, exon_dict, start_dict, strand_dict)
    else:
        raise NotImplementedError(f"Unknown run_type: {run_type}")

    # Map RefSeq NC identifiers to chromosome names
    refseq_to_chr = {
        f"NC_{i:06d}": f"chr{i}" for i in range(1, 23)
    } | {"NC_000023": "chrX", "NC_000024": "chrY", "NC_012920": "chrMT"}

    for tid, chrom in chr_dict.items():
        base = chrom.split(".")[0]
        if base in refseq_to_chr:
            chr_dict[tid] = refseq_to_chr[base]

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


    

def extract_3utr_region(transcript_lists, exon_dict, cds_end_dict, strand_dict, additiona_nuc = 100):

    utrs = {}

    for key in transcript_lists:
        exons = exon_dict.get(key, [])
        if not exons or key not in cds_end_dict or not cds_end_dict[key]:
            utrs[key] = []
            continue

        strand = strand_dict[key]
        exons_sorted = sorted([list(e) for e in exons], key=lambda x: x[0])

        if strand == '+':
            tx_exons = exons_sorted
            cds_end = cds_end_dict[key][1]
            slices = []

            remaining = additiona_nuc
            for s, e in reversed(tx_exons):
                if e < cds_end:  
                    continue
                take_right = min(e, cds_end)
                left_bound = max(s, take_right - remaining + 1)
                if left_bound <= take_right and remaining > 0:
                    slices.append([left_bound, take_right])
                    remaining -= (take_right - left_bound + 1)
                if remaining <= 0:
                    break

            if remaining > 0:
                for s, e in reversed(tx_exons):
                    if e >= cds_end: 
                        continue
                    take = min(e - s + 1, remaining)
                    if take > 0:
                        slices.append([e - take + 1, e])
                        remaining -= take
                    if remaining <= 0:
                        break

            for s, e in tx_exons:
                if e <= cds_end:
                    continue
                left = max(s, cds_end + 1)
                if left <= e:
                    slices.append([left, e])

            slices.sort(key=lambda x: x[0])

            if slices:
                slices[-1][1] -= 0
                if slices[-1][1] < slices[-1][0]:
                    slices.pop()

            utrs[key] = slices


        elif strand == '-':
            tx = exons_sorted 
            cds_end = cds_end_dict[key][0] 
            slices = []

            remaining = additiona_nuc
            for i, (s, e) in enumerate(tx):
                if s <= cds_end <= e:
                    right = min(e, cds_end + remaining - 1)
                    slices.append([cds_end, right])
                    remaining -= (right - cds_end + 1)
                    j = i + 1
                    while remaining > 0 and j < len(tx):
                        s2, e2 = tx[j]
                        take = min(e2 - s2 + 1, remaining)
                        slices.append([s2, s2 + take - 1])
                        remaining -= take
                        j += 1
                    break

            for s, e in tx:
                if e < cds_end:
                    slices.append([s, e])
                elif s < cds_end:
                    slices.append([s, cds_end - 1])

            slices.sort(key=lambda x: x[0])

            if slices:
                slices[0][1] += 0
                if slices[0][1] < slices[0][0]:
                    slices.pop(0)

            utrs[key] = slices

        else:
            utrs[key] = []

    return utrs
    
    
def extract_utr_region(transcript_lists, exon_dict, start_dict, strand_dict, additiona_nuc = 100):

    utrs = {}

    for key in transcript_lists:
        exons = exon_dict.get(key, [])
        if not exons or not start_dict.get(key):
            continue

        strand = strand_dict[key]
        exons_sorted = sorted([list(e) for e in exons], key=lambda x: x[0])

        if strand == '+':
            tx_exons = exons_sorted
            cds_start = start_dict[key][0]

            acc = 0
            slices = []
            done = False

            for s, e in tx_exons:
                if e < cds_start:
                    acc += (e - s + 1)
                elif s <= cds_start <= e:
                    utr5_last_len = max(0, (cds_start - 1) - s + 1)
                    acc += utr5_last_len
                    need = acc + additiona_nuc
                    remaining = need
                    for s2, e2 in tx_exons:
                        if remaining <= 0: break
                        seg_len = e2 - s2 + 1
                        take = min(seg_len, remaining)
                        slices.append([s2, s2 + take - 1])
                        remaining -= take
                    done = True
                    break
                else:
                    break

            if not done:
                slices = []

            if slices:
                slices[-1][1] -= 0
                if slices[-1][1] < slices[-1][0]:
                    slices.pop()

            utrs[key] = slices

        elif strand == '-':
            tx_exons = list(reversed(exons_sorted))
            cds_start = start_dict[key][1]

            acc = 0
            slices_tx = []
            found = False

            for s, e in tx_exons:
                if s > cds_start:
                    acc += (e - s + 1)
                elif s <= cds_start <= e:
                    utr5_last_len = max(0, e - (cds_start + 1) + 1)
                    acc += utr5_last_len
                    need = acc + additiona_nuc
                    remaining = need
                    for s2, e2 in tx_exons:
                        if remaining <= 0: break
                        seg_len = e2 - s2 + 1
                        take = min(seg_len, remaining)
                        slices_tx.append([e2 - take + 1, e2])
                        remaining -= take
                    found = True
                    break
                else:
                    break

            if not found:
                slices = []
            else:
                slices = sorted(slices_tx, key=lambda x: x[0])

            if slices:
                slices[0][1] -= 1
                if slices[0][1] < slices[0][0]:
                    slices.pop(0)
            slices = sorted(slices, key=lambda x: x[0], reverse = True)
            utrs[key] = slices

        else:
            utrs[key] = []

    return utrs


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="loads potential DMS candidates")
    
    # Define arguments
    parser.add_argument("--gff_path", type=str, required=True, help="Path to gff file")
    parser.add_argument("--bamfile", type=str, required=True, help="Path to sample bam file")
    parser.add_argument("--target_folder", type=str, required=True, help="Path to output folder")
    parser.add_argument("--run_type", type=str, required=True, help="run type")
    
    args = parser.parse_args()



    # 1l beta
    utr_regions, gene_name_dict, strand_dict, chr_dict = etract_refseq_utr(args.gff_path, args.run_type)

    return_chr_type = {"NC_000001":"chr1","NC_000002":"chr2","NC_000003":"chr3","NC_000004":"chr4","NC_000005":"chr5","NC_000006":"chr6","NC_000007":"chr7","NC_000008":"chr8","NC_000009":"chr9"
    ,"NC_000010":"chr10","NC_000011":"chr11","NC_000012":"chr12","NC_000013":"chr13","NC_000014":"chr14","NC_000015":"chr15","NC_000016":"chr16","NC_000017":"chr17","NC_000018":"chr18","NC_000019":"chr19",
    "NC_000020":"chr20","NC_000021":"chr21","NC_000022":"chr22","NC_000023":"chrX","NC_000024":"chrY","NC_012920":"chrMT" }
    
    
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
   
    
    
