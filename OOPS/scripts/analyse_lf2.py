import glob
import pandas as pd
import numpy as np
import subprocess as sp
#import matplotlib.pyplot as plt
#import seaborn as sns
from sklearn.mixture import GaussianMixture
import numpy as np
import argparse
import warnings
warnings.filterwarnings("ignore")

def sort_lists_by_two_indices(list1, index1, index2):
    # Get the sorted indices based on index1 and then index2
    sorted_indices = sorted(range(len(list1)),
                            key=lambda i: (list1[i][index1], int(list1[i][index2])))
    
    # Sort both lists using the same indices
    #sorted_list1 = [list1[i] for i in sorted_indices]
    #reordered_list2 = [list2[i] for i in sorted_indices]
    
    return sorted_indices

def write_table_data(data, name):
    
    df = pd.DataFrame(data, columns = ["chr", "window_start", "window_end", "transcript_id", "gene", "baseMean", "padj", "lf"])
    df.to_csv(name)

    return



def get_isoform_expression(tids, path = "/media/sven/Intenso/projects/24/oops/giulia_output/output_full_files_split_il/"):
    regions = []
    mean_list = []
    
    for tid in tids:
        if tid == "":continue
        try:
            df = pd.read_csv(path + tid + "_stat_results.csv")

            region = df['Unnamed: 0']
            mean_ = df['baseMean']
            for enum, reg in enumerate(region):
                if reg not in regions:
                    regions.append(reg)
                    mean_list.append(mean_[enum])
        except:
            continue
    
    if len(mean_list) < 10: return 0
    q1 = np.quantile(mean_list, 0.2)

    mean_list = np.mean(mean_list)*0.3
    

    return mean_list



def collect_isoforms(path):
    
    isoform_dict = {}
    tid_dict = {}
    
    file_ = open(path)
    direction_dict = {}
    
    for line in file_:
    
        
        if len(line.split("\t")) < 5: continue
        
        gene_id = line.split('gene_id "')[1].split('";')[0]
        t_id = line.split('transcript_id "')[1].split('";')[0]
        dir_ = line.split("\t")[6]

        
        if gene_id not in isoform_dict.keys():
            isoform_dict[gene_id] = []
        if t_id not in isoform_dict[gene_id]: 
            isoform_dict[gene_id].append(t_id)
            direction_dict[gene_id] = dir_

        if t_id not in tid_dict.keys():
            tid_dict[t_id] = gene_id
            
    file_.close()

    return isoform_dict, tid_dict, direction_dict

def gmm_outlier_detection(baseline, test_vector, threshold=0.01):
    """
    Identifies outliers in test_vector based on the Gaussian Mixture Model (GMM)
    using the baseline vector.

    Parameters:
    - baseline (array-like): Baseline vector to fit the GMM.
    - test_vector (array-like): Vector to test for outliers.
    - threshold (float): Probability threshold below which values are considered outliers.

    Returns:
    - outliers (list): Values in test_vector that are outliers.
    """
    baseline = np.array(baseline).reshape(-1, 1)
    test_vector = np.array(test_vector).reshape(-1, 1)

    # Fit a Gaussian Mixture Model with one or two components
    gmm = GaussianMixture(n_components=1, random_state=42)
    gmm.fit(baseline)

    # Compute probability density of the test vector
    log_probs = gmm.score_samples(test_vector)  # Log-likelihood of test samples
    probs = np.exp(log_probs)  # Convert log likelihood to probability

    # Find outliers based on low probability threshold
    outliers = test_vector[probs < threshold].flatten()

    return outliers, probs < threshold
    
def modified_z_score_outliers(baseline, test_vector, threshold=2):
    """
    Identifies outliers in test_vector based on the modified z-score method 
    using the baseline vector.

    Parameters:
    - baseline (array-like): Baseline vector to compute the median and MAD.
    - test_vector (array-like): Vector to test for outliers.
    - threshold (float): Modified z-score threshold (default = 3.5).

    Returns:
    - outliers (list): Values in test_vector that are outliers.
    """
    baseline = np.array(baseline)
    test_vector = np.array(test_vector)

    median = np.median(baseline)
    mad = np.median(np.abs(baseline - median))  # Median Absolute Deviation

    modified_z_scores = 0.6745 * (test_vector - median) / (mad + 1e-9)  # Avoid div by zero
    outliers = test_vector[np.abs(modified_z_scores) > threshold]
    
    

    return outliers, np.abs(modified_z_scores) > threshold
    
    
    
def plot_two_violins(data1, data2, labels, title="violin plot", xlabel="positive/negative lf", ylabel="lf"):
    """
    Creates a violin plot with two violins.

    Parameters:
        data1 (list or array): Data for the first violin.
        data2 (list or array): Data for the second violin.
        labels (tuple): Labels for the two violins (e.g., ("Category 1", "Category 2")).
        title (str): Title of the plot.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
    """
    if len(labels) != 2:
        raise ValueError("Labels must be a tuple with two elements.")
    
    # Prepare data for seaborn
    data = pd.DataFrame({
        "correlation": data1 + data2,
        "finetuning": [labels[0]] * len(data1) + [labels[1]] * len(data2)
    })

    # Create the violin plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(x="finetuning", y="correlation", data=data)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(axis='y')
    plt.show()




def get_exons(file_):

    file_ = open(file_).readlines()
    
    dict_ =  {}
    exon_dict = {}
    
    for line_ in file_:
    
        id_ = line_.split("\t")[-1][1:-3]
        id_ = line_.split("\t")[3]

        
        if id_ not in dict_.keys():
            exon_dict[id_] = []
            dict_[id_] = line_.split("\t")[0]

        exon_dict[id_].append(line_.split("\t")[1] + "-" + line_.split("\t")[2])


    return dict_, exon_dict

def calculate_coverage_ratio(unique_regions, region, mean, print_,id_):

    unique_coverage_dict= {}
    base_coverage = []
    
    if print_: print(region)
    if print_: print(mean)
    
    if print_: print(id_)
    
    for ex in unique_regions:
        unique_coverage_dict[ex] = []

    for enum, reg in enumerate(region):
        found = False
        for ex in unique_regions:
            if int(ex.split("-")[1]) - int(ex.split("-")[0]) < 10: continue

            
            if reg >= int(ex.split("-")[0]) and reg +40 <= int(ex.split("-")[1]):
                found = True
                if ex not in unique_coverage_dict.keys(): unique_coverage_dict[ex] = []
                unique_coverage_dict[ex].append(mean[enum])
                
            elif reg < int(ex.split("-")[0]) and reg +40 >= int(ex.split("-")[1]):
                found = True
                if ex not in unique_coverage_dict.keys(): unique_coverage_dict[ex] = []
                unique_coverage_dict[ex].append(mean[enum])
        else:
            base_coverage.append(mean[enum])
    
    if print_:print("INFO")
    if print_:print(unique_regions)
    if print_:print(unique_coverage_dict)
    if print_:print(base_coverage)
    #base_mean = np.mean(base_coverage)
    base_mean = np.quantile(base_coverage, 0.3)
    print("base mean")
    if print_:print(base_mean)
    invalid = False
    
    for key in unique_coverage_dict.keys():
        if print_: print("comparsion")
        if print_: print(key)
        if print_: print(np.mean(unique_coverage_dict[key]))
        if print_: print(base_mean )
        
        if np.mean(unique_coverage_dict[key]) < base_mean:
            if print_:print("basemean to large")
            if print_:print(key)
            if print_:print(np.mean(unique_coverage_dict[key]))
            if print_:print(base_mean )
            invalid = True
            
        if len(unique_coverage_dict[key]) == 0 and (int(key.split("-")[1]) - int(key.split("-")[0]) < 50): 
            if print_:print("empty list")
            if print_:print(key)
            if print_:print(unique_coverage_dict[key])
            invalid = True
            
            
    if print_: print(invalid)
    return invalid

def calculate_isoform_coverage(iso_dict, tid_dict, exon_dict):
    

    
    unique_regions = {}
    
    counter = 0


    for key in iso_dict.keys():
         key_val = key
         if key_val == "GAPDH": print(key)

         tids = iso_dict[key]
         if key_val == "GAPDH": print(tids)
         
         
         for key in tids:
            if key == "": continue
            if key not in exon_dict.keys(): continue
            if key_val == "GAPDH":print(key)
            if key_val == "GAPDH":print(exon_dict[key])
            counter = counter + 1

            #if counter > 10000: continue
    
            for ex in exon_dict[key]:
                if key_val == "GAPDH":print(ex)
                overlap_exists = False
                for key2 in tids:
                    if key2 == "": continue
                    if key2 not in exon_dict.keys(): continue
                    for ex2 in exon_dict[key2]:
                        if key2 == key: continue

                        if int(ex.split("-")[1]) < int(ex2.split("-")[1]) and int(ex.split("-")[1]) > int(ex2.split("-")[0]):
                    
                            unique_region = str(ex.split("-")[1])+ "-" + str(ex2.split("-")[1])
                            
                            if key2 not in unique_regions.keys(): unique_regions[key2] = []
                            if key == "NM_001256799.3": print(ex)
                            if key == "NM_001256799.3": print(ex2)
                            if key == "NM_001256799.3": print(unique_region)
                            unique_regions[key2].append(unique_region)
                            overlap_exists = True
                            
                        if int(ex.split("-")[0]) < int(ex2.split("-")[0]) and int(ex.split("-")[1]) <= int(ex2.split("-")[1]) and int(ex.split("-")[1]) > int(ex2.split("-")[0]):

                            unique_region = str(ex.split("-")[0])+ "-" + str(ex2.split("-")[0])
                            if key not in unique_regions.keys(): unique_regions[key] = []
                            if key == "NM_001256799.3": print(ex)
                            if key == "NM_001256799.3": print(ex2)
                            if key == "NM_001256799.3": print(unique_region)
                            unique_regions[key].append(unique_region)
                            overlap_exists = True
                            
                        if int(ex.split("-")[0]) == int(ex2.split("-")[0]) and int(ex.split("-")[1]) == int(ex2.split("-")[1]):
                        
                        
                            overlap_exists = True
                       ###### completely unique exon
                if key_val == "GAPDH":print("overlap_exists")
                if key_val == "GAPDH": print(overlap_exists)
                if overlap_exists == False:
                    if key not in unique_regions.keys(): unique_regions[key] = []
                    unique_regions[key].append(ex)
                
            if key_val == "GAPDH":print("here")
            if key_val == "GAPDH":print(unique_regions[key])
    
    return unique_regions
    
    
if __name__ == "__main__":


    cmdline_parser = argparse.ArgumentParser('creates track file for OOPS-seq')

    cmdline_parser.add_argument('-f', '--deseq_output',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-z', '--track_output',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-c', '--track_output_p',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-d', '--track_output_p2',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-e', '--track_output_n',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-i', '--track_table',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-j', '--track_table_p',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-k', '--track_table_p_z_score',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-l', '--track_table_p_gm',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-a', '--positive_tracks',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-r', '--refseq_file',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-b', '--exons',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-y', '--track_output2',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-m', '--file_trackb2',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-n', '--file_track2_full',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-x', '--file_track_fp2',
                                default="",
                                help='output html',
                                required = True,
                                type=str)
                                
    args, unknowns = cmdline_parser.parse_known_args()
    
    



    files = glob.glob(args.deseq_output + "/*")
    file_a = open(args.track_output, "w")
    file_a2 = open(args.track_output2, "w")
    file_track = open(args.track_output_p, "w")
    file_trackb = open(args.track_output_p2, "w")
    file_track2 = open(args.track_output_n, "w")
    file_track2b = open(args.track_output_n, "w")
    file_track_fp = open(args.positive_tracks,"w")
    file_trackb2 = open(args.file_trackb2,"w")
    file_track2_full = open(args.file_track2_full,"w")
    file_track_fp2 = open(args.file_track_fp2,"w")
    
    
    # file_trackb2 file_track2_full file_track_fp2
    
    
    
    header = 'track name=window description="window" color=0,255,0,\n#chrom chromStart chromEnd\n'

    write_list = []
    regions_save = []
    regions_save_full = []

    iso_dict, tid_dict, direction_dict = collect_isoforms(args.refseq_file)
    print(args.exons)
    chr_dict, exon_dict = get_exons(args.exons)


    unique_regions = calculate_isoform_coverage(iso_dict, tid_dict, exon_dict)#



    pos_lf = []
    neg_lf = []
    pos_lf_save = []
    neg_lf_save = []
    pos_lf_save_data = []
    pos_lf_save_full = []
    pos_lf_exons = []
    neg_lf_exons = []
    all_exons = []
    region_save_data = []
    invalid_coutner = 0
    valid_coutner = 0

    failure_count = 0
    for enum, file_ in enumerate(files):
    
    
        print(file_)
        try:
            df = pd.read_csv(file_)
        except:
            failure_count = failure_count + 1
            continue

        region = df['Unnamed: 0']
        padj = df['padj']
        lf = df['log2FoldChange']
        mean = df['baseMean']
    
        id_ = file_.split("/")[-1].split("_")[0] +"_" + file_.split("/")[-1].split("_")[1]
        id_ = id_.split("|")[0]
        

    
        if id_ not in tid_dict.keys(): continue
        gene = tid_dict[id_]
        all_tids = iso_dict[gene]
        direction = direction_dict[gene]
        

    
        if tid_dict[id_] == "GAPDH": print_ = True
        else: print_ = False


        q1 = get_isoform_expression(all_tids)

        if id_ not in unique_regions.keys(): continue
    
    
    
        #try:
        chr_ = chr_dict[id_]
        exons = exon_dict[id_]
        exons = sorted(exons, key=lambda e: int(e.split("-")[0]))
        invalid = calculate_coverage_ratio(unique_regions[id_], region, mean, print_, id_)
        if invalid: invalid_coutner = invalid_coutner + 1
        else: valid_coutner = valid_coutner + 1
        
        if invalid: continue
        
        #except:
        #    continue
        
        exon_cov_dict = {}

        for enum, reg in enumerate(region):
    
            if padj[enum] > 0 and padj[enum] < 0.05 and mean[enum] > 20:
                found = False

                correct_exon = None
            
                for enum, ex in enumerate(exons): 


                    if reg >= int(ex.split("-")[0]) and reg <= int(ex.split("-")[1]):
                        found = True
                        correct_exon = int(ex.split("-")[1])
                        if enum < len(exons): next_exon = exons[enum+1:]
                        full_correct_exon = ex
                        if ex not in exon_cov_dict.keys(): exon_cov_dict[ex] = []
                        exon_cov_dict[ex].append(mean[enum])
                    
        
                assert found == True
                
                save_tmp_list = []
                save_full_tmp = []
                save_data_tmp = []
            
                if reg+40 <= correct_exon:
                    save = chr_ + "\t" + str(reg-1) + "\t" + str(min(reg+40, correct_exon)) + "\n"
                    save_full = chr_ + "\t" + str(reg) + "\t" + str(min(reg+40, correct_exon)) + "\t"  + id_ + "_" + str(enum) + "\t" + "bla" + "\t" + direction +"\n"
                    save_data = [chr_, reg, min(reg+40, correct_exon), id_, gene, mean[enum], padj[enum], lf[enum]]
                    
                    save_tmp_list.append(save)
                    save_full_tmp.append(save_full)
                    save_data_tmp.append(save_data)
                    
                    
                    
                else: 
                    save = chr_ + "\t" + str(reg-1) + "\t" + str(min(reg+40, correct_exon)) + "\n"
                    save_full = chr_ + "\t" + str(reg) + "\t" + str(min(reg+40, correct_exon)) + "\t"  + id_ + "_" + str(enum) + "\t" + "bla" + "\t" + direction +"\n"
                    save_data = [chr_, reg, min(reg+40, correct_exon), id_, gene, mean[enum], padj[enum], lf[enum]]
                    rest_ = reg+40 - correct_exon
                    
                    save_tmp_list.append(save)
                    save_full_tmp.append(save_full)
                    save_data_tmp.append(save_data)
                    
                    count_exons = 0
                    while rest_ > 0:
                        if count_exons >= len(next_exon): break
                        current_exon_save = next_exon[count_exons]
                        current_exon_save_start = int(current_exon_save.split("-")[0])
                        current_exon_save_end = int(current_exon_save.split("-")[1])
                        rest_ = current_exon_save_start+rest_ - current_exon_save_end
                        
                        save = chr_ + "\t" + str(current_exon_save_start-1) + "\t" + str(min(current_exon_save_start+rest_, current_exon_save_end)) + "\n"
                        save_full = chr_ + "\t" + str(current_exon_save_start) + "\t" + str(min(current_exon_save_start+rest_, current_exon_save_end)) + "\t"  + id_ + "_" + str(enum) + "\t" + "bla" + "\t" + direction +"\n"
                        save_data = [chr_, reg, min(current_exon_save_start+rest_, current_exon_save_end), id_, gene, mean[enum], padj[enum], lf[enum]]
                        
                        
                        save_tmp_list.append(save)
                        save_full_tmp.append(save_full)
                        save_data_tmp.append(save_data)
                        count_exons = count_exons + 1
                        
                        
                        
                for save_count, save in enumerate(save_tmp_list): 
    
                    if save not in regions_save:
                        regions_save.append(save)
                        regions_save_full.append(save_full_tmp[save_count])
                        region_save_data.append(save_data_tmp[save_count])
                        all_exons.append(full_correct_exon)
                        
                        
                        """
                        here
                        """
                        if lf[enum] >0: 
                            pos_lf.append(lf[enum])
                            pos_lf_save.append(save)
                            pos_lf_save_data.append(save_data_tmp[save_count])
                            pos_lf_save_full.append(save_full_tmp[save_count])
                            pos_lf_exons.append(full_correct_exon)
                    
                    
                        if lf[enum] <0: 
                            neg_lf.append(lf[enum])
                            neg_lf_save.append(save)
                            neg_lf_exons.append(full_correct_exon)

        """
        for exon in exon_cov_dict.keys():
            if np.mean(exon_cov_dict[exon]) < q1:
        
                regions_save = [p for enum, p in enumerate(regions_save) if all_exons[enum] != exon]
                regions_save_full =[p for enum, p in enumerate(regions_save_full) if all_exons[enum] != exon]
                region_save_data = [p for enum, p in enumerate(region_save_data) if all_exons[enum] != exon]
        
                pos_lf = [p for enum, p in enumerate(pos_lf) if pos_lf_exons[enum] != exon]
                pos_lf_save = [p for enum, p in enumerate(pos_lf_save) if pos_lf_exons[enum] != exon]
                pos_lf_save_full = [p for enum, p in enumerate(pos_lf_save_full) if pos_lf_exons[enum] != exon]
                pos_lf_save_data =[p for enum, p in enumerate(pos_lf_save_data) if pos_lf_exons[enum] != exon]
                neg_lf = [p for enum, p in enumerate(neg_lf) if neg_lf_exons[enum] != exon]
                neg_lf_save = [p for enum, p in enumerate(neg_lf_save) if neg_lf_exons[enum] != exon]
            
                pos_lf_exons = [p for enum, p in enumerate(pos_lf_exons) if pos_lf_exons[enum] != exon]
                neg_lf_exons = [p for enum, p in enumerate(neg_lf_exons) if neg_lf_exons[enum] != exon]
    
        """
    print(len(regions_save))
    print(valid_coutner)
    print(invalid_coutner)




    """
    unify peaks
    """
    def collapse_regions(pos_lf, table_data=[]):
        print("table_data")
        print(table_data)
        pos_lf_collapsed = []
        current_start = None
        current_end = None
        current_chr = None
        current_tid = ""
        
        for enum, reg in enumerate(pos_lf):
            start = int(reg.split("\t")[1])
            end = int(reg.split("\t")[2]) 
            chr_ = reg.split("\t")[0]
    
            if current_start == None:
    
                current_start = start
                current_end = end
                current_chr = chr_
                print("enum")
                print(enum)
                if table_data != []:print(table_data[enum])
                if table_data != []:print(table_data[enum][3])
                if table_data != []:print(table_data[enum][4])
                if table_data != []: current_tid = table_data[enum][3] + "|" + table_data[enum][4]
    
            if current_end >= start and chr_ == current_chr:
                current_end = end
    
            else:
                pos_lf_collapsed.append(current_chr + "\t" + str(current_start) + "\t" + str(current_end) + "\t" + current_tid + "\n")
        
                current_start = start
                current_end = end
                current_chr = chr_
                if table_data != []: current_tid = table_data[enum][3] + "|" + table_data[enum][4]
        
        if current_chr != None:
            pos_lf_collapsed.append(str(current_chr) + "\t" + str(current_start) + "\t" + str(current_end) + "\t" + current_tid + "\n")    


        return pos_lf_collapsed


    ind = sort_lists_by_two_indices([s.split("\t") for s in regions_save], 0, 1)
    regions_save = [regions_save[i] for i in ind]
    print("region_save_data")
    print(region_save_data)
    region_save_data_sorted = [region_save_data[i] for i in ind]
    regions_save = collapse_regions(regions_save, region_save_data_sorted)


    file_a.write(header)
    file_track.write(header)
    file_track2.write(header)
    file_track_fp.write(header)

    for reg in regions_save:
         file_a.write(reg)

    for reg in regions_save_full:
         file_a2.write(reg)

    print(region_save_data)
    write_table_data(np.array(region_save_data), args.track_table)


    print("positive")
    print(np.mean(pos_lf))
    print(np.quantile(pos_lf, 0.1))
    print(np.quantile(pos_lf, 0.2))
    print(np.quantile(pos_lf, 0.3))
    print(np.quantile(pos_lf, 0.4))
    print(np.quantile(pos_lf, 0.5))
    print(np.quantile(pos_lf, 0.6))
    print(np.quantile(pos_lf, 0.7))
    print(np.quantile(pos_lf, 0.8))
    print(np.quantile(pos_lf, 0.9))

    print("negative")
    print(np.mean(neg_lf))
    print(np.quantile(neg_lf, 0.1))
    print(np.quantile(neg_lf, 0.2))
    print(np.quantile(neg_lf, 0.3))
    print(np.quantile(neg_lf, 0.4))
    print(np.quantile(neg_lf, 0.5))
    print(np.quantile(neg_lf, 0.6))
    print(np.quantile(neg_lf, 0.7))
    print(np.quantile(neg_lf, 0.8))
    print(np.quantile(neg_lf, 0.9))

    ## violiin plot
    print("data")
    print(len(pos_lf_save))
    print(len(neg_lf_save))


    ## save positive full

    ind = sort_lists_by_two_indices([s.split("\t") for s in pos_lf_save], 0, 1)
    pos_lf_save_order = [pos_lf_save[i] for i in ind]
    pos_lf_save_data_sorted = [pos_lf_save_data[i] for i in ind]
    print("here")
    print(pos_lf_save_data)
    pos_lf_save2 = collapse_regions(pos_lf_save_order, pos_lf_save_data_sorted)


    for d in pos_lf_save2:
        file_track_fp.write(d)
    
    for d in pos_lf_save_full:
        file_track_fp2.write(d)
    

    write_table_data(np.array(pos_lf_save_data), args.track_table_p)

    outliers, index = modified_z_score_outliers(np.abs(neg_lf), pos_lf)
    print(len(outliers))
    print(index)
    data1 = np.array(pos_lf_save)[index]
    data1b = np.array(pos_lf_save_full)[index]


    ind = sort_lists_by_two_indices([s.split("\t") for s in data1], 0, 1)
    data1 = [data1[i] for i in ind]

    pos_lf_save_data_sorted = np.array(pos_lf_save_data)[index]
    pos_lf_save_data_sorted = [pos_lf_save_data_sorted[i] for i in ind]


    data1 =  collapse_regions(data1, pos_lf_save_data_sorted)


    write_table_data(np.array(pos_lf_save_data)[index], args.track_table_p_z_score)



    for d in data1:
        file_track.write(d)

    for d in data1b:
        file_track2_full.write(d)


    outliers, index = modified_z_score_outliers(pos_lf, np.abs(neg_lf))
    print(len(outliers))
    print(index)

    data1 = np.array(neg_lf_save)[index]

    ind = sort_lists_by_two_indices([s.split("\t") for s in data1], 0, 1)
    data1 = [data1[i] for i in ind]

    data1 =  collapse_regions(data1)
    for d in data1:
        file_track2.write(d)

    
    outliers, index = gmm_outlier_detection(np.abs(neg_lf), pos_lf)
    print(len(outliers))
    print(index)
    print(pos_lf_save)
    data1 = np.array(pos_lf_save)[index]
    data1b = np.array(pos_lf_save_full)[index]

    ind = sort_lists_by_two_indices([s.split("\t") for s in data1], 0, 1)
    data1 = [data1[i] for i in ind]

    pos_lf_save_data_sorted = np.array(pos_lf_save_data)[index]
    pos_lf_save_data_sorted = [pos_lf_save_data_sorted[i] for i in ind]

    print("data1")
    print(data1)

    data1 =  collapse_regions(data1, pos_lf_save_data_sorted)

    #ind = sort_lists_by_two_indices([s.split("\t") for s in data1], 0, 1)
    #data1 = [data1[i] for i in ind]


    for d in data1:
        file_trackb.write(d)
    
    for d in data1b:
        file_trackb2.write(d)

    
    write_table_data(np.array(pos_lf_save_data)[index], args.track_table_p_gm)

    outliers, index = gmm_outlier_detection(pos_lf, np.abs(neg_lf))
    print(len(outliers))
    print(index)
    data1 = np.array(neg_lf_save)[index]
    print(len(data1))
    data1 =  collapse_regions(data1)

    #ind = sort_lists_by_two_indices([s.split("\t") for s in data1], 0, 1)
    #data1 = [data1[i] for i in ind]


    for d in data1:
        file_track2b.write(d)
        
        
        ## turn off deseq controll
