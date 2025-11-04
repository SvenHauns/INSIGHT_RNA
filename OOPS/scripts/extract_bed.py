import argparse
import numpy as np

refseq_to_chr = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9":  "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9":  "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
    "NC_012920.1":  "chrM"
}


def gtf_to_exon_bed_with_chr(gtf_file, output_bed):
    with open(gtf_file) as infile, open(output_bed, 'w') as out:
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != "exon":
                continue

            raw_chrom = fields[0]
            chrom = refseq_to_chr.get(raw_chrom)
            if chrom is None:
                continue 

            start = int(fields[3]) - 1  # convert to 0-based for BED
            end = int(fields[4])
            strand = fields[6]
            attr_field = fields[8]

            attrs = {}
            for attr in attr_field.strip().split(';'):
                if attr.strip():
                    key, value = attr.strip().split(' ', 1)
                    attrs[key] = value.strip('"')

            transcript_id = attrs.get("transcript_id", "NA")
            gene_name = attrs.get("gene_name", attrs.get("gene_id", "NA"))

            out.write(f"{chrom}\t{start}\t{end}\t{transcript_id}\t{gene_name}\t{strand}\n")




from collections import defaultdict

def generate_windows_from_exons(bed_file, window_size=40, step_size=20):
    exon_dict = defaultdict(list)

    # Read exons grouped by transcript ID
    with open(bed_file) as f:
        for line in f:
            chrom, start, end, tx_id, gene, strand = line.strip().split("\t")
            exon_dict[tx_id].append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "gene": gene,
                "strand": strand
            })

    windows = []

    for tx_id, exons in exon_dict.items():
        # Sort by start position
        exons.sort(key=lambda x: x["start"])
        

        # Intra-exon windows
        start_offsett = []
        start_offsett2 = []
        added_exon = np.inf
        

        
        for exon_num, exon in enumerate(exons):
                        
            pos = exon["start"]
            counter_left = 0
            counter_right = 0
            counter_middle = 0
            
            while pos <=  exon["end"]:

                if pos + window_size <= exon["end"]:
                
                    window_size_offset = window_size
                    
                    add = "junction_left_" + str(exon_num) + "_" + str(counter_left) if start_offsett != [] else "intra"
                    
                    if exon_num > added_exon:
                        add = "junction_left_" + str(exon_num) + "_" + str(counter_left) if start_offsett != [] else "intra"
                        if start_offsett != []: counter_left = counter_left + 1
                    
                        if start_offsett != []:  window_size_offset = window_size - start_offsett[0] 
                        else: window_size_offset = window_size
                        if start_offsett != []: start_offsett.pop(0)
                
                    windows.append((
                    exon["chrom"],
                    pos,
                    pos + window_size_offset,
                    tx_id,
                    exon["gene"],
                    exon["strand"],
                    add
                    ))
                    

                    if pos + window_size_offset - 20 < exon["start"]: pos = exon["start"] 
                    else: pos = pos + window_size_offset - 20

                    #pos += step_size
                    if start_offsett == []:
                        added_exon = np.inf
                        
                elif pos <= exon["end"] and pos + window_size > exon["end"] and exon_num > added_exon:
                
                    #if start_offsett != []: counter_middle = counter_middle + 1
                    if start_offsett != []:  
                        window_size_offset = window_size - start_offsett[0] 
                        current_start_offset = start_offsett[0] 
                    else: window_size_offset = window_size
                    if start_offsett != []: end = min(pos + window_size_offset, exon["end"])
                    else: end = exon["end"]
                    
                    
                    if start_offsett != []: start_offsett.pop(0)
                    
                    if end < exon["end"] or ((end - pos) + current_start_offset) == 40:
                        add = "junction_left_" + str(exon_num) + "_" + str(counter_left)
                        counter_left = counter_left + 1
                    else:
                        add =  "junction_middle_" + str(exon_num) + "_" + str(counter_left) + "_" + str(counter_right)
                        counter_left = counter_left + 1
                        counter_right = counter_right + 1
                        
                    
                    windows.append((
                    exon["chrom"],
                    pos,
                    end,
                    tx_id,
                    exon["gene"],
                    exon["strand"],
                    add
                    ))
                    

                    if (end - pos) + current_start_offset != 40:
                        start_offsett2.append((end - pos) + current_start_offset)

                    if start_offsett == []:
                        start_offsett = start_offsett2
                        added_exon = exon_num
                        start_offsett2 = []
                        
                    if start_offsett == []:
                        added_exon = np.inf
                        
                    
                    assert (end - pos) + current_start_offset <= 40

                    
                    if pos + window_size_offset - 20 < exon["start"]: pos = exon["start"] 
                    else: pos = pos + window_size_offset - 20

                elif pos <= exon["end"] and pos + window_size > exon["end"] and (exon_num > added_exon) == False:
                

                    windows.append((
                    exon["chrom"],
                    pos,
                    exon["end"],
                    tx_id,
                    exon["gene"],
                    exon["strand"],
                    "junction_right_" + str(exon_num) + "_" + str(counter_right)
                    ))

                    
                    counter_right = counter_right + 1
                    start_offsett.append(exon["end"] - pos)
                    pos += step_size
                    added_exon = exon_num

            
    return windows
    
    
def write_windows_to_bed(windows, out_file):
    with open(out_file, "w") as f:
        for w in windows:
            chrom, start, end, tx_id, gene, strand, wtype = w
            name = f"{tx_id}|{gene}|{wtype}"
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")


if __name__ == '__main__':
    cmdline_parser = argparse.ArgumentParser('extract bed file windows')

    cmdline_parser.add_argument('-f', '--input_file',
                                default="",
                                help='input file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-o', '--output_file',
                                default="",
                                help='output file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-w', '--output_window_file',
                                default="",
                                help='output window file',
                                required = True,
                                type=str)

    args, unknowns = cmdline_parser.parse_known_args()
    gtf_to_exon_bed_with_chr(args.input_file, args.output_file)

    windows = generate_windows_from_exons(args.output_file, window_size=40, step_size=20)
    write_windows_to_bed(windows, args.output_window_file)

    
    
    

