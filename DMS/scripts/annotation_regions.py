"""GTF transcript/UTR extraction.

Slicing helpers use GTF's 1-based inclusive convention. extract_refseq_utr
converts once to BED0 and orders intervals from transcript 5-prime to 3-prime.
The established 100-nucleotide coding flank is preserved for both UTR modes.
"""


def extract_refseq_utr(gff_path, run_type="full"):

    exon_dict, start_dict, strand_dict, gene_name_dict, chr_dict = {}, {}, {}, {}, {}
    transcripts = []

    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            if "transcript_id" not in line:
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

            if transcript_id not in exon_dict:
                transcripts.append(transcript_id)
                exon_dict[transcript_id] = []
                start_dict[transcript_id] = []
                gene_name_dict[transcript_id] = ""

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

    # GTF and the UTR slicing helpers use 1-based inclusive coordinates.
    # Convert once at the annotation boundary; all downstream data is BED0.
    utr_regions = {
        tid: sorted([[start - 1, end] for start, end in regions],
                    key=lambda exon: exon[0], reverse=(strand_dict.get(tid) == "-"))
        for tid, regions in utr_regions.items()
    }
    return utr_regions, gene_name_dict, strand_dict, chr_dict


def extract_full_exon(transcript_lists, exon_dict, start_dict, strand_dict):
    """Retain all annotated exons, including transcripts without a start codon."""
    return {key: list(exon_dict[key]) for key in transcript_lists if exon_dict[key]}


def extract_3utr_region(transcript_lists, exon_dict, cds_end_dict, strand_dict, coding_flank=100):

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

            remaining = coding_flank
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
                if slices[-1][1] < slices[-1][0]:
                    slices.pop()

            utrs[key] = slices


        elif strand == '-':
            tx = exons_sorted
            cds_end = cds_end_dict[key][0]
            slices = []

            remaining = coding_flank
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
                if slices[0][1] < slices[0][0]:
                    slices.pop(0)

            utrs[key] = slices

        else:
            utrs[key] = []

    return utrs


def extract_utr_region(transcript_lists, exon_dict, start_dict, strand_dict, coding_flank=100):

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
                    need = acc + coding_flank
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
                    need = acc + coding_flank
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
                slices = slices_tx

            slices = sorted(slices, key=lambda x: x[0], reverse=True)
            utrs[key] = slices

        else:
            utrs[key] = []

    return utrs
