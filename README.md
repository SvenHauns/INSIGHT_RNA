# InsightRNA
This snakemake pipeline processes DMS and OOPS-seq files and integrates information about secondary structure into a a single structural profile.

# OOPS
provides the snakemake pipeline to process OOPS-seq data

# DMS
provides the snakemake pipeline to process DMS based data

# Integration
provides the final integration fo the datasets


## requirements
- refseq dataset
- hg382bit
- ribonanzanet model
- matching eClip dataset

For convenience we also provide a refseq file and hg38.2bit here:
https://drive.google.com/file/d/1FMUAAIMRE5GKaikTGa1R3B7qQsnHWs4T/view?usp=sharing

An example for the eClip dataset can be found here:
https://drive.google.com/file/d/1Qn9BWyIP4R4-jcJx0D-dekL5dFhQYdrm/view?usp=sharing


## Example

```
python main.py --samples {sample-for-DMS} --samples_oops {sample-for-oops} --nc_samples_oops {nc-samples-oops} --cl_samples_oops {cl-samples-oops} --use-conda
```
