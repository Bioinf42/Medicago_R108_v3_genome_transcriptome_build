markdown

# *Medicago truncatula* R108 genome V3 T2T Nodule Transcriptome Pipeline 

---
## Overview
This pipeline processes unstranded paired-end RNA-seq data to build a transcriptome for *Medicago truncatula* R108 root nodules, also reffered to as *Medicago littoralis*. 
It performs quality control, trimming, alignment to the R108 V3 T2T reference genome, transcript assembly, merging, and comparison to the reference annotation, producing a comprehensive nodule transcriptome for downstream analysis. It provides it in a docker container
to be used on any sytsem. 

---
## Table of Contents
- [Overview](#overview)
- [Prerequisites](#prerequisites)  
- [Repository Contents](#repository-contents)  
- [Directory structures](#directory-structures)
- [Installation](#installation)  
- [Usage](#usage)  
- [Outputs](#outputs)  
- [Transcript Filtering](#transcript-filtering) 
- [License](#license)  
- [Citation](#citation) 

---
## Prerequisites
- System: Linux (e.g., Ubuntu)

### Software included (`env.yaml`)
- 'Conda-forge' : Used to manage and download other software 
- 'Python-3.11' : Custom merge and filtering script
- 'Java (OpenJDK 11) : Required for Trimmomatic
- 'fastqc_v0.12.1' : For quality control 
- 'Trimmomatic-0.39' : For trimming and adaptor removal 
- 'hisat2-2.2.1' : For alignment
- 'samtools-1.18' : For BAM conversion and sorting
- 'stringtie-3.0.0' : For transcript assembly
- 'gffread-0.12.7' : For extracting transcript sequence
- 'gffcompare-0.11.2' : To compare transcriptome gtf to original genome gtf

### Input Files:
- Paired-end RNA FASTQ files (e.g., `SAMPLE_R1_001.fq.gz`, `SAMPLE_R2_001.fq.gz`)
- Reference genome - "R108_T2T.v3.0.fa"
(https://figshare.com/ndownloader/files/56630867)
- Reference annotation - "R108_T2T.v3.0.gff"
(https://figshare.com/ndownloader/articles/29665022/versions/3?folder_path=Medicago_gene_anno_newVersion)

### Hardware Recommendations 
- Multi-core CPU (12+ threads recommended)
- RAM at least 32GB
- Sufficient disk space for FASTQ and output files.

---
## Repository Contents
- `run_docker.sh`: The docker container launcher script.

- `transcriptome_build_pipe.sh`: Main end-to-end pipeline script, makes directories and runs tools (trimming, alignment, assembly, merging, comparison). 

- `env.yaml` : packages and dependencies. 

- `Dockerfile` : instructions docker uses to build the image and runs micromamba. 

---

## Directory structures

The repository has the following structure upon download:

```bash
Base/
├── docker/
│   ├── work/
│   ├── Dockerfile
│   └── env.yaml
├── filtering_scripts/
│   ├── filter_tx_tpm_matrix.py 
│   └── make_tx_tpm_matrix.py
├── raw_data/
├── reference/
│    ├── R108_T2T.v3.0_fa.zip
│    └── R108_T2T.v3.0_gtf.zip
├── run_docker.sh
├── transcriptome_build_pipe.sh
├── README.md        
└── LICENSE   
```

The repository after running the pipeline:

```bash
Base/
├── docker/
│   ├── work/
│   ├── Dockerfile
│   └── env.yaml
├── filtering_scripts/
│   ├── filter_tx_tpm_matrix.py 
│   └── make_tx_tpm_matrix.py
├── raw_data/
├── reference/
│    ├── R108_T2T.v3.0.fa
│    ├── R108_T2T.v3.0.gtf
│    ├── R108_T2T.v3.0.exon.txt
│    └── R108_T2T.v3.0.splicesite.txt
├── run_docker.sh
├── transcriptome_build_pipe.sh
├── README.md  
├── LICENSE     
├── logs
└── results 
    ├── fastqc_untrimmed/
    ├── trimmed/
    ├── fastqc_trimmed/
    ├── index/
    ├── aligned_bam/
    ├── stringtie_counts
    ├── stringtie_transcripts/
    ├── stringtie_merged/
    ├── gtf_compare/
    ├── stringtie_filtered/
    └── final_gtf/
        └── gff_final_compare/
```


---
## Installation
Follow these steps to install all required tools and to set up your environment before running the pipeline.

1. **Clone the repository**
    
```bash
git clone https://github.com/yourusername/R108_Nodule_tx_Git.git
```

2. **Navigate to the base directory of the repository and make the scripts executable**
     
```bash
chmod +x run_docker.sh
chmod +x rtranscriptome_build_pipe.sh
```

3. **Unzip the Genome `.fa` and annotation `.gtf` files**

These are provided in the reference directory. Or download from the links provided above, using gffread tool`.gtf` was produced from `.gff3` file provided by paper authors.
Only a gtf file will work with this pipeline.  



1. **Run the respository**

Stay in the base directory and run the two commands seperately. Check your system and run the second command with the appropriate number of threads. Reccomend to use 1-2 threads less than 
the system has. 
    
```bash
docker build --no-cache -t txpipe:0.3 -f docker/Dockerfile .
STEP=all THREADS=14 IMAGE=txpipe:0.3 ./run_docker.sh

```

---
## Outputs

There will be 11 output folders in results directory produced for a successful run.

1. **fastqc_untrimmed**: HTML & `.zip` reports from the raw FASTQ files (one pair of reports per library). Open the HTML file for each sample in a web browser and eyeball: 
per-base quality, over-represented sequences, adapter contamination. A quick primer is on the [FastQC project page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

2. **trimmed**: four FASTQ files per sample, all gzip-compressed: 

- `*_R1_paired_001.fq.gz` forward reads that kept their mate 
- `*_R1_unpaired_001.fq.gz` forward reads whose mate was discarded
- `*_R2_paired_001.fq.gz` reverse reads that kept their mate    
- `*_R2_unpaired_001.fq.gz` reverse reads whose mate was discarded

  Paired vs. unpaired? Trimmomatic processes each read independently; if trimming pushes one mate below the `MINLEN` threshold (50 nt here) that read is dropped while its partner is kept. 
  ‘*Paired*’ files therefore contain read pairs that both survived; ‘*unpaired*’ files hold the few orphan reads that survived alone. Down-stream aligners only use the *paired* files.


3. **fastqc_trimmed**: HTML & `.zip` reports of the trimmed data, sequence quality should be higher and adaptors should now be removed. 

4. **index**: Contains the HISAT2 genome index files (*.ht2) built from the reference FASTA (and splice/exon information from the GTF). These files are required for fast RNA-seq alignment
and are reused across runs.

5. **aligned_bam**:  One `.sorted.bam` and `.sorted.bam.bai` per sample, this is a compressed `.sam` file converted into binary coordinate-sorted and indexed. These are the files that are used by StringTie; you can view them in IGV.
6. **stringtie_counts**: Per-sample abundance tables (`*_counts.tab`). Expression estimates (FPKM, TPM, read counts) for every transcript feature in the matching GTF
7. **stringtie_transcript**: a GTF is produced per sample. Each file lists all transcript models assembled from that sample alone.
8. **stringtie_merged**: individual GTF files were merged producing a initial transcriptome outputting the file `R108_merged_nodule.gtf`.
9. **gtf_compare**: Has the output of the comparison of the generated transcriptome annotation `R108_merged_nodule.gtf` to the original genome annoation `R108_T2T.v3.0.fa`
10. **stringtie_filtered**: BAM files were realigned using the generated transcriptome annoation `R108_merged_nodule.gtf`. Python scripts provided create a matrix from the expresson values for each sample and then produce `pass_tx.txt` with genes that pass the low expression filter. 
11. **final_gtf**: The generated transcriptome annotation `R108_merged_nodule.gtf` has its low level transcripts and isoforms removed generating `R108_merged_nodule_filtered.gtf`. It also contains the folder `gff_final_compare/` which compares this final annoation to the original genome annoation. 

```
