markdown
# *Medicago truncatula* V3 T2T R108 nodule transcriptome: tutorial & pipeline

---
## Overview
This pipeline processes unstranded paired-end RNA-seq data to build a transcriptome for *Medicago truncatula* R108 root nodules, also referred to as *Medicago littoralis*. 
It performs quality control, trimming, alignment to the R108 V3 T2T reference genome, transcript assembly, merging, and comparison to the reference annotation, producing a comprehensive nodule transcriptome for downstream analysis. A filtering step is also included to remove lowly expressed transcripts. The pipeline is packaged in a Docker image for reproducibility and deployment ease. 

---
## Table of Contents
- [Overview](#overview)
- [Prerequisites](#prerequisites)  
- [Repository Contents](#repository-contents)  
- [Directory structures](#directory-structures)
- [Installation](#installation)  
- [Usage](#usage)  
- [Outputs](#outputs)  
- [License](#license)  
- [Citation](#citation) 

---
## Prerequisites
- System: Linux (e.g., Ubuntu)
- Docker 

### Software included (`env.yaml`)
- `Python-3.11` : Custom merge and filtering script
- `Java (OpenJDK 11)` : Required for Trimmomatic
- `fastqc_v0.12.1` : For quality control 
- `Trimmomatic-0.39` : For trimming and adaptor removal 
- `hisat2-2.2.1` : For alignment
- `samtools-1.18` : For BAM conversion and sorting
- `stringtie-3.0.0` : For transcript assembly
- `gffread-0.12.7` : For extracting transcript sequence
- `gffcompare-0.11.2` : To compare transcriptome gtf to original genome gtf

### Input Files:
- Paired-end RNA FASTQ files (e.g., `SAMPLE_R1_001.fq.gz`, `SAMPLE_R2_001.fq.gz`)
- Reference genome - "R108_T2T.v3.0.fa"
(https://figshare.com/ndownloader/files/56630867)
- Reference annotation - Downloaded "R108_T2T.v3.0.gff" then converted to "R108_T2T.v3.0.gtf" using gffread. The `.gtf` is provided which the pipeline uses.
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
V3_T2T_R108_nodule_transcriptome/
├── docker/
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
V3_T2T_R108_nodule_transcriptome/
├── docker/
│   ├── Dockerfile
│   └── env.yaml
├── filtering_scripts/
│   ├── filter_tx_tpm_matrix.py 
│   └── make_tx_tpm_matrix.py
├── raw_data/
├── reference/
│    ├── R108_T2T.v3.0.fa
│    └──  R108_T2T.v3.0.gtf
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
git clone https://github.com/Bioinf42/V3_T2T_R108_nodule_transcriptome.git
```

2. **Navigate to the base directory of the repository and make the scripts executable**

The directory name `V3_T2T_R108_nodule_transcriptome/`
     
```bash
cd V3_T2T_R108_nodule_transcriptome
chmod +x run_docker.sh
chmod +x transcriptome_build_pipe.sh
```

3. **Download the reference genome and annotation with test RNA data**

Stay in the base directory and use wget command or download from browser. It is also provided in releases on the github repo page.

```bash
wget https://github.com/Bioinf42/V3_T2T_R108_nodule_transcriptome/releases/download/v0.3/Test_input.zip
```

Unzip the `Test_input.zip` 

```bash
unzip Test_input.zip
```

Move into `Test_input`, in the reference directory unzip the `.gtf.zip` and `.fa.zip` files.

```bash
cd Test_input
unzip reference/R108_T2T.v3.0_fa.zip
unzip reference/R108_T2T.v3.0_gtf.zip 
```

Do not change directory. Move the reference and raw_data folders into `V3_T2T_R108_nodule_transcriptome/`

```bash
mv raw_data ../
mv reference ../
```
---

## Usage
1. **Run the repository**
Stay in the base directory `V3_T2T_R108_nodule_transcriptome/` and run the two commands separately. Check your system and run the second command with the appropriate number of threads. Recommend to use 1-2 threads less than 
the system has. 

Build the docker image:
```bash
docker build --no-cache -t txpipe:0.3 -f docker/Dockerfile .
```
Run the pipeline:
```bash
STEP=all THREADS=14 IMAGE=txpipe:0.3 ./run_docker.sh
```

The directory and file with the final output gtf is `results/final_gtf/R108_merged_nodule_filtered.gtf`

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
7. **stringtie_transcripts**: a GTF is produced per sample. Each file lists all transcript models assembled from that sample alone.
8. **stringtie_merged**: individual GTF files were merged producing a initial transcriptome outputting the file `R108_merged_nodule.gtf`.
9. **gtf_compare**: Has the output of the comparison of the generated transcriptome annotation `R108_merged_nodule.gtf` to the original genome annotation `R108_T2T.v3.0.gtf`
10. **stringtie_filtered**: BAM files were reused to estimate expression using the generated transcriptome annotation `R108_merged_nodule.gtf`. Python scripts provided create a matrix from the expression values for each sample and then produce `pass_tx.txt` with genes that pass the low expression filter. 
11. **final_gtf**: The generated transcriptome annotation `R108_merged_nodule.gtf` has its low level transcripts and isoforms removed generating `R108_merged_nodule_filtered.gtf`. It also contains the folder `gff_final_compare/` which compares this final annotation to the original genome annotation. 

## License

- **Code**: MIT License © 2025 <Matthew Jolly>. See `LICENSE`.
- **Text/docs/figures in this repo**: Creative Commons Attribution 4.0 International (**CC BY 4.0**).
- **Third-party tools and external genomes/annotations/data** remain under their own licenses/terms; follow the providers’ terms.

## Citation

**This repository**
- Matthew Jolly (Gifford Lab, University of Warwick). *Medicago truncatula V3 T2T R108 nodule transcriptome: tutorial & pipeline*. GitHub, `v0.1.0` (2025)

**Software**

- 'fastqc_v0.12.1' : Andrews S. (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data
- 'Trimmomatic-0.39' : Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30(15):2114-2120. doi:10.1093/bioinformatics/btu170
- 'hisat2-2.2.1' : Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol. 2019;37(8):907-915. doi:10.1038/s41587-019-0201-4
- 'samtools-1.18' : Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-2079. doi:10.1093/bioinformatics/btp352
- 'stringtie-3.0.0' : Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol. 2015;33(3):290-295. doi:10.1038/nbt.3122
- 'gffread-0.12.7' : Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. F1000Res. 2020;9:ISCB Comm J-304. Published 2020 Apr 28. doi:10.12688/f1000research.23297.2

**Genome and Annotation**

- Shen L, Yi C, Liu Y, Han F, Feng J. Two complete telomere-to-telomere genome assemblies of Medicago reveal the landscape and evolution of its centromeres. Molecular Plant. 2025;18(9):1409-12. https://doi.org/10.1016/j.molp.2025.07.016
