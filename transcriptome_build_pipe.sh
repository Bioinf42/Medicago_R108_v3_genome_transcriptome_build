#!/usr/bin/env bash
set -euo pipefail

#Set the number of threads. Can be altered.
THREADS="${THREADS:-14}"

#Set the directories used in the pipeline. 
ROOT="${ROOT:-/work}"
RAW="$ROOT/raw_data"
REF="$ROOT/reference"
RES="$ROOT/results"
FQU="$ROOT/results/fastq_untrimmed"
FQT="$ROOT/results/fastq_trimmed"
TRM="$ROOT/results/trimmed"
IDX="$ROOT/results/index"
ALB="$ROOT/results/aligned_bam"
STT="$ROOT/results/stringtie_transcripts"
STC="$ROOT/results/stringtie_counts"
STF="$ROOT/results/stringtie_filtered"
STFGTF="$ROOT/results/stringtie_filtered/expr_gtf"
STFTAB="$ROOT/results/stringtie_filtered/expr_tab"
STM="$ROOT/results/stringtie_merged"
COM="$ROOT/results/gtf_compare"
FGTF="$ROOT/results/final_gtf"
FGTFFC="$ROOT/results/final_gtf/gff_final_compare"
LOG="$ROOT/logs"



# Create required folders (safe to run even if they already exist)
mkdir -p \
  "$RAW" "$REF" "$RES" "$LOG" \
  "$FQU" "$FQT" "$TRM" "$IDX" "$ALB" \
  "$STT" "$STC" "$STF" "$STM" "$COM" "$STFGTF" \
  "$STFTAB" "$FGTF" "$FGTFFC"

STEP="${STEP:-all}"

if [[ "$STEP" == "all" || "$STEP" == "fastqc_untrimmed" ]]; then
  shopt -s nullglob
  #Run FastQC on all *.fq.gz files in $RAW in parallel (up to $THREADS at once)
  printf '%s\0' "$RAW"/*.fq.gz | xargs -0 -r -n 1 -P "$THREADS" fastqc -o "$FQU"
else
   echo "Skipping FastQC (untrimmed) step (STEP=$STEP)"
fi

if [[ "$STEP" == "all" || "$STEP" == "trim" ]]; then
  shopt -s nullglob
  #Run trimmomatic on untrimmed .fq files
  ADAPTERS="${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE-2.fa"

  for r1 in "$RAW"/*_R1_001.fq.gz; do
    sample="$(basename "$r1" _R1_001.fq.gz)"
    r2="$RAW/${sample}_R2_001.fq.gz"
   
    trimmomatic PE -threads "$THREADS" -phred33 \
      "$r1" "$r2" \
      "${TRM}/${sample}_R1_paired_001.fq.gz" \
      "${TRM}/${sample}_R1_unpaired_001.fq.gz" \
      "${TRM}/${sample}_R2_paired_001.fq.gz" \
      "${TRM}/${sample}_R2_unpaired_001.fq.gz" \
      "ILLUMINACLIP:${ADAPTERS}:2:30:10:2:True" \
      SLIDINGWINDOW:4:20 MINLEN:50
  done
else
   echo "Skipping Trimmomatic trimming step (STEP=$STEP)"
fi

if [[ "$STEP" == "all" || "$STEP" == "fastqc_trimmed" ]]; then
shopt -s nullglob
  #Run FastQC on all trimmed *.fq.gz files in $RAW in parallel (up to $THREADS at once)
  printf '%s\0' "$TRM"/*.fq.gz | xargs -0 -r -n 1 -P "$THREADS" fastqc -o "$FQT"
else
   echo "Skipping FastQC (trimmed) step (STEP=$STEP)"
fi


if [[ "$STEP" != "all" && "$STEP" != "hisat_index" ]]; then
  echo "Skipping HISAT index build step (STEP=$STEP)"
else
  #Make exon and splice files
  hisat2_extract_exons.py "${REF}"/R108_T2T.v3.0.gtf > "${REF}"/R108_V3_exon.txt
  hisat2_extract_splice_sites.py "${REF}"/R108_T2T.v3.0.gtf > "${REF}"/R108_V3_splicesite.txt

  #Make the hisat index using the exon and splice site files. 
  hisat2-build -p "$THREADS" \
    --ss "${REF}"/R108_V3_splicesite.txt \
    --exon "${REF}"/R108_V3_exon.txt \
    "${REF}"/R108_T2T.v3.0.fa \
    "${IDX}"/R108_V3_index

fi


if [[ "$STEP" != "all" && "$STEP" != "hisat_align" ]]; then
  echo "Skipping HISAT align step (STEP=$STEP)"
else
#Run the hisat alignment, use dta setting for splice awareness. Also pipe 
#directly into samtools to sort and index, this prevents SAM file buildup which 
#takes up a lot of memory.
  for r1 in "$TRM"/*_R1_paired_001.fq.gz; do 
      sample="$(basename "$r1" _R1_paired_001.fq.gz)"
      r2="$TRM/${sample}_R2_paired_001.fq.gz"

      hisat2 -p "$THREADS" --dta \
      -x "${IDX}"/R108_V3_index \
      -1 "$r1" \
      -2 "$r2" \
    | samtools sort -@ "$THREADS" -o "$ALB/${sample}.sorted.bam"

    samtools index "$ALB/${sample}.sorted.bam"

  done

fi


if [[ "$STEP" != "all" && "$STEP" != "stringtie_assemble" ]]; then
  echo "Skipping StringTie assemble step (STEP=$STEP)"
else
#Run stringtie with strict filters to remove false positives
  for r1 in "$ALB"/*.sorted.bam; do 
      sample="$(basename "$r1" .sorted.bam)"
      stringtie "$r1" -p "$THREADS" \
        -G "${REF}"/R108_T2T.v3.0.gtf \
        -o "${STT}"/"${sample}_transcripts.gtf" \
        -A "${STC}"/"${sample}_counts.tab" \
        -c 3 -j 5 -f 0.20 -m 300
  done

fi

if [[ "$STEP" != "all" && "$STEP" != "stringtie_merge" ]]; then
  echo "Skipping StringTie merge step (STEP=$STEP)"
else

#Merge created gtf files to create a single transcriptome gtf file. 
  stringtie --merge -p "$THREADS" -l R108v3_nod \
      -G "${REF}"/R108_T2T.v3.0.gtf \
      -o "${STM}/R108_merged_nodule.gtf" \
      "${STT}"/*.gtf
fi


if [[ "$STEP" != "all" && "$STEP" != "gff_compare" ]]; then
  echo "Skipping gffcompare (merged vs reference) step (STEP=$STEP)"
else
#Compare the gtf file to the original 
   gffcompare -r "${REF}"/R108_T2T.v3.0.gtf \
     -o "${COM}/compare" \
     "${STM}/R108_merged_nodule.gtf"

fi


if [[ "$STEP" != "all" && "$STEP" != "filtering_counts_1" ]]; then
  echo "Skipping StringTie quantify-for-filtering step (STEP=$STEP)"
else
# For each aligned BAM, run StringTie in estimate mode (-e) against the merged annotation
# to produce per-sample expression GTFs and abundance tables.
  for r1 in "$ALB"/*.sorted.bam; do 
        sample="$(basename "$r1" .sorted.bam)"
        stringtie "$r1" -p "$THREADS" -e -B \
           -G "${STM}"/R108_merged_nodule.gtf \
           -o "${STFGTF}"/"${sample}_transcripts.gtf" \
           -A "${STFTAB}"/"${sample}_counts.tab" 
  done

fi


if [[ "$STEP" != "all" && "$STEP" != "python_filtering_1" ]]; then
  echo "Skipping TPM matrix build step (STEP=$STEP)"
else
#Merge the gtf files produced using the make_tx_tpm_matrix.py python script.
   cd "${STF}"
   python /work/filtering_scripts/make_tx_tpm_matrix.py

fi

if [[ "$STEP" != "all" && "$STEP" != "python_filtering_2" ]]; then
  echo "Skipping TPM filter step (STEP=$STEP)"
else
#filter the matrix to remove low expression transcripts
   cd "${STF}"
   python /work/filtering_scripts/filter_tx_tpm_matrix.py

fi


if [[ "$STEP" != "all" && "$STEP" != "final_filtering" ]]; then
  echo "Skipping final GTF filtering step (STEP=$STEP)"
else
#use the pass list to filter the earlier merged gtf file. 
   gffread "${STM}/R108_merged_nodule.gtf" \
           --ids "${STF}/pass_tx.txt" \
           -T \
           -o "${FGTF}/R108_merged_nodule_filtered.gtf"

fi


if [[ "$STEP" != "all" && "$STEP" != "final_compare" ]]; then
  echo "Skipping final gffcompare step (STEP=$STEP)"
else
#Compare the final gtf file to original gtf file
   gffcompare -r "${REF}"/R108_T2T.v3.0.gtf \
              -o "${FGTFFC}/gff_compare" \
              "${FGTF}/R108_merged_nodule_filtered.gtf"

fi