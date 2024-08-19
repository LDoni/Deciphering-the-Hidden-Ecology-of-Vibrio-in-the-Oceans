# Deciphering-the-Hidden-Ecology-of-Vibrio-in-the-Oceans

# Metagenomic Data Processing Pipeline

This pipeline processes metagenomic data by performing the following steps: trimming reads, classifying them with Kraken2 using the RefSeq database, extracting Vibrio reads, classifying these reads with Enterobase, and finally refining the classifications using Bracken.

## Steps Overview

1. **Trimming**: Pre-process the raw reads to remove low-quality bases and adapters.
2. **Kraken2 Classification with RefSeq**: Classify the trimmed reads at the genus level using the RefSeq database.
3. **Bracken Genus-Level Refinement (RefSeq)**: Refine the genus-level classifications using Bracken.
4. **Extract Vibrio Reads (RefSeq)**: Extract reads classified as Vibrio from the RefSeq output.
5. **Kraken2 Classification with Enterobase**: Classify the extracted Vibrio reads at the species level using the Enterobase database.
6. **Bracken Species-Level Refinement (Enterobase)**: Refine the species-level classifications using Bracken.
7. **Extract Vibrio Reads (Enterobase)**: Extract Vibrio reads from the Enterobase output.

## Script

```bash
#!/bin/bash

# Set the current path variable
currpath=$(pwd)

# Set the sample name file
SAMPLE="samples_TARA.txt"

# RefSeq Kraken2 database directory
REFSEQ_KRA_DB="/home/LAPO/krakenDB_bacteria_refseq"

# Enterobase Kraken2 database directory
ENT_KRA_DB="/home/LAPO/Desktop/TARA/enterobase_vibrio_KRAKEN2_db"

# Output directories
TRIM_OUT="${currpath}/trimming"
REFSEQ_OUT="${TRIM_OUT}/kraken_refseq"
ENT_OUT="${TRIM_OUT}/kraken_entero/entero_extraction"

#######################
# 1. Trimming of reads
#######################

cat $SAMPLE | parallel -j 10 'trim_galore \
--cores 5 --paired {}_1.fastq.gz {}_2.fastq.gz \
--trim-n --illumina --fastqc_args "-outdir ${TRIM_OUT}/quality_control" -o ${TRIM_OUT}'

###############################################################
# 2. Kraken2 classification with RefSeq database after trimming
###############################################################

cat $SAMPLE | parallel -j 10 \
"kraken2 --db ${REFSEQ_KRA_DB} \
--paired ${TRIM_OUT}/{}_1_val_1.fq.gz ${TRIM_OUT}/{}_2_val_2.fq.gz \
--threads 20 \
--use-names \
--gzip-compressed \
--report ${REFSEQ_OUT}/{}_report-kraken_refseq.txt \
--output ${REFSEQ_OUT}/{}_output_refseq.kraken"

#########################################################
# 3. Bracken genus-level refinement after Kraken2 (RefSeq)
#########################################################

cat $SAMPLE | parallel -j 10 \
"bracken -d ${REFSEQ_KRA_DB} \
-i ${REFSEQ_OUT}/{}_report-kraken_refseq.txt \
-o ${REFSEQ_OUT}/{}_bracken_genus_refseq.txt \
-l G"

######################################################################
# 4. Extract reads classified as Vibrio using RefSeq classified reads
######################################################################

cat $SAMPLE | parallel -j 10 \
"python3 ~/extract_kraken_reads.py \
-k ${REFSEQ_OUT}/{}_output_refseq.kraken \
-s1 ${TRIM_OUT}/{}_1_val_1.fq.gz \
-s2 ${TRIM_OUT}/{}_2_val_2.fq.gz \
--taxid 662 --include-children \
-o ${REFSEQ_OUT}/{}_extracted_reads_vibrio_1.fa \
-o2 ${REFSEQ_OUT}/{}_extracted_reads_vibrio_2.fa"

#####################################################################################
# 5. Kraken2 classification with Enterobase database using Vibrio extracted reads
#####################################################################################

cat $SAMPLE | parallel -j 10 \
"kraken2 --db ${ENT_KRA_DB} \
--paired ${REFSEQ_OUT}/{}_extracted_reads_vibrio_1.fa ${REFSEQ_OUT}/{}_extracted_reads_vibrio_2.fa \
--threads 20 \
--use-names \
--report ${ENT_OUT}/{}_report-kraken_entero.txt \
--output ${ENT_OUT}/{}_output_entero.kraken"

###############################################################
# 6. Bracken species-level refinement after Kraken2 (Enterobase)
###############################################################

cat $SAMPLE | parallel -j 10 \
"bracken -d ${ENT_KRA_DB} \
-i ${ENT_OUT}/{}_report-kraken_entero.txt \
-o ${ENT_OUT}/{}_bracken_species_entero.txt \
-l S"

####################################################################
# 7. Extract reads classified as Vibrio using Enterobase classified reads
####################################################################

cat $SAMPLE | parallel -j 10 \
"python3 ~/extract_kraken_reads.py \
-k ${ENT_OUT}/{}_output_entero.kraken \
-s1 ${REFSEQ_OUT}/{}_extracted_reads_vibrio_1.fa \
-s2 ${REFSEQ_OUT}/{}_extracted_reads_vibrio_2.fa \
--taxid 662 --include-children \
-o ${ENT_OUT}/{}_extracted_reads_vibrio_entero_1.fa \
-o2 ${ENT_OUT}/{}_extracted_reads_vibrio_entero_2.fa"
