#!/bin/bash

rawreads_folder=/workdir

sample_name="F1"

name_r1="file_R1.fastq.gz"
name_r2="file_R2.fastq.gz"

threads_deborah=5

# -------------------------------------------

blackhole="/blackhole_db"

r1="${rawreads_folder}/${name_r1}"
r2="${rawreads_folder}/${name_r2}"

# Creo la cartella di output per il campione che sto analizzando
output_folder="/outdir/${sample_name}"
mkdir -p $output_folder


# _of sta per "output_folder"

#######################
######## FASTP ########
#######################

fastp_of="${output_folder}/fastp"
mkdir -p $fastp_of

fastp_r1="${fastp_of}/${sample_name}_R1.fastq.gz"
fastp_r2="${fastp_of}/${sample_name}_R2.fastq.gz"


fastp -i ${r1} -I ${r2} -o ${fastp_r1} -O ${fastp_r2} \
    --unpaired1 ${fastp_of}/unpaired_R1.fastq.gz \
    --unpaired2 ${fastp_of}/unpaired_R2.fastq.gz \
    --json ${fastp_of}/report.json \
    --html ${fastp_of}/report.html \
    --failed_out ${fastp_of}/failed.fastq.gz \
    --thread ${threads_deborah} -q 20 u 40


#######################
######## KAIJU ########
#######################

kaiju_of="${output_folder}/kaiju"
mkdir -p $kaiju_of

kaiju_db="${blackhole}/kaiju"

kaiju -v -t ${kaiju_db}/nodes.dmp -f ${kaiju_db}/kaiju_db_refseq.fmi \
    -z ${threads_deborah} \
    -i ${fastp_r1} -j ${fastp_r2} \
    -o ${kaiju_of}
