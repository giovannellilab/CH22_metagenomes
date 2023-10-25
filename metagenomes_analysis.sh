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


#######################
####### ASSEMBLY ######
#######################


metaspades_of="${output_folder}/metaspades"
mkdir -p $metaspades_of


spades.py --meta --only-assembler -t ${threads_deborah}  -1 ${fastp_r1} -2 ${fastp_r2} -o ${metaspades_of}
    

#######################
##### READMAPPING #####
#######################


readmapping_of="${output_folder}/readmapping"
mkdir -p "$readmapping_of/stats"

sam_file="${readmapping_of}/read_mapping.sam"
bam_file="${readmapping_of}/read_mapping.bam"
sorted_bam="${readmapping_of}/read_mapping_sorted.bam"
indexed_bam="${readmapping_of}/read_mapping_sorted.bam.bai"

bbmap.sh ${threads_deborah} ref=${metaspades_of}/contigs.fasta \
 in=${fastp_r1} in2=${fastp_r2} out=${sam_file} covstats=${readmapping_of}/covstats.tsv nodisk \
 bhist=${readmapping_of}/stats/base_composition_by_pos_hist.txt \
 qhist=${readmapping_of}/stats/quality_by_pos_hist.txt \
 aqhist=${readmapping_of}/stats/average_read_quality_hist.txt \
 lhist=${readmapping_of}/stats/read_length_hist.txt \
 ihist=${readmapping_of}/stats/inserted_size_hist.txt \
 ehist=${readmapping_of}/stats/error_per_read_hist.txt \
 qahist=${readmapping_of}/stats/quality_accuracy_hist.txt \
 indelhist=${readmapping_of}/stats/indel_length_hist.txt \
 mhist=${readmapping_of}/stats/match_sub_del_ins_rates_by_location.txt \
 gchist=${readmapping_of}/stats/gc_content_hist.txt \
 idhist=${readmapping_of}/stats/read_count_vs_perc_identity_hist.txt \
 scafstats=${readmapping_of}/stats/reads_mapped_to_scaffold.txt

samtools view -S -b ${sam_file} > ${bam_file}
samtools sort ${bam_file} -o ${sorted_bam}
samtools index ${sorted_bam} ${indexed_bam}
jgi_summarize_bam_contig_depths --outputDepth ${readmapping_of}/metabat2_depth.txt ${sorted_bam}
awk '{print $1\"\t\"$2}' ${readmapping_of}/covstats.tsv | grep -v '^#' > ${readmapping_of}/maxbin2_abundance.txt



#######################
#### MULTIBINNERS #####
#######################  

multibinners_of="${output_folder}/multibinners"

concoct_folder="${multibinners_of}/concoct"
maxbin_folder="${multibinners_of}/maxbin2"
metabat_folder="${multibinners_of}/metabat2"

log_concoct="${concoct_folder}/log.out"
log_maxbin="${maxbin_folder}/log.out"
log_metabat="${metabat_folder}/log.out"

mkdir -p ${concoct_folder} ${maxbin_folder} ${metabat_folder}

mkdir -p "${concoct_folder}/bins"
cut_up_fasta.py ${metaspades_of}/contigs.fasta -c 10000 -o 0 --merge_last -b ${concoct_folder}/contigs_10K.bed > ${concoct_folder}/contigs_10K.fa
concoct_coverage_table.py ${concoct_folder}/contigs_10K.bed ${sorted_bam} > ${concoct_folder}/coverage_table.tsv
concoct --composition_file ${concoct_folder}/contigs_10K.fa --coverage_file ${concoct_folder}/coverage_table.tsv -b {concoct_folder}/ >> ${log_concoct} 2>&1
merge_cutup_clustering.py ${concoct_folder}/clustering_gt1000.csv > ${concoct_folder}/clustering_merged.csv
extract_fasta_bins.py ${metaspades_of}/contigs.fasta ${concoct_folder}/clustering_merged.csv --output_path ${concoct_folder}/bins


 
mkdir -p ${maxbin_folder}/bins
run_MaxBin.pl -contig ${metaspades_of}/contigs.fasta -abund ${readmapping_of}/maxbin2_abundance.txt -thread ${threads_deborah} -out ${maxbin_folder}/output >> ${log_maxbin} 2>&1
mv ${maxbin_folder}/*.fasta {maxbin_folder}/bins/



mkdir -p "${metabat_folder}/bins"
metabat2 --inFile ${metaspades_of}/contigs.fasta --abdFile ${readmapping_of}/metabat2_depth.txt -o ${metabat_folder}/output >> ${log_metabat} 2>&1
mv ${metabat_folder}/*.fa ${metabat_folder}/bins/


#######################
###### DASTOOL ########
#######################


dastool_of="${output_folder}/dastool"
mkdir -p "$dastool_of"


Fasta_to_Contig2Bin.sh -i ${multibinners_of}/concoct/bins -e fasta > ${dastool_of}/concoct_dastool.tsv
Fasta_to_Contig2Bin.sh -i ${multibinners_of}/maxbin2/bins -e fasta > ${dastool_of}/maxbin2_dastool.tsv
Fasta_to_Contig2Bin.sh -i ${multibinners_of}/metabat2/bins -e fasta > ${dastool_of}/metabat2_dastool.tsv
-i {output.folder}/concoct_dastool.tsv,{output.folder}/maxbin2_dastool.tsv,${dastool_of}/metabat2_dastool.tsv \
            -l concoct,maxbin2,metabat2 \
            -c {input.contig_path}/contigs.fasta \
            -o {output.folder}/das_tool \
            --write_bin_evals \
            --write_bins \
            --threads ${threads_deborah} \
            --score_threshold 0.6 \
            --search_engine diamond
        cd ${dastool_of} && mv das_tool_DASTool_bins bins





