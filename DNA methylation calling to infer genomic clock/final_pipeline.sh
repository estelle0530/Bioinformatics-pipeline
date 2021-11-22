#!/bin/bash

#SBATCH -n 12 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 20:00:00 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=10G # Memory in GB (see also --mem-per-cpu)
#SBATCH -o pipeline.out # Standard out goes to this file
#SBATCH -e pipeline.err # Standard err goes to this file

# for every SRR in the list of SRRs file
# for srr in $(cat SRR_list.txt)
# do
# # call the bash script that does the fastq dump, passing it the SRR number next in file
# sbatch single_script.sh $srr
# done

cat *_1.fastq > RRBS_WT_rep1_1.fastq
cat *_2.fastq > RRBS_WT_rep1_2.fastq

module load python/2.7.14-fasrc01
module load cutadapt/1.8.1-fasrc01
module load fastqc/0.11.8-fasrc01
module load TrimGalore/0.5.0-fasrc01

trim_galore --paired --fastqc RRBS_WT_rep1_1.fastq RRBS_WT_rep1_2.fastq

module load bowtie2/2.3.2-fasrc02
module load samtools/1.10-fasrc01

perl ./Bismark-0.23.0/bismark --genome ./rrbs-fastq -1 RRBS_WT_rep1_1_val_1.fq -2 RRBS_WT_rep1_2_val_2.fq
perl ./Bismark-0.23.0/bismark_methylation_extractor -p --no_overlap --CX --bedGraph --gzip --cytosine_report --genome_folder ./rrbs-fastq RRBS_WT_rep1_1_val_1_bismark_bt2_pe.bam
