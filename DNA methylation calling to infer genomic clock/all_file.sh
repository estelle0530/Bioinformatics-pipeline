#!/bin/bash

#SBATCH -n 12 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 10:00:00 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=10G # Memory in GB (see also --mem-per-cpu)
#SBATCH -o srr_read.out # Standard out goes to this file
#SBATCH -e srr_read.err # Standard err goes to this file

# for every SRR in the list of SRRs file
for srr in $(cat SRR_list.txt)
do
# call the bash script that does the fastq dump, passing it the SRR number next in file
sbatch single_script.sh $srr
done