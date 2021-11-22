#!/bin/bash

#SBATCH -n 12 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 10:00:00 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=10G # Memory in GB (see also --mem-per-cpu)
#SBATCH -o single_read.out # Standard out goes to this file
#SBATCH -e single_read.err # Standard err goes to this file

module load sratoolkit/2.8.0-fasrc01 

# for single end reads only
fastq-dump --split-3 $1 

# for paired end reads only
# fastq-dump --split-3  $1