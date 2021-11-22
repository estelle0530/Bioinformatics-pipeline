#!/bin/bash

#SBATCH -n 12 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 20:00:00 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=100G # Memory in GB (see also --mem-per-cpu)
#SBATCH -o pipeline.out # Standard out goes to this file
#SBATCH -e pipeline.err # Standard err goes to this file

module load bowtie2/2.3.2-fasrc02
module load samtools/1.10-fasrc01

#perl ./Bismark-0.23.0/bismark_genome_preparation ../mc38/
perl ./Bismark-0.23.0/bismark_genome_preparation ../rdna-edit-lemos.fasta