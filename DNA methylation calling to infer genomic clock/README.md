# DNA methylation calling to infer genomic clock

## Data processing pipeline 

[bismark_genome.sh](https://github.com/estelleyao0530/Bioinformatics-pipeline/blob/main/DNA%20methylation%20calling%20to%20infer%20genomic%20clock/bismark_genome.sh) : Bismark prepare genome step run only once for all samples

[single_script.sh](https://github.com/estelleyao0530/Bioinformatics-pipeline/blob/main/DNA%20methylation%20calling%20to%20infer%20genomic%20clock/single_script.sh) : fastq dump for single SRR file

[all_file.sh](https://github.com/estelleyao0530/Bioinformatics-pipeline/blob/main/DNA%20methylation%20calling%20to%20infer%20genomic%20clock/all_file.sh) : call single_script.sh for parallel downloads of a list of SRR files

[final_pipeline.sh](https://github.com/estelleyao0530/Bioinformatics-pipeline/blob/main/DNA%20methylation%20calling%20to%20infer%20genomic%20clock/final_pipeline.sh) : scripted functions for all data processing packages required to generate coverage file

**The last scripts must run in the the same directory where Bismark 0.23.0 is downloaded**

## rDNA clock script 
[calc_age.py](https://github.com/estelleyao0530/Bioinformatics-pipeline/blob/main/DNA%20methylation%20calling%20to%20infer%20genomic%20clock/calc_age.py) : takes in coverage file generated from the data pipeline as input and output rDNA estimated age

## RNA-seq analysis 
[rna_seq_analysis.R](https://github.com/estelleyao0530/Bioinformatics-pipeline/blob/main/DNA%20methylation%20calling%20to%20infer%20genomic%20clock/rna_seq_analysis.R): analyze gene expression, remake Bismark M-bias plot
