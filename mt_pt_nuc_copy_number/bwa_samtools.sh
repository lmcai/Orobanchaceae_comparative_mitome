#!/bin/bash
#SBATCH -J bwa            # job name
#sbatch -o mt.o%j             # output and error file name (%j expands to slurm jobid)
#sbatch -e mt.e%j             # output and error file name (%j expands to slurm jobid)
#SBATCH -N 1                        # number of nodes requested
#SBATCH -n 16                       # total number of tasks to run in parallel
#SBATCH -p normal              # queue (partition) 
#SBATCH -t 24:00:00                 # run time (hh:mm:ss) 
source activate getorg
gunzip < ../../cleaned_reads/$1/$2_1.fq.gz >../../cleaned_reads/$1.R1.fq
gunzip < ../../cleaned_reads/$1/$3_1.fq.gz >>../../cleaned_reads/$1.R1.fq
gunzip < ../../cleaned_reads/$1/$2_2.fq.gz >../../cleaned_reads/$1.R2.fq
gunzip < ../../cleaned_reads/$1/$3_2.fq.gz >>../../cleaned_reads/$1.R2.fq

bwa index $1.fasta
bwa mem -t16 -C $1.fasta ../../cleaned_reads/$1.R1.fq ../../cleaned_reads/$1.R2.fq | samtools sort -@16 -o $1.sort.bam
bwa mem -t16 -C $1.fasta ../../cleaned_reads/$2 ../../cleaned_reads/$3 | samtools sort -@16 -o $1.sort.bam
samtools depth $1.sort.bam > $1.cov.tsv
