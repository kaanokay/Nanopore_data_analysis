#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=kao25@hi.is # for example uname@hi.is
#SBATCH --partition=mimir # request node from a specific partition
#SBATCH --nodes=1 # number of nodes [mimir 66 ve 72 arasında bir node atayacak rastgele]
#SBATCH --ntasks-per-node=32 # 32 cores per node
#SBATCH --mem-per-cpu=3900 # MB RAM per cpu core
#SBATCH --time=14-00:00:00 # run for 4 hours maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

# To able to run SeqKit from conda envrionment first thing that should be done:

. ~/.bashrc

ml purge

# Load Anaconda3
ml load Anaconda3/2022.05

# Activate SeqKit environment in conda
conda activate /hpcdata/Mimir/kao25/SeqKit

# Define how many cores will be used
cores=128

# Define the list of fastq files to be deduplicated
trimmed_fastq_files=(trimmed_Ctr6_1.fq trimmed_Ctr6_2.fq trimmed_Ctr6_3.fq trimmed_Kmt2a_1.fq trimmed_Kmt2a_2.fq trimmed_Kmt2a_3.fq)

# Loop through the list of fastq files
for fq in "${trimmed_fastq_files[@]}"
do
    # Remove duplicates from the trimmed fastq file using seqkit
    seqkit rmdup \
            -j $cores \
            -n \
            -D "${fq%.fq}_duplicate_statistics.csv" \
            -o trimmed_and_duplicate_removed_$fq \
            $fq
done

# ----------------------------------------

# If you get compressed fastq files then perform command below:

# 1. remove duplicates

# Define how many cores will be used
cores=128

# Define the list of fastq files to be deduplicated
trimmed_fastq_files=(trimmed_Ctr6_1.fq trimmed_Ctr6_2.fq trimmed_Ctr6_3.fq trimmed_Kmt2a_1.fq trimmed_Kmt2a_2.fq trimmed_Kmt2a_3.fq)

# Loop through the list of fastq files
for fq in "${trimmed_fastq_files[@]}"
do
    # Remove duplicates from the trimmed fastq file using seqkit
    seqkit rmdup \
            -j $cores \
            -n \
            -D "${fq%.fq}_duplicate_statistics.csv" \
            -o - \
            $fq | gzip > trimmed_and_duplicate_removed_${fq%.fq}.fq.gz
done
