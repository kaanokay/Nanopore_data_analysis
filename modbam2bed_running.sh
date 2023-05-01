#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=kao25@hi.is # for example uname@hi.is
#SBATCH --partition=mimir # request node from a specific partition
#SBATCH --nodes=1 # number of nodes
#SBATCH --ntasks-per-node=32 # 32 cores per node
#SBATCH --mem-per-cpu=3900 # MB RAM per cpu core
#SBATCH --time=14-00:00:00 # run for 14 days maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

# Load modules

module load SAMtools/1.15-GCC-11.2.0
module load modbam2bed/0.6.2

# Run modbam2bed

modbam2bed --aggregate -e -m 5mC --cpg -p Crebbp_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/UCSC_Mus_musculus.GRCm38_genome/UCSC_Mus_musculus.GRCm38_genome.fa.gz Crebbp_merged_and_aligned_sorted.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Ctr3_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/UCSC_Mus_musculus.GRCm38_genome/UCSC_Mus_musculus.GRCm38_genome.fa.gz Ctr3_merged_and_aligned_sorted.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Dnmt1_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/UCSC_Mus_musculus.GRCm38_genome/UCSC_Mus_musculus.GRCm38_genome.fa.gz Dnmt1_merged_and_aligned_sorted.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Ezh2_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/UCSC_Mus_musculus.GRCm38_genome/UCSC_Mus_musculus.GRCm38_genome.fa.gz Ezh2_merged_and_aligned_sorted.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Hdac6_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/UCSC_Mus_musculus.GRCm38_genome/UCSC_Mus_musculus.GRCm38_genome.fa.gz Hdac6_merged_and_aligned_sorted.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Kdm6a_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/UCSC_Mus_musculus.GRCm38_genome/UCSC_Mus_musculus.GRCm38_genome.fa.gz Kdm6a_merged_and_aligned_sorted.bam
exit
