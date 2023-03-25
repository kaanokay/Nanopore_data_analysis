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

# Run modbam2bed from conda envrionment

modbam2bed --aggregate -e -m 5mC --cpg -p Ctr6_1_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/Mus_musculus.GRCm38_genome/Mus_musculus.GRCm38_refence_genome_without_unconventional_chromosomes/Mus_musculus.GRCm38_reference_genome.fa.bgz Ctr6_1_aligned_modified_and_mapping_QC_filtered.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Ctr6_2_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/Mus_musculus.GRCm38_genome/Mus_musculus.GRCm38_refence_genome_without_unconventional_chromosomes/Mus_musculus.GRCm38_reference_genome.fa.bgz Ctr6_2_aligned_modified_and_mapping_QC_filtered.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Ctr6_3_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/Mus_musculus.GRCm38_genome/Mus_musculus.GRCm38_refence_genome_without_unconventional_chromosomes/Mus_musculus.GRCm38_reference_genome.fa.bgz Ctr6_3_aligned_modified_and_mapping_QC_filtered.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Kmt2a_1_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/Mus_musculus.GRCm38_genome/Mus_musculus.GRCm38_refence_genome_without_unconventional_chromosomes/Mus_musculus.GRCm38_reference_genome.fa.bgz Kmt2a_1_aligned_modified_and_mapping_QC_filtered.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Kmt2a_2_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/Mus_musculus.GRCm38_genome/Mus_musculus.GRCm38_refence_genome_without_unconventional_chromosomes/Mus_musculus.GRCm38_reference_genome.fa.bgz Kmt2a_2_aligned_modified_and_mapping_QC_filtered.bam
modbam2bed --aggregate -e -m 5mC --cpg -p Kmt2a_3_CpGs --threads 128 -c /hpcdata/Mimir/kao25/Juan_Nanopore_data/Mus_musculus.GRCm38_genome/Mus_musculus.GRCm38_refence_genome_without_unconventional_chromosomes/Mus_musculus.GRCm38_reference_genome.fa.bgz Kmt2a_3_aligned_modified_and_mapping_QC_filtered.bam
exit
