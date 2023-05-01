#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kao25@hi.is # for example uname@hi.is
#SBATCH --partition=mimir  # request node from a specific partition
#SBATCH --nodes=1                 # number of nodes
#SBATCH --ntasks-per-node=20      # 48 cores per node (96 in total)
#SBATCH --mem-per-cpu=3900        # MB RAM per cpu core
#SBATCH --time=14-00:00:00         # run for 4 hours maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

ml load minimap2/2.24-GCCcore-11.2.0
ml load SAMtools/1.15-GCC-11.2.0

minimap2 -k17 -ax map-ont --secondary=yes -t 128 -y /hpcdata/Mimir/kao25/Juan_Nanopore_data/UCSC_Mus_musculus.GRCm38_genome/UCSC_Mus_musculus.GRCm38_genome.fa.gz Chd1_low_quality_reads.fq | samtools view -Sb - --threads 128 | samtools sort - --threads 128 > Chd1_low_quality_reads_seconday_alignments.bam
minimap2 -k17 -ax map-ont --secondary=yes -t 128 -y /hpcdata/Mimir/kao25/Juan_Nanopore_data/UCSC_Mus_musculus.GRCm38_genome/UCSC_Mus_musculus.GRCm38_genome.fa.gz Chd1_high_quality_reads.fq | samtools view -Sb - --threads 128 | samtools sort - --threads 128 > Chd1_high_quality_reads_seconday_alignments.bam
exit
