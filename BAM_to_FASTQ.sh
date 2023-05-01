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

module load SAMtools/1.15-GCC-11.2.0

# Important note: To save storage size, we need to obtain fastq files as fq.gz extension instead of .fq extension.
# So, for this, command below can be used:

samtools fastq -TMM,ML --threads 128 Chd1_merged.bam | gzip > Chd1_merged.fq.gz
exit

# MM and ML tags for methylation are kept in fastq file after conversion of BAM file to fastq file.

# Its also possible to get uncompressed fastq files through command below:

# samtools fastq -TMM,ML --threads 128 Chd1_merged.bam > Chd1_merged.fq

