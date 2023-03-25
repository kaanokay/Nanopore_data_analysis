#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kao25@hi.is # for example uname@hi.is
#SBATCH --partition=mimir # request node from a specific partition
#SBATCH --nodes=1 # number of nodes [mimir 66 ve 72 arasında bir node atayacak rastgele]
#SBATCH --ntasks-per-node=32 # 32 cores per node
#SBATCH --mem-per-cpu=3900 # MB RAM per cpu core
#SBATCH --time=0-01:00:00 # run for 4 hours maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes
lscpu
echo $HOSTNAME
lsmem
module load HISAT2/2.2.1-gompi-2021b
hisat2 -p 128 --dta -x /hpcdata/Mimir/kao25/AliRıza_Bey_HPC_training/Human_genome/Homo_sapiens.GRCh38.dna.primary_assembly_index --rna-strandness FR -1 /hpcdata/Mimir/kao25/AliRıza_Bey_HPC_training/SRR11192598_1.fast -2 SRR11192598_2.fast -S /hpcdata/Mimir/kao25/AliRıza_Bey_HPC_training/SRR11192598___new.sam
exit
