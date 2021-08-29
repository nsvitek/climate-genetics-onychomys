#!/bin/bash
#
#SBATCH --job-name=R_fst_r60
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=1:30:00
#SBATCH -p short-40core
#SBATCH --output %x_%j.o
#SBATCH --error %x_%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natasha.vitek@stonybrook.edu   # Your email address

module load shared
module load R/4.0.2

cd /gpfs/scratch/nvitek/onychomys_leucogaster

Rscript oleu_fst_r60.R
