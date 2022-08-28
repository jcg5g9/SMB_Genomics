
#! /bin/bash

#SBATCH --account=biosci
#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J admix # give the job a custom name
#SBATCH -o admix_results-%j.out  # give the job output a custom name
#SBATCH -t 0-03:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 7  # number of cores (AKA tasks)

# Commands here run only on the first core

# Modules to load

module load admixture_linux/admixture_linux-1.3.

## Commands with srun will run on neo cores in the neoocation

admixture --seed=2324 --cv=10 ../../data/processed_bed/all_samples/01_all.bed 1

