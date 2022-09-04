#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J SMB_Genomics_BayeScan  # give the job a custom name
#SBATCH -o bayescan_results_smb_lineages-%j.out  # give the job output a custom name
#SBATCH -t 1-00:00:00  # two day time limit
#SBATCH --account=biosci
#SBATCH -N 1  # number of nodes
#SBATCH -n 21  # number of cores (AKA tasks)

# Load modules needed
module load rss/rss-2020
module load bayescan/bayescan-2.1

# Commands here run only on the first core

# Commands with srun will run on all cores in the allocation

## analysis 1: all black bass samples
# bayescan 01_all -od ../../data/bayescan_output -o 01_all

## analysis 2: species complex samples (Neosho Bass and Smallmouth Bass)
# bayescan 02_species_complex -od ../../data/bayescan_output -o 02_species_complex

## analysis 3: Neosho Bass samples (Neosho Bass only)
# bayescan 03_neosho -od ../../data/bayescan_output -o 03_neosho

## analysis 1: all black bass samples
# bayescan 04_smallmouth -od ../../data/bayescan_output -o 04_smallmouth
