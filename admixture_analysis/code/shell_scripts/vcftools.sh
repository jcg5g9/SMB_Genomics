#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J filter_admixed  # give the job a custom name
#SBATCH -o filtering_missing.out  # give the job output a custom name
#SBATCH -t 0-02:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Load VCFtools module

module load rss/rss-2020
module load vcftools/vcftools-v0.1.14

# Commands with srun will run on all cores in the allocation

## filter 01: omit Spotted Bass X Neosho Bass hybrid (BFC10) 
vcftools --vcf ../../data/processed_vcf/popgenvcf --remove ../../filtering_data/bfc10.txt --out ../../data/filtering_data/01_popgen_spb_hybrid

