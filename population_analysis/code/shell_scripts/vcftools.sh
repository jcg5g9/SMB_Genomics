#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J vcftools  # give the job a custom name
#SBATCH -o filter_spotted.out  # give the job output a custom name
#SBATCH -t 0-02:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Load VCFtools module

module load rss/rss-2020
module load vcftools/vcftools-v0.1.14

# Commands with srun will run on all cores in the allocation
## filter 01: omit Spotted Bass
vcftools --vcf ../../data/processed_vcf/finerad.vcf --remove ../../data/filtering_data/bfc10.txt --recode --recode-INFO-all --out ../../data/processed_vcf/01_finerad_spb_hybrid
