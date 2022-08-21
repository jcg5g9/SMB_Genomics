#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J get_missing  # give the job a custom name
#SBATCH -o vcf_tools_missing-%j.out  # give the job output a custom name
#SBATCH -t 0-02:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Load VCFtools module

module load rss/rss-2020
module load vcftools/vcftools-v0.1.14

# Commands with srun will run on all cores in the allocation

#vcftools --vcf ../../../raw_data/smb_genomics.vcf --minDP 15 --recode --recode-INFO-all --out ../../data/processed_vcf/01_filter_depth
vcftools --vcf ../../data/processed_vcf/01_filter_depth.vcf --missing-indv --out ../../data/filtering_data/smb_genomnics_missing_genotype
