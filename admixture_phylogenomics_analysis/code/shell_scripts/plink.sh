#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J plink_conversion # give the job a custom name
#SBATCH -o get_plink_bed-%j.out  # give the job output a custom name
#SBATCH -t 0-01:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Commands here run only on the first core

# Modules to load

module load plink/plink-high-contig-1.90p

# Commands with srun will run on all cores in the allocation

## Conversion 01: Convert full vcf data (popgen.vcf) to .bed format
plink --vcf ../../data/processed_vcf/popgen.vcf --allow-extra-chr --vcf-half-call m --make-bed --out ../../data/processed_bed/all_samples/01_all

## Conversion 02: Convert full vcf data (02_popgen_pure.vcf) to .bed format
plink --vcf ../../data/processed_vcf/02_popgen_pure.vcf --allow-extra-chr --vcf-half-call m --make-bed --out ../../data/processed_bed/pure_samples/02_pure
