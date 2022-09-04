#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J plink # give the job a custom name
#SBATCH -o plink.out  # give the job output a custom name
#SBATCH -t 0-01:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

## Commands here run only on the first core

# (UNCOMMENT THIS LINE TO RUN THE CODE)

## Modules to load

module load rss/rss-2020
module load plink/plink-high-contig-1.90p

## Conversion 01: Convert VCF data for all black bass samples (excluding Spotted Bass hybrid, BFC10) to .bed format
#plink --vcf ../../data/processed_vcf/01_popgen_spb_hybrid.vcf --allow-extra-chr --vcf-half-call m --make-bed --out ../../data/processed_bed/01_all_output/01_all

## Conversion 02: Convert VCF data for species complex samples (Neosho Bass and Smallmouth Bass) to .bed format
plink --vcf ../../data/processed_vcf/02_popgen_species_complex.vcf --allow-extra-chr --vcf-half-call m --make-bed --out ../../data/processed_bed/02_species_complex_output/02_species_complex

## Conversion 03: Convert VCF data for Neosho Bass only to .bed format
plink --vcf ../../data/processed_vcf/03_popgen_neosho.vcf --allow-extra-chr --vcf-half-call m --make-bed --out ../../data/processed_bed/03_neosho_output/03_neosho

## Conversion 04: Convert VCF data for Smallmouth Bass only to .bed format
plink --vcf ../../data/processed_vcf/04_popgen_smallmouth.vcf --allow-extra-chr --vcf-half-call m --make-bed --out ../../data/processed_bed/04_smallmouth_output/04_smallmouth
