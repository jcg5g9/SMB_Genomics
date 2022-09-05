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

## filter 01: omit Spotted Bass individuals to get samples at the species complex level (Smallmouth Bass and Neosho Bass only)
#vcftools --vcf ../../data/processed_vcf/01_popgen_spb_hybrid.vcf --remove ../../data/filtering_data/spotted_bass.txt --recode --recode-INFO-all --out ../../data/processed_vcf/02_popgen_species_complex

## filter 02: omit Smallmouth Bass samples to get only Neosho samples for the Neosho Bass population level
#vcftools --vcf ../../data/processed_vcf/02_popgen_species_complex.vcf --remove ../../data/filtering_data/smallmouth_bass.txt --recode --recode-INFO-all --out ../../data/processed_vcf/03_popgen_neosho

## filter 03: omit Neosho Bass samples to get only Smallmouth Bass samples for the Smallmouth Bass population level
#vcftools --vcf ../../data/processed_vcf/02_popgen_species_complex.vcf --remove ../../data/filtering_data/neosho_bass.txt --recode --recode-INFO-all --out ../../data/processed_vcf/04_popgen_smallmouth

## filter 04A: Get outlier SNPs for all black bass samples
#vcftools --vcf ../../data/processed_vcf/01_popgen_spb_hybrid.vcf --positions ../../data/filtering_data/all_outlier_snps.txt --recode --recode-INFO-all --out ../../data/processed_vcf/01_popgen_all_outlier_snps

## filter 04B: Get neutral SNPs for all black bass samples
#vcftools --vcf ../../data/processed_vcf/01_popgen_spb_hybrid.vcf --positions ../../data/filtering_data/all_neutral_snps.txt --recode --recode-INFO-all --out ../../data/processed_vcf/01_popgen_all_neutral_snps

## filter 05A: Get outlier SNPs for species complex samples
#vcftools --vcf ../../data/processed_vcf/02_popgen_species_complex.vcf --positions ../../data/filtering_data/species_complex_outlier_snps.txt --recode --recode-INFO-all --out ../../data/processed_vcf/02_popgen_species_complex_outlier_snps

## filter 05B: Get neutral SNPs for species complex samples
#vcftools --vcf ../../data/processed_vcf/02_popgen_species_complex.vcf --positions ../../data/filtering_data/species_complex_neutral_snps.txt --recode --recode-INFO-all --out ../../data/processed_vcf/02_popgen_species_complex_neutral_snps

## filter 06A: Get outlier SNPs for Neosho Bass samples
#vcftools --vcf ../../data/processed_vcf/03_popgen_neosho.vcf --positions ../../data/filtering_data/neosho_outlier_snps.txt --recode --recode-INFO-all --out ../../data/processed_vcf/03_popgen_neosho_outlier_snps

## filter 06B: Get neutral SNPs for Neosho Bass samples
#vcftools --vcf ../../data/processed_vcf/03_popgen_neosho.vcf --positions ../../data/filtering_data/neosho_neutral_snps.txt --recode --recode-INFO-all --out ../../data/processed_vcf/03_popgen_neosho_neutral_snps

## filter 07A: Get outlier SNPs for Smallmouth Bass samples
vcftools --vcf ../../data/processed_vcf/04_popgen_smallmouth.vcf --positions ../../data/filtering_data/smallmouth_outlier_snps.txt --recode --recode-INFO-all --out ../../data/processed_vcf/04_popgen_smallmouth_outlier_snps

## filter 07B: Get neutral SNPs for Smallmouth Bass samples
vcftools --vcf ../../data/processed_vcf/04_popgen_smallmouth.vcf --positions ../../data/filtering_data/smallmouth_neutral_snps.txt --recode --recode-INFO-all --out ../../data/processed_vcf/04_popgen_smallmouth_neutral_snps
