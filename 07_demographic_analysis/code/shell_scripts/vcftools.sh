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

## filter 01: keep BAYOU and WHITE samples
vcftools --vcf ../../data/processed_vcf/02_popgen_species_complex_neutral_snps.vcf --keep ../../data/filtering_data/bayou_white.txt --recode --recode-INFO-all --out ../../data/processed_vcf/01_popgen_bayou_white

## filter 02: keep ELK and WHITE population samples
vcftools --vcf ../../data/processed_vcf/02_popgen_species_complex_neutral_snps.vcf --keep ../../data/filtering_data/elk_white.txt --recode --recode-INFO-all --out ../../data/processed_vcf/02_popgen_elk_white

## filter 03: keep ILLI and SKIA population samples
vcftools --vcf ../../data/processed_vcf/02_popgen_species_complex_neutral_snps.vcf --keep ../../data/filtering_data/illi_skia.txt --recode --recode-INFO-all --out ../../data/processed_vcf/03_popgen_illi_skia

## filter 04: keep UPPARK and WHITE population samples
vcftools --vcf ../../data/processed_vcf/02_popgen_species_complex_neutral_snps.vcf --keep ../../data/filtering_data/uppark_white.txt --recode --recode-INFO-all --out ../../data/processed_vcf/04_popgen_uppark_white
