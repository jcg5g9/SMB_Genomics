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

## missing file generation: generate a file of missing genotype call rates for all samples
#vcftools --vcf ../../../raw_data/smb_genomics_edited.vcf --missing-indv --out ../../data/filtering_data/smb_genomics_missing

## filter 01: omit poor quality individuals from dataset (BFORK32, GRSPB34, ER05)
#vcftools --vcf ../../../raw_data/smb_genomics_edited.vcf --remove ../../data/filtering_data/smb_genomics_badsamples.txt --recode --recode-INFO-all --out ../../data/processed_vcf/01_filter_badsamples

## filter 02: filter out SNPs with less than 15X depth of coverage
#vcftools --vcf ../../data/processed_vcf/01_filter_badsamples.vcf --minDP 15 --recode --recode-INFO-all --out ../../data/processed_vcf/02_filter_depth

## filter 03: filter out SNPs with phred quality score less than or equal to 20
#vcftools --vcf ../../data/processed_vcf/02_filter_depth.vcf --minQ 20 --recode --recode-INFO-all --out ../../data/processed_vcf/03_filter_qual

## filter 04: filter out SNPs with less than 80% missing genotype calls across samples
#vcftools --vcf ../../data/processed_vcf/03_filter_qual.vcf --max-missing 0.8 --recode --recode-INFO-all --out ../../data/processed_vcf/04_filter_genotype

### Important Note:
### After this step, we split into two distinct datasets:
### 1) Dataset A: popgen: the full, thinned dataset, keeping only one SNP per RAD tag (to avoid tight linkage) to be used for population genomic analysis, and 
### 2) Dataset B: finerad: the full, non-thinned dataset, to be used for haplotype inference in fineradstructure

## filter A_01: filter to retain only one SNP per RAD tag in the pop gen dataset (dataset A)
#vcftools --vcf ../../data/processed_vcf/04_filter_genotype.vcf --thin 100 --recode --recode-INFO-all --out ../../data/processed_vcf/A_01_popgen_filter_thin

## filter A_02: filter pop gen dataset (dataset A) to omit SNPs with minor allele count less than or equal to 2
#vcftools --vcf ../../data/processed_vcf/A_01_popgen_filter_thin.vcf --mac 2 --recode --recode-INFO-all --out ../../data/processed_vcf/A_02_popgen_filter_mac

## filter B_01:	filter out SNPs with minor allele count less than or equal to 2 (minor allele frequency ~ 0.011) in the finerad dataset (dataset B)
#vcftools --vcf ../../data/processed_vcf/04_filter_genotype.vcf --mac 2 --recode --recode-INFO-all --out ../../data/processed_vcf/B_01_finerad_filter_mac


