#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J get_snp_table  # give the job a custom name
#SBATCH -o snp_table_vcf_allfilters-%j.out # give the job output a custom name
#SBATCH -t 0-02:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Load GATK module

module load rss/rss-2020
module load gatk/gatk-4.0.1.1

# Commands with srun will run on all cores in the allocation

## snp table generation: generate a snp	table for the popgen dataset (dataset A)
#gatk VariantsToTable \
#	-V ../../data/processed_vcf/A_02_popgen_filter_mac.vcf \
#	-F CHROM \
#	-F POS \
#	-F ID \
#	-F QUAL \
#	-F AC \
#	-F HET \
#	-F HOM-REF \
#	-F HOM-VAR \
#	-F NO-CALL \
#	-F VAR \
#	-F NCALLED \
#	-O ../../data/filtering_data/A_02_popgen_snps.txt

## snp table generation: generate a snp table for the finerad dataset (dataset B)
gatk VariantsToTable \
        -V ../../data/processed_vcf/B_01_finerad_filter_mac.vcf \
        -F CHROM \
        -F POS \
        -F ID \
        -F QUAL \
        -F AC \
        -F HET \
        -F HOM-REF \
        -F HOM-VAR \
        -F NO-CALL \
        -F VAR \
        -F NCALLED \
        -O ../../data/filtering_data/B_01_finerad_snps.txt
