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

gatk VariantsToTable \
	-V smb_genomics_test_qual.vcf \
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
	-O smb_genomics_test_qual.txt
