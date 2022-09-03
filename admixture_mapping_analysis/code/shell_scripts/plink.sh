#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J get_ped  # give the job a custom name
#SBATCH -o plink.out  # give the job output a custom name
#SBATCH -t 0-01:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Commands here run only on the first core

# Modules to load

module load plink/plink-high-contig-1.90p

# Commands with srun will run on all cores in the allocation

plink --vcf ../../data/processed_vcf/01_popgen_spb_hybrid.vcf --allow-extra-chr --vcf-half-call m --recode --out ../../data/processed_ped/popgen
