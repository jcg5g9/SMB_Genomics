#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J snphylo  # give the job a custom name
#SBATCH -o snphylo_pure.out  # give the job output a custom name
#SBATCH -t 0-24:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000

## UNCOMMENT STARTING HERE:

## Commands here run only on the first core
module load rss/rss-2020
module load snphylo/snphylo-2016-02-04-python-2.7.14-tk

snphylo.sh -v ../../data/processed_vcf/02_popgen_pure_snphylo.vcf -r -a 288336 -b -B 1000 -o GLVR11_1

