#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J ILLI_SKIA  # give the job a custom name
#SBATCH -o ILLI_SKIA-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00  # two day time limit
#SBATCH --mem 100G

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)

# Commands here run only on the first core
module load miniconda3

# Commands with srun will run on all cores in the allocation

## DADI run 01: Run DADI on BAYOU and WHITE populations
#python3 ../analysis_scripts/bayou_white.py

## DADI run 02:	Run DADI on ELK and WHITE populations
#python3 ../analysis_scripts/elk_white.py

## DADI run 03:	Run DADI on ILLI and SKIA populations
#python3 ../analysis_scripts/illi_skia.py

## DADI run 04:	Run DADI on UPPARK and WHITE populations
#python3 ../analysis_scripts/uppark_white.py
