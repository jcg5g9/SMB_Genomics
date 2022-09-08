#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J dadi  # give the job a custom name
#SBATCH -o dadi.out  # give the job output a custom name
#SBATCH -t 2-00:00  # two day time limit
#SBATCH --mem 100G

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)

# Commands here run only on the first core
module load miniconda3

### Run DADI

## DADI run 01: Run DADI on BAYOU and WHITE populations
#python3 ../analysis_scripts/bayou_white.py

## DADI run 02:	Run DADI on ELK and WHITE populations
#python3 ../analysis_scripts/elk_white.py

## DADI run 03:	Run DADI on ILLI and SKIA populations
#python3 ../analysis_scripts/illi_skia.py

## DADI run 04:	Run DADI on UPPARK and WHITE populations
#python3 ../analysis_scripts/uppark_white.py

### Summarize DADI output

## Summarize DADI output 01: Summarize DADI results for BAYOU and WHITE populations
#python3 ../base_scripts/Summarize_Outputs.py ../../data/dadi_output/bayou_white/

## Summarize DADI output 02: Summarize DADI results for ELK and WHITE populations
#python3 ../base_scripts/Summarize_Outputs.py ../../data/dadi_output/elk_white/

## Summarize DADI output 03: Summarize DADI results for ILLI and SKIA populations
#python3 ../base_scripts/Summarize_Outputs.py ../../data/dadi_output/illi_skia/

## Summarize DADI output 04:	Summarize DADI results for UPPARK and WHITE populations
#python3 ../base_scripts/Summarize_Outputs.py ../../data/dadi_output/uppark_white/

### Plot DADI

## Plot DADI 01:	Plot DADI on BAYOU and WHITE populations
#python3 ../plotting_scripts/bayou_white.py

## Plot DADI 02:	Plot DADI on ELK and WHITE populations
#python3 ../plotting_scripts/elk_white.py

## Plot DADI 03:	Plot DADI on ILLI and SKIA populations
#python3 ../plotting_scripts/illi_skia.py

## Plot DADI 04:	Plot DADI on UPPARK and WHITE populations
#python3 ../plotting_scripts/uppark_white.py

