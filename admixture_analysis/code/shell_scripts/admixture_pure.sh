#! /bin/bash

#SBATCH --account=biosci
#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J admix # give the job a custom name
#SBATCH -o 02_pure.out  # give the job output a custom name (CHANGE THIS OUTFILE NAME TO MATCH OUTPUT Q AND P FILES)
#SBATCH -t 0-03:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 7  # number of cores (AKA tasks)

# Commands here run only on the first core

## UNCOMMENT LINES STARTING HERE 

echo "### Starting at: $(date) ###"


## load packages (LEAVE THIS COMMENTED)
module load rss/rss-2020
module load admixture_linux/admixture_linux-1.3.0

COMMANDA=`head -n ${SLURM_ARRAY_TASK_ID} ../batch_cmd_lists/02_pure_batch_cmd_list.txt | tail -n 1`
eval $COMMANDA

echo "### Ending at: $(date) ###"
