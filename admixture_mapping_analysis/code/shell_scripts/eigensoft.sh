#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J eigensoft # give the job a custom name
#SBATCH -o eigensoft.out  # give the job output a custom name
#SBATCH -t 0-01:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Commands here run only on the first core

# Modules to load

module load rss/rss-2020
module load eigensoft/eigensoft-7.2.1

# Commands with srun will run on all cores in the allocation

convertf -p ../batch_cmd_lists/par.PED.EIGENSTRAT
