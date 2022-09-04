#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J pgd_coversion  # give the job a custom name
#SBATCH -o pgd_results-%j.out  # give the job output a custom name
#SBATCH -t 0-02:00  # two hour time limit
#SBATCH --mem-per-cpu=8G

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Commands here run only on the first core
module load rss/rss-2020
module load pgdspider/pgdspider-2.1.1.5 

# Commands with srun will run on all cores in the allocation
java -Xmx8g -Xms8g -jar ${PGDSPIDER_ROOT}/PGDSpider2-cli.jar -inputfile ../../data/processed_vcf/03_popgen_spb.vcf -inputformat VCF -outputfile 01_all -outputformat GESTE_BAYE_SCAN -spid ../spid_files/spid_files/01_all.spid
