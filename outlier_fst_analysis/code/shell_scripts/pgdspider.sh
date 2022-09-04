#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J pgdspider_coversion  # give the job a custom name
#SBATCH -o pgdspider.out  # give the job output a custom name
#SBATCH -t 0-02:00  # two hour time limit
#SBATCH --mem-per-cpu=8G

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Commands here run only on the first core
module load rss/rss-2020
module load pgdspider/pgdspider-2.1.1.5 

# Commands with srun will run on all cores in the allocation

## conversion 01: all black bass samples (excluding the Spotted Bass hybrid, BFC10)
#java -Xmx8g -Xms8g -jar ${PGDSPIDER_ROOT}/PGDSpider2-cli.jar -inputfile ../../data/processed_vcf/01_popgen_spb_hybrid.vcf -inputformat VCF -outputfile ../../data/processed_bayescan/01_all -outputformat GESTE_BAYE_SCAN -spid ../spid_files/01_all.spid

## conversion 02: species complex samples (Neosho Bass and Smallmouth Bass, excluding Spotted Bass)
#java -Xmx8g -Xms8g -jar ${PGDSPIDER_ROOT}/PGDSpider2-cli.jar -inputfile ../../data/processed_vcf/02_popgen_species_complex.vcf -inputformat VCF -outputfile ../../data/processed_bayescan/02_species_complex -outputformat GESTE_BAYE_SCAN -spid ../spid_files/02_species_complex.spid

## conversion 03: Neosho Bass samples
#java -Xmx8g -Xms8g -jar ${PGDSPIDER_ROOT}/PGDSpider2-cli.jar -inputfile ../../data/processed_vcf/03_popgen_neosho.vcf -inputformat VCF -outputfile ../../data/processed_bayescan/03_neosho -outputformat GESTE_BAYE_SCAN -spid ../spid_files/03_neosho.spid

## conversion 04: Smallmouth Bass samples
java -Xmx8g -Xms8g -jar ${PGDSPIDER_ROOT}/PGDSpider2-cli.jar -inputfile ../../data/processed_vcf/04_popgen_smallmouth.vcf -inputformat VCF -outputfile ../../data/processed_bayescan/04_smallmouth -outputformat GESTE_BAYE_SCAN -spid ../spid_files/04_smallmouth.spid
