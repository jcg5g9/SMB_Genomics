#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J finerad  # give the job a custom name
#SBATCH -o finerad.out  # give the job output a custom name
#SBATCH -t 0-02:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Commands here run only on the first core
module load rss/rss-2020
module load fineradstructure/fineradstructure-0.3.2

# Commands with srun will run on all cores in the allocation

## file conversion 01: convert VCF into haplotype chunks for all Smallmouth Bass and Neosho Bass
#RADpainter hapsFromVCF ../../data/processed_vcf/finerad.vcf > ../../data/hap_maps/01_smb_neosho.txt
#RADpainter paint ../../data/hap_maps/01_smb_neosho.txt

## file conversion 02: convert VCF into haplotype chunks for all pure samples
#RADpainter hapsFromVCF ../../data/processed_vcf/A_03_finerad_pure.vcf > ../../data/hap_maps/02_pure.txt
#RADpainter paint ../../data/hap_maps/02_pure.txt

## file conversion 03: convert VCF into haplotype chunks for all admixed samples
#RADpainter hapsFromVCF ../../data/processed_vcf/B_03_finerad_admixed.vcf > ../../data/hap_maps/03_admixed.txt
#RADpainter paint ../../data/hap_maps/03_admixed.txt

## Coancestry 01: Smallmouth Bass and Neosho Bass
#finestructure -x 100000 -y 100000 -z 1000 ../../hap_maps/01_smb_neosho_chunks.out ../../mcmc_files/01_smb_neosho_chunks.mcmc.xml
#finestructure -m T -x 10000 ../../hap_maps/01_smb_neosho_chunks.out ../../mcmc_files/01_smb_neosho_chunks.mcmc.xml 01_smb_neosho_chunks.mcmcTree.xml

## Coancestry 02: Pure samples
#finestructure -x 100000 -y 100000 -z 1000 ../../hap_maps/02_pure_chunks.out ../../mcmc_files/02_pure_chunks.mcmc.xml
#finestructure -m T -x 10000 ../../hap_maps/02_pure_chunks.out ../../mcmc_files/02_pure_chunks.mcmc.xml 02_pure_chunks.mcmcTree.xml

## Coancestry 03:  Smallmouth Bass and Neosho Bass
#finestructure -x 100000 -y 100000 -z 1000 ../../hap_maps/03_admixed_chunks.out ../../mcmc_files/03_admixed_chunks.mcmc.xml
#finestructure -m T -x 10000 ../../hap_maps/03_admixed_chunks.out ../../mcmc_files/03_admixed_chunks.mcmc.xml 03_admixed_chunks.mcmcTree.xml
