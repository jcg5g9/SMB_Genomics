#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J SMB_fineRAD  # give the job a custom name
#SBATCH -o finerad_structure_withfilters-%j.out  # give the job output a custom name
#SBATCH -t 0-02:00  # two hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)

# Commands here run only on the first core
module load fineradstructure/fineradstructure-0.3.2

# Commands with srun will run on all cores in the allocation

## file conversion 01: convert VCF into haplotype chunks for all Smallmouth Bass and Neosho Bass
RADpainter hapsFromVCF ../../vcfdata/vcf_pass2_badsamples_qual_missing_mac2_het45_admixed.vcf > ../fineradstructure/hap_withfilters_admixed.txt
# RADpainter paint ../fineradstructure/hap_withfilters_admixed.txt
# finestructure -x 100000 -y 100000 -z 1000 ../fineradstructure/hap_withfilters_admixed_chunks.out ../fineradstructure/hap_withfilters_admixed_chunks.mcmc.xml
# finestructure -m T -x 10000 ../fineradstructure/hap_withfilters_admixed_chunks.out ../fineradstructure/hap_withfilters_admixed_chunks.mcmc.xml hap_withfilters_admixed_chunks.mcmcTree.xml
