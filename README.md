# Reproduce analyses for Gunn et al. (<i>in revision</i>)
<font size="+1">Follow the steps listed below in the <b><i>Analyses</i></b> section to reproduce analyses for Gunn et al. (<i>in revision</i>). Each step below gives a summary of the analysis and directs you to a general code file which then works through the analysis step-by-step. This general file will usually point you to other Rmd code, bash shell scripts, or python scripts. Each analysis is contained within subdirectories of the same name in the main R project directory.</font>

## Project: Population genomic analysis of Smallmouth Bass and Neosho Bass in the Central Interior Highlands
We investigated the extent of genomic divergence, local directional selection, and admixture between the Smallmouth Bass (<i>Micropterus dolomieu</i>) and the Neosho Bass (<i>M. velox</i>) in the Central Interior Highlands (CIH) ecoregion of central north America. Specifically, we used ddRADseq data to assessed the phylogenomic relationship between and within species, characterizing inter- and intraspecific diversity and SNPs potentially under local directional selection at the population level. Additionally, we inferred the relative timing of admxiture in Neosho Bass streams where there is known introgressive hybridization with Smallmouth Bass to understand the influence of natural, historic geographic factors on mixing (stream capture or transient flooding) vs. anthropogenic factors (i.e., non-native introductions through stocking), which is known to have occurred widely in these economically valuable species. We ultimately hoped to provide novel insights into the diversity of endemic, ecologically important and popular sport fish in the CIH.

## Analyses

### Analysis 1: Generating Species Native Range Maps
In this analysis, we generated easily readible maps displaying the native distributions of the two species of interest, Smallmouth Bass and Neosho Bass. We generated two types of maps: 1) a full range map, in which the full native range of each species is displayed, and 2) a close-up map of the Central Interior Highlands (CIH), where the paraptry of the species' ranges is shown in detail. In R, we generated only georeferenced outlines of these maps. Shapes representing stream sites and/or populations were superimposed a posteriori on the maps in PowerPoint.

#### Run the code: `map_analysis/smb_genomics_map_analysis.Rmd`

### Analysis 2: SNP Filtering, Data Processing, and Preliminary Calculations
In this analysis, we performed further quality filtering on the processed and genotyped SNPs generated at Floragenex, Inc. for Smallmouth Bass and Neosho Bass. Specifically, we screened the processed data for SNPs with greater than 15X read depth; fish samples with less than 20% genotype calls across all SNPs ("badsamples"); SNPs with a phred quality score less than 20 ("qual"); and SNPs with greater than 20% missing genotype calls across fish individuals ('missing').

#### Run the code: `filtering_processing_analysis/smb_genomics_filtering_processing_analysis.Rmd`

### Analysis 3: Admixture and phylogenomics
In this analysis, we used the popgen.vcf data generated in Analysis 2 (SNP Filtering...) to assess population genomic structure of Spotted Bass, Smallmouth Bass, and Neosho Bass in the CIH. Specifically, we conducted an initial screen of hybridization and gene flow by running maximum likelihood clustering on the full, filtered dataset and identified any individuals of interspecific origin between Spotted Bass and all other Interior Highlands fish (Smallmouth Bass and Neosho Bass) and individuals of interspecific origin between Smallmouth Bass and Neosho Bass. After removing hybrids, we conducted a separate admixture analysis and complementary phylogenomic analysis on the "pure" genomic samples to estimate genomic divergence between species.

#### Run the code: `filtering_processing_analysis/smb_genomics_admixture_phylogenomics_analysis.Rmd`

### Analysis 4: Population Inference
In this analysis, we used the finerad.vcf data generated in Analysis 2 (SNP Filtering...) to assess fine-scale coancestry between Smallmouth Bass, and Neosho Bass in the CIH using haplotype inference (excluding Spotted Bass). Specifically, we estimated coancestry in 1) the full dataset, with all pure and admixed individuals, excluding the Spotted Bass X Smallmouth Bass hybrid (BFC10) inferred from population genomic analysis in Analysis 3, 2) the pure dataset, with only pure individuals of Smallmouth Bass and Neosho Bass, and 3) the admixed dataset, with only admixed individuals of Neosho Bass (no admixed Smallmouth Bass were detected).

#### Run the code: `population_analysis/smb_genomics_population_analysis.Rmd`

### Analysis 5: Admixture Mapping
In this analysis, we used the popgen.vcf data generated in Analysis 2 (SNP Filtering...) to assess assess the relative timing of admixture events between Smallmouth Bass and Neosho Bass. Specifically, we used moment statistics in MatLab with the software program MIXMAPPER to build a scaffold phylogeny with significantly pure (non-admixed) populations (based on our <i>a posteriori</i>) discovered populations in Analysis 4) and to map significantly admixed populations (also based on our discovered populations in Analysis 4) onto the tree. 

#### Run the code: `population_analysis/smb_genomics_admixture_mapping_analysis.Rmd`
