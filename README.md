# Reproduce analyses for Gunn et al. (DATE)
<font size="+1">Follow the steps listed below in the <b><i>Analyses</i></b> section to reproduce analyses for Gunn et al. (DATE). Each step below gives a summary of the analysis and directs you to a general code file which then works through the analysis step-by-step. This general file will usually point you to other Rmd code, bash shell scripts, or python scripts. Each analysis is contained within subdirectories of the same name in the main R project directory.</font>

# Project: Population genomic analysis of Smallmouth Bass and Neosho Bass in the Central Interior Highlands
We investigated the extent of genomic divergence, local directional selection, and admixture between the Smallmouth Bass (<i>Micropterus dolomieu</i>) and the Neosho Bass (<i>M. velox</i>) in the Central Interior Highlands (CIH) ecoregion of central north America. Specifically, we used ddRADseq data to assessed the phylogenomic relationship between and within species, characterizing inter- and intraspecific diversity and SNPs potentially under local directional selection at the population level. Additionally, we inferred the relative timing of admxiture in Neosho Bass streams where there is known introgressive hybridization with Smallmouth Bass to understand the influence of natural, historic geographic factors on mixing (stream capture or transient flooding) vs. anthropogenic factors (i.e., non-native introductions through stocking), which is known to have occurred widely in these economically valuable species. We ultimately hoped to provide novel insights into the diversity of endemic, ecologically important and popular sport fish in the CIH.

## Analyses

### Analysis 1: Generating Species Native Range Maps
In this analysis, we generated easily readible maps displaying the native distributions of the two species of interest, Smallmouth Bass and Neosho Bass. We generated two types of maps: 1) a full range map, in which the full native range of each species is displayed, and 2) a close-up map of the Central Interior Highlands (CIH), where the paraptry of the species' ranges is shown in detail. In R, we generated only georeferenced outlines of these maps. Shapes representing stream sites and/or populations were superimposed a posteriori on the maps in PowerPoint.

#### Run the code: `map_analysis/smb_genomics_map_analysis.Rmd`

### Analysis 2: SNP Filtering, Data Processing, and Preliminary Calculations


#### Run the code: `allele_frequenc_analysis/smb_genomics_allele_frequency_analysis.Rmd`
