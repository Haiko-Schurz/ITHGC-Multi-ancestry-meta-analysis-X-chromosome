# ITHGC-Multi-ancestry-meta-analysis-X-chromosome

These scripts were used to do the X-chromosome specific analysis of the International Tuberculosis host genetics consortium X chromosome specific GWAS meta-analysis. Submission of this manuscript is currently in progress, but it is a companion paper to the ITHGC autosome meta-analysis paper (https://doi.org/10.7554/eLife.84394). 

**ITHGC_concordance_effect_direction_Xlinked.R**: Script used to do the concordance in direction of effects analysis between the different datasets, sex-stratified for males and females and combined.

**SexBiasAdmixture.R**: Script to assess significant differences in ancestral components between males and females. Ancestral proportions of the X chromosome and autosome were calculated using the Admixture version software https://dalexander.github.io/admixture/download.html.  

**Modeling X chromosome inactivation. **
These scripts were used to model X chromosome inactivation in the X-linked genetic association testing. THe scripts were translated from MatLab code and altered. Original script and methods were obtainerd from: Wang J and Shete S, 2014, Genetic Epidemiology
**fxchrom.R:** Script to performing association test for X-chromosome, while modelling different variations of X chromosome inactivation for females. 
**fLLR.R:** Log-likelihood ratio testing function for the fxchrom.R script to determine the most likely X chromosome inactivation model. 
**fperm.R:** Optional permutation function to adjust for multiple test correction when assessing multiple X-linked genetic variants. This function was replaced by a faster false discovry rate adjustment function, but can be implemnted if desired. 
