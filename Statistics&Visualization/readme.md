### R scripts (RStudio version 1.2.5042) used to run Statistics + visualization in Madi et al. 2023 (Phage predation is a biomarker for disease severity and shapes pathogen genetic diversity in cholera patients).

#### Content :

* Post-inStrain filtering to remove false positive SNVs (adapted from Garud and Good et al. 2019)
* Figures: Code to generate figures and supplementary figures.
  
* GAMs : Code for the generalized additive mixed models with :
   -  Vc relative abundance from metagenomics reads as a function of phage, antibiotics and their interactions.
   -  Vc absolute abundance from qPCR data as a function of phage, antibiotics and their interactions.
   -  Average frequency of NS mutations in Vc as a function of phage, antibiotics and the anti-phage resistance profile (Vc and phage abundance from metagenomics reads).

* GLMM : the generalized linear mixed models with SNVs count in Vc as a function of phage and intibiotics (Vc and phage abundance from metagenomics reads).
  
* Ordinations :
   -  Principal component analysis (PCA) with all species annotated with Kraken2/Bracken.
   -  Redundancy analysis (RDA) with the most dominant species as a function of phages, antibiotics and patient's metadata.
   -  Multiple Factor Analysis (MFA) between the most dominant species and all the antibiotics resistance genes predicted with deepARG.
     
* Indicator species analysis to associate species with different sickness stages.
