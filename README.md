# The role of the microbiota in respiratory virus-bacterial pathobiont relationships in the upper respiratory tract

Author: Matthew Kelly <a href="https://orcid.org/0000-0001-8819-2315" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a>  
Last update: October 22, 2024

This repository contains the files and code necessary to replicate the analyses presented in the manuscript '_The role of the microbiota in respiratory virus-bacterial pathobiont relationships in the upper respiratory tract_', which has been submitted to medRxiv as a preprint (MEDRXIV/2024/315478). The overall objective of this manuscript was to investigate relationships between respiratory viruses, the URT bacterial microbiota, and bacterial respiratory pathobionts during infancy. 

## Overview

This repository contains the data files ([`Data_Files/`](Data_Files/)), scripts ([`Scripts/`](Scripts/)), and outputs ([`Outputs/`](Outputs/)) for each of the analyses presented in this manuscript. The raw sequencing files are publicly available in the Sequence Read Archive ([PRJNA698366](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA698366)). 

## File Description

- [`Data_Files`](Data_Files/)/

  - [`phy.bots.nps.rds`](Data_Files/phy.bots.nps.rds): phyloseq object containing raw read counts and metadata for all sequenced samples
  - [`phy.inf.np.16s.rds`](Data_Files/phy.inf.np.16s.rds): phyloseq object containing processed sequencing data from infant nasopharyngeal samples (n=2235) after contaminant removal
  - [`phy.inf.np.picrust.rds`](Data_Files/phy.inf.np.picrust.rds): phyloseq object containing microbial pathway abundances for infant nasopharyngeal samples (n=2235) 
  - [`metadata_inf_np.csv`](Data_Files/metadata_inf_np.csv): metadata file with raw data associated with all infant study visit (n=2409)
  - [`metadata_inf_np_RV.csv`](Data_Files/metadata_inf_np_RV.csv): metadata file with processed data associated with all infant study visits (n=2409)
  - [`infant.png`](Data_Files/infant.png): infant image file for Figure 1

- [`Scripts`](Scripts/)/

  - [`RV_Data_Preprocessing.R`](Scripts/RV_Data_Preprocessing.R): processing of raw sequencing data and metadata; inputs: [`phy.bots.nps.rds`](Data_Files/phy.bots.nps.rds), [`metadata_inf_np.csv`](Data_Files/metadata_inf_np.csv); outputs: [`phy.inf.np.16s.rds`](Data_Files/phy.inf.np.16s.rds), [`phy.inf.np.picrust.rds`](Data_Files/phy.inf.np.picrust.rds), [`metadata_inf_np_RV.csv`](Data_Files/metadata_inf_np_RV.csv)
  - [`RV_PCR_Analyses.R`](Scripts/RV_PCR_Analyses.R): analyses of PCR data on respiratory viruses and bacterial pathobionts; inputs: [`metadata_inf_np_RV.csv`](Data_Files/metadata_inf_np_RV.csv)
  - [`RV_16S_Analyses.R`](Scripts/RV_16S_Analyses.R): analyses of upper respiratory microbiota data; inputs: [`phy.inf.np.16s.rds`](Data_Files/phy.inf.np.16s.rds), [`phy.inf.np.picrust.rds`](Data_Files/phy.inf.np.picrust.rds); outputs: output files contained in [`Output_Files`](Output_Files/)
  - [`RV_Figure1.R`](Scripts/RV_Figure1.R): Figure 1; inputs: [`phy.inf.np.16s.rds`](Data_Files/phy.inf.np.16s.rds), [`metadata_inf_np_RV.csv`](Data_Files/metadata_inf_np_RV.csv), [`infant.png`](Data_Files/infant.png)
  - [`RV_Figure2.R`](Scripts/RV_Figure2.R): Figure 2; input: [`metadata_inf_np_RV.csv`](Data_Files/metadata_inf_np_RV.csv)
  - [`RV_Figure3.R`](Scripts/RV_Figure3.R): Figure 3; inputs: [`phy.inf.np.16s.rds`](Data_Files/phy.inf.np.16s.rds), model output files contained in [`Output_Files/Maaslin2/`](Output_Files/Maaslin2/)
  - [`RV_Figure4.R`](Scripts/RV_Figure4.R): Figure 4; inputs: output files contained in [`Output_Files/Random_Forest/`](Output_Files/Random_Forest/), [`Output_Files/Mixed_Effect_Logistic_Regression/`](Output_Files/Mixed_Effect_Logistic_Regression/)
  - [`RV_FigureS1.R`](Scripts/RV_FigureS1.R): Supplementary Figure 1; input: [`phy.inf.np.16s.rds`](Data_Files/phy.inf.np.16s.rds)
  - [`RV_FigureS2.R`](Scripts/RV_FigureS2.R): Supplementary Figure 2; input: [`phy.inf.np.16s.rds`](Data_Files/phy.inf.np.16s.rds)
  - [`RV_FigureS3.R`](Scripts/RV_FigureS3.R): Supplementary Figure 3; output files contained in [`Output_Files/`](Output_Files/)
