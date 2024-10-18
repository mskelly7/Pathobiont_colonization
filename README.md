# The dynamics of bacterial respiratory pathobiont colonization among infants in Botswana

Author: Matthew Kelly <a href="https://orcid.org/0000-0001-8819-2315" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a>  
Last update: October 14, 2024

This repository contains the files and code necessary to replicate the analyses presented in the manuscript '_The role of the microbiota in respiratory virus-bacterial pathobiont relationships in the upper respiratory tract_', which has been submitted to medRxiv as a preprint (https://www.medrxiv.org/content/______________________). The overall objective of this manuscript was ____________________

## Overview

This repository contains the data files ([`Data_Files/`](Data_Files/)), scripts ([`Scripts/`](Scripts/)), and outputs ([`Outputs/`](Outputs/)) for each of the analyses presented in this manuscript. The raw sequencing files are publicly available in the Sequence Read Archive ([PRJNA1024980](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1024980)). 

## File Description

- [`Data_Files`](Data_Files/)/

  - [`phy.bots.nps.rds`](Data_Files/phy.bots.nps.rds): phyloseq object containing raw read counts and metadata for all sequenced samples
  - [`phy.inf.np.16s.rds`](Data_Files/phy.inf.np.16s.rds): phyloseq object containing processed sequencing data from infant nasopharyngeal samples (n=2235) after contaminant removal
  - [`phy.inf.np.picrust.rds`](Data_Files/phy.inf.np.picrust.rds): phyloseq object containing microbial pathway abundances for infant nasopharyngeal samples (n=2235) 
  - [`metadata_inf_np.csv`](Data_Files/metadata_inf_np.csv): metadata file with raw data associated with all infant study visit (n=2409)
  - [`metadata_inf_np_RV.csv`](Data_Files/metadata_inf_np_RV.csv): metadata file with processed data associated with all infant study visits (n=2409)

- [`Scripts`](Scripts/)/

  - [`RV_Data_Preprocessing.R`](Scripts/RV_Data_Preprocessing.R): processing of raw sequencing data and metadata; inputs: [`phy.bots.nps.rds`](Data_Files/phy.bots.nps.rds), [`metadata_inf_np.csv`](Data_Files/metadata_inf_np.csv); outputs: [`phy.inf.np.16s.rds`](Data_Files/phy.inf.np.16s.rds), [`phy.inf.np.picrust.rds`](Data_Files/phy.inf.np.picrust.rds), [`metadata_inf_np_RV.csv`](Data_Files/metadata_inf_np_RV.csv)
  
  - [`RNASeq_Figure1.R`](Scripts/RNASeq_Figure1.R): Figure 1; inputs: [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds), [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds), FGSEA ("modules") output files in [`1_COVID_Neg_by_Age`](Output_Files/1_COVID_Neg_by_Age/)

- [`Output_Files`](Output_Files/)/ 
