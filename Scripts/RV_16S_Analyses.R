# Botswana Infant Microbiome Study - RV-Bacterial Pathogen Analyses
# Matthew Kelly, MD, MPH 
# Analyses of upper respiratory microbiome data
# Last update: March 9, 2025

remove(list=ls())
setwd("____________________") 
set.seed(1234)

version
library(phyloseq)
packageVersion("phyloseq")
library(tidyverse)
library(dplyr)
library(plyr)
library(data.table)
library(gridExtra)
library(httr)
library(reshape2)
library(DT)
library(ggplot2)
library(cowplot)
library(DataCombine)
library(vegan)
packageVersion("vegan")
library(ggpubr)
library(Maaslin2)
packageVersion("Maaslin2")
library(nlme)
packageVersion("nlme")
library(lme4)
packageVersion("lme4")
library(survminer)
library(gplots)
library(RColorBrewer)
library(splines)
library(pROC)
packageVersion("pROC")
library(randomForest)
packageVersion("randomForest")
library(rJava)
library(xlsx)
library(pheatmap)
library(knitr)

phy.inf.np.16s <- readRDS("phy.inf.np.16s.rds")
metadata_16s <- data.frame(sample_data(phy.inf.np.16s))
nsamples(phy.inf.np.16s)
ntaxa(phy.inf.np.16s)
# Create metadata file excluding 0-month samples
metadata_16s_noM0 <- subset(metadata_16s, month!="0")
phy.picrust <- readRDS("phy.inf.np.picrust.rds")
phy.picrust
table(sample_names(phy.picrust)==sample_names(phy.inf.np.16s))

# Create phyloseq object with only samples with complete data for respiratory viruses and bacterial pathobionts
phy.complete <- subset_samples(phy.inf.np.16s, !is.na(inf_rv_yn) & !is.na(inf_multi_hi) & !is.na(inf_multi_mc) & !is.na(inf_multi_sa) & !is.na(inf_multi_sp))
nsamples(phy.complete)
# Remove ASVs that are not present in any remaining samples
ntaxa(phy.complete)
phy.complete <- prune_taxa(taxa_sums(phy.complete) > 0, phy.complete)
ntaxa(phy.complete)

# Create metadata file for samples with complete data for respiratory viruses and bacterial pathobionts
metadata_complete <- data.frame(sample_data(phy.complete))
nrow(metadata_complete)
table(metadata_complete$inf_rv_yn, useNA="always")
table(metadata_complete$inf_multi_hi, useNA="always")
table(metadata_complete$inf_multi_mc, useNA="always")
table(metadata_complete$inf_multi_sa, useNA="always")
table(metadata_complete$inf_multi_sp, useNA="always")

# ******************************************************************************************************
# How does the prevalence of specific nasopharyngeal microbiota clusters change with age during infancy?
# ******************************************************************************************************

table(metadata_16s$month, metadata_16s$cluster)

# MIXED EFFECT LOGISTIC REGRESSION MODELS evaluating changes in the prevalences of nasopharyngeal microbiota profiles during infancy
# These analyses exclude samples collected at birth

metadata_16s_noM0$OTH[metadata_16s_noM0$cluster!="OTH"] <- 0
metadata_16s_noM0$OTH[metadata_16s_noM0$cluster=="OTH"] <- 1
metadata_16s_noM0$OTH <- as.numeric(metadata_16s_noM0$OTH)
logit_OTH_mo <- glmer(OTH ~ month + (1 | study_id), data = metadata_16s_noM0, family = binomial, 
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_OTH_mo)))
exp(tab <- cbind(Est = fixef(logit_OTH_mo), LL = fixef(logit_OTH_mo) - 1.96 * se, UL = fixef(logit_OTH_mo) + 1.96*se))
summary(logit_OTH_mo)

metadata_16s_noM0$STA[metadata_16s_noM0$cluster!="STA"] <- 0
metadata_16s_noM0$STA[metadata_16s_noM0$cluster=="STA"] <- 1
metadata_16s_noM0$STA <- as.numeric(metadata_16s_noM0$STA)
logit_STA_mo <- glmer(STA ~ month + (1 | study_id), data = metadata_16s_noM0, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_STA_mo)))
exp(tab <- cbind(Est = fixef(logit_STA_mo), LL = fixef(logit_STA_mo) - 1.96 * se, UL = fixef(logit_STA_mo) + 1.96*se))
summary(logit_STA_mo)

metadata_16s_noM0$COR[metadata_16s_noM0$cluster!="COR"] <- 0
metadata_16s_noM0$COR[metadata_16s_noM0$cluster=="COR"] <- 1
metadata_16s_noM0$COR <- as.numeric(metadata_16s_noM0$COR)
logit_COR_mo <- glmer(COR ~ month + (1 | study_id), data = metadata_16s_noM0, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_COR_mo)))
exp(tab <- cbind(Est = fixef(logit_COR_mo), LL = fixef(logit_COR_mo) - 1.96 * se, UL = fixef(logit_COR_mo) + 1.96*se))
summary(logit_COR_mo)

metadata_16s_noM0$STR[metadata_16s_noM0$cluster!="STR"] <- 0
metadata_16s_noM0$STR[metadata_16s_noM0$cluster=="STR"] <- 1
metadata_16s_noM0$STR <- as.numeric(metadata_16s_noM0$STR)
logit_STR_mo <- glmer(STR ~ month + (1 | study_id), data = metadata_16s_noM0, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_STR_mo)))
exp(tab <- cbind(Est = fixef(logit_STR_mo), LL = fixef(logit_STR_mo) - 1.96 * se, UL = fixef(logit_STR_mo) + 1.96*se))
summary(logit_STR_mo)

metadata_16s_noM0$CD[metadata_16s_noM0$cluster!="CD"] <- 0
metadata_16s_noM0$CD[metadata_16s_noM0$cluster=="CD"] <- 1
metadata_16s_noM0$CD <- as.numeric(metadata_16s_noM0$CD)
logit_CD_mo <- glmer(CD ~ month + (1 | study_id), data = metadata_16s_noM0, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_CD_mo)))
exp(tab <- cbind(Est = fixef(logit_CD_mo), LL = fixef(logit_CD_mo) - 1.96 * se, UL = fixef(logit_CD_mo) + 1.96*se))
summary(logit_CD_mo)

metadata_16s_noM0$HAE[metadata_16s_noM0$cluster!="HAE"] <- 0
metadata_16s_noM0$HAE[metadata_16s_noM0$cluster=="HAE"] <- 1
metadata_16s_noM0$HAE <- as.numeric(metadata_16s_noM0$HAE)
logit_HAE_mo <- glmer(HAE ~ month + (1 | study_id), data = metadata_16s_noM0, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_HAE_mo)))
exp(tab <- cbind(Est = fixef(logit_HAE_mo), LL = fixef(logit_HAE_mo) - 1.96 * se, UL = fixef(logit_HAE_mo) + 1.96*se))
summary(logit_HAE_mo)

metadata_16s_noM0$CDM[metadata_16s_noM0$cluster!="CDM"] <- 0
metadata_16s_noM0$CDM[metadata_16s_noM0$cluster=="CDM"] <- 1
metadata_16s_noM0$CDM <- as.numeric(metadata_16s_noM0$CDM)
logit_CDM_mo <- glmer(CDM ~ month + (1 | study_id), data = metadata_16s_noM0, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_CDM_mo)))
exp(tab <- cbind(Est = fixef(logit_CD_mo), LL = fixef(logit_CDM_mo) - 1.96 * se, UL = fixef(logit_CDM_mo) + 1.96*se))
summary(logit_CDM_mo)

metadata_16s_noM0$MOR[metadata_16s_noM0$cluster!="MOR"] <- 0
metadata_16s_noM0$MOR[metadata_16s_noM0$cluster=="MOR"] <- 1
metadata_16s_noM0$MOR <- as.numeric(metadata_16s_noM0$MOR)
logit_MOR_mo <- glmer(MOR ~ month + (1 | study_id), data = metadata_16s_noM0, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_MOR_mo)))
exp(tab <- cbind(Est = fixef(logit_MOR_mo), LL = fixef(logit_MOR_mo) - 1.96 * se, UL = fixef(logit_MOR_mo) + 1.96*se))
summary(logit_MOR_mo)

# ***************************************************************************************
# What is the impact of respiratory virus detection on nasopharyngeal microbiota cluster? 
# ***************************************************************************************

metadata_cluster <- metadata_complete[order(metadata_complete$study_id, metadata_complete$month),]
table(metadata_cluster$cluster)
metadata_cluster$cluster_original[metadata_cluster$cluster=="COR"] <- "COR"
metadata_cluster$cluster_original[metadata_cluster$cluster=="STR"] <- "STR"
metadata_cluster$cluster_original[metadata_cluster$cluster=="STA"] <- "STA"
metadata_cluster$cluster_original[metadata_cluster$cluster=="CD"] <- "CD"
metadata_cluster$cluster_original[metadata_cluster$cluster=="MOR"] <- "MOR"
metadata_cluster$cluster_original[metadata_cluster$cluster=="CDM"] <- "CDM"
metadata_cluster$cluster_original[metadata_cluster$cluster=="HAE"] <- "HAE"
metadata_cluster$cluster_original[metadata_cluster$cluster=="OTH"] <- "OTH"
metadata_cluster <- suppressWarnings(slide(metadata_cluster, "cluster_original", TimeVar="month", GroupVar="study_id", NewVar="cluster_lag", slideBy = -1))
metadata_cluster <- metadata_cluster[,c("SampleID","study_id","month","inf_rv_yn","inf_multi_hi",
                        "inf_multi_mc","inf_multi_sa","inf_multi_sp", "sex", "lbw", "mat_hiv", "residence", "wood", "num_kids", "season", 
                        "breastmilk", "inf_abx_any", "pcv", "cluster_original", "cluster_lag")]
metadata_cluster <- subset(metadata_cluster, !is.na(cluster_original) & !is.na(cluster_lag))
metadata_cluster$unstable[metadata_cluster$cluster_original=="OTH" & metadata_cluster$cluster_lag=="OTH"] <- 0
metadata_cluster$unstable[metadata_cluster$cluster_original=="OTH" & metadata_cluster$cluster_lag!="OTH"] <- 1
metadata_cluster$unstable[metadata_cluster$cluster_original=="STA" & metadata_cluster$cluster_lag=="STA"] <- 0
metadata_cluster$unstable[metadata_cluster$cluster_original=="STA" & metadata_cluster$cluster_lag!="STA"] <- 1
metadata_cluster$unstable[metadata_cluster$cluster_original=="COR" & metadata_cluster$cluster_lag=="COR"] <- 0
metadata_cluster$unstable[metadata_cluster$cluster_original=="COR" & metadata_cluster$cluster_lag!="COR"] <- 1
metadata_cluster$unstable[metadata_cluster$cluster_original=="STR" & metadata_cluster$cluster_lag=="STR"] <- 0
metadata_cluster$unstable[metadata_cluster$cluster_original=="STR" & metadata_cluster$cluster_lag!="STR"] <- 1
metadata_cluster$unstable[metadata_cluster$cluster_original=="CD" & metadata_cluster$cluster_lag=="CD"] <- 0
metadata_cluster$unstable[metadata_cluster$cluster_original=="CD" & metadata_cluster$cluster_lag!="CD"] <- 1
metadata_cluster$unstable[metadata_cluster$cluster_original=="HAE" & metadata_cluster$cluster_lag=="HAE"] <- 0
metadata_cluster$unstable[metadata_cluster$cluster_original=="HAE" & metadata_cluster$cluster_lag!="HAE"] <- 1
metadata_cluster$unstable[metadata_cluster$cluster_original=="CDM" & metadata_cluster$cluster_lag=="CDM"] <- 0
metadata_cluster$unstable[metadata_cluster$cluster_original=="CDM" & metadata_cluster$cluster_lag!="CDM"] <- 1
metadata_cluster$unstable[metadata_cluster$cluster_original=="MOR" & metadata_cluster$cluster_lag=="MOR"] <- 0
metadata_cluster$unstable[metadata_cluster$cluster_original=="MOR" & metadata_cluster$cluster_lag!="MOR"] <- 1
metadata_cluster$cluster_original <- as.factor(metadata_cluster$cluster_original)
metadata_cluster$unstable <- as.factor(metadata_cluster$unstable)
table(metadata_cluster$cluster_original, metadata_cluster$unstable)
prop.table(table(metadata_cluster$cluster_original, metadata_cluster$unstable),1)

logit_stability <- glmer(unstable ~ inf_rv_yn + month + 
                           (1 | study_id), data = metadata_cluster, family = binomial, 
                         control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_stability)))
exp(tab <- cbind(Est = fixef(logit_stability), LL = fixef(logit_stability) - 1.96 * se, UL = fixef(logit_stability) + 1.96*se))
summary(logit_stability)
# Respiratory virus infection tends to be associated with increased microbiome profile instability (p=0.046)

# MIXED EFFECT LOGISTIC REGRESSION MODELS evaluating for shifts in microbiota profiles with respiratory viruses 

metadata_complete$OTH[metadata_complete$cluster!="OTH"] <- 0
metadata_complete$OTH[metadata_complete$cluster=="OTH"] <- 1
metadata_complete$OTH <- as.numeric(metadata_complete$OTH)
logit_OTH <- glmer(OTH ~ inf_rv_yn + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any +
                     (1 | study_id), data = metadata_complete, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_OTH)))
exp(tab <- cbind(Est = fixef(logit_OTH), LL = fixef(logit_OTH) - 1.96 * se, UL = fixef(logit_OTH) + 1.96*se))
summary(logit_OTH)

metadata_complete$STA[metadata_complete$cluster!="STA"] <- 0
metadata_complete$STA[metadata_complete$cluster=="STA"] <- 1
metadata_complete$STA <- as.numeric(metadata_complete$STA)
logit_STA <- glmer(STA ~ inf_rv_yn + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any +
                     (1 | study_id), data = metadata_complete, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_STA)))
exp(tab <- cbind(Est = fixef(logit_STA), LL = fixef(logit_STA) - 1.96 * se, UL = fixef(logit_STA) + 1.96*se))
summary(logit_STA)

metadata_complete$COR[metadata_complete$cluster!="COR"] <- 0
metadata_complete$COR[metadata_complete$cluster=="COR"] <- 1
metadata_complete$COR <- as.numeric(metadata_complete$COR)
logit_COR <- glmer(COR ~ inf_rv_yn + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any +
                     (1 | study_id), data = metadata_complete, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_COR)))
exp(tab <- cbind(Est = fixef(logit_COR), LL = fixef(logit_COR) - 1.96 * se, UL = fixef(logit_COR) + 1.96*se))
summary(logit_COR)

metadata_complete$STR[metadata_complete$cluster!="STR"] <- 0
metadata_complete$STR[metadata_complete$cluster=="STR"] <- 1
metadata_complete$STR <- as.numeric(metadata_complete$STR)
logit_STR <- glmer(STR ~ inf_rv_yn + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any +
                     (1 | study_id), data = metadata_complete, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_STR)))
exp(tab <- cbind(Est = fixef(logit_STR), LL = fixef(logit_STR) - 1.96 * se, UL = fixef(logit_STR) + 1.96*se))
summary(logit_STR)

metadata_complete$CD[metadata_complete$cluster!="CD"] <- 0
metadata_complete$CD[metadata_complete$cluster=="CD"] <- 1
metadata_complete$CD <- as.numeric(metadata_complete$CD)
logit_CD <- glmer(CD ~ inf_rv_yn + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any +
                     (1 | study_id), data = metadata_complete, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_CD)))
exp(tab <- cbind(Est = fixef(logit_CD), LL = fixef(logit_CD) - 1.96 * se, UL = fixef(logit_CD) + 1.96*se))
summary(logit_CD)

metadata_complete$HAE[metadata_complete$cluster!="HAE"] <- 0
metadata_complete$HAE[metadata_complete$cluster=="HAE"] <- 1
metadata_complete$HAE <- as.numeric(metadata_complete$HAE)
logit_HAE <- glmer(HAE ~ inf_rv_yn + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any +
                    (1 | study_id), data = metadata_complete, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_HAE)))
exp(tab <- cbind(Est = fixef(logit_HAE), LL = fixef(logit_HAE) - 1.96 * se, UL = fixef(logit_HAE) + 1.96*se))
summary(logit_HAE)

metadata_complete$CDM[metadata_complete$cluster!="CDM"] <- 0
metadata_complete$CDM[metadata_complete$cluster=="CDM"] <- 1
metadata_complete$CDM <- as.numeric(metadata_complete$CDM)
logit_CDM <- glmer(CDM ~ inf_rv_yn + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any +
                     (1 | study_id), data = metadata_complete, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_CDM)))
exp(tab <- cbind(Est = fixef(logit_CDM), LL = fixef(logit_CDM) - 1.96 * se, UL = fixef(logit_CDM) + 1.96*se))
summary(logit_CDM)

metadata_complete$MOR[metadata_complete$cluster!="MOR"] <- 0
metadata_complete$MOR[metadata_complete$cluster=="MOR"] <- 1
metadata_complete$MOR <- as.numeric(metadata_complete$MOR)
logit_MOR <- glmer(MOR ~ inf_rv_yn + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any +
                    (1 | study_id), data = metadata_complete, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
se <- sqrt(diag(vcov(logit_MOR)))
exp(tab <- cbind(Est = fixef(logit_MOR), LL = fixef(logit_MOR) - 1.96 * se, UL = fixef(logit_MOR) + 1.96*se))
summary(logit_MOR)

rm(list=grep("logit",ls(),value=TRUE,invert=FALSE), tab, se)

# Impact of CLUSTER on pathobiont colonization

# Create datasets for analyses of pathobiont acquisition
acquisition <- metadata_16s
acquisition <- acquisition[order(acquisition$study_id, acquisition$month),]
acquisition$month <- as.numeric(acquisition$month)
acquisition$study_id <- as.factor(acquisition$study_id)
acquisition$inf_multi_hi[acquisition$inf_multi_hi=="Y"] <- 1
acquisition$inf_multi_hi[acquisition$inf_multi_hi=="N"] <- 0
acquisition$inf_multi_hi <- as.numeric(acquisition$inf_multi_hi)
acquisition$inf_multi_mc[acquisition$inf_multi_mc=="Y"] <- 1
acquisition$inf_multi_mc[acquisition$inf_multi_mc=="N"] <- 0
acquisition$inf_multi_mc <- as.numeric(acquisition$inf_multi_mc)
acquisition$inf_multi_sa[acquisition$inf_multi_sa=="Y"] <- 1
acquisition$inf_multi_sa[acquisition$inf_multi_sa=="N"] <- 0
acquisition$inf_multi_sa <- as.numeric(acquisition$inf_multi_sa)
acquisition$inf_multi_sp[acquisition$inf_multi_sp=="Y"] <- 1
acquisition$inf_multi_sp[acquisition$inf_multi_sp=="N"] <- 0
acquisition$inf_multi_sp <- as.numeric(acquisition$inf_multi_sp)
acquisition <- suppressWarnings(slide(acquisition, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="hi_lag", slideBy = -1))
acquisition <- suppressWarnings(slide(acquisition, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="mc_lag", slideBy = -1))
acquisition <- suppressWarnings(slide(acquisition, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="sa_lag", slideBy = -1))
acquisition <- suppressWarnings(slide(acquisition, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="sp_lag", slideBy = -1))
acquisition$cluster_original[acquisition$cluster=="COR"] <- "COR"
acquisition$cluster_original[acquisition$cluster=="STR"] <- "STR"
acquisition$cluster_original[acquisition$cluster=="STA"] <- "STA"
acquisition$cluster_original[acquisition$cluster=="CD"] <- "CD"
acquisition$cluster_original[acquisition$cluster=="MOR"] <- "MOR"
acquisition$cluster_original[acquisition$cluster=="CDM"] <- "CDM"
acquisition$cluster_original[acquisition$cluster=="HAE"] <- "HAE"
acquisition$cluster_original[acquisition$cluster=="OTH"] <- "OTH"
acquisition <- suppressWarnings(slide(acquisition, "cluster_original", TimeVar="month", GroupVar="study_id", NewVar="cluster_lag", slideBy = -1))
nrow(acquisition)
acquisition <- acquisition[,c("SampleID","study_id","month","inf_rv_yn","inf_rv_cat2","inf_multi_hi","hi_lag",
                        "inf_multi_mc","mc_lag","inf_multi_sa","sa_lag","inf_multi_sp","sp_lag", "sex", "lbw", "mat_hiv", "residence", "wood", "num_kids", "season", 
                        "breastmilk", "inf_abx_any", "pcv", "cluster_original", "cluster_lag")]
# Remove 0-month timepoint (do not predict outcome at 0 months)
acquisition <- subset(acquisition, month!="0")
nrow(acquisition)
# Remove intervals with baseline or persistent colonization
acquisition_hi <- subset(acquisition, hi_lag!=1)
nrow(acquisition_hi)
acquisition_mc <- subset(acquisition, mc_lag!=1)
nrow(acquisition_mc)
acquisition_sa <- subset(acquisition, sa_lag!=1)
nrow(acquisition_sa)
acquisition_sp <- subset(acquisition, sp_lag!=1)
nrow(acquisition_sp)
remove(acquisition)

# H. influenzae
acquisition_hi <- within(acquisition_hi, cluster_lag <- relevel(as.factor(cluster_lag), ref = "COR"))
logit_hi <- glmer(inf_multi_hi ~ cluster_lag + month + inf_rv_yn + 
                    sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                  data = acquisition_hi, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_hi)
logit_hi_se <- sqrt(diag(vcov(logit_hi)))
exp(cbind(Est = fixef(logit_hi), LL = fixef(logit_hi) - 1.96 * logit_hi_se, UL = fixef(logit_hi) + 1.96*logit_hi_se))

# M. catarrhalis
acquisition_mc <- within(acquisition_mc, cluster_lag <- relevel(as.factor(cluster_lag), ref = "COR"))
logit_mc <- glmer(inf_multi_mc ~ cluster_lag + month + inf_rv_yn + 
                    sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                  data = acquisition_mc, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_mc)
logit_mc_se <- sqrt(diag(vcov(logit_mc)))
exp(cbind(Est = fixef(logit_mc), LL = fixef(logit_mc) - 1.96 * logit_mc_se, UL = fixef(logit_mc) + 1.96*logit_mc_se))

# S. aureus
acquisition_sa <- within(acquisition_sa, cluster_lag <- relevel(as.factor(cluster_lag), ref = "COR"))
logit_sa <- glmer(inf_multi_sa ~ cluster_lag + month + inf_rv_yn + 
                    sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                  data = acquisition_sa, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_sa)
logit_sa_se <- sqrt(diag(vcov(logit_sa)))
exp(cbind(Est = fixef(logit_sa), LL = fixef(logit_sa) - 1.96 * logit_sa_se, UL = fixef(logit_sa) + 1.96*logit_sa_se))

# S. pneumoniae
acquisition_sp <- within(acquisition_sp, cluster_lag <- relevel(as.factor(cluster_lag), ref = "COR"))
logit_sp <- glmer(inf_multi_sp ~ cluster_lag + month + inf_rv_yn + 
                    sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + pcv + (1 | study_id), 
                  data = acquisition_sp, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_sp)
logit_sp_se <- sqrt(diag(vcov(logit_sp)))
exp(cbind(Est = fixef(logit_sp), LL = fixef(logit_sp) - 1.96 * logit_sp_se, UL = fixef(logit_sp) + 1.96*logit_sp_se))

# *************************************************************
# Preprocessing of data for analyses of URT microbiota features 
# *************************************************************

rf_palette <- c('#377eb8', '#e41a1c', '#984ea3', '#4daf4a', '#a65628', 'gray30')
pathogen_palette <- c("#483d8b", "#b22222", "#228b21", "#7a388a")
barplot_palette <- c("darkslateblue","seagreen")

theme_barplot <- theme(axis.text.y = element_blank(), 
                       axis.ticks.y = element_blank(), 
                       plot.title = element_text(size=8, hjust = 0.5),
                       axis.text.x = element_text(size=7.5, color="black"), 
                       axis.title.x = element_blank(),
                       legend.position = "none", 
                       axis.line.x = element_line(linewidth = 0.1, linetype = "solid", colour = "black"),
                       panel.grid.major.y = element_blank(), 
                       panel.border = element_blank(),
                       plot.margin = unit(c(0.3, 0.5, 0.2, 0), "cm"))

table(sample_names(phy.picrust)==sample_names(phy.inf.np.16s))

# Limit analyses to samples with complete infection status
tb <- table(sample_data(phy.complete)$study_id)
table(tb)
# remove subjects with only < 2 time point
sub_keep <- names(tb[tb>=2])
phy.complete <- subset_samples(phy.complete, study_id %in% sub_keep)
# check distribution of number of samples per subject
table(table(sample_data(phy.complete)$study_id))
phy.complete
# rename samples 
sample_names(phy.complete) <- gsub("-","_", sample_names(phy.complete))
metadata_complete <- data.frame(sample_data(phy.complete))
metadata_complete$SampleID <- gsub("-","_", metadata_complete$SampleID)
dim(metadata_complete)
# rename picrust samples 
sample_names(phy.picrust) <- gsub("-","_", sample_names(phy.picrust))

# Check important variables
# Infection status
var_infect <- c("inf_rv_yn", "inf_multi_hi", "inf_multi_mc", "inf_multi_sa", "inf_multi_sp")
metadata_complete[,var_infect] <- lapply(metadata_complete[,var_infect], as.factor)
summary(metadata_complete[,var_infect])
# Clinical variables
# Check if any of the clinical variables are time-dependent
var_clinical <- c("pcv", "num_kids", "mat_hiv", "breastmilk", "inf_abx_any", "season")
metadata_complete[,var_clinical[-c(1,2)]] <- lapply(metadata_complete[,var_clinical[-c(1,2)]], as.factor)
for (ii in 1:length(var_clinical)){
  meta_uni <- unique(metadata_complete[,c("study_id", var_clinical[ii])])
  print(dim(meta_uni))
}
summary(metadata_complete[,var_clinical])

# Truncate variable for number of children
table(unique(metadata_complete[, c("study_id", "num_kids")])$num_kids)
metadata_complete$num_kids_trunc <- pmin(metadata_complete$num_kids, 4)

# Save processed data
sample_data(phy.complete) <- metadata_complete
summary(metadata_complete)
saveRDS(phy.complete, file="Pixu_Analyses/phyloseq_proc.RDS")
write.csv(data.frame(sample_data(phy.complete)), file="Pixu_Analyses/metadata_proc.csv")
write.csv(otu_table(phy.complete), file="Pixu_Analyses/countdata_proc.csv")
write.csv(tax_table(phy.complete), file="Pixu_Analyses/taxdata_proc.csv")
phy.picrust <- prune_samples(sample_names(phy.complete), phy.picrust)
tab.picrust <- t(otu_table(phy.picrust))
write.csv(tab.picrust, file="Pixu_Analyses/funcdata_proc.csv")

# ********************************************
# Analyses of data for URT microbiota features 
# ********************************************

metadata <- read.csv("Pixu_Analyses/metadata_proc.csv", header=T, row.names=1)
countdata <- read.csv("Pixu_Analyses/countdata_proc.csv", header=T, row.names=1)
taxdata <- read.csv("Pixu_Analyses/taxdata_proc.csv", header=T, row.names=1)
taxdata$feature <- taxdata$ASV
metadata$study_id <- as.factor(metadata$study_id)
func_reladata <- read.csv("Pixu_Analyses/funcdata_proc.csv", header=T, row.names=1)
table(rownames(metadata)==rownames(countdata))

# Construct spline variables for month, default has degree=3 and no knots 
month_spline <- bs(metadata$month)
colnames(month_spline) <- paste0("spline", 1:ncol(month_spline))
metadata <- cbind(metadata, month_spline)

par(mfrow=c(1,3))
plot(metadata$month, month_spline[,1])
plot(metadata$month, month_spline[,2])
plot(metadata$month, month_spline[,3])

# Filter ASV data
cut_ra <- 0.001
cut_prev <- 0.05
propdata <- countdata/rowSums(countdata)
tax_keep <- apply(propdata, 2, function(x){mean(x>=cut_ra)>=cut_prev})
table(tax_keep)
count_filter <- countdata[rownames(metadata), tax_keep]
count_filter$ASV_other <- rowSums(countdata[rownames(metadata), !tax_keep])
dim(count_filter)
write.csv(count_filter, file="Pixu_Analyses/count_filter.csv")

# Filter pathway data
cut_ra <- 0.01
cut_prev <- 0.05
ind_keep_func <- colMeans(func_reladata>=cut_ra)>=cut_prev
table(ind_keep_func)
func_rela_filter <- func_reladata[rownames(metadata), ind_keep_func]
func_rela_filter$pathway_other <- rowSums(func_reladata[rownames(metadata), !ind_keep_func])
dim(func_rela_filter)
# Some pathways are 100% correlated
cormat <- cor(func_rela_filter)
ind <- which(upper.tri(cormat, diag=F) , arr.ind = TRUE )
data.frame(col = dimnames(cormat)[[2]][ind[,2]] ,
           row = dimnames(cormat)[[1]][ind[,1]] ,
           correlation = cormat[ind]) %>%
  filter(correlation==1)
ft.rm <- c("PWY0.1319", "VALSYN.PWY")
func_rela_filter <- func_rela_filter[,!colnames(func_rela_filter)%in%ft.rm]
dim(func_rela_filter)

# Log transform ASV composition 
log_comp <- count_filter
log_comp[log_comp==0] <- 0.5
log_comp <- log(log_comp/rowSums(log_comp))
# log transform functional data composition,
func_log_comp <- func_rela_filter
func_rela_zero <- func_rela_filter==0
for (ii in 1:nrow(func_log_comp)){
  func_log_comp[ii,func_rela_zero[ii,]] <-
    min(func_log_comp[ii,!func_rela_zero[ii,]])/2
}
func_log_comp <- log(func_log_comp)

# Check variables of interest
var_infect0 <- c("inf_rv_yn", "inf_multi_hi", "inf_multi_mc", "inf_multi_sa", "inf_multi_sp")
pathogen_name <- c("Respiratory virus", "H. influenzae", "M. catarrhalis", "S. aureus", "S. pneumoniae")
var_clinical <- c("num_kids_trunc", "inf_rv_yn", "mat_hiv", "breastmilk", "inf_abx_any", "season", "lbw", "wood", "residence", "sex", "pcv")

var_month <- paste0("spline", 1:3)
metadata[,var_infect0] <- lapply(metadata[,var_infect0], as.factor)
metadata[,var_clinical[3:11]] <- lapply(metadata[,var_clinical[3:11]], as.factor)
summary(metadata[,var_infect0])
summary(metadata[,var_clinical])

# Find microbiome sampleID for the preceding timepoint

# Variables as response
var_infect <- c("inf_multi_hi", "inf_multi_mc", "inf_multi_sa", "inf_multi_sp")
var_infect_prev <- paste0(c("inf_multi_hi", "inf_multi_mc", "inf_multi_sa", "inf_multi_sp"), "_prev")
var_infect_new <- c("inf_hi_new", "inf_mc_new", "inf_sa_new", "inf_sp_new")

# Find the previous time point
meta_prev <- NULL
vec_studyid <- unique(metadata$study_id)
for (ii in 1:length(vec_studyid)){
  meta_sub <- metadata %>%
    filter(study_id==vec_studyid[ii]) %>%
    dplyr::select(SampleID, study_id, month, all_of(var_infect_new), all_of(var_infect)) %>%
    arrange(month)
  tab_prev <- meta_sub[,c("SampleID", "month", var_infect)]
  colnames(tab_prev) <- paste0(colnames(tab_prev), "_prev")
  tab_prev <- tab_prev[-nrow(tab_prev),]
  rownames(tab_prev) <- meta_sub$SampleID[-1]
  meta_prev <- rbind(meta_prev, tab_prev)
}
meta_prev <- cbind(metadata[rownames(meta_prev),], meta_prev)

# Format data for random forest models
# For this analysis, we use age(t), RV(t), clinical(t) and microbiome(t-1) predicting infection(t)

var_ASV <- colnames(count_filter)
var_pathway <- colnames(func_rela_filter)

# construct data.frame for random forest prediction
data_forpred <- meta_prev[rownames(meta_prev), c("SampleID", "SampleID_prev", "study_id", "month", "month_prev", colnames(month_spline), var_clinical, var_infect, var_infect_new, var_infect_prev)]
data_forpred <- cbind(data_forpred, 
                      log_comp[data_forpred$SampleID_prev,],
                      func_log_comp[data_forpred$SampleID_prev,])

for (ii in 1:length(var_infect)){
  # only use new infections and no infections
  data_forpred_ii <- data_forpred %>% 
    filter(!!rlang::sym(var_infect_prev[ii])=="N")
  print(var_infect[ii])
  # tabulate the number of Y and N
  print(table(data_forpred_ii[,var_infect[ii]]))
  # tabulate the months
  print(table(data_forpred_ii$month))
  # tabulate the length of time between predictors and response variables
  print(table(data_forpred_ii$month-data_forpred_ii$month_prev))
}
write.csv(data_forpred, file="Pixu_Analyses/data_forpred.csv")

# ***************************************************************************************************************
# Random forest models evaluating the utility of the URT microbiota for the prediction of pathobiont colonization
# ***************************************************************************************************************

# Random forest models with all time points
# Leave-one-subject-out for prediction

# RF model 1: age ONLY

res_rf_month <- vector(length(var_infect), mode="list")
names(res_rf_month) <- var_infect
# run randomForest
for (ii in 1:length(var_infect)){
  # don't use samples that are not new infection
  data_forpred_ii <- data_forpred %>% 
    filter(!!rlang::sym(var_infect_prev[ii])=="N") %>%
    filter(month-month_prev<=2)
  fmlr <- paste(var_infect[ii], "~month")
  print(fmlr)
  # leave-one-subject-out for AUC
  set.seed(1234)
  vec_id <- unique(data_forpred_ii$study_id)
  predprob_ii <- rep(NA, nrow(data_forpred_ii))
  for (jj in 1:length(vec_id)){
    test_id <- data_forpred_ii$study_id==vec_id[jj]
    data_train <- data_forpred_ii[!test_id,]
    data_test <- data_forpred_ii[test_id,]
    res <- randomForest(as.formula(fmlr), data=data_train, strata=data_train[,var_infect[ii]],
                        sampsize=c(100,100))
    predprob_ii[test_id] <- predict(res, newdata=data_test, type = "prob")[,"Y"]
  }
  res_rf_month[[ii]] <- list(predprob=data.frame(study_id=data_forpred_ii$study_id,
                                                 SampleID=rownames(data_forpred_ii),
                                                 observed=data_forpred_ii[,var_infect[ii]],
                                                 predprob=predprob_ii))
}
saveRDS(res_rf_month, file="Pixu_Analyses/results_v2/res_rf_month_only.RDS")

# Summarize AUC and feature importance
res_rf_month <- readRDS("Pixu_Analyses/results_v2/res_rf_month_only.RDS")
roc_auc_month_rf <- matrix(NA, length(var_infect), 3)
colnames(roc_auc_month_rf) <- c("CI_lower", "AUC", "CI_upper")
rownames(roc_auc_month_rf) <- var_infect
for (ii in 1:length(var_infect)){
  res <- res_rf_month[[ii]]$predprob
  predprob <- res$predprob
  vec_id <- unique(res$study_id)
  ind_Y <- res$observed=="Y"
  set.seed(1234)
  roc_auc_month_rf[ii,] <- ci.auc(roc(cases=predprob[ind_Y],
                                      controls=predprob[!ind_Y]), 
                                  method="bootstrap")
}

roc_auc_month_rf
write.xlsx(roc_auc_month_rf, file="Pixu_Analyses/figure_table_v2/randomForest_month_only.xlsx", 
           sheetN="AUC", append=F)

# RF model 2: age + clinical covariates (includes respiratory virus infection)

res_rf_clinical <- vector(length(var_infect), mode="list")
names(res_rf_clinical) <- var_infect
# run randomForest
for (ii in 1:length(var_infect)){
  # don't use samples that are not new infection
  data_forpred_ii <- data_forpred %>% 
    filter(!!rlang::sym(var_infect_prev[ii])=="N") %>%
    filter(month-month_prev<=2)
  var_adj <- c("month", var_clinical[-length(var_clinical)])
  if(var_infect[ii]=="inf_multi_sp") var_adj <- c("month", var_clinical)
  fmlr <- paste(var_infect[ii], "~", paste(var_adj, collapse=" + "))
  print(fmlr)
  # leave-one-subject-out for AUC
  set.seed(1234)
  vec_id <- unique(data_forpred_ii$study_id)
  importance_ii <- matrix(NA, length(vec_id), length(var_adj))
  colnames(importance_ii) <- var_adj
  rownames(importance_ii) <- vec_id
  predprob_ii <- rep(NA, nrow(data_forpred_ii))
  for (jj in 1:length(vec_id)){
    test_id <- data_forpred_ii$study_id==vec_id[jj]
    data_train <- data_forpred_ii[!test_id,]
    data_test <- data_forpred_ii[test_id,]
    res <- randomForest(as.formula(fmlr), data=data_train, strata=data_train[,var_infect[ii]],
                        sampsize=c(100,100))
    predprob_ii[test_id] <- predict(res, newdata=data_test, type = "prob")[,"Y"]
    importance_ii[jj,] <- res$importance
  }
  res_rf_clinical[[ii]] <- list(predprob=data.frame(study_id=data_forpred_ii$study_id,
                                                    SampleID=rownames(data_forpred_ii),
                                                    observed=data_forpred_ii[,var_infect[ii]],
                                                    predprob=predprob_ii),
                                importance=importance_ii)
}
saveRDS(res_rf_clinical, file="Pixu_Analyses/results_v2/res_rf_clinical.RDS")

# Summarize AUC and feature importance

res_rf_clinical <- readRDS("Pixu_Analyses/results_v2/res_rf_clinical.RDS")
roc_auc_clinical <- matrix(NA, length(var_infect), 3)
colnames(roc_auc_clinical) <- c("CI_lower", "AUC", "CI_upper")
rownames(roc_auc_clinical) <- var_infect
importance_clinical <- matrix(NA, 1+length(var_clinical), length(var_infect))
rownames(importance_clinical) <- c("month", var_clinical)
colnames(importance_clinical) <- var_infect
for (ii in 1:length(var_infect)){
  res <- res_rf_clinical[[ii]]$predprob
  predprob <- res$predprob
  vec_id <- unique(res$study_id)
  ind_Y <- res$observed=="Y"
  set.seed(1234)
  roc_auc_clinical[ii,] <- ci.auc(roc(cases=predprob[ind_Y],
                                      controls=predprob[!ind_Y]), 
                                  method="bootstrap")
  if (var_infect[ii]=="inf_multi_sp"){
    importance_clinical[,ii] <- apply(res_rf_clinical[[ii]]$importance, 2, mean)
  }else{
    importance_clinical[-nrow(importance_clinical),ii] <- apply(res_rf_clinical[[ii]]$importance, 2, mean)
  }
}

roc_auc_clinical
pheatmap(importance_clinical)

write.xlsx(roc_auc_clinical, file="Pixu_Analyses/figure_table_v2/randomForest_clinical.xlsx", 
           sheetN="AUC", append=F)
write.xlsx(importance_clinical, file="Pixu_Analyses/figure_table_v2/randomForest_clinical.xlsx", 
           sheetN="Average Importance", append=T)

# RF model 3: microbiome ONLY

res_rf_asv <- vector(length(var_infect), mode="list")
names(res_rf_asv) <- var_infect
# run randomForest
for (ii in 1:length(var_infect)){
  # don't use samples that are not new infection
  data_forpred_ii <- data_forpred %>% 
    filter(!!rlang::sym(var_infect_prev[ii])=="N") %>%
    filter(month-month_prev<=2)
  fmlr <- paste(var_infect[ii], "~", paste(var_ASV, collapse=" + "))
  print(fmlr)
  # leave-one-subject-out for AUC
  set.seed(1234)
  vec_id <- unique(data_forpred_ii$study_id)
  importance_ii <- matrix(NA, length(vec_id), length(var_ASV))
  colnames(importance_ii) <- var_ASV
  rownames(importance_ii) <- vec_id
  predprob_ii <- rep(NA, nrow(data_forpred_ii))
  for (jj in 1:length(vec_id)){
    test_id <- data_forpred_ii$study_id==vec_id[jj]
    data_train <- data_forpred_ii[!test_id,]
    data_test <- data_forpred_ii[test_id,]
    res <- randomForest(as.formula(fmlr), data=data_train, strata=data_train[,var_infect[ii]],
                        sampsize=c(100,100))
    predprob_ii[test_id] <- predict(res, newdata=data_test, type = "prob")[,"Y"]
    importance_ii[jj,] <- res$importance
  }
  res_rf_asv[[ii]] <- list(predprob=data.frame(study_id=data_forpred_ii$study_id,
                                               SampleID=rownames(data_forpred_ii),
                                               observed=data_forpred_ii[,var_infect[ii]],
                                               predprob=predprob_ii),
                           importance=importance_ii)
}
saveRDS(res_rf_asv, file="Pixu_Analyses/results_v2/res_rf_asv.RDS")

# Summarize AUC and feature importance

res_rf_asv <- readRDS("Pixu_Analyses/results_v2/res_rf_asv.RDS")
roc_auc_asv <- matrix(NA, length(var_infect), 3)
colnames(roc_auc_asv) <- c("CI_lower", "AUC", "CI_upper")
rownames(roc_auc_asv) <- var_infect
importance_asv <- matrix(NA, length(var_ASV), length(var_infect))
rownames(importance_asv) <- var_ASV
colnames(importance_asv) <- var_infect
for (ii in 1:length(var_infect)){
  res <- res_rf_asv[[ii]]$predprob
  predprob <- res$predprob
  vec_id <- unique(res$study_id)
  ind_Y <- res$observed=="Y"
  set.seed(1234)
  roc_auc_asv[ii,] <- ci.auc(roc(cases=predprob[ind_Y],
                                 controls=predprob[!ind_Y]), 
                             method="bootstrap")
  importance_asv[,ii] <- apply(res_rf_asv[[ii]]$importance, 2, mean)
}

roc_auc_asv
pheatmap(importance_asv)

write.xlsx(roc_auc_asv, file="Pixu_Analyses/figure_table_v2/randomForest_asv.xlsx", 
           sheetN="AUC", append=F)
write.xlsx(importance_asv, file="Pixu_Analyses/figure_table_v2/randomForest_asv.xlsx", 
           sheetN="Average Importance", append=T)

# RF model 4: microbiome + functional pathway data

res_rf_asv_pathway <- vector(length(var_infect), mode="list")
names(res_rf_asv_pathway) <- var_infect
# run randomForest
for (ii in 1:length(var_infect)){
  # don't use samples that are not new infection
  data_forpred_ii <- data_forpred %>% 
    filter(!!rlang::sym(var_infect_prev[ii])=="N") %>%
    filter(month-month_prev<=2)
  fmlr <- paste(var_infect[ii], "~", paste(c(var_ASV, var_pathway), collapse=" + "))
  print(fmlr)
  # leave-one-subject-out for AUC
  set.seed(1234)
  vec_id <- unique(data_forpred_ii$study_id)
  importance_ii <- matrix(NA, length(vec_id), length(var_ASV)+length(var_pathway))
  colnames(importance_ii) <- c(var_ASV, var_pathway)
  rownames(importance_ii) <- vec_id
  predprob_ii <- rep(NA, nrow(data_forpred_ii))
  for (jj in 1:length(vec_id)){
    test_id <- data_forpred_ii$study_id==vec_id[jj]
    data_train <- data_forpred_ii[!test_id,]
    data_test <- data_forpred_ii[test_id,]
    res <- randomForest(as.formula(fmlr), data=data_train, strata=data_train[,var_infect[ii]],
                        sampsize=c(100,100))
    predprob_ii[test_id] <- predict(res, newdata=data_test, type = "prob")[,"Y"]
    importance_ii[jj,] <- res$importance
  }
  res_rf_asv_pathway[[ii]] <- list(predprob=data.frame(study_id=data_forpred_ii$study_id,
                                                       SampleID=rownames(data_forpred_ii),
                                                       observed=data_forpred_ii[,var_infect[ii]],
                                                       predprob=predprob_ii),
                                   importance=importance_ii)
}
saveRDS(res_rf_asv_pathway, file="Pixu_Analyses/results_v2/res_rf_asv_pathway.RDS")

# Summarize AUC and feature importance

res_rf_asv_pathway <- readRDS("Pixu_Analyses/results_v2/res_rf_asv_pathway.RDS")
roc_auc_asv_pathway <- matrix(NA, length(var_infect), 3)
colnames(roc_auc_asv_pathway) <- c("CI_lower", "AUC", "CI_upper")
rownames(roc_auc_asv_pathway) <- var_infect
importance_asv_pathway <- matrix(NA, length(var_ASV)+length(var_pathway), length(var_infect))
rownames(importance_asv_pathway) <- c(var_ASV, var_pathway)
colnames(importance_asv_pathway) <- var_infect
for (ii in 1:length(var_infect)){
  res <- res_rf_asv_pathway[[ii]]$predprob
  predprob <- res$predprob
  vec_id <- unique(res$study_id)
  ind_Y <- res$observed=="Y"
  set.seed(1234)
  roc_auc_asv_pathway[ii,] <- ci.auc(roc(cases=predprob[ind_Y],
                                         controls=predprob[!ind_Y]), 
                                     method="bootstrap")
  importance_asv_pathway[,ii] <-
    apply(res_rf_asv_pathway[[ii]]$importance, 2, mean)
}

roc_auc_asv_pathway
pheatmap(importance_asv_pathway)

write.xlsx(roc_auc_asv_pathway, file="Pixu_Analyses/figure_table_v2/randomForest_asv_pathway.xlsx", 
           sheetN="AUC", append=F)
write.xlsx(importance_asv_pathway, file="Pixu_Analyses/figure_table_v2/randomForest_asv_pathway.xlsx", 
           sheetN="Average Importance", append=T)

# RF model 5: microbiome + clinical variables (including age and respiratory virus infection)

res_rf_asv_clinical <- vector(length(var_infect), mode="list")
names(res_rf_asv_clinical) <- var_infect
# run randomForest
for (ii in 1:length(var_infect)){
  # don't use samples that are not new infection
  data_forpred_ii <- data_forpred %>% 
    filter(!!rlang::sym(var_infect_prev[ii])=="N") %>%
    filter(month-month_prev<=2)
  var_adj <- c("month", var_clinical[-length(var_clinical)])
  if(var_infect[ii]=="inf_multi_sp") var_adj <- c("month", var_clinical)
  fmlr <- paste(var_infect[ii], "~", paste(c(var_adj, var_ASV), collapse=" + "))
  print(fmlr)
  # leave-one-subject-out for AUC
  set.seed(1234)
  vec_id <- unique(data_forpred_ii$study_id)
  importance_ii <- matrix(NA, length(vec_id), length(var_adj)+length(var_ASV))
  colnames(importance_ii) <- c(var_adj, var_ASV)
  rownames(importance_ii) <- vec_id
  predprob_ii <- rep(NA, nrow(data_forpred_ii))
  for (jj in 1:length(vec_id)){
    test_id <- data_forpred_ii$study_id==vec_id[jj]
    data_train <- data_forpred_ii[!test_id,]
    data_test <- data_forpred_ii[test_id,]
    res <- randomForest(as.formula(fmlr), data=data_train, strata=data_train[,var_infect[ii]],
                        sampsize=c(100,100))
    predprob_ii[test_id] <- predict(res, newdata=data_test, type = "prob")[,"Y"]
    importance_ii[jj,] <- res$importance
  }
  res_rf_asv_clinical[[ii]] <- list(predprob=data.frame(study_id=data_forpred_ii$study_id,
                                                        SampleID=rownames(data_forpred_ii),
                                                        observed=data_forpred_ii[,var_infect[ii]],
                                                        predprob=predprob_ii),
                                    importance=importance_ii)
}
saveRDS(res_rf_asv_clinical, file="Pixu_Analyses/results_v2/res_rf_asv_clinical.RDS")

# Summarize AUC and feature importance

res_rf_asv_clinical <- readRDS("Pixu_Analyses/results_v2/res_rf_asv_clinical.RDS")
roc_auc_asv_clinical <- matrix(NA, length(var_infect), 3)
colnames(roc_auc_asv_clinical) <- c("CI_lower", "AUC", "CI_upper")
rownames(roc_auc_asv_clinical) <- var_infect
importance_asv_clinical <- matrix(NA, 1+length(var_clinical)+length(var_ASV), length(var_infect))
rownames(importance_asv_clinical) <- c("month", var_clinical, var_ASV)
colnames(importance_asv_clinical) <- var_infect
for (ii in 1:length(var_infect)){
  res <- res_rf_asv_clinical[[ii]]$predprob
  predprob <- res$predprob
  vec_id <- unique(res$study_id)
  ind_Y <- res$observed=="Y"
  set.seed(1234)
  roc_auc_asv_clinical[ii,] <- ci.auc(roc(cases=predprob[ind_Y],
                                          controls=predprob[!ind_Y]), 
                                      method="bootstrap")
  if (var_infect[ii]=="inf_multi_sp"){
    importance_asv_clinical[,ii] <- apply(res_rf_asv_clinical[[ii]]$importance, 2, mean)
  }else{
    importance_asv_clinical[-length(var_clinical)-1,ii] <- apply(res_rf_asv_clinical[[ii]]$importance, 2, mean)
  }
}

roc_auc_asv_clinical
pheatmap(importance_asv_clinical)

write.xlsx(roc_auc_asv_clinical, file="Pixu_Analyses/figure_table_v2/randomForest_asv_clinical.xlsx", 
           sheetN="AUC", append=F)
write.xlsx(importance_asv_clinical, file="Pixu_Analyses/figure_table_v2/randomForest_asv_clinical.xlsx", 
           sheetN="Average Importance", append=T)

# RF model 6: clinical variables (including age and respiratory virus infection), microbiome, and pathway data (kitchen sink model)

res_rf_all <- vector(length(var_infect), mode="list")
names(res_rf_all) <- var_infect
# run randomForest
for (ii in 1:length(var_infect)){
  # don't use samples that are not new infection
  data_forpred_ii <- data_forpred %>% 
    filter(!!rlang::sym(var_infect_prev[ii])=="N") %>%
    filter(month-month_prev<=2)
  var_adj <- c("month", var_clinical[-length(var_clinical)])
  if(var_infect[ii]=="inf_multi_sp") var_adj <- c("month", var_clinical)
  fmlr <- paste(var_infect[ii], "~", paste(c(var_adj, var_ASV, var_pathway), collapse=" + "))
  print(fmlr)
  # leave-one-subject-out for AUC
  set.seed(1234)
  vec_id <- unique(data_forpred_ii$study_id)
  importance_ii <- matrix(NA, length(vec_id), length(var_adj)+length(var_ASV)+length(var_pathway))
  colnames(importance_ii) <- c(var_adj, var_ASV, var_pathway)
  rownames(importance_ii) <- vec_id
  predprob_ii <- rep(NA, nrow(data_forpred_ii))
  for (jj in 1:length(vec_id)){
    test_id <- data_forpred_ii$study_id==vec_id[jj]
    data_train <- data_forpred_ii[!test_id,]
    data_test <- data_forpred_ii[test_id,]
    res <- randomForest(as.formula(fmlr), data=data_train, strata=data_train[,var_infect[ii]],
                        sampsize=c(100,100))
    predprob_ii[test_id] <- predict(res, newdata=data_test, type = "prob")[,"Y"]
    importance_ii[jj,] <- res$importance
  }
  res_rf_all[[ii]] <- list(predprob=data.frame(study_id=data_forpred_ii$study_id,
                                               SampleID=rownames(data_forpred_ii),
                                               observed=data_forpred_ii[,var_infect[ii]],
                                               predprob=predprob_ii),
                           importance=importance_ii)
}
saveRDS(res_rf_all, file="Pixu_Analyses/results_v2/res_rf_all.RDS")

# Summarize AUC and feature importance

res_rf_all <- readRDS("Pixu_Analyses/results_v2/res_rf_all.RDS")
roc_auc_all <- matrix(NA, length(var_infect), 3)
colnames(roc_auc_all) <- c("CI_lower", "AUC", "CI_upper")
rownames(roc_auc_all) <- var_infect
importance_all <- matrix(NA,
                         1+length(var_clinical)+length(var_ASV)+length(var_pathway),
                         length(var_infect))
rownames(importance_all) <- c("month", var_clinical, var_ASV, var_pathway)
colnames(importance_all) <- var_infect
for (ii in 1:length(var_infect)){
  res <- res_rf_all[[ii]]$predprob
  predprob <- res$predprob
  vec_id <- unique(res$study_id)
  ind_Y <- res$observed=="Y"
  set.seed(1234)
  roc_auc_all[ii,] <- ci.auc(roc(cases=predprob[ind_Y],
                                 controls=predprob[!ind_Y]), 
                             method="bootstrap")
  if (var_infect[ii]=="inf_multi_sp"){
    importance_all[,ii] <- apply(res_rf_all[[ii]]$importance, 2, mean)
  }else{
    importance_all[-length(var_clinical)-1,ii] <- apply(res_rf_all[[ii]]$importance, 2, mean)
  }
}

roc_auc_all
pheatmap(importance_all)

write.xlsx(roc_auc_all, file="Pixu_Analyses/figure_table_v2/randomForest_all.xlsx", 
           sheetN="AUC", append=F)
write.xlsx(importance_all, file="Pixu_Analyses/figure_table_v2/randomForest_all.xlsx", 
           sheetN="Average Importance", append=T)

# *************************************************************************************
# Association between individual ASVs and acquisition of specific bacterial pathobionts
# Mixed-effect logistic regression models adjusting for month variables only
# *************************************************************************************

# infect(t) as outcome, and age(t) and ASV(t-1) as predictor (this is a supplement to the random forest analysis)

# Centralize the microbial features to improve convergence
data_forpred_center <- data_forpred
data_forpred_center[,var_ASV] <- scale(data_forpred[,var_ASV], center=TRUE, scale=FALSE)
res_all <- vector(mode='list', length=4)
names(res_all) <- var_infect
for (ii in 1:length(var_infect)){
  # don't use samples that are not new infection
  data_forpred_ii <- data_forpred_center %>% 
    filter(!!rlang::sym(var_infect_prev[ii])=="N") %>%
    filter(month-month_prev<=2)
  var_adj <- var_month
  res_ASV <- matrix(NA, length(var_ASV), 5)
  colnames(res_ASV) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "convergence")
  rownames(res_ASV) <- var_ASV
  for (jj in 1:length(var_ASV)){
    fmlr <- paste(var_infect[ii], "~", paste(c(var_ASV[jj], var_adj), collapse=" + "), "+ (1|study_id)")
    print(fmlr)
    res <- summary(glmer(fmlr, data=data_forpred_ii, family="binomial", control=glmerControl(optimizer="nlminbwrap")))
    res_ASV[jj,-5] <- res$coefficients[2,]
  }
  res_ASV <- as.data.frame(res_ASV)
  # FDR adjusted p-value
  res_ASV$qvalue <- p.adjust(res_ASV$`Pr(>|z|)`, method="fdr")
  res_all[[ii]] <- res_ASV
}
saveRDS(res_all, file="Pixu_Analyses/results_v2/glmer_ASV_acquisition_adj_month.RDS")

# Save result in spreadsheet and print significant variables

res_all <- readRDS("Pixu_Analyses/results_v2/glmer_ASV_acquisition_adj_month.RDS")
for (ii in 1:length(var_infect)){
  write.xlsx(res_all[[ii]], file="Pixu_Analyses/figure_table_v2/glmer_ASV_acquisition_adj_month.xlsx", 
             sheetN=var_infect[ii], append=ifelse(ii==1,F,T))
}

for (ii in 1:length(var_infect)){
  res <- res_all[[ii]] %>%
    dplyr::arrange(qvalue) %>%
    dplyr::filter(qvalue<=0.2)
  res <- cbind(taxdata[rownames(res),], res)
  print(kable(res[,-c(1:3)], caption=paste("For variable", var_infect[ii])))
  cat("\n")
}

# ***************************************************************************************************
# Analyses of the impact of respiratory virus infection or pathogen acquisition on the URT microbiota
# ***************************************************************************************************

# What is the impact of respiratory virus detection or acquisition of bacterial respiratory pathogens on ASVs within the nasopharyngeal microbiota?
# Maaslin2: ASV(t) as outcome, infect(t) as predictor

set.seed(1234)
for (ii in 1:length(var_infect0)){
  outdir <- file.path("Pixu_Analyses/results_v2/maaslin_impact", var_infect0[ii])
  if (!dir.exists(outdir)){
    dir.create(outdir)
  }else{
    file.remove(file.path(outdir, list.files(outdir)), recursive=TRUE)
  }
  var_adj <- var_clinical[-11]
  if(var_infect0[ii]=="inf_multi_sp") var_adj <- var_clinical
  fit_data <- Maaslin2(count_filter, 
                       metadata, 
                       output=outdir, 
                       min_abundance = 0,
                       min_prevalence = 0,
                       fixed_effects = c(var_infect0[ii], 
                                         paste0("spline", 1:3), 
                                         var_adj),
                       random_effects = c("study_id"),
                       max_significance = 0.20,
                       plot_scatter = FALSE)
  res <- fit_data$results %>% 
    mutate(coef_dir=ifelse(coef>0, "Pos", "Neg"))
  res <- merge(taxdata, res,  by="feature") %>%
    dplyr::filter(qval<=0.20)
  write.csv(res, file=file.path("Pixu_Analyses/results_v2/maaslin_impact", var_infect0[ii], "res_w_taxa.csv"))
}
for (ii in 1:length(var_infect0)){
  res <- read.csv(file.path("Pixu_Analyses/results_v2/maaslin_impact", var_infect0[ii], "res_w_taxa.csv"), header=TRUE, row.names=1)
  res <- res %>%
    dplyr::filter(metadata==var_infect0[ii], qval<=0.20) %>%
    dplyr::arrange(qval)
  print(kable(res[,c(5:9, 15, 19, 12:14, 16:18)], caption=paste("For variable", var_infect0[ii])))
  cat("\n")
}