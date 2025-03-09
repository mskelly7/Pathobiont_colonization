# Botswana Infant Microbiome Study - RV-Bacterial Pathobiont Analyses
# Matthew Kelly, MD, MPH 
# Analyses of respiratory virus and bacterial pathobiont PCR data
# Last update: March 9, 2025

remove(list=ls())
setwd("_______________")
set.seed(1234)

version
library(tidyverse)
library(dplyr)
library(plyr)
library(data.table)
library(gridExtra)
library(reshape2)
library(survminer)
library(DataCombine)
library(ggsignif)
library(grid)
library(lme4)
packageVersion("lme4")

metadata_inf_np <- read.csv("metadata_inf_np_RV.csv")
nrow(metadata_inf_np)
metadata_inf_np$inf_rv_cat3[metadata_inf_np$inf_rv_yn=="Y" & metadata_inf_np$current_uri=="Y"] <- "Symptomatic"
metadata_inf_np$inf_rv_cat3[metadata_inf_np$inf_rv_yn=="Y" & metadata_inf_np$current_uri=="N"] <- "Asymptomatic"
metadata_inf_np$inf_rv_cat3[metadata_inf_np$inf_rv_yn=="N"] <- "None"
metadata_inf_np$inf_rv_cat3 <- as.factor(metadata_inf_np$inf_rv_cat3)
table(metadata_inf_np$inf_rv_cat3, useNA="always")

# Create dataset for summarizing infant characteristics
infants <- metadata_inf_np %>% group_by(study_id) %>% filter(row_number(month) == 1)
nrow(infants)

# Create datasets for analyses of viral and bacterial PCR data
metadata_rv <- subset(metadata_inf_np, !is.na(inf_rv_yn))
nrow(metadata_rv)
metadata_bact <- subset(metadata_inf_np, !is.na(inf_multi_hi) & !is.na(inf_multi_mc) & !is.na(inf_multi_sa) & !is.na(inf_multi_sp))
nrow(metadata_bact)
metadata_bact_noM0 <- subset(metadata_bact, month!="0")
nrow(metadata_bact_noM0)
metadata_rv_bact <- subset(metadata_inf_np, !is.na(inf_rv_yn) & !is.na(inf_multi_hi) & !is.na(inf_multi_mc) & 
                             !is.na(inf_multi_sa) & !is.na(inf_multi_sp))
nrow(metadata_rv_bact)

# Create datasets for analyses of pathobiont acquisition
acquisition <- metadata_rv_bact
acquisition <- acquisition[order(acquisition$study_id, acquisition$month),]
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
nrow(acquisition)
acquisition <- acquisition[,c("SampleID","study_id","age_days","month","inf_rv_yn","inf_rv_cat2","inf_rv_cat3","inf_multi_hi","hi_lag","inf_multi_mc","mc_lag",
                              "inf_multi_sa","sa_lag","inf_multi_sp","sp_lag", "sex", "lbw", "mat_hiv", "residence", "wood", "num_kids", "season", 
                              "breastmilk", "inf_abx_any", "pcv", "hib")]
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

# ********************************************************************
# GENERATE SUMMARY DATA ON RESPIRATORY VIRUS INFECTIONS IN THIS COHORT
# ********************************************************************

table(metadata_rv$inf_rv_yn)
prop.table(table(metadata_rv$inf_rv_yn))
table(metadata_rv$inf_rv_cat)
prop.table(table(metadata_rv$inf_rv_cat))
table(metadata_rv$inf_rv_cat, metadata_rv$month)
prop.table(table(metadata_rv$inf_rv_cat, metadata_rv$month),2)
table(metadata_rv$inf_rv_yn, metadata_rv$month)
prop.table(table(metadata_rv$inf_rv_yn, metadata_rv$month),2)

metadata_rv$inf_rv_pos[metadata_rv$inf_rv_yn=="N"] <- 0
metadata_rv$inf_rv_pos[metadata_rv$inf_rv_yn=="Y"] <- 1
logit_rv_age <- glm(inf_rv_pos ~ month, data=metadata_rv, family="binomial")
summary(logit_rv_age)
confint(logit_rv_age)
# Prevalence of respiratory virus infection increases with increasing age (p<0.0001)

# What proportion of respiratory virus infections were associated with current RTI symptoms? 
table(metadata_rv$current_uri, metadata_rv$inf_rv_cat)
prop.table(table(metadata_rv$inf_rv_yn, metadata_rv$current_uri),1)
chisq.test(metadata_rv$inf_rv_yn, metadata_rv$current_uri)
prop.table(table(metadata_rv$inf_rv_cat, metadata_rv$current_uri),1)
prop.table(table(metadata_rv$hmpv, metadata_rv$current_uri),1)
prop.table(table(metadata_rv$flu, metadata_rv$current_uri),1)
prop.table(table(metadata_rv$rsv, metadata_rv$current_uri),1)
# Current RTI symptoms in >70% of flu, HMPV, and RSV infections
prop.table(table(metadata_rv$rhino, metadata_rv$current_uri),1)
prop.table(table(metadata_rv$adeno, metadata_rv$current_uri),1)
# Current RTI symptoms in <50% of rhino and adeno infections
prop.table(table(metadata_rv$paraflu, metadata_rv$current_uri),1)
prop.table(table(metadata_rv$sars2, metadata_rv$current_uri),1)

# *************************************************************
# GENERATE SUMMARY DATA ON BACTERIAL PATHOBIONTS IN THIS COHORT
# *************************************************************

table(metadata_bact$inf_multi_sp)
prop.table(table(metadata_bact$inf_multi_sp))
table(metadata_bact$inf_multi_hi)
prop.table(table(metadata_bact$inf_multi_hi))
table(metadata_bact$inf_multi_mc)
prop.table(table(metadata_bact$inf_multi_mc))
table(metadata_bact$inf_multi_sa)
prop.table(table(metadata_bact$inf_multi_sa))

# Logistic regression models to evaluate associations between infant age and the prevalence of bacterial pathobiont colonization

metadata_bact$inf_multi_hi_pos[metadata_bact$inf_multi_hi=="N"] <- 0
metadata_bact$inf_multi_hi_pos[metadata_bact$inf_multi_hi=="Y"] <- 1
logit_hi_age <- glm(inf_multi_hi_pos ~ month, data=metadata_bact, family="binomial")
summary(logit_hi_age)
confint(logit_hi_age)

metadata_bact$inf_multi_mc_pos[metadata_bact$inf_multi_mc=="N"] <- 0
metadata_bact$inf_multi_mc_pos[metadata_bact$inf_multi_mc=="Y"] <- 1
logit_mc_age <- glm(inf_multi_mc_pos ~ month, data=metadata_bact, family="binomial")
summary(logit_mc_age)
confint(logit_mc_age)

metadata_bact$inf_multi_sp_pos[metadata_bact$inf_multi_sp=="N"] <- 0
metadata_bact$inf_multi_sp_pos[metadata_bact$inf_multi_sp=="Y"] <- 1
logit_sp_age <- glm(inf_multi_sp_pos ~ month, data=metadata_bact, family="binomial")
summary(logit_sp_age)
confint(logit_sp_age)

# For Staph aureus, focusing on decline in prevalence following peak at 1 month of age
metadata_bact$inf_multi_sa_pos[metadata_bact$inf_multi_sa=="N"] <- 0
metadata_bact$inf_multi_sa_pos[metadata_bact$inf_multi_sa=="Y"] <- 1
metadata_bact_noM0$inf_multi_sa_pos[metadata_bact_noM0$inf_multi_sa=="N"] <- 0
metadata_bact_noM0$inf_multi_sa_pos[metadata_bact_noM0$inf_multi_sa=="Y"] <- 1
logit_sa_age <- glm(inf_multi_sa_pos ~ month, data=metadata_bact_noM0, family="binomial")
summary(logit_sa_age)
confint(logit_sa_age)

# Mixed effect logistic regression models to evaluate associations between respiratory virus infection and bacterial pathobiont ACQUISITION and RTI symptoms
metadata_rv_bact$current_rti[metadata_rv_bact$current_uri=="N"] <- 0
metadata_rv_bact$current_rti[metadata_rv_bact$current_uri=="Y"] <- 1
logit_rti_acq <- glmer(current_rti ~ inf_rv_yn + inf_hi_new + inf_mc_new + inf_sa_new + inf_sp_new + month + (1 | study_id), 
                   data = metadata_rv_bact, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_rti_acq)
logit_rti_acq_se <- sqrt(diag(vcov(logit_rti_acq)))
exp(cbind(Est = fixef(logit_rti_acq), LL = fixef(logit_rti_acq) - 1.96 * logit_rti_acq_se, UL = fixef(logit_rti_acq) + 1.96*logit_rti_acq_se))

# Mixed effect logistic regression model to evaluate for an association with pneumococcal colonization density and RTI symptoms
metadata_rv_bact_sp <- subset(metadata_rv_bact, inf_sp_new=="Y")
logit_rti_sp <- glmer(current_rti ~ inf_rv_yn + inf_hi_new + inf_mc_new + inf_sa_new + inf_sp + month + (1 | study_id), 
                   data = metadata_rv_bact_sp, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_rti_sp)
logit_rti_sp_se <- sqrt(diag(vcov(logit_rti_sp)))
exp(cbind(Est = fixef(logit_rti_sp), LL = fixef(logit_rti_sp) - 1.96 * logit_rti_sp_se, UL = fixef(logit_rti_sp) + 1.96*logit_rti_sp_se))

# ***************************************************************************************
# ANALYSES EVALUATING FOR ASSOCIATIONS WITH THE ODDS OF BACTERIAL PATHOBIONT COLONIZATION
# ***************************************************************************************

# Mixed effect logistic regression models evaluating associations between infant characteristics, virus detection, and pathobiont colonization on odds of pathobiont acquisition
# Infant characteristics evaluated were age (visit month), season, maternal HIV infection, number of kids in the household, breastfeeding, antibiotics
# Analyses for S. pneumoniae acquisition were additionally adjusted for PCV-13 doses 

logit_hi <- glmer(inf_multi_hi ~ inf_rv_yn + mc_lag + sa_lag + sp_lag + month + 
                    sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                      data = acquisition_hi, family = binomial, 
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_hi)
logit_hi_se <- sqrt(diag(vcov(logit_hi)))
exp(cbind(Est = fixef(logit_hi), LL = fixef(logit_hi) - 1.96 * logit_hi_se, UL = fixef(logit_hi) + 1.96*logit_hi_se))
# Respiratory virus infection and preceding M. catarrhalis/S. pneumoniae colonization are associated with higher odds of H. influenzae acquisition
# A higher number of children in the household is associated with an increased odds of H. influenzae acquisition

logit_mc <- glmer(inf_multi_mc ~ inf_rv_yn + hi_lag + sa_lag + sp_lag + month + 
                    sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                  data = acquisition_mc, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_mc)
logit_mc_se <- sqrt(diag(vcov(logit_mc)))
exp(cbind(Est = fixef(logit_mc), LL = fixef(logit_mc) - 1.96 * logit_mc_se, UL = fixef(logit_mc) + 1.96*logit_mc_se))
# Trend toward increased odds of M. catarrhalis acquisition with respiratory virus infection
# No association between preceding pathobiont colonization and the odds of M. catarrhalis acquisition
# Lower risk of M. catarrhalis acquisition during the rainy (warm) season

logit_sa <- glmer(inf_multi_sa ~ inf_rv_yn + hi_lag + mc_lag + sp_lag + month + 
                    sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                  data = acquisition_sa, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_sa)
logit_sa_se <- sqrt(diag(vcov(logit_sa)))
exp(cbind(Est = fixef(logit_sa), LL = fixef(logit_sa) - 1.96 * logit_sa_se, UL = fixef(logit_sa) + 1.96*logit_sa_se))
# Respiratory virus infection is associated with a lower odds of S. aureus acquisition
# Trend toward lower odds of S. aureus acquisition with preceding S. pneumoniae acquisition
# Higher odds of S. aureus acquisition during the rainy season
# Lower odds of S. aureus acquisition with breastfeeding, receipt antibiotic treatment 

logit_sp <- glmer(inf_multi_sp ~ inf_rv_yn + hi_lag + mc_lag + sa_lag + month + 
                    sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + pcv + (1 | study_id), 
                  data = acquisition_sp, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_sp)
logit_sp_se <- sqrt(diag(vcov(logit_sp)))
exp(cbind(Est = fixef(logit_sp), LL = fixef(logit_sp) - 1.96 * logit_sp_se, UL = fixef(logit_sp) + 1.96*logit_sp_se))
# Respiratory virus infection and M. catarrhalis colonization associated with higher odds of S. pneumoniae acquisition
# Trend toward higher odds of S. pneumoniae acquisition with preceding H. influenzae acquisition
# A higher number of children in the household and urban residence is associated with higher odds of S. pneumoniae acquisition
# No association between PCV-13 vaccination status and the odds of S. pneumoniae acquisition

# Pneumococcal colonization density among infants with S. pneumoniae carriage
# Linear regression model evaluating associations between infant characteristics, viral infections, and other pathobiont colonization on pneumococcal colonization density
metadata_rv_sp <- subset(metadata_rv_bact, inf_sp>0 & inf_multi_sp=="Y")
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_rv_yn, summary)
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_multi_hi, summary)
summary(lm(inf_sp ~ inf_rv_yn + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + 
             pcv, data=metadata_rv_sp))
# Respiratory virus detection is associated with higher pneumococcal colonization density in samples with Sp colonization 
# Higher pneumococcal colonization density with H. influenzae colonization, decreasing age, a higher number of children in the household, antibiotic exposure
# No association with PCV-13 vaccination status

# Analyses categorizing respiratory viruses as involving rhino/entero only vs. other respiratory virus infections

# H. influenzae - respiratory virus infections (rhino/entero vs. other) 
logit_hi2 <- glmer(inf_multi_hi ~ inf_rv_cat2 + mc_lag + sa_lag + sp_lag + month + 
                    sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                  data = acquisition_hi, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_hi2)
logit_hi2_se <- sqrt(diag(vcov(logit_hi2)))
exp(cbind(Est = fixef(logit_hi2), LL = fixef(logit_hi2) - 1.96 * logit_hi2_se, UL = fixef(logit_hi2) + 1.96*logit_hi2_se))
# Higher odds of H. influenzae acquisition with "other" respiratory virus infections 

# M. catarrhalis - respiratory virus infections (rhino/entero vs. other) 
logit_mc2 <- glmer(inf_multi_mc ~ inf_rv_cat2 + hi_lag + sa_lag + sp_lag + month + 
                     sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                   data = acquisition_mc, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_mc2)
logit_mc2_se <- sqrt(diag(vcov(logit_mc2)))
exp(cbind(Est = fixef(logit_mc2), LL = fixef(logit_mc2) - 1.96 * logit_mc2_se, UL = fixef(logit_mc2) + 1.96*logit_mc2_se))
# Higher odds of M. catarrhalis acquisition with "other" respiratory virus infections 

# S. aureus - respiratory virus infections (rhino/entero vs. other) 
logit_sa2 <- glmer(inf_multi_sa ~ inf_rv_cat2 + hi_lag + mc_lag + sp_lag + month + 
                     sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                   data = acquisition_sa, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_sa2)
logit_sa2_se <- sqrt(diag(vcov(logit_sa2)))
exp(cbind(Est = fixef(logit_sa2), LL = fixef(logit_sa2) - 1.96 * logit_sa2_se, UL = fixef(logit_sa2) + 1.96*logit_sa2_se))
# Both rhino/entero and other respiratory viruses associated with lower odds of S. aureus acquisition

# S. pneumoniae - respiratory virus infections (rhino/entero vs. other) 
logit_sp2 <- glmer(inf_multi_sp ~ inf_rv_cat2 + hi_lag + mc_lag + sa_lag + month + 
                     sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + pcv + (1 | study_id), 
                   data = acquisition_sp, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_sp2)
logit_sp2_se <- sqrt(diag(vcov(logit_sp2)))
exp(cbind(Est = fixef(logit_sp2), LL = fixef(logit_sp2) - 1.96 * logit_sp2_se, UL = fixef(logit_sp2) + 1.96*logit_sp2_se))
# Both rhino/entero and other respiratory virus infections associated with an increased odds of S. pneumoniae colonization

# PNEUMOCOCCAL COLONIZATION DENSITY
table(metadata_rv_sp$inf_rv_cat2)
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_rv_cat2, summary)
metadata_rv_sp$inf_rv_cat2 <- as.factor(metadata_rv_sp$inf_rv_cat2)
summary(lm(inf_sp ~ relevel(inf_rv_cat2, ref="Other") + inf_multi_hi + inf_multi_mc + inf_multi_sa + 
             month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))
# No difference in pneumococcal colonization density comparing infections involving rhino/entero only vs. infections with other viruses

# Analyses categorizing respiratory viruses based on the presence of symptoms

table(acquisition$inf_rv_cat3, useNA="always")

# H. influenzae - respiratory virus infections (symptomatic vs. asymptomatic) 
logit_hi3 <- glmer(inf_multi_hi ~ relevel(inf_rv_cat3, ref="None") + mc_lag + sa_lag + sp_lag + month + 
                     sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                   data = acquisition_hi, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_hi3)
logit_hi3_se <- sqrt(diag(vcov(logit_hi3)))
exp(cbind(Est = fixef(logit_hi3), LL = fixef(logit_hi3) - 1.96 * logit_hi3_se, UL = fixef(logit_hi3) + 1.96*logit_hi3_se))
# Only symptomatic respiratory virus infection associated with increased odds of H. influenzae acquisition

# M. catarrhalis - respiratory virus infections (symptomatic vs. asymptomatic)  
logit_mc3 <- glmer(inf_multi_mc ~ relevel(inf_rv_cat3, ref="None") + hi_lag + sa_lag + sp_lag + month + 
                     sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                   data = acquisition_mc, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_mc3)
logit_mc3_se <- sqrt(diag(vcov(logit_mc3)))
exp(cbind(Est = fixef(logit_mc3), LL = fixef(logit_mc3) - 1.96 * logit_mc3_se, UL = fixef(logit_mc3) + 1.96*logit_mc3_se))
# Only symptomatic respiratory virus infection associated with increased odds of M. catarrhalis acquisition

# S. aureus - respiratory virus infections (symptomatic vs. asymptomatic) 
logit_sa3 <- glmer(inf_multi_sa ~ relevel(inf_rv_cat3, ref="None") + hi_lag + mc_lag + sp_lag + month + 
                     sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + (1 | study_id), 
                   data = acquisition_sa, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_sa3)
logit_sa3_se <- sqrt(diag(vcov(logit_sa3)))
exp(cbind(Est = fixef(logit_sa3), LL = fixef(logit_sa3) - 1.96 * logit_sa3_se, UL = fixef(logit_sa3) + 1.96*logit_sa3_se))
# Both symptomatic and asymptomatic respiratory virus infection associated with reduced odds of S. aureus acquisition

# S. pneumoniae - respiratory virus infections (symptomatic vs. asymptomatic)  
logit_sp3 <- glmer(inf_multi_sp ~ relevel(inf_rv_cat3, ref="None") + hi_lag + mc_lag + sa_lag + month + 
                     sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + pcv + (1 | study_id), 
                   data = acquisition_sp, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(logit_sp3)
logit_sp3_se <- sqrt(diag(vcov(logit_sp3)))
exp(cbind(Est = fixef(logit_sp3), LL = fixef(logit_sp3) - 1.96 * logit_sp3_se, UL = fixef(logit_sp3) + 1.96*logit_sp3_se))
# Only symptomatic respiratory virus infection associated with increased odds of S. pneumoniae acquisition

# PNEUMOCOCCAL COLONIZATION DENSITY
table(metadata_rv_sp$inf_rv_cat3, useNA="always")
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_rv_cat3, summary)
metadata_rv_sp$inf_rv_cat3 <- as.factor(metadata_rv_sp$inf_rv_cat3)
summary(lm(inf_sp ~ relevel(inf_rv_cat3, ref="Asymptomatic") + inf_multi_hi + inf_multi_mc + inf_multi_sa + 
             month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))
# No significant difference in pneumococcal colonization density in symptomatic vs. asymptomatic detections

# **************
# ABSTRACT/INTRO
# **************

nrow(metadata_inf_np)
suppressWarnings(nrow(infants))

# *******
# RESULTS
# *******

# Respiratory virus infections and bacterial pathobiont colonization occur frequently during infancy

nrow(metadata_inf_np)
suppressWarnings(nrow(infants))
summary(infants$bw)
nrow(metadata_rv)
table(metadata_rv$inf_rv_yn)
prop.table(table(metadata_rv$inf_rv_yn))
prop.table(table(metadata_rv$rhino))
prop.table(table(metadata_rv$adeno))
prop.table(table(metadata_rv$rsv))

summary(logit_rv_age)
prop.table(table(metadata_rv$month, metadata_rv$inf_rv_yn),1)
prop.table(table(metadata_bact$month, metadata_bact$inf_multi_hi),1)
summary(logit_hi_age)
prop.table(table(metadata_bact$month, metadata_bact$inf_multi_mc),1)
summary(logit_mc_age)
prop.table(table(metadata_bact$month, metadata_bact$inf_multi_sp),1)
summary(logit_sp_age)
# The model for Staph aureus excludes the birth time point
prop.table(table(metadata_bact$month, metadata_bact$inf_multi_sa),1)
summary(logit_sa_age)

nrow(metadata_inf_np)
table(metadata_inf_np$current_uri, useNA="always")
prop.table(table(metadata_inf_np$current_uri))
summary(logit_rti_acq)
exp(cbind(Est = fixef(logit_rti_acq), LL = fixef(logit_rti_acq) - 1.96*logit_rti_acq_se, UL = fixef(logit_rti_acq) + 1.96*logit_rti_acq_se))

prop.table(table(metadata_rv$rsv, metadata_rv$current_uri),1)
prop.table(table(metadata_rv$flu, metadata_rv$current_uri),1)
prop.table(table(metadata_rv$hmpv, metadata_rv$current_uri),1)
prop.table(table(metadata_rv$adeno, metadata_rv$current_uri),1)
prop.table(table(metadata_rv$rhino, metadata_rv$current_uri),1)

summary(logit_rti_acq)
exp(cbind(Est = fixef(logit_rti_acq), LL = fixef(logit_rti_acq) - 1.96*logit_rti_acq_se, UL = fixef(logit_rti_acq) + 1.96*logit_rti_acq_se))
summary(logit_rti_sp)
exp(cbind(Est = fixef(logit_rti_sp), LL = fixef(logit_rti_sp) - 1.96 * logit_rti_sp_se, UL = fixef(logit_rti_sp) + 1.96*logit_rti_sp_se))
prop.table(table(metadata_bact$inf_hi_new, metadata_bact$inf_abx_any, useNA="always"),1)
prop.table(table(metadata_bact$inf_mc_new, metadata_bact$inf_abx_any, useNA="always"),1)
prop.table(table(metadata_bact$inf_sp_new, metadata_bact$inf_abx_any, useNA="always"),1)
metadata_bact$inf_any_new[metadata_bact$inf_hi_new=="N" & metadata_bact$inf_mc_new=="N" & metadata_bact$inf_sp_new=="N"] <- "N"
metadata_bact$inf_any_new[metadata_bact$inf_hi_new=="Y" | metadata_bact$inf_mc_new=="Y" | metadata_bact$inf_sp_new=="Y"] <- "Y"
prop.table(table(metadata_bact$inf_any_new, metadata_bact$inf_abx_any, useNA="always"),1)

# Respiratory viruses and environmental exposures associated with bacterial pathobiont colonization 

summary(logit_hi)
exp(cbind(Est = fixef(logit_hi), LL = fixef(logit_hi) - 1.96 * logit_hi_se, UL = fixef(logit_hi) + 1.96*logit_hi_se))
summary(logit_mc)
exp(cbind(Est = fixef(logit_mc), LL = fixef(logit_mc) - 1.96 * logit_mc_se, UL = fixef(logit_mc) + 1.96*logit_mc_se))
summary(logit_sp)
exp(cbind(Est = fixef(logit_sp), LL = fixef(logit_sp) - 1.96 * logit_sp_se, UL = fixef(logit_sp) + 1.96*logit_sp_se))
summary(logit_sa)
exp(cbind(Est = fixef(logit_sa), LL = fixef(logit_sa) - 1.96 * logit_sa_se, UL = fixef(logit_sa) + 1.96*logit_sa_se))

table(metadata_rv$inf_rv_cat3)
table(metadata_rv$inf_rv_cat2)
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_rv_yn, summary)
summary(lm(inf_sp ~ inf_rv_yn + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + season + mat_hiv + num_kids + breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_rv_cat3, summary)
summary(lm(inf_sp ~ relevel(inf_rv_cat3, ref="Asymptomatic") + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + season + mat_hiv + num_kids + 
             breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_rv_cat2, summary)
summary(lm(inf_sp ~ relevel(inf_rv_cat2, ref="Other") + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + season + mat_hiv + num_kids + 
             breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))

summary(logit_hi)
exp(cbind(Est = fixef(logit_hi), LL = fixef(logit_hi) - 1.96 * logit_hi_se, UL = fixef(logit_hi) + 1.96*logit_hi_se))
summary(logit_mc)
exp(cbind(Est = fixef(logit_mc), LL = fixef(logit_mc) - 1.96 * logit_mc_se, UL = fixef(logit_mc) + 1.96*logit_mc_se))
summary(logit_sa)
exp(cbind(Est = fixef(logit_sa), LL = fixef(logit_sa) - 1.96 * logit_sa_se, UL = fixef(logit_sa) + 1.96*logit_sa_se))
summary(logit_sp)
exp(cbind(Est = fixef(logit_sp), LL = fixef(logit_sp) - 1.96 * logit_sp_se, UL = fixef(logit_sp) + 1.96*logit_sp_se))

table(metadata_rv_sp$inf_multi_sp, metadata_rv_sp$inf_sp_new)
summary(lm(inf_sp ~ inf_rv_yn + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + season + mat_hiv + num_kids + breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$pcv, summary)

# Patterns of pathobiont colonization during infancy support synergistic interspecies relationships

fisher.test(metadata_bact$inf_multi_hi_pos, metadata_bact$inf_multi_mc_pos)
fisher.test(metadata_bact$inf_multi_mc_pos, metadata_bact$inf_multi_sp_pos)
fisher.test(metadata_bact$inf_multi_hi_pos, metadata_bact$inf_multi_sp_pos)
fisher.test(metadata_bact$inf_multi_sa_pos, metadata_bact$inf_multi_hi_pos)
fisher.test(metadata_bact$inf_multi_sa_pos, metadata_bact$inf_multi_mc_pos)
fisher.test(metadata_bact$inf_multi_sa_pos, metadata_bact$inf_multi_sp_pos)

summary(logit_hi)
exp(cbind(Est = fixef(logit_hi), LL = fixef(logit_hi) - 1.96 * logit_hi_se, UL = fixef(logit_hi) + 1.96*logit_hi_se))
summary(logit_mc)
exp(cbind(Est = fixef(logit_mc), LL = fixef(logit_mc) - 1.96 * logit_mc_se, UL = fixef(logit_mc) + 1.96*logit_mc_se))
summary(logit_sp)
exp(cbind(Est = fixef(logit_sp), LL = fixef(logit_sp) - 1.96 * logit_sp_se, UL = fixef(logit_sp) + 1.96*logit_sp_se))

table(metadata_rv_sp$inf_multi_sp, metadata_rv_sp$inf_sp_new)
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_multi_hi, summary)
summary(lm(inf_sp ~ inf_rv_yn + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + season + mat_hiv + num_kids + breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))

summary(logit_sa)
exp(cbind(Est = fixef(logit_sa), LL = fixef(logit_sa) - 1.96 * logit_sa_se, UL = fixef(logit_sa) + 1.96*logit_sa_se))

# *******
# TABLE 1
# *******

table(infants$sex)
prop.table(table(infants$sex))
summary(infants$bw)
table(infants$mat_hiv, infants$inf_hiv)
prop.table(table(infants$mat_hiv, infants$inf_hiv))
table(infants$residence)
prop.table(table(infants$residence))
table(infants$electric)
prop.table(table(infants$electric))
table(infants$wood)
prop.table(table(infants$wood))
table(infants$mat_educ)
prop.table(table(infants$mat_educ))
summary(infants$num_kids)
summary(infants$num_adol)
summary(infants$num_adults)
table(infants$season)
prop.table(table(infants$season))

pcv <- metadata_inf_np
pcv$pcv[!is.na(pcv$pcv1) | !is.na(pcv$pcv2) | !is.na(pcv$pcv3)] <- 1
pcv$pcv[is.na(pcv$pcv1) & !is.na(pcv$pcv2) & !is.na(pcv$pcv3)] <- 0
pcv$pcv <- as.numeric(pcv$pcv)
pcv_any <- pcv %>% group_by(study_id) %>% tally(pcv)
pcv_any$n <- as.numeric(pcv_any$n)
pcv_any$pcv_any <- NA
pcv_any$pcv_any[pcv_any$n>0] <- "Y"
pcv_any$pcv_any[pcv_any$n==0] <- "N"
table(pcv_any$pcv_any, useNA="always")
prop.table(table(pcv_any$pcv_any))
remove(pcv, pcv_any)

bm <- metadata_inf_np
bm$bm[bm$breastmilk=="Y"] <- 1
bm$bm[bm$breastmilk=="N"] <- 0
bm$bm <- as.numeric(bm$bm)
bm_any <- bm %>% group_by(study_id) %>% tally(bm)
bm_any$n <- as.numeric(bm_any$n)
bm_any$bm_any <- NA
bm_any$bm_any[bm_any$n>0] <- "Y"
bm_any$bm_any[bm_any$n==0] <- "N"
table(bm_any$bm_any, useNA="always")
prop.table(table(bm_any$bm_any))
remove(bm, bm_any)

abx <- metadata_inf_np
abx$abx[abx$inf_abx_any=="Y"] <- 1
abx$abx[abx$inf_abx_any=="N"] <- 0
abx$abx <- as.numeric(abx$abx)
abx_any <- abx %>% group_by(study_id) %>% tally(abx)
abx_any$n <- as.numeric(abx_any$n)
abx_any$abx_any <- NA
abx_any$abx_any[abx_any$n>0] <- "Y"
abx_any$abx_any[abx_any$n==0] <- "N"
table(abx_any$abx_any, useNA="always")
prop.table(table(abx_any$abx_any))
abx$amox[abx$inf_abx_amox=="Y"] <- 1
abx$amox[abx$inf_abx_amox=="N"] <- 0
abx$amox <- as.numeric(abx$amox)
amox_any <- abx %>% group_by(study_id) %>% tally(amox)
amox_any$n <- as.numeric(amox_any$n)
amox_any$amox_any <- NA
amox_any$amox_any[amox_any$n>0] <- "Y"
amox_any$amox_any[amox_any$n==0] <- "N"
table(amox_any$amox_any, useNA="always")
prop.table(table(amox_any$amox_any))
abx$metro[abx$inf_abx_metro=="Y"] <- 1
abx$metro[abx$inf_abx_metro=="N"] <- 0
abx$metro <- as.numeric(abx$metro)
metro_any <- abx %>% group_by(study_id) %>% tally(metro)
metro_any$n <- as.numeric(metro_any$n)
metro_any$metro_any <- NA
metro_any$metro_any[metro_any$n>0] <- "Y"
metro_any$metro_any[metro_any$n==0] <- "N"
table(metro_any$metro_any, useNA="always")
prop.table(table(metro_any$metro_any))
abx$cotrim[abx$inf_abx_cotrim=="Y"] <- 1
abx$cotrim[abx$inf_abx_cotrim=="N"] <- 0
abx$cotrim <- as.numeric(abx$cotrim)
cotrim_any <- abx %>% group_by(study_id) %>% tally(cotrim)
cotrim_any$n <- as.numeric(cotrim_any$n)
cotrim_any$cotrim_any <- NA
cotrim_any$cotrim_any[cotrim_any$n>0] <- "Y"
cotrim_any$cotrim_any[cotrim_any$n==0] <- "N"
table(cotrim_any$cotrim_any, useNA="always")
prop.table(table(cotrim_any$cotrim_any))
remove(abx, abx_any, amox_any, metro_any, cotrim_any)

# *******
# TABLE 2
# *******

summary(logit_hi)
exp(cbind(Est = fixef(logit_hi), LL = fixef(logit_hi) - 1.96 * logit_hi_se, UL = fixef(logit_hi) + 1.96*logit_hi_se))
summary(logit_mc)
exp(cbind(Est = fixef(logit_mc), LL = fixef(logit_mc) - 1.96 * logit_mc_se, UL = fixef(logit_mc) + 1.96*logit_mc_se))
summary(logit_sa)
exp(cbind(Est = fixef(logit_sa), LL = fixef(logit_sa) - 1.96 * logit_sa_se, UL = fixef(logit_sa) + 1.96*logit_sa_se))
summary(logit_sp)
exp(cbind(Est = fixef(logit_sp), LL = fixef(logit_sp) - 1.96 * logit_sp_se, UL = fixef(logit_sp) + 1.96*logit_sp_se))

# ********************
# SUPPLEMENTAL TABLE 1
# ********************

table(metadata_rv$adeno)
prop.table(table(metadata_rv$adeno))
table(metadata_rv$adeno, metadata_rv$current_uri)
prop.table(table(metadata_rv$adeno, metadata_rv$current_uri),1)

table(metadata_rv$hmpv)
prop.table(table(metadata_rv$hmpv))
table(metadata_rv$hmpv, metadata_rv$current_uri)
prop.table(table(metadata_rv$hmpv, metadata_rv$current_uri),1)

table(metadata_rv$flu)
prop.table(table(metadata_rv$flu))
table(metadata_rv$flu, metadata_rv$current_uri)
prop.table(table(metadata_rv$flu, metadata_rv$current_uri),1)
table(metadata_rv$flu, metadata_rv$inf_rv)
13/nrow(metadata_rv) # Flu A sample prevalence
1/nrow(metadata_rv) # Flu B sample prevalence
flu_a <- metadata_rv[grep("FLU A", metadata_rv$inf_rv),]
table(flu_a$current_uri)
prop.table(table(flu_a$current_uri))
flu_b <- metadata_rv[grep("FLU B", metadata_rv$inf_rv),]
table(flu_b$current_uri)
prop.table(table(flu_b$current_uri))

table(metadata_rv$paraflu)
prop.table(table(metadata_rv$paraflu))
table(metadata_rv$paraflu, metadata_rv$current_uri)
prop.table(table(metadata_rv$paraflu, metadata_rv$current_uri),1)
table(metadata_rv$paraflu, metadata_rv$inf_rv)

table(metadata_rv$paraflu)
prop.table(table(metadata_rv$paraflu))
table(metadata_rv$paraflu, metadata_rv$current_uri)
prop.table(table(metadata_rv$paraflu, metadata_rv$current_uri),1)
table(metadata_rv$paraflu, metadata_rv$inf_rv)
8/nrow(metadata_rv) # Paraflu 1 sample prevalence
4/nrow(metadata_rv) # Paraflu 2 sample prevalence
27/nrow(metadata_rv) # Paraflu 3 sample prevalence
p1 <- metadata_rv[grep("P1", metadata_rv$inf_rv),]
table(p1$current_uri)
prop.table(table(p1$current_uri))
p2 <- metadata_rv[grep("P2", metadata_rv$inf_rv),]
table(p2$current_uri)
prop.table(table(p2$current_uri))
p3 <- metadata_rv[grep("P3", metadata_rv$inf_rv),]
table(p3$current_uri)
prop.table(table(p3$current_uri))

table(metadata_rv$rhino)
prop.table(table(metadata_rv$rhino))
table(metadata_rv$rhino, metadata_rv$current_uri)
prop.table(table(metadata_rv$rhino, metadata_rv$current_uri),1)

table(metadata_rv$rsv)
prop.table(table(metadata_rv$rsv))
table(metadata_rv$rsv, metadata_rv$current_uri)
prop.table(table(metadata_rv$rsv, metadata_rv$current_uri),1)

table(metadata_rv$sars2)
prop.table(table(metadata_rv$sars2))
table(metadata_rv$sars2, metadata_rv$current_uri)
prop.table(table(metadata_rv$sars2, metadata_rv$current_uri),1)

table(metadata_rv$inf_rv_cat)
prop.table(table(metadata_rv$inf_rv_cat))
table(metadata_rv$inf_rv_cat, metadata_rv$current_uri)
prop.table(table(metadata_rv$inf_rv_cat, metadata_rv$current_uri),1)

# ********************
# SUPPLEMENTAL TABLE 2
# ********************

table(acquisition_hi$inf_rv_cat3)
table(acquisition_hi$inf_rv_cat3, acquisition_hi$inf_multi_hi)
prop.table(table(acquisition_hi$inf_rv_cat3, acquisition_hi$inf_multi_hi),1)
summary(logit_hi3)
exp(cbind(Est = fixef(logit_hi3), LL = fixef(logit_hi3) - 1.96 * logit_hi3_se, UL = fixef(logit_hi3) + 1.96*logit_hi3_se))

table(acquisition_mc$inf_rv_cat3)
table(acquisition_mc$inf_rv_cat3, acquisition_mc$inf_multi_mc)
prop.table(table(acquisition_mc$inf_rv_cat3, acquisition_mc$inf_multi_mc),1)
summary(logit_mc3)
exp(cbind(Est = fixef(logit_mc3), LL = fixef(logit_mc3) - 1.96 * logit_mc3_se, UL = fixef(logit_mc3) + 1.96*logit_mc3_se))

table(acquisition_sa$inf_rv_cat3)
table(acquisition_sa$inf_rv_cat3, acquisition_sa$inf_multi_sa)
prop.table(table(acquisition_sa$inf_rv_cat3, acquisition_sa$inf_multi_sa),1)
summary(logit_sa3)
exp(cbind(Est = fixef(logit_sa3), LL = fixef(logit_sa3) - 1.96 * logit_sa3_se, UL = fixef(logit_sa3) + 1.96*logit_sa3_se))

table(acquisition_sp$inf_rv_cat3)
table(acquisition_sp$inf_rv_cat3, acquisition_sp$inf_multi_sp)
prop.table(table(acquisition_sp$inf_rv_cat3, acquisition_sp$inf_multi_sp),1)
summary(logit_sp3)
exp(cbind(Est = fixef(logit_sp3), LL = fixef(logit_sp3) - 1.96 * logit_sp3_se, UL = fixef(logit_sp3) + 1.96*logit_sp3_se))

# ********************
# SUPPLEMENTAL TABLE 3
# ********************

table(acquisition_hi$inf_rv_cat2)
table(acquisition_hi$inf_rv_cat2, acquisition_hi$inf_multi_hi)
prop.table(table(acquisition_hi$inf_rv_cat2, acquisition_hi$inf_multi_hi),1)
summary(logit_hi2)
exp(cbind(Est = fixef(logit_hi2), LL = fixef(logit_hi2) - 1.96 * logit_hi2_se, UL = fixef(logit_hi2) + 1.96*logit_hi2_se))

table(acquisition_mc$inf_rv_cat2)
table(acquisition_mc$inf_rv_cat2, acquisition_mc$inf_multi_mc)
prop.table(table(acquisition_mc$inf_rv_cat2, acquisition_mc$inf_multi_mc),1)
summary(logit_mc2)
exp(cbind(Est = fixef(logit_mc2), LL = fixef(logit_mc2) - 1.96 * logit_mc2_se, UL = fixef(logit_mc2) + 1.96*logit_mc2_se))

table(acquisition_sa$inf_rv_cat2)
table(acquisition_sa$inf_rv_cat2, acquisition_sa$inf_multi_sa)
prop.table(table(acquisition_sa$inf_rv_cat2, acquisition_sa$inf_multi_sa),1)
summary(logit_sa2)
exp(cbind(Est = fixef(logit_sa2), LL = fixef(logit_sa2) - 1.96 * logit_sa2_se, UL = fixef(logit_sa2) + 1.96*logit_sa2_se))

table(acquisition_sp$inf_rv_cat2)
table(acquisition_sp$inf_rv_cat2, acquisition_sp$inf_multi_sp)
prop.table(table(acquisition_sp$inf_rv_cat2, acquisition_sp$inf_multi_sp),1)
summary(logit_sp2)
exp(cbind(Est = fixef(logit_sp2), LL = fixef(logit_sp2) - 1.96 * logit_sp2_se, UL = fixef(logit_sp2) + 1.96*logit_sp2_se))