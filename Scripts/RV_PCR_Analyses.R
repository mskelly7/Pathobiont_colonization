# Botswana Infant Microbiome Study - RV-Bacterial Pathobiont Analyses
# Matthew Kelly, MD, MPH 
# Analyses of respiratory virus and bacterial pathobiont PCR data
# Last update: October 18, 2024

remove(list=ls())
setwd("__________________") 
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
library(survival)
packageVersion("survival")
library(ggsignif)
library(grid)
library(lme4)
packageVersion("lme4")

metadata_inf_np <- read.csv("metadata_inf_np_RV.csv")
nrow(metadata_inf_np)

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

# Create datasets for survival analyses
survival <- metadata_rv_bact
survival <- survival[order(survival$study_id, survival$month),]
survival$month <- as.numeric(survival$month)
survival$study_id <- as.factor(survival$study_id)
survival <- suppressWarnings(slide(survival, "age_days", TimeVar="month", GroupVar="study_id", NewVar="age_lag", slideBy = -1))
names(survival)[names(survival) == "age_days"] <- "stop_age"
survival$stop <- as.numeric(survival$stop)
names(survival)[names(survival) == "month"] <- "stop_month"
names(survival)[names(survival) == "age_lag"] <- "start_age"
survival$inf_multi_hi[survival$inf_multi_hi=="Y"] <- 1
survival$inf_multi_hi[survival$inf_multi_hi=="N"] <- 0
survival$inf_multi_hi <- as.numeric(survival$inf_multi_hi)
survival$inf_multi_mc[survival$inf_multi_mc=="Y"] <- 1
survival$inf_multi_mc[survival$inf_multi_mc=="N"] <- 0
survival$inf_multi_mc <- as.numeric(survival$inf_multi_mc)
survival$inf_multi_sa[survival$inf_multi_sa=="Y"] <- 1
survival$inf_multi_sa[survival$inf_multi_sa=="N"] <- 0
survival$inf_multi_sa <- as.numeric(survival$inf_multi_sa)
survival$inf_multi_sp[survival$inf_multi_sp=="Y"] <- 1
survival$inf_multi_sp[survival$inf_multi_sp=="N"] <- 0
survival$inf_multi_sp <- as.numeric(survival$inf_multi_sp)
survival <- suppressWarnings(slide(survival, "inf_multi_hi", TimeVar="stop_month", GroupVar="study_id", NewVar="hi_lag", slideBy = -1))
survival <- suppressWarnings(slide(survival, "inf_multi_mc", TimeVar="stop_month", GroupVar="study_id", NewVar="mc_lag", slideBy = -1))
survival <- suppressWarnings(slide(survival, "inf_multi_sa", TimeVar="stop_month", GroupVar="study_id", NewVar="sa_lag", slideBy = -1))
survival <- suppressWarnings(slide(survival, "inf_multi_sp", TimeVar="stop_month", GroupVar="study_id", NewVar="sp_lag", slideBy = -1))
nrow(survival)
survival <- survival[,c("SampleID","study_id","start_age","stop_age","stop_month","inf_rv_yn","inf_rv_cat2","inf_multi_hi","hi_lag",
                        "inf_multi_mc","mc_lag","inf_multi_sa","sa_lag","inf_multi_sp","sp_lag", "sex", "lbw", "mat_hiv", "residence", "wood", "num_kids", "season", 
                        "breastmilk", "inf_abx_any", "pcv", "hib")]
# Remove 0-month timepoint (do not predict outcome at 0 months)
survival <- subset(survival, stop_month!="0")
nrow(survival)
# Remove intervals with baseline or persistent colonization
survival_hi <- subset(survival, hi_lag!=1)
nrow(survival_hi)
survival_mc <- subset(survival, mc_lag!=1)
nrow(survival_mc)
survival_sa <- subset(survival, sa_lag!=1)
nrow(survival_sa)
survival_sp <- subset(survival, sp_lag!=1)
nrow(survival_sp)

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
# ANALYSES EVALUATING FOR ASSOCIATIONS WITH THE RISK OF BACTERIAL PATHOBIONT COLONIZATION
# ***************************************************************************************

# Multivariable recurrent event Cox PH models evaluating associations between infant characteristics, virus detection, and pathobiont colonization on risk of specific pathobionts
# Infant characteristics evaluated were age (study timepoint), season, maternal HIV infection, number of kids in the household, breastfeeding, antibiotics
# Analyses for S. pneumoniae acquisition were additionally adjusted for PCV-13 doses 

# H. influenzae
cox_hi <- coxph(Surv(start_age, stop_age, inf_multi_hi) ~ inf_rv_yn + mc_lag + sa_lag + sp_lag + stop_month + 
                  sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any +  
                  cluster(study_id), id=study_id, data=survival_hi)
summary(cox_hi)
cox_hi_ph <- cox.zph(cox_hi) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)
cox_hi_ph
suppressWarnings(ggcoxzph(cox_hi_ph)) # generate plot of scaled Schoenfeld residuals
remove(cox_hi_ph)
# M. catarrhalis colonization is associated with an increased hazard of H. influenzae colonization
# A higher number of children in the household is associated with an increased hazard of H. influenzae colonization

# M. catarrhalis 
cox_mc <- coxph(Surv(start_age, stop_age, inf_multi_mc) ~ inf_rv_yn + hi_lag + sa_lag + sp_lag + stop_month +
                  sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + 
                  cluster(study_id), id=study_id, data=survival_mc)
summary(cox_mc)
cox_mc_ph <- cox.zph(cox_mc) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)
cox_mc_ph
suppressWarnings(ggcoxzph(cox_mc_ph)) # generate plot of scaled Schoenfeld residuals
remove(cox_mc_ph)
# RV detection and S. pneumoniae colonization are associated with an increased hazard of M. catarrhalis colonization
# Maternal HIV infection, higher number of children in the household associated with an increased hazard of M. catarrhalis colonization
# Breastfeeding associated with a lower hazard of M. catarrhalis colonization

# S. aureus
cox_sa <- coxph(Surv(start_age, stop_age, inf_multi_sa) ~ inf_rv_yn + hi_lag + mc_lag + sp_lag + stop_month + 
                  sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + 
                  cluster(study_id), id=study_id, data=survival_sa)
summary(cox_sa)
cox_sa_ph <- cox.zph(cox_sa) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)
cox_sa_ph
suppressWarnings(ggcoxzph(cox_sa_ph)) # generate plot of scaled Schoenfeld residuals
remove(cox_sa_ph)
# RV detection is associated with a lower risk of S. aureus colonization
# Higher hazard of S. aureus colonization during the rainy season

# S. pneumoniae - additionally adjust for PCV-13 doses
cox_sp <- coxph(Surv(start_age, stop_age, inf_multi_sp) ~ inf_rv_yn + hi_lag + mc_lag + sa_lag + stop_month + 
                  sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + pcv + 
                  cluster(study_id), id=study_id, data=survival_sp)
summary(cox_sp)
cox_sp_ph <- cox.zph(cox_sp) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)
cox_sp_ph
suppressWarnings(ggcoxzph(cox_sp_ph)) # generate plot of scaled Schoenfeld residuals
remove(cox_sp_ph)
# RV detection and M. catarrhalis colonization associated with higher hazard of S. pneumoniae colonization
# A higher number of children in the household and urban residence is associated with an increased hazard of S. pneumoniae colonization
# PCV-13 vaccination associated with a lower hazard of S. pneumoniae colonization

# Pneumococcal colonization density among infants with S. pneumoniae carriage
# Linear regression model evaluating associations between infant characteristics, viral infections, and other pathobiont colonization on pneumococcal colonization density
metadata_rv_sp <- subset(metadata_rv_bact, inf_sp>0 & inf_multi_sp=="Y")
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_rv_yn, summary)
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_multi_hi, summary)
summary(lm(inf_sp ~ inf_rv_yn + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + 
             pcv, data=metadata_rv_sp))
# RV detection is associated with higher pneumococcal colonization density in samples with Sp colonization (linear regression, p=0.006)
# Higher pneumococcal colonization density with H. influenzae colonization (p<0.0001), decreasing age (p=0.001), a higher number of children in the household (p=0.009),
# and antibiotic exposure since the last study visit (p=0.006)
# No association with number of PCV-13 doses

# Analyses categorizing respiratory viruses as involving rhino/entero only vs. other respiratory virus infections

# H. influenzae - respiratory virus infections (rhino/entero vs. other) 
cox_hi2 <- coxph(Surv(start_age, stop_age, inf_multi_hi) ~ inf_rv_cat2 + mc_lag + sa_lag + sp_lag + stop_month + 
                   sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + 
                   cluster(study_id), id=study_id, data=survival_hi)
summary(cox_hi2)
cox_hi2_ph <- cox.zph(cox_hi2) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)
cox_hi2_ph
suppressWarnings(ggcoxzph(cox_hi2_ph)) # generate plot of scaled Schoenfeld residuals
remove(cox_hi2_ph)
# Higher hazard of H. influenzae colonization with "other" respiratory virus infections 

# M. catarrhalis - respiratory virus infections (rhino/entero vs. other) 
cox_mc2 <- coxph(Surv(start_age, stop_age, inf_multi_mc) ~ inf_rv_cat2 + hi_lag + sa_lag + sp_lag + stop_month + 
                   sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + 
                   cluster(study_id), id=study_id, data=survival_mc)
summary(cox_mc2)
cox_mc2_ph <- cox.zph(cox_mc2) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)
cox_mc2_ph
suppressWarnings(ggcoxzph(cox_mc2_ph)) # generate plot of scaled Schoenfeld residuals
remove(cox_mc2_ph)
# Higher hazard of H. influenzae colonization with "other" respiratory virus infections (p=0.003)

# S. aureus - respiratory virus infections (rhino/entero vs. other) 
cox_sa2 <- coxph(Surv(start_age, stop_age, inf_multi_sa) ~ inf_rv_cat2 + hi_lag + mc_lag + sp_lag + stop_month + 
                   sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + 
                   cluster(study_id), id=study_id, data=survival_sa)
summary(cox_sa2)
cox_sa2_ph <- cox.zph(cox_sa2) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)
cox_sa2_ph
suppressWarnings(ggcoxzph(cox_sa2_ph)) # generate plot of scaled Schoenfeld residuals
remove(cox_sa2_ph)
# Both rhino/entero and other respiratory viruses associated with lower hazard of S. aureus colonization

# S. pneumoniae - respiratory virus infections (rhino/entero vs. other) 
cox_sp2 <- coxph(Surv(start_age, stop_age, inf_multi_sp) ~ inf_rv_cat2 + hi_lag + mc_lag + sa_lag + stop_month + 
                   sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + pcv + 
                   cluster(study_id), id=study_id, data=survival_sp)
summary(cox_sp2)
cox_sp2_ph <- cox.zph(cox_sp2) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)
cox_sp2_ph
suppressWarnings(ggcoxzph(cox_sp2_ph)) # generate plot of scaled Schoenfeld residuals
remove(cox_sp2_ph)
# Both rhino/entero and other respiratory virus infections associated with an increased risk of S. pneumoniae colonization

# PNEUMOCOCCAL COLONIZATION DENSITY
table(metadata_rv_sp$inf_rv_cat2)
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_rv_cat2, summary)
metadata_rv_sp$inf_rv_cat2 <- as.factor(metadata_rv_sp$inf_rv_cat2)
summary(lm(inf_sp ~ relevel(inf_rv_cat2, ref="Other") + inf_multi_hi + inf_multi_mc + inf_multi_sa + 
             month + sex + lbw + mat_hiv + residence + wood + num_kids + season + breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))
# No difference in pneumococcal colonization density comparing infections involving rhino/entero only vs. infections with other viruses

# ********
# ABSTRACT
# ********

suppressWarnings(nrow(infants))
nrow(metadata_inf_np)
table(metadata_rv$inf_rv_yn)



# *******
# RESULTS
# *******

# Respiratory virus infections and bacterial pathobiont colonization occur frequently during infancy

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

# Infant and sample characteristics associated with the dynamics of bacterial pathobiont colonization 

summary(cox_hi)
summary(cox_mc)
summary(cox_sp)
summary(cox_sa)

table(metadata_rv$inf_rv_cat2)
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_rv_yn, summary)
summary(lm(inf_sp ~ inf_rv_yn + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + season + mat_hiv + num_kids + breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))
summary(lm(inf_sp ~ relevel(inf_rv_cat2, ref="Other") + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + season + mat_hiv + num_kids + 
             breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))

summary(cox_hi)
summary(cox_mc)
summary(cox_sp)
summary(cox_sa)

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

summary(cox_hi)
summary(cox_mc)
summary(cox_sp)
summary(cox_sa)

table(metadata_rv_sp$inf_multi_sp, metadata_rv_sp$inf_sp_new)
tapply(metadata_rv_sp$inf_sp, metadata_rv_sp$inf_multi_hi, summary)
summary(lm(inf_sp ~ inf_rv_yn + inf_multi_hi + inf_multi_mc + inf_multi_sa + month + season + mat_hiv + num_kids + breastmilk + inf_abx_any + pcv, data=metadata_rv_sp))

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

summary(cox_hi)
summary(cox_mc)
summary(cox_sa)
summary(cox_sp)

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

table(survival_hi$inf_rv_yn)
table(survival_hi$inf_rv_cat2)
table(survival_hi$inf_rv_yn, survival_hi$inf_multi_hi)
table(survival_hi$inf_rv_cat2, survival_hi$inf_multi_hi)
prop.table(table(survival_hi$inf_rv_yn, survival_hi$inf_multi_hi),1)
prop.table(table(survival_hi$inf_rv_cat2, survival_hi$inf_multi_hi),1)
summary(cox_hi)
summary(cox_hi2)

table(survival_mc$inf_rv_yn)
table(survival_mc$inf_rv_cat2)
table(survival_mc$inf_rv_yn, survival_mc$inf_multi_mc)
table(survival_mc$inf_rv_cat2, survival_mc$inf_multi_mc)
prop.table(table(survival_mc$inf_rv_yn, survival_mc$inf_multi_mc),1)
prop.table(table(survival_mc$inf_rv_cat2, survival_mc$inf_multi_mc),1)
summary(cox_mc)
summary(cox_mc2)

table(survival_sa$inf_rv_yn)
table(survival_sa$inf_rv_cat2)
table(survival_sa$inf_rv_yn, survival_sa$inf_multi_sa)
table(survival_sa$inf_rv_cat2, survival_sa$inf_multi_sa)
prop.table(table(survival_sa$inf_rv_yn, survival_sa$inf_multi_sa),1)
prop.table(table(survival_sa$inf_rv_cat2, survival_sa$inf_multi_sa),1)
summary(cox_sa)
summary(cox_sa2)

table(survival_sp$inf_rv_yn)
table(survival_sp$inf_rv_cat2)
table(survival_sp$inf_rv_yn, survival_sp$inf_multi_sp)
table(survival_sp$inf_rv_cat2, survival_sp$inf_multi_sp)
prop.table(table(survival_sp$inf_rv_yn, survival_sp$inf_multi_sp),1)
prop.table(table(survival_sp$inf_rv_cat2, survival_sp$inf_multi_sp),1)
summary(cox_sp)
summary(cox_sp2)