# Botswana Infant Microbiome Study - RV-Bacterial Analyses
# Matthew Kelly, MD, MPH 
# Supplementary Figure 2
# Last updated: October 21, 2024

remove(list=ls())
setwd("_____________________") 
set.seed(1234)

version
library(phyloseq)
packageVersion("phyloseq")
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(cowplot)
library(ggtext)
library(DataCombine)
theme_set(theme_light()) 

# Create color palettes
cluster_cols_9 <- c("gray20", # OTH
                            "mediumorchid4", # STA
                            "darkslateblue", # COR
                            "coral3", # STR
                            "indianred4", # CD
                            "peru", # HAE
                            "dodgerblue3", # CDM
                            "darkolivegreen4") # MOR
                            
# Load required datasets
phy.inf.np.16s <- read_rds("phy.inf.np.16s.rds")
metadata_16s <- data.frame(sample_data(phy.inf.np.16s))

library(ggplot2)
library(ggalluvial)
require(dplyr)
library(formattable)

theme_barplot <-   theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                          axis.text.y = element_text(size=8, color="black"), 
                          axis.title.y = element_text(angle=90, size=9, color="black"), 
                          axis.text.x = element_text(size=8.5, color="black"),
                          axis.title.x = element_text(size=9, color="black"),
                          legend.position = "right", legend.text=element_text(size=8, color="black"), legend.title=element_blank(),
                          strip.background =element_rect(fill="grey20"), 
                          strip.text = element_text(size=8, colour = 'white', margin = ggplot2::margin(-0.01,0,-0.01,0, "cm")), 
                          plot.title = element_text(size = 9, hjust = 0.5, margin = ggplot2::margin(0,0,-1.5,0)), plot.margin = unit(c(0.3, 0, 0, 0), "cm"),
                          legend.key.size = unit(0.4, 'cm'), legend.box.spacing = unit(0, "pt"), panel.spacing = unit(0,'lines'))

# H. influenzae colonization

cluster_hi <- metadata_16s[,c("study_id", "month", "cluster", "inf_multi_hi")]
cluster_01 <- subset(cluster_hi, (month=="0" | month=="1") & !is.na(inf_multi_hi))
cluster_01 <- cluster_01[order(cluster_01$study_id, cluster_01$month),]
cluster_01 <- subset(cluster_01, month!="0" | inf_multi_hi=="N")
cluster_01 <- cluster_01 %>% group_by(study_id) %>% filter(n() == 2)
cluster_01$timing <- ""
cluster_01$timing[cluster_01$month=="0"] <- "Pre"
cluster_01$timing[cluster_01$month=="1"] <- "Post"
cluster_01 <- suppressWarnings(slide(cluster_01, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_hi_lag", slideBy = 1))
cluster_01$hi_yn[cluster_01$inf_multi_hi=="Y" | cluster_01$inf_hi_lag=="Y"] <- "Yes"
cluster_01$hi_yn[cluster_01$inf_multi_hi=="N" & (cluster_01$inf_hi_lag=="N" | is.na(cluster_01$inf_hi_lag))] <- "No"
cluster_01$study_id <- paste0(cluster_01$study_id, "_01")
cluster_01 <- cluster_01[,c("study_id", "timing", "cluster", "hi_yn")]

cluster_12 <- subset(cluster_hi, (month=="1" | month=="2") & !is.na(inf_multi_hi))
cluster_12 <- cluster_12[order(cluster_12$study_id, cluster_12$month),]
cluster_12 <- subset(cluster_12, month!="1" | inf_multi_hi=="N")
cluster_12 <- cluster_12 %>% group_by(study_id) %>% filter(n() == 2)
cluster_12$timing <- ""
cluster_12$timing[cluster_12$month=="1"] <- "Pre"
cluster_12$timing[cluster_12$month=="2"] <- "Post"
cluster_12 <- suppressWarnings(slide(cluster_12, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_hi_lag", slideBy = 1))
cluster_12$hi_yn[cluster_12$inf_multi_hi=="Y" | cluster_12$inf_hi_lag=="Y"] <- "Yes"
cluster_12$hi_yn[cluster_12$inf_multi_hi=="N" & (cluster_12$inf_hi_lag=="N" | is.na(cluster_12$inf_hi_lag))] <- "No"
cluster_12$study_id <- paste0(cluster_12$study_id, "_12")
cluster_12 <- cluster_12[,c("study_id", "timing", "cluster", "hi_yn")]

cluster_23 <- subset(cluster_hi, (month=="2" | month=="3") & !is.na(inf_multi_hi))
cluster_23 <- cluster_23[order(cluster_23$study_id, cluster_23$month),]
cluster_23 <- subset(cluster_23, month!="2" | inf_multi_hi=="N")
cluster_23 <- cluster_23 %>% group_by(study_id) %>% filter(n() == 2)
cluster_23$timing <- ""
cluster_23$timing[cluster_23$month=="2"] <- "Pre"
cluster_23$timing[cluster_23$month=="3"] <- "Post"
cluster_23 <- suppressWarnings(slide(cluster_23, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_hi_lag", slideBy = 1))
cluster_23$hi_yn[cluster_23$inf_multi_hi=="Y" | cluster_23$inf_hi_lag=="Y"] <- "Yes"
cluster_23$hi_yn[cluster_23$inf_multi_hi=="N" & (cluster_23$inf_hi_lag=="N" | is.na(cluster_23$inf_hi_lag))] <- "No"
cluster_23$study_id <- paste0(cluster_23$study_id, "_23")
cluster_23 <- cluster_23[,c("study_id", "timing", "cluster", "hi_yn")]

cluster_34 <- subset(cluster_hi, (month=="3" | month=="4") & !is.na(inf_multi_hi))
cluster_34 <- cluster_34[order(cluster_34$study_id, cluster_34$month),]
cluster_34 <- subset(cluster_34, month!="3" | inf_multi_hi=="N")
cluster_34 <- cluster_34 %>% group_by(study_id) %>% filter(n() == 2)
cluster_34$timing <- ""
cluster_34$timing[cluster_34$month=="3"] <- "Pre"
cluster_34$timing[cluster_34$month=="4"] <- "Post"
cluster_34 <- suppressWarnings(slide(cluster_34, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_hi_lag", slideBy = 1))
cluster_34$hi_yn[cluster_34$inf_multi_hi=="Y" | cluster_34$inf_hi_lag=="Y"] <- "Yes"
cluster_34$hi_yn[cluster_34$inf_multi_hi=="N" & (cluster_34$inf_hi_lag=="N" | is.na(cluster_34$inf_hi_lag))] <- "No"
cluster_34$study_id <- paste0(cluster_34$study_id, "_34")
cluster_34 <- cluster_34[,c("study_id", "timing", "cluster", "hi_yn")]

cluster_45 <- subset(cluster_hi, (month=="4" | month=="5") & !is.na(inf_multi_hi))
cluster_45 <- cluster_45[order(cluster_45$study_id, cluster_45$month),]
cluster_45 <- subset(cluster_45, month!="4" | inf_multi_hi=="N")
cluster_45 <- cluster_45 %>% group_by(study_id) %>% filter(n() == 2)
cluster_45$timing <- ""
cluster_45$timing[cluster_45$month=="4"] <- "Pre"
cluster_45$timing[cluster_45$month=="5"] <- "Post"
cluster_45 <- suppressWarnings(slide(cluster_45, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_hi_lag", slideBy = 1))
cluster_45$hi_yn[cluster_45$inf_multi_hi=="Y" | cluster_45$inf_hi_lag=="Y"] <- "Yes"
cluster_45$hi_yn[cluster_45$inf_multi_hi=="N" & (cluster_45$inf_hi_lag=="N" | is.na(cluster_45$inf_hi_lag))] <- "No"
cluster_45$study_id <- paste0(cluster_45$study_id, "_45")
cluster_45 <- cluster_45[,c("study_id", "timing", "cluster", "hi_yn")]

cluster_56 <- subset(cluster_hi, (month=="5" | month=="6") & !is.na(inf_multi_hi))
cluster_56 <- cluster_56[order(cluster_56$study_id, cluster_56$month),]
cluster_56 <- subset(cluster_56, month!="5" | inf_multi_hi=="N")
cluster_56 <- cluster_56 %>% group_by(study_id) %>% filter(n() == 2)
cluster_56$timing <- ""
cluster_56$timing[cluster_56$month=="5"] <- "Pre"
cluster_56$timing[cluster_56$month=="6"] <- "Post"
cluster_56 <- suppressWarnings(slide(cluster_56, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_hi_lag", slideBy = 1))
cluster_56$hi_yn[cluster_56$inf_multi_hi=="Y" | cluster_56$inf_hi_lag=="Y"] <- "Yes"
cluster_56$hi_yn[cluster_56$inf_multi_hi=="N" & (cluster_56$inf_hi_lag=="N" | is.na(cluster_56$inf_hi_lag))] <- "No"
cluster_56$study_id <- paste0(cluster_56$study_id, "_56")
cluster_56 <- cluster_56[,c("study_id", "timing", "cluster", "hi_yn")]

cluster_68 <- subset(cluster_hi, (month=="6" | month=="8") & !is.na(inf_multi_hi))
cluster_68 <- cluster_68[order(cluster_68$study_id, cluster_68$month),]
cluster_68 <- subset(cluster_68, month!="6" | inf_multi_hi=="N")
cluster_68 <- cluster_68 %>% group_by(study_id) %>% filter(n() == 2)
cluster_68$timing <- ""
cluster_68$timing[cluster_68$month=="6"] <- "Pre"
cluster_68$timing[cluster_68$month=="8"] <- "Post"
cluster_68 <- suppressWarnings(slide(cluster_68, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_hi_lag", slideBy = 1))
cluster_68$hi_yn[cluster_68$inf_multi_hi=="Y" | cluster_68$inf_hi_lag=="Y"] <- "Yes"
cluster_68$hi_yn[cluster_68$inf_multi_hi=="N" & (cluster_68$inf_hi_lag=="N" | is.na(cluster_68$inf_hi_lag))] <- "No"
cluster_68$study_id <- paste0(cluster_68$study_id, "_68")
cluster_68 <- cluster_68[,c("study_id", "timing", "cluster", "hi_yn")]

cluster_810 <- subset(cluster_hi, (month=="8" | month=="10") & !is.na(inf_multi_hi))
cluster_810 <- cluster_810[order(cluster_810$study_id, cluster_810$month),]
cluster_810 <- subset(cluster_810, month!="8" | inf_multi_hi=="N")
cluster_810 <- cluster_810 %>% group_by(study_id) %>% filter(n() == 2)
cluster_810$timing <- ""
cluster_810$timing[cluster_810$month=="8"] <- "Pre"
cluster_810$timing[cluster_810$month=="10"] <- "Post"
cluster_810 <- suppressWarnings(slide(cluster_810, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_hi_lag", slideBy = 1))
cluster_810$hi_yn[cluster_810$inf_multi_hi=="Y" | cluster_810$inf_hi_lag=="Y"] <- "Yes"
cluster_810$hi_yn[cluster_810$inf_multi_hi=="N" & (cluster_810$inf_hi_lag=="N" | is.na(cluster_810$inf_hi_lag))] <- "No"
cluster_810$study_id <- paste0(cluster_810$study_id, "_810")
cluster_810 <- cluster_810[,c("study_id", "timing", "cluster", "hi_yn")]

cluster_1012 <- subset(cluster_hi, (month=="10" | month=="12") & !is.na(inf_multi_hi))
cluster_1012 <- cluster_1012[order(cluster_1012$study_id, cluster_1012$month),]
cluster_1012 <- subset(cluster_1012, month!="10" | inf_multi_hi=="N")
cluster_1012 <- cluster_1012 %>% group_by(study_id) %>% filter(n() == 2)
cluster_1012$timing <- ""
cluster_1012$timing[cluster_1012$month=="10"] <- "Pre"
cluster_1012$timing[cluster_1012$month=="12"] <- "Post"
cluster_1012 <- suppressWarnings(slide(cluster_1012, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_hi_lag", slideBy = 1))
cluster_1012$hi_yn[cluster_1012$inf_multi_hi=="Y" | cluster_1012$inf_hi_lag=="Y"] <- "Yes"
cluster_1012$hi_yn[cluster_1012$inf_multi_hi=="N" & (cluster_1012$inf_hi_lag=="N" | is.na(cluster_1012$inf_hi_lag))] <- "No"
cluster_1012$study_id <- paste0(cluster_1012$study_id, "_1012")
cluster_1012 <- cluster_1012[,c("study_id", "timing", "cluster", "hi_yn")]

cluster_all <- rbind(cluster_01, cluster_12, cluster_23, cluster_34, cluster_45, cluster_56, cluster_68, cluster_810, cluster_1012)

counts <- cluster_all %>% dplyr::count(hi_yn) 
cluster_all <- merge(cluster_all, counts, by="hi_yn")
cluster_all$pct <- (1/(cluster_all$n))*2
cluster_all$timing <- factor(cluster_all$timing, levels=c("Pre","Post"))

alpha_hi <- c(0.3,0.3,1.0,0.3,0.3,0.3,0.3,0.3,0.3,0.3,1.0,0.3,0.3,0.3,0.3,0.3,
              0.3,0.3,1.0,0.3,0.3,0.3,0.3,0.3,0.3,0.3,1.0,0.3,0.3,0.3,0.3)
# There are only 31 arguments because there are no OTH samples in the post-Hi acquisition group
title_hi <- expression(paste(italic("H. influenzae"), " acquisition"))
hi_alluvium <- ggplot(cluster_all, aes(x = timing, stratum = cluster, alluvium = study_id, y = pct, fill = cluster)) +
  scale_fill_manual(values=cluster_cols_9) +
  geom_flow(alpha=0.3) +
  geom_stratum(width=0.38, color="gray20", alpha=alpha_hi) + ylab("Proportion of samples") + xlab("") +
  theme_barplot + facet_wrap(~hi_yn) + ggtitle(title_hi) + 
  guides(fill=guide_legend(override.aes = list(alpha = c(0.3,0.3,0.3,0.3,0.3,1.0,0.3,0.3))))

png(file="R_Plots/Figure_S2/Figure_S2_hi.png", width = 3.75, height = 2.75, units = 'in', res = 300)
plot(hi_alluvium)
dev.off()

# M. catarrhalis colonization

cluster_mc <- metadata_16s[,c("study_id", "month", "cluster", "inf_multi_mc")]
cluster_01 <- subset(cluster_mc, (month=="0" | month=="1") & !is.na(inf_multi_mc))
cluster_01 <- cluster_01[order(cluster_01$study_id, cluster_01$month),]
cluster_01 <- subset(cluster_01, month!="0" | inf_multi_mc=="N")
cluster_01 <- cluster_01 %>% group_by(study_id) %>% filter(n() == 2)
cluster_01$timing <- ""
cluster_01$timing[cluster_01$month=="0"] <- "Pre"
cluster_01$timing[cluster_01$month=="1"] <- "Post"
cluster_01 <- suppressWarnings(slide(cluster_01, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_mc_lag", slideBy = 1))
cluster_01$mc_yn[cluster_01$inf_multi_mc=="Y" | cluster_01$inf_mc_lag=="Y"] <- "Yes"
cluster_01$mc_yn[cluster_01$inf_multi_mc=="N" & (cluster_01$inf_mc_lag=="N" | is.na(cluster_01$inf_mc_lag))] <- "No"
cluster_01$study_id <- paste0(cluster_01$study_id, "_01")
cluster_01 <- cluster_01[,c("study_id", "timing", "cluster", "mc_yn")]

cluster_12 <- subset(cluster_mc, (month=="1" | month=="2") & !is.na(inf_multi_mc))
cluster_12 <- cluster_12[order(cluster_12$study_id, cluster_12$month),]
cluster_12 <- subset(cluster_12, month!="1" | inf_multi_mc=="N")
cluster_12 <- cluster_12 %>% group_by(study_id) %>% filter(n() == 2)
cluster_12$timing <- ""
cluster_12$timing[cluster_12$month=="1"] <- "Pre"
cluster_12$timing[cluster_12$month=="2"] <- "Post"
cluster_12 <- suppressWarnings(slide(cluster_12, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_mc_lag", slideBy = 1))
cluster_12$mc_yn[cluster_12$inf_multi_mc=="Y" | cluster_12$inf_mc_lag=="Y"] <- "Yes"
cluster_12$mc_yn[cluster_12$inf_multi_mc=="N" & (cluster_12$inf_mc_lag=="N" | is.na(cluster_12$inf_mc_lag))] <- "No"
cluster_12$study_id <- paste0(cluster_12$study_id, "_12")
cluster_12 <- cluster_12[,c("study_id", "timing", "cluster", "mc_yn")]

cluster_23 <- subset(cluster_mc, (month=="2" | month=="3") & !is.na(inf_multi_mc))
cluster_23 <- cluster_23[order(cluster_23$study_id, cluster_23$month),]
cluster_23 <- subset(cluster_23, month!="2" | inf_multi_mc=="N")
cluster_23 <- cluster_23 %>% group_by(study_id) %>% filter(n() == 2)
cluster_23$timing <- ""
cluster_23$timing[cluster_23$month=="2"] <- "Pre"
cluster_23$timing[cluster_23$month=="3"] <- "Post"
cluster_23 <- suppressWarnings(slide(cluster_23, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_mc_lag", slideBy = 1))
cluster_23$mc_yn[cluster_23$inf_multi_mc=="Y" | cluster_23$inf_mc_lag=="Y"] <- "Yes"
cluster_23$mc_yn[cluster_23$inf_multi_mc=="N" & (cluster_23$inf_mc_lag=="N" | is.na(cluster_23$inf_mc_lag))] <- "No"
cluster_23$study_id <- paste0(cluster_23$study_id, "_23")
cluster_23 <- cluster_23[,c("study_id", "timing", "cluster", "mc_yn")]

cluster_34 <- subset(cluster_mc, (month=="3" | month=="4") & !is.na(inf_multi_mc))
cluster_34 <- cluster_34[order(cluster_34$study_id, cluster_34$month),]
cluster_34 <- subset(cluster_34, month!="3" | inf_multi_mc=="N")
cluster_34 <- cluster_34 %>% group_by(study_id) %>% filter(n() == 2)
cluster_34$timing <- ""
cluster_34$timing[cluster_34$month=="3"] <- "Pre"
cluster_34$timing[cluster_34$month=="4"] <- "Post"
cluster_34 <- suppressWarnings(slide(cluster_34, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_mc_lag", slideBy = 1))
cluster_34$mc_yn[cluster_34$inf_multi_mc=="Y" | cluster_34$inf_mc_lag=="Y"] <- "Yes"
cluster_34$mc_yn[cluster_34$inf_multi_mc=="N" & (cluster_34$inf_mc_lag=="N" | is.na(cluster_34$inf_mc_lag))] <- "No"
cluster_34$study_id <- paste0(cluster_34$study_id, "_34")
cluster_34 <- cluster_34[,c("study_id", "timing", "cluster", "mc_yn")]

cluster_45 <- subset(cluster_mc, (month=="4" | month=="5") & !is.na(inf_multi_mc))
cluster_45 <- cluster_45[order(cluster_45$study_id, cluster_45$month),]
cluster_45 <- subset(cluster_45, month!="4" | inf_multi_mc=="N")
cluster_45 <- cluster_45 %>% group_by(study_id) %>% filter(n() == 2)
cluster_45$timing <- ""
cluster_45$timing[cluster_45$month=="4"] <- "Pre"
cluster_45$timing[cluster_45$month=="5"] <- "Post"
cluster_45 <- suppressWarnings(slide(cluster_45, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_mc_lag", slideBy = 1))
cluster_45$mc_yn[cluster_45$inf_multi_mc=="Y" | cluster_45$inf_mc_lag=="Y"] <- "Yes"
cluster_45$mc_yn[cluster_45$inf_multi_mc=="N" & (cluster_45$inf_mc_lag=="N" | is.na(cluster_45$inf_mc_lag))] <- "No"
cluster_45$study_id <- paste0(cluster_45$study_id, "_45")
cluster_45 <- cluster_45[,c("study_id", "timing", "cluster", "mc_yn")]

cluster_56 <- subset(cluster_mc, (month=="5" | month=="6") & !is.na(inf_multi_mc))
cluster_56 <- cluster_56[order(cluster_56$study_id, cluster_56$month),]
cluster_56 <- subset(cluster_56, month!="5" | inf_multi_mc=="N")
cluster_56 <- cluster_56 %>% group_by(study_id) %>% filter(n() == 2)
cluster_56$timing <- ""
cluster_56$timing[cluster_56$month=="5"] <- "Pre"
cluster_56$timing[cluster_56$month=="6"] <- "Post"
cluster_56 <- suppressWarnings(slide(cluster_56, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_mc_lag", slideBy = 1))
cluster_56$mc_yn[cluster_56$inf_multi_mc=="Y" | cluster_56$inf_mc_lag=="Y"] <- "Yes"
cluster_56$mc_yn[cluster_56$inf_multi_mc=="N" & (cluster_56$inf_mc_lag=="N" | is.na(cluster_56$inf_mc_lag))] <- "No"
cluster_56$study_id <- paste0(cluster_56$study_id, "_56")
cluster_56 <- cluster_56[,c("study_id", "timing", "cluster", "mc_yn")]

cluster_68 <- subset(cluster_mc, (month=="6" | month=="8") & !is.na(inf_multi_mc))
cluster_68 <- cluster_68[order(cluster_68$study_id, cluster_68$month),]
cluster_68 <- subset(cluster_68, month!="6" | inf_multi_mc=="N")
cluster_68 <- cluster_68 %>% group_by(study_id) %>% filter(n() == 2)
cluster_68$timing <- ""
cluster_68$timing[cluster_68$month=="6"] <- "Pre"
cluster_68$timing[cluster_68$month=="8"] <- "Post"
cluster_68 <- suppressWarnings(slide(cluster_68, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_mc_lag", slideBy = 1))
cluster_68$mc_yn[cluster_68$inf_multi_mc=="Y" | cluster_68$inf_mc_lag=="Y"] <- "Yes"
cluster_68$mc_yn[cluster_68$inf_multi_mc=="N" & (cluster_68$inf_mc_lag=="N" | is.na(cluster_68$inf_mc_lag))] <- "No"
cluster_68$study_id <- paste0(cluster_68$study_id, "_68")
cluster_68 <- cluster_68[,c("study_id", "timing", "cluster", "mc_yn")]

cluster_810 <- subset(cluster_mc, (month=="8" | month=="10") & !is.na(inf_multi_mc))
cluster_810 <- cluster_810[order(cluster_810$study_id, cluster_810$month),]
cluster_810 <- subset(cluster_810, month!="8" | inf_multi_mc=="N")
cluster_810 <- cluster_810 %>% group_by(study_id) %>% filter(n() == 2)
cluster_810$timing <- ""
cluster_810$timing[cluster_810$month=="8"] <- "Pre"
cluster_810$timing[cluster_810$month=="10"] <- "Post"
cluster_810 <- suppressWarnings(slide(cluster_810, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_mc_lag", slideBy = 1))
cluster_810$mc_yn[cluster_810$inf_multi_mc=="Y" | cluster_810$inf_mc_lag=="Y"] <- "Yes"
cluster_810$mc_yn[cluster_810$inf_multi_mc=="N" & (cluster_810$inf_mc_lag=="N" | is.na(cluster_810$inf_mc_lag))] <- "No"
cluster_810$study_id <- paste0(cluster_810$study_id, "_810")
cluster_810 <- cluster_810[,c("study_id", "timing", "cluster", "mc_yn")]

cluster_1012 <- subset(cluster_mc, (month=="10" | month=="12") & !is.na(inf_multi_mc))
cluster_1012 <- cluster_1012[order(cluster_1012$study_id, cluster_1012$month),]
cluster_1012 <- subset(cluster_1012, month!="10" | inf_multi_mc=="N")
cluster_1012 <- cluster_1012 %>% group_by(study_id) %>% filter(n() == 2)
cluster_1012$timing <- ""
cluster_1012$timing[cluster_1012$month=="10"] <- "Pre"
cluster_1012$timing[cluster_1012$month=="12"] <- "Post"
cluster_1012 <- suppressWarnings(slide(cluster_1012, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_mc_lag", slideBy = 1))
cluster_1012$mc_yn[cluster_1012$inf_multi_mc=="Y" | cluster_1012$inf_mc_lag=="Y"] <- "Yes"
cluster_1012$mc_yn[cluster_1012$inf_multi_mc=="N" & (cluster_1012$inf_mc_lag=="N" | is.na(cluster_1012$inf_mc_lag))] <- "No"
cluster_1012$study_id <- paste0(cluster_1012$study_id, "_1012")
cluster_1012 <- cluster_1012[,c("study_id", "timing", "cluster", "mc_yn")]

cluster_all <- rbind(cluster_01, cluster_12, cluster_23, cluster_34, cluster_45, cluster_56, cluster_68, cluster_810, cluster_1012)

counts <- cluster_all %>% dplyr::count(mc_yn) 
cluster_all <- merge(cluster_all, counts, by="mc_yn")
cluster_all$pct <- (1/(cluster_all$n))*2
cluster_all$timing <- factor(cluster_all$timing, levels=c("Pre","Post"))

alpha_mc <- c(1.0,1.0,0.3,0.3,0.3,0.3,0.3,0.3,1.0,1.0,0.3,0.3,0.3,0.3,0.3,0.3,
              1.0,1.0,0.3,0.3,0.3,0.3,0.3,0.3,1.0,1.0,0.3,0.3,0.3,0.3,0.3)
# There are only 31 arguments because there are no OTH samples in the post-Mcat acquisition group
title_mc <- expression(paste(italic("M. catarrhalis"), " acquisition"))
mc_alluvium <- ggplot(cluster_all, aes(x = timing, stratum = cluster, alluvium = study_id, y = pct, fill = cluster)) +
  scale_fill_manual(values=cluster_cols_9) +
  geom_flow(alpha=0.3) +
  geom_stratum(width=0.38, color="gray20", alpha=alpha_mc) + ylab("Proportion of samples") + xlab("") +
  theme_barplot + facet_wrap(~mc_yn) + ggtitle(title_mc) +
  guides(fill=guide_legend(override.aes = list(alpha = c(0.3,0.3,0.3,0.3,0.3,0.3,1.0,1.0))))

png(file="R_Plots/Figure_S2/Figure_S2_mc.png", width = 3.75, height = 2.75, units = 'in', res = 300)
plot(mc_alluvium)
dev.off()

# S. aureus colonization

cluster_sa <- metadata_16s[,c("study_id", "month", "cluster", "inf_multi_sa")]
cluster_01 <- subset(cluster_sa, (month=="0" | month=="1") & !is.na(inf_multi_sa))
cluster_01 <- cluster_01[order(cluster_01$study_id, cluster_01$month),]
cluster_01 <- subset(cluster_01, month!="0" | inf_multi_sa=="N")
cluster_01 <- cluster_01 %>% group_by(study_id) %>% filter(n() == 2)
cluster_01$timing <- ""
cluster_01$timing[cluster_01$month=="0"] <- "Pre"
cluster_01$timing[cluster_01$month=="1"] <- "Post"
cluster_01 <- suppressWarnings(slide(cluster_01, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_sa_lag", slideBy = 1))
cluster_01$sa_yn[cluster_01$inf_multi_sa=="Y" | cluster_01$inf_sa_lag=="Y"] <- "Yes"
cluster_01$sa_yn[cluster_01$inf_multi_sa=="N" & (cluster_01$inf_sa_lag=="N" | is.na(cluster_01$inf_sa_lag))] <- "No"
cluster_01$study_id <- paste0(cluster_01$study_id, "_01")
cluster_01 <- cluster_01[,c("study_id", "timing", "cluster", "sa_yn")]

cluster_12 <- subset(cluster_sa, (month=="1" | month=="2") & !is.na(inf_multi_sa))
cluster_12 <- cluster_12[order(cluster_12$study_id, cluster_12$month),]
cluster_12 <- subset(cluster_12, month!="1" | inf_multi_sa=="N")
cluster_12 <- cluster_12 %>% group_by(study_id) %>% filter(n() == 2)
cluster_12$timing <- ""
cluster_12$timing[cluster_12$month=="1"] <- "Pre"
cluster_12$timing[cluster_12$month=="2"] <- "Post"
cluster_12 <- suppressWarnings(slide(cluster_12, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_sa_lag", slideBy = 1))
cluster_12$sa_yn[cluster_12$inf_multi_sa=="Y" | cluster_12$inf_sa_lag=="Y"] <- "Yes"
cluster_12$sa_yn[cluster_12$inf_multi_sa=="N" & (cluster_12$inf_sa_lag=="N" | is.na(cluster_12$inf_sa_lag))] <- "No"
cluster_12$study_id <- paste0(cluster_12$study_id, "_12")
cluster_12 <- cluster_12[,c("study_id", "timing", "cluster", "sa_yn")]

cluster_23 <- subset(cluster_sa, (month=="2" | month=="3") & !is.na(inf_multi_sa))
cluster_23 <- cluster_23[order(cluster_23$study_id, cluster_23$month),]
cluster_23 <- subset(cluster_23, month!="2" | inf_multi_sa=="N")
cluster_23 <- cluster_23 %>% group_by(study_id) %>% filter(n() == 2)
cluster_23$timing <- ""
cluster_23$timing[cluster_23$month=="2"] <- "Pre"
cluster_23$timing[cluster_23$month=="3"] <- "Post"
cluster_23 <- suppressWarnings(slide(cluster_23, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_sa_lag", slideBy = 1))
cluster_23$sa_yn[cluster_23$inf_multi_sa=="Y" | cluster_23$inf_sa_lag=="Y"] <- "Yes"
cluster_23$sa_yn[cluster_23$inf_multi_sa=="N" & (cluster_23$inf_sa_lag=="N" | is.na(cluster_23$inf_sa_lag))] <- "No"
cluster_23$study_id <- paste0(cluster_23$study_id, "_23")
cluster_23 <- cluster_23[,c("study_id", "timing", "cluster", "sa_yn")]

cluster_34 <- subset(cluster_sa, (month=="3" | month=="4") & !is.na(inf_multi_sa))
cluster_34 <- cluster_34[order(cluster_34$study_id, cluster_34$month),]
cluster_34 <- subset(cluster_34, month!="3" | inf_multi_sa=="N")
cluster_34 <- cluster_34 %>% group_by(study_id) %>% filter(n() == 2)
cluster_34$timing <- ""
cluster_34$timing[cluster_34$month=="3"] <- "Pre"
cluster_34$timing[cluster_34$month=="4"] <- "Post"
cluster_34 <- suppressWarnings(slide(cluster_34, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_sa_lag", slideBy = 1))
cluster_34$sa_yn[cluster_34$inf_multi_sa=="Y" | cluster_34$inf_sa_lag=="Y"] <- "Yes"
cluster_34$sa_yn[cluster_34$inf_multi_sa=="N" & (cluster_34$inf_sa_lag=="N" | is.na(cluster_34$inf_sa_lag))] <- "No"
cluster_34$study_id <- paste0(cluster_34$study_id, "_34")
cluster_34 <- cluster_34[,c("study_id", "timing", "cluster", "sa_yn")]

cluster_45 <- subset(cluster_sa, (month=="4" | month=="5") & !is.na(inf_multi_sa))
cluster_45 <- cluster_45[order(cluster_45$study_id, cluster_45$month),]
cluster_45 <- subset(cluster_45, month!="4" | inf_multi_sa=="N")
cluster_45 <- cluster_45 %>% group_by(study_id) %>% filter(n() == 2)
cluster_45$timing <- ""
cluster_45$timing[cluster_45$month=="4"] <- "Pre"
cluster_45$timing[cluster_45$month=="5"] <- "Post"
cluster_45 <- suppressWarnings(slide(cluster_45, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_sa_lag", slideBy = 1))
cluster_45$sa_yn[cluster_45$inf_multi_sa=="Y" | cluster_45$inf_sa_lag=="Y"] <- "Yes"
cluster_45$sa_yn[cluster_45$inf_multi_sa=="N" & (cluster_45$inf_sa_lag=="N" | is.na(cluster_45$inf_sa_lag))] <- "No"
cluster_45$study_id <- paste0(cluster_45$study_id, "_45")
cluster_45 <- cluster_45[,c("study_id", "timing", "cluster", "sa_yn")]

cluster_56 <- subset(cluster_sa, (month=="5" | month=="6") & !is.na(inf_multi_sa))
cluster_56 <- cluster_56[order(cluster_56$study_id, cluster_56$month),]
cluster_56 <- subset(cluster_56, month!="5" | inf_multi_sa=="N")
cluster_56 <- cluster_56 %>% group_by(study_id) %>% filter(n() == 2)
cluster_56$timing <- ""
cluster_56$timing[cluster_56$month=="5"] <- "Pre"
cluster_56$timing[cluster_56$month=="6"] <- "Post"
cluster_56 <- suppressWarnings(slide(cluster_56, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_sa_lag", slideBy = 1))
cluster_56$sa_yn[cluster_56$inf_multi_sa=="Y" | cluster_56$inf_sa_lag=="Y"] <- "Yes"
cluster_56$sa_yn[cluster_56$inf_multi_sa=="N" & (cluster_56$inf_sa_lag=="N" | is.na(cluster_56$inf_sa_lag))] <- "No"
cluster_56$study_id <- paste0(cluster_56$study_id, "_56")
cluster_56 <- cluster_56[,c("study_id", "timing", "cluster", "sa_yn")]

cluster_68 <- subset(cluster_sa, (month=="6" | month=="8") & !is.na(inf_multi_sa))
cluster_68 <- cluster_68[order(cluster_68$study_id, cluster_68$month),]
cluster_68 <- subset(cluster_68, month!="6" | inf_multi_sa=="N")
cluster_68 <- cluster_68 %>% group_by(study_id) %>% filter(n() == 2)
cluster_68$timing <- ""
cluster_68$timing[cluster_68$month=="6"] <- "Pre"
cluster_68$timing[cluster_68$month=="8"] <- "Post"
cluster_68 <- suppressWarnings(slide(cluster_68, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_sa_lag", slideBy = 1))
cluster_68$sa_yn[cluster_68$inf_multi_sa=="Y" | cluster_68$inf_sa_lag=="Y"] <- "Yes"
cluster_68$sa_yn[cluster_68$inf_multi_sa=="N" & (cluster_68$inf_sa_lag=="N" | is.na(cluster_68$inf_sa_lag))] <- "No"
cluster_68$study_id <- paste0(cluster_68$study_id, "_68")
cluster_68 <- cluster_68[,c("study_id", "timing", "cluster", "sa_yn")]

cluster_810 <- subset(cluster_sa, (month=="8" | month=="10") & !is.na(inf_multi_sa))
cluster_810 <- cluster_810[order(cluster_810$study_id, cluster_810$month),]
cluster_810 <- subset(cluster_810, month!="8" | inf_multi_sa=="N")
cluster_810 <- cluster_810 %>% group_by(study_id) %>% filter(n() == 2)
cluster_810$timing <- ""
cluster_810$timing[cluster_810$month=="8"] <- "Pre"
cluster_810$timing[cluster_810$month=="10"] <- "Post"
cluster_810 <- suppressWarnings(slide(cluster_810, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_sa_lag", slideBy = 1))
cluster_810$sa_yn[cluster_810$inf_multi_sa=="Y" | cluster_810$inf_sa_lag=="Y"] <- "Yes"
cluster_810$sa_yn[cluster_810$inf_multi_sa=="N" & (cluster_810$inf_sa_lag=="N" | is.na(cluster_810$inf_sa_lag))] <- "No"
cluster_810$study_id <- paste0(cluster_810$study_id, "_810")
cluster_810 <- cluster_810[,c("study_id", "timing", "cluster", "sa_yn")]

cluster_1012 <- subset(cluster_sa, (month=="10" | month=="12") & !is.na(inf_multi_sa))
cluster_1012 <- cluster_1012[order(cluster_1012$study_id, cluster_1012$month),]
cluster_1012 <- subset(cluster_1012, month!="10" | inf_multi_sa=="N")
cluster_1012 <- cluster_1012 %>% group_by(study_id) %>% filter(n() == 2)
cluster_1012$timing <- ""
cluster_1012$timing[cluster_1012$month=="10"] <- "Pre"
cluster_1012$timing[cluster_1012$month=="12"] <- "Post"
cluster_1012 <- suppressWarnings(slide(cluster_1012, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_sa_lag", slideBy = 1))
cluster_1012$sa_yn[cluster_1012$inf_multi_sa=="Y" | cluster_1012$inf_sa_lag=="Y"] <- "Yes"
cluster_1012$sa_yn[cluster_1012$inf_multi_sa=="N" & (cluster_1012$inf_sa_lag=="N" | is.na(cluster_1012$inf_sa_lag))] <- "No"
cluster_1012$study_id <- paste0(cluster_1012$study_id, "_1012")
cluster_1012 <- cluster_1012[,c("study_id", "timing", "cluster", "sa_yn")]

cluster_all <- rbind(cluster_01, cluster_12, cluster_23, cluster_34, cluster_45, cluster_56, cluster_68, cluster_810, cluster_1012)

counts <- cluster_all %>% dplyr::count(sa_yn) 
cluster_all <- merge(cluster_all, counts, by="sa_yn")
cluster_all$pct <- (1/(cluster_all$n))*2
cluster_all$timing <- factor(cluster_all$timing, levels=c("Pre","Post"))

alpha_sa <- c(0.3,0.3,0.3,0.3,0.3,0.3,1.0,0.3,0.3,0.3,0.3,0.3,0.3,0.3,1.0,0.3,
              0.3,0.3,0.3,0.3,0.3,0.3,1.0,0.3,0.3,0.3,0.3,0.3,0.3,0.3,1.0)
# There are only 31 arguments because there are no OTH samples in the post-SA acquisition group
title_sa <- expression(paste(italic("S. aureus"), " acquisition"))
sa_alluvium <- ggplot(cluster_all, aes(x = timing, stratum = cluster, alluvium = study_id, y = pct, fill = cluster)) +
  scale_fill_manual(values=cluster_cols_9) +
  geom_flow(alpha=0.3) +
  geom_stratum(width=0.38, color="gray20", alpha=alpha_sa) + ylab("Proportion of samples") + xlab("") +
  theme_barplot + facet_wrap(~sa_yn) + ggtitle(title_sa) +
guides(fill=guide_legend(override.aes = list(alpha = c(0.3,1.0,0.3,0.3,0.3,0.3,0.3,0.3))))

png(file="R_Plots/Figure_S2/Figure_S2_sa.png", width = 3.75, height = 2.75, units = 'in', res = 300)
plot(sa_alluvium)
dev.off()

# S. pneumoniae colonization

cluster_sp <- metadata_16s[,c("study_id", "month", "cluster", "inf_multi_sp")]
cluster_01 <- subset(cluster_sp, (month=="0" | month=="1") & !is.na(inf_multi_sp))
cluster_01 <- cluster_01[order(cluster_01$study_id, cluster_01$month),]
cluster_01 <- subset(cluster_01, month!="0" | inf_multi_sp=="N")
cluster_01 <- cluster_01 %>% group_by(study_id) %>% filter(n() == 2)
cluster_01$timing <- ""
cluster_01$timing[cluster_01$month=="0"] <- "Pre"
cluster_01$timing[cluster_01$month=="1"] <- "Post"
cluster_01 <- suppressWarnings(slide(cluster_01, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_sp_lag", slideBy = 1))
cluster_01$sp_yn[cluster_01$inf_multi_sp=="Y" | cluster_01$inf_sp_lag=="Y"] <- "Yes"
cluster_01$sp_yn[cluster_01$inf_multi_sp=="N" & (cluster_01$inf_sp_lag=="N" | is.na(cluster_01$inf_sp_lag))] <- "No"
cluster_01$study_id <- paste0(cluster_01$study_id, "_01")
cluster_01 <- cluster_01[,c("study_id", "timing", "cluster", "sp_yn")]

cluster_12 <- subset(cluster_sp, (month=="1" | month=="2") & !is.na(inf_multi_sp))
cluster_12 <- cluster_12[order(cluster_12$study_id, cluster_12$month),]
cluster_12 <- subset(cluster_12, month!="1" | inf_multi_sp=="N")
cluster_12 <- cluster_12 %>% group_by(study_id) %>% filter(n() == 2)
cluster_12$timing <- ""
cluster_12$timing[cluster_12$month=="1"] <- "Pre"
cluster_12$timing[cluster_12$month=="2"] <- "Post"
cluster_12 <- suppressWarnings(slide(cluster_12, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_sp_lag", slideBy = 1))
cluster_12$sp_yn[cluster_12$inf_multi_sp=="Y" | cluster_12$inf_sp_lag=="Y"] <- "Yes"
cluster_12$sp_yn[cluster_12$inf_multi_sp=="N" & (cluster_12$inf_sp_lag=="N" | is.na(cluster_12$inf_sp_lag))] <- "No"
cluster_12$study_id <- paste0(cluster_12$study_id, "_12")
cluster_12 <- cluster_12[,c("study_id", "timing", "cluster", "sp_yn")]

cluster_23 <- subset(cluster_sp, (month=="2" | month=="3") & !is.na(inf_multi_sp))
cluster_23 <- cluster_23[order(cluster_23$study_id, cluster_23$month),]
cluster_23 <- subset(cluster_23, month!="2" | inf_multi_sp=="N")
cluster_23 <- cluster_23 %>% group_by(study_id) %>% filter(n() == 2)
cluster_23$timing <- ""
cluster_23$timing[cluster_23$month=="2"] <- "Pre"
cluster_23$timing[cluster_23$month=="3"] <- "Post"
cluster_23 <- suppressWarnings(slide(cluster_23, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_sp_lag", slideBy = 1))
cluster_23$sp_yn[cluster_23$inf_multi_sp=="Y" | cluster_23$inf_sp_lag=="Y"] <- "Yes"
cluster_23$sp_yn[cluster_23$inf_multi_sp=="N" & (cluster_23$inf_sp_lag=="N" | is.na(cluster_23$inf_sp_lag))] <- "No"
cluster_23$study_id <- paste0(cluster_23$study_id, "_23")
cluster_23 <- cluster_23[,c("study_id", "timing", "cluster", "sp_yn")]

cluster_34 <- subset(cluster_sp, (month=="3" | month=="4") & !is.na(inf_multi_sp))
cluster_34 <- cluster_34[order(cluster_34$study_id, cluster_34$month),]
cluster_34 <- subset(cluster_34, month!="3" | inf_multi_sp=="N")
cluster_34 <- cluster_34 %>% group_by(study_id) %>% filter(n() == 2)
cluster_34$timing <- ""
cluster_34$timing[cluster_34$month=="3"] <- "Pre"
cluster_34$timing[cluster_34$month=="4"] <- "Post"
cluster_34 <- suppressWarnings(slide(cluster_34, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_sp_lag", slideBy = 1))
cluster_34$sp_yn[cluster_34$inf_multi_sp=="Y" | cluster_34$inf_sp_lag=="Y"] <- "Yes"
cluster_34$sp_yn[cluster_34$inf_multi_sp=="N" & (cluster_34$inf_sp_lag=="N" | is.na(cluster_34$inf_sp_lag))] <- "No"
cluster_34$study_id <- paste0(cluster_34$study_id, "_34")
cluster_34 <- cluster_34[,c("study_id", "timing", "cluster", "sp_yn")]

cluster_45 <- subset(cluster_sp, (month=="4" | month=="5") & !is.na(inf_multi_sp))
cluster_45 <- cluster_45[order(cluster_45$study_id, cluster_45$month),]
cluster_45 <- subset(cluster_45, month!="4" | inf_multi_sp=="N")
cluster_45 <- cluster_45 %>% group_by(study_id) %>% filter(n() == 2)
cluster_45$timing <- ""
cluster_45$timing[cluster_45$month=="4"] <- "Pre"
cluster_45$timing[cluster_45$month=="5"] <- "Post"
cluster_45 <- suppressWarnings(slide(cluster_45, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_sp_lag", slideBy = 1))
cluster_45$sp_yn[cluster_45$inf_multi_sp=="Y" | cluster_45$inf_sp_lag=="Y"] <- "Yes"
cluster_45$sp_yn[cluster_45$inf_multi_sp=="N" & (cluster_45$inf_sp_lag=="N" | is.na(cluster_45$inf_sp_lag))] <- "No"
cluster_45$study_id <- paste0(cluster_45$study_id, "_45")
cluster_45 <- cluster_45[,c("study_id", "timing", "cluster", "sp_yn")]

cluster_56 <- subset(cluster_sp, (month=="5" | month=="6") & !is.na(inf_multi_sp))
cluster_56 <- cluster_56[order(cluster_56$study_id, cluster_56$month),]
cluster_56 <- subset(cluster_56, month!="5" | inf_multi_sp=="N")
cluster_56 <- cluster_56 %>% group_by(study_id) %>% filter(n() == 2)
cluster_56$timing <- ""
cluster_56$timing[cluster_56$month=="5"] <- "Pre"
cluster_56$timing[cluster_56$month=="6"] <- "Post"
cluster_56 <- suppressWarnings(slide(cluster_56, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_sp_lag", slideBy = 1))
cluster_56$sp_yn[cluster_56$inf_multi_sp=="Y" | cluster_56$inf_sp_lag=="Y"] <- "Yes"
cluster_56$sp_yn[cluster_56$inf_multi_sp=="N" & (cluster_56$inf_sp_lag=="N" | is.na(cluster_56$inf_sp_lag))] <- "No"
cluster_56$study_id <- paste0(cluster_56$study_id, "_56")
cluster_56 <- cluster_56[,c("study_id", "timing", "cluster", "sp_yn")]

cluster_68 <- subset(cluster_sp, (month=="6" | month=="8") & !is.na(inf_multi_sp))
cluster_68 <- cluster_68[order(cluster_68$study_id, cluster_68$month),]
cluster_68 <- subset(cluster_68, month!="6" | inf_multi_sp=="N")
cluster_68 <- cluster_68 %>% group_by(study_id) %>% filter(n() == 2)
cluster_68$timing <- ""
cluster_68$timing[cluster_68$month=="6"] <- "Pre"
cluster_68$timing[cluster_68$month=="8"] <- "Post"
cluster_68 <- suppressWarnings(slide(cluster_68, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_sp_lag", slideBy = 1))
cluster_68$sp_yn[cluster_68$inf_multi_sp=="Y" | cluster_68$inf_sp_lag=="Y"] <- "Yes"
cluster_68$sp_yn[cluster_68$inf_multi_sp=="N" & (cluster_68$inf_sp_lag=="N" | is.na(cluster_68$inf_sp_lag))] <- "No"
cluster_68$study_id <- paste0(cluster_68$study_id, "_68")
cluster_68 <- cluster_68[,c("study_id", "timing", "cluster", "sp_yn")]

cluster_810 <- subset(cluster_sp, (month=="8" | month=="10") & !is.na(inf_multi_sp))
cluster_810 <- cluster_810[order(cluster_810$study_id, cluster_810$month),]
cluster_810 <- subset(cluster_810, month!="8" | inf_multi_sp=="N")
cluster_810 <- cluster_810 %>% group_by(study_id) %>% filter(n() == 2)
cluster_810$timing <- ""
cluster_810$timing[cluster_810$month=="8"] <- "Pre"
cluster_810$timing[cluster_810$month=="10"] <- "Post"
cluster_810 <- suppressWarnings(slide(cluster_810, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_sp_lag", slideBy = 1))
cluster_810$sp_yn[cluster_810$inf_multi_sp=="Y" | cluster_810$inf_sp_lag=="Y"] <- "Yes"
cluster_810$sp_yn[cluster_810$inf_multi_sp=="N" & (cluster_810$inf_sp_lag=="N" | is.na(cluster_810$inf_sp_lag))] <- "No"
cluster_810$study_id <- paste0(cluster_810$study_id, "_810")
cluster_810 <- cluster_810[,c("study_id", "timing", "cluster", "sp_yn")]

cluster_1012 <- subset(cluster_sp, (month=="10" | month=="12") & !is.na(inf_multi_sp))
cluster_1012 <- cluster_1012[order(cluster_1012$study_id, cluster_1012$month),]
cluster_1012 <- subset(cluster_1012, month!="10" | inf_multi_sp=="N")
cluster_1012 <- cluster_1012 %>% group_by(study_id) %>% filter(n() == 2)
cluster_1012$timing <- ""
cluster_1012$timing[cluster_1012$month=="10"] <- "Pre"
cluster_1012$timing[cluster_1012$month=="12"] <- "Post"
cluster_1012 <- suppressWarnings(slide(cluster_1012, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_sp_lag", slideBy = 1))
cluster_1012$sp_yn[cluster_1012$inf_multi_sp=="Y" | cluster_1012$inf_sp_lag=="Y"] <- "Yes"
cluster_1012$sp_yn[cluster_1012$inf_multi_sp=="N" & (cluster_1012$inf_sp_lag=="N" | is.na(cluster_1012$inf_sp_lag))] <- "No"
cluster_1012$study_id <- paste0(cluster_1012$study_id, "_1012")
cluster_1012 <- cluster_1012[,c("study_id", "timing", "cluster", "sp_yn")]

cluster_all <- rbind(cluster_01, cluster_12, cluster_23, cluster_34, cluster_45, cluster_56, cluster_68, cluster_810, cluster_1012)

counts <- cluster_all %>% dplyr::count(sp_yn) 
cluster_all <- merge(cluster_all, counts, by="sp_yn")
cluster_all$pct <- (1/(cluster_all$n))*2
cluster_all$timing <- factor(cluster_all$timing, levels=c("Pre","Post"))

alpha_sp <- c(0.3,0.3,0.3,0.3,1.0,0.3,0.3,0.3,0.3,0.3,0.3,0.3,1.0,0.3,0.3,0.3,
              0.3,0.3,0.3,0.3,1.0,0.3,0.3,0.3,0.3,0.3,0.3,0.3,1.0,0.3,0.3)
# There are only 31 arguments because there are no OTH samples in the post-S. pneumoniae acquisition group
title_sp <- expression(paste(italic("S. pneumoniae"), " acquisition"))
sp_alluvium <- ggplot(cluster_all, aes(x = timing, stratum = cluster, alluvium = study_id, y = pct, fill = cluster)) +
  scale_fill_manual(values=cluster_cols_9) +
  geom_flow(alpha=0.3) +
  geom_stratum(width=0.38, color="gray20", alpha=alpha_sp) + ylab("Proportion of samples") + xlab("") +
  theme_barplot + facet_wrap(~sp_yn) + ggtitle(title_sp) +
  guides(fill=guide_legend(override.aes = list(alpha = c(0.3,0.3,0.3,1.0,0.3,0.3,0.3,0.3))))

png(file="R_Plots/Figure_S2/Figure_S2_sp.png", width = 3.75, height = 2.75, units = 'in', res = 300)
plot(sp_alluvium)
dev.off()

png(file="R_Plots/Figure_S2.png", width = 7.5, height = 5.5, units = 'in', res = 1200)
plot_grid(hi_alluvium, NULL, mc_alluvium, NULL, NULL, NULL, sa_alluvium, NULL, sp_alluvium, 
          labels=c("a","","b","","","","c","","d"), nrow=3, rel_heights=c(1,-0.05,1), ncol=3, rel_widths=c(1,0.05,1), label_size=9) 
dev.off()