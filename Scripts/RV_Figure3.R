# Botswana Infant Microbiome Study - RV-Bacterial Analyses
# Matthew Kelly, MD, MPH 
# Figure 3
# Last updated: March 9, 2025

remove(list=ls())
setwd("___________________") 
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
cluster_cols_8 <- c("gray10", # OTH
                    "mediumorchid4", # STA
                    "darkslateblue", # COR
                    "coral3", # STR
                    "indianred4", # CD
                    "peru", # HAE
                    "dodgerblue3", # CDM
                    "darkolivegreen4") # MOR

# Load required datasets
phy.inf.np.16s <- read_rds("phy.inf.np.16s.rds")
metadata_inf_np_16s <- data.frame(sample_data(phy.inf.np.16s))
phy.agglom <- tax_glom(phy.inf.np.16s, taxrank = 'Genus')
ntaxa(phy.agglom)
phy.relative <- transform_sample_counts(phy.agglom, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.relative))

# Fig. 3a - create an NMDS plot on Bray-Curtis distances

ord <- ordinate(phy.relative, method="NMDS", distance="bray")

ord_theme   <-  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                      axis.text.y = element_text(size=8, color="black"), 
                      axis.title.y = element_text(size=9, color="black"),
                      axis.text.x  = element_text(size=8, color="black"),
                      axis.title.x = element_text(size=9, color="black"),
                      legend.position = "right", legend.text = element_text(size=8, color="black"), 
                      legend.title = element_blank(), legend.key = element_blank(), plot.margin = unit(c(0.5, -0.2, 0, 0), "cm"),
                      legend.box.spacing = unit(0, "pt"), legend.key.size = unit(0.4, 'cm'))

cluster_nmds <- plot_ordination(phy.relative, ord, color="cluster") +
  geom_point(shape=16, size=1.5) + ord_theme + scale_color_manual(values=cluster_cols_8) +
  stat_ellipse(geom="polygon", alpha=0, type="t", level=0.9, linewidth=0.5)  

png(file="R_Plots/Figure_3/Figure_3a.png", width = 3.375, height = 2.8, units = 'in', res = 1200)
plot(cluster_nmds)
dev.off()

# Fig. 3b - cluster barplot

# Create dataframes with overall relative abundances of phyla and genera
relative_inf <- psmelt(phy.relative)
relative_inf$Phylum <- as.character(relative_inf$Phylum)
phyla_abundances <- aggregate(relative_inf$Abundance, by=list(Phylum=relative_inf$Phylum), FUN=sum)
phyla_abundances$x <- as.numeric((phyla_abundances$x)/(nsamples(phy.relative)))
names(phyla_abundances)[names(phyla_abundances) == "x"] <- "phyla_Ab"
phyla_abundances <- phyla_abundances[order(-phyla_abundances$phyla_Ab),] 
head(phyla_abundances, 5)
sum(phyla_abundances$phyla_Ab)    # Should sum to 1 (sum of relative abundances of phyla)
nrow(phyla_abundances)            # Corresponds to # of unique phyla
relative_inf$Genus <- as.character(relative_inf$Genus)
genera_abundances <- aggregate(relative_inf$Abundance, by=list(Phylum=relative_inf$Phylum, Genus=relative_inf$Genus,
                                                                      OTU=relative_inf$OTU), FUN=mean)
genera_abundances <- aggregate(genera_abundances$x, by=list(Phylum=genera_abundances$Phylum,
                                                            Genus=genera_abundances$Genus), FUN=sum)
names(genera_abundances)[names(genera_abundances) == "x"] <- "genus_Ab"
sum(genera_abundances$genus_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(genera_abundances)            # Corresponds to # of unique genera
genera_abundances <- genera_abundances[order(-genera_abundances$genus_Ab),]  
TOPGenera <- unique(genera_abundances$Genus[1:18])
genus_df <- genera_abundances[genera_abundances$Genus %in% TOPGenera,]
genus_df$Genus <- factor(genus_df$Genus, levels = genus_df$Genus[order(-genus_df$genus_Ab)])
head(genera_abundances, 25)

# Rename genera other than top genera as "Other" in creating dataframe relative_inf 
relative_inf$Genus[!(relative_inf$Genus %in% TOPGenera)] <- "Other"
TOPGenera
relative_inf$Genus[relative_inf$Genus=="Bacteria;p;c;o;f;g"] <- "Other"
relative_inf$Genus[relative_inf$Genus=="Enterobacteriaceae;g"] <- "Other"
relative_inf$Genus[relative_inf$Genus=="Alphaproteobacteria;o;f;g"] <- "Other"
relative_inf$Genus[relative_inf$Genus=="Actinobacteria;o;f;g"] <- "Other"
relative_inf$Genus[relative_inf$Genus=="Neisseriaceae [G-1]"] <- "Other"
unique(relative_inf$Genus)
sum(relative_inf$Abundance)
relative_inf$month <- as.factor(relative_inf$month)
relative_inf$month2 <- ""
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="0"] <- "I0"
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="1"] <- "I1"
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="2"] <- "I2"
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="3"] <- "I3"
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="4"] <- "I4"
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="5"] <- "I5"
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="6"] <- "I6"
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="8"] <- "I8"
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="10"] <- "I10"
relative_inf$month2[relative_inf$Source=="INF" & relative_inf$month=="12"] <- "I12"
relative_inf$month2 <- factor(relative_inf$month2, levels=c("I0", "I1", "I2", "I3", "I4", "I5", "I6", "I8", "I10", "I12"))
relative_inf$Genus <- as.factor(relative_inf$Genus)
table(relative_inf$Genus)
relative_inf$Genus <- factor(relative_inf$Genus, levels=c("Acinetobacter", "Corynebacterium", "Dolosigranulum", "Gemella", 
                                                          "Haemophilus", "Lactobacillus", "Micrococcus", "Moraxella", "Neisseria", "Prevotella", 
                                                          "Pseudomonas", "Staphylococcus", "Streptococcus", "Veillonella", "Other"))

theme_barplot <-   theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                         axis.title.x = element_text(size=9, color="black"), 
                         axis.text.y = element_text(size=8, color="black"), 
                         axis.title.y = element_text(angle=90, size=9, color="black"), 
                         axis.text.x = element_text(size=8, color="black"),
                         legend.position = "right", legend.text=element_text(size=8, color="black"), legend.text.align = 0,
                         legend.box.spacing = unit(0, "pt"),
                         legend.title=element_blank(), legend.key.size = unit(0.35, 'cm'), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) 

earthy_cols_15 <- c("rosybrown2", # Acinetobacter
                    "darkslateblue", # Corynebacterium
                    "indianred4", # Dolosigranulum
                    "cornflowerblue", # Gemella
                    "peru", # Haemophilus
                    "indianred1", # Lactobacillus
                    "gray50", # Micrococcus                 
                    "darkolivegreen4", # Moraxella
                    "dodgerblue3", # Neisseria
                    "plum3", # Prevotella
                    "navajowhite3", # Pseudomonas 
                    "mediumorchid4", # Staphylococcus 
                    "coral3", # Streptococcus
                    "forestgreen", # Veillonella
                    "gray10") # Other

cluster_barplot <- ggplot(arrange(relative_inf, Genus), aes(x=cluster, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="fill") + theme_barplot + 
  scale_fill_manual(values=earthy_cols_15, 
                    labels=c(expression(italic("Acinetobacter")), expression(italic("Corynebacterium")), expression(italic("Dolosigranulum")), expression(italic("Gemella")), 
                             expression(italic("Haemophilus")), expression(italic("Lactobacillus")), expression(italic("Micrococcus")), expression(italic("Moraxella")), 
                             expression(italic("Neisseria")), expression(italic("Prevotella")), expression(italic("Pseudomonas")), expression(italic("Staphylococcus")), 
                             expression(italic("Streptococcus")), expression(italic("Veillonella")), "Other")) + 
  xlab("Microbiota type") + ylab("Relative abundance") 

png(file="R_Plots/Figure_3/Figure_3b.png", width = 4.125, height = 2.8, units = 'in', res = 1200)
plot(cluster_barplot)
dev.off()

# Fig. 3c - alluvial diagram for microbiome cluster transitions by month

library(ggplot2)
library(ggalluvial)
require(dplyr)
library(formattable)

cluster_month <- metadata_inf_np_16s[,c("study_id", "month", "cluster", "inf_rv_yn")]
counts <- cluster_month %>% dplyr::count(month) 
cluster_month <- merge(cluster_month, counts, by="month")
cluster_month$pct <- (1/(cluster_month$n))
cluster_month$month2[cluster_month$month=="0"] <- "m0"
cluster_month$month2[cluster_month$month=="1"] <- "m1"
cluster_month$month2[cluster_month$month=="2"] <- "m2"
cluster_month$month2[cluster_month$month=="3"] <- "m3"
cluster_month$month2[cluster_month$month=="4"] <- "m4"
cluster_month$month2[cluster_month$month=="5"] <- "m5"
cluster_month$month2[cluster_month$month=="6"] <- "m6"
cluster_month$month2[cluster_month$month=="8"] <- "m8"
cluster_month$month2[cluster_month$month=="10"] <- "m10"
cluster_month$month2[cluster_month$month=="12"] <- "m12"
cluster_month$month2 <- factor(cluster_month$month2, levels=c("m0","m1","m2","m3","m4","m5","m6","m8","m10","m12"))

theme_barplot2 <-   theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                          axis.text.y = element_text(size=8, color="black"), 
                          axis.title.y = element_text(angle=90, size=9, color="black"),
                          axis.text.x = element_text(size=8.5, color="black"),
                          axis.title.x = element_text(size=9, color="black"), 
                          legend.position = "right", legend.text=element_text(size=8, color="black"), legend.title=element_blank(),
                          legend.box.spacing = unit(0, "pt"), plot.margin = unit(c(0.5, 0.2, 0, 0), "cm"),
                          legend.key.size = unit(0.4, 'cm'))

cluster_alluvium <- ggplot(cluster_month, aes(x = month2, stratum = cluster, alluvium = study_id, y = pct, fill = cluster)) +
  scale_fill_manual(values=cluster_cols_8) +
  geom_flow() +
  geom_stratum(width=0.38, color="gray20") + ylab("Proportion of samples") + xlab("") +
  theme_barplot2

png(file="R_Plots/Figure_3/Figure_3c.png", width = 4.125, height = 3.2, units = 'in', res = 1200)
plot(cluster_alluvium)
dev.off()

# Fig. 3d - alluvial diagram for microbiome cluster transitions by respiratory virus infection

library(ggplot2)
library(ggalluvial)
require(dplyr)
library(formattable)

cluster_rv <- metadata_inf_np_16s[,c("study_id", "month", "cluster", "inf_rv_yn")]
cluster_01 <- subset(cluster_rv, (month=="0" | month=="1") & !is.na(inf_rv_yn))
cluster_01 <- cluster_01[order(cluster_01$study_id, cluster_01$month),]
cluster_01 <- subset(cluster_01, month!="0" | inf_rv_yn=="N")
cluster_01 <- cluster_01 %>% group_by(study_id) %>% filter(n() == 2)
cluster_01$timing <- ""
cluster_01$timing[cluster_01$month=="0"] <- "Pre"
cluster_01$timing[cluster_01$month=="1"] <- "Post"
cluster_01 <- suppressWarnings(slide(cluster_01, "inf_rv_yn", TimeVar="month", GroupVar="study_id", NewVar="inf_rv_lag", slideBy = 1))
cluster_01$rv_yn[cluster_01$inf_rv_yn=="Y" | cluster_01$inf_rv_lag=="Y"] <- "Yes"
cluster_01$rv_yn[cluster_01$inf_rv_yn=="N" & (cluster_01$inf_rv_lag=="N" | is.na(cluster_01$inf_rv_lag))] <- "No"
cluster_01$study_id <- paste0(cluster_01$study_id, "_01")
cluster_01 <- cluster_01[,c("study_id", "timing", "cluster", "rv_yn")]

cluster_12 <- subset(cluster_rv, (month=="1" | month=="2") & !is.na(inf_rv_yn))
cluster_12 <- cluster_12[order(cluster_12$study_id, cluster_12$month),]
cluster_12 <- subset(cluster_12, month!="1" | inf_rv_yn=="N")
cluster_12 <- cluster_12 %>% group_by(study_id) %>% filter(n() == 2)
cluster_12$timing <- ""
cluster_12$timing[cluster_12$month=="1"] <- "Pre"
cluster_12$timing[cluster_12$month=="2"] <- "Post"
cluster_12 <- suppressWarnings(slide(cluster_12, "inf_rv_yn", TimeVar="month", GroupVar="study_id", NewVar="inf_rv_lag", slideBy = 1))
cluster_12$rv_yn[cluster_12$inf_rv_yn=="Y" | cluster_12$inf_rv_lag=="Y"] <- "Yes"
cluster_12$rv_yn[cluster_12$inf_rv_yn=="N" & (cluster_12$inf_rv_lag=="N" | is.na(cluster_12$inf_rv_lag))] <- "No"
cluster_12$study_id <- paste0(cluster_12$study_id, "_12")
cluster_12 <- cluster_12[,c("study_id", "timing", "cluster", "rv_yn")]

cluster_23 <- subset(cluster_rv, (month=="2" | month=="3") & !is.na(inf_rv_yn))
cluster_23 <- cluster_23[order(cluster_23$study_id, cluster_23$month),]
cluster_23 <- subset(cluster_23, month!="2" | inf_rv_yn=="N")
cluster_23 <- cluster_23 %>% group_by(study_id) %>% filter(n() == 2)
cluster_23$timing <- ""
cluster_23$timing[cluster_23$month=="2"] <- "Pre"
cluster_23$timing[cluster_23$month=="3"] <- "Post"
cluster_23 <- suppressWarnings(slide(cluster_23, "inf_rv_yn", TimeVar="month", GroupVar="study_id", NewVar="inf_rv_lag", slideBy = 1))
cluster_23$rv_yn[cluster_23$inf_rv_yn=="Y" | cluster_23$inf_rv_lag=="Y"] <- "Yes"
cluster_23$rv_yn[cluster_23$inf_rv_yn=="N" & (cluster_23$inf_rv_lag=="N" | is.na(cluster_23$inf_rv_lag))] <- "No"
cluster_23$study_id <- paste0(cluster_23$study_id, "_23")
cluster_23 <- cluster_23[,c("study_id", "timing", "cluster", "rv_yn")]

cluster_34 <- subset(cluster_rv, (month=="3" | month=="4") & !is.na(inf_rv_yn))
cluster_34 <- cluster_34[order(cluster_34$study_id, cluster_34$month),]
cluster_34 <- subset(cluster_34, month!="3" | inf_rv_yn=="N")
cluster_34 <- cluster_34 %>% group_by(study_id) %>% filter(n() == 2)
cluster_34$timing <- ""
cluster_34$timing[cluster_34$month=="3"] <- "Pre"
cluster_34$timing[cluster_34$month=="4"] <- "Post"
cluster_34 <- suppressWarnings(slide(cluster_34, "inf_rv_yn", TimeVar="month", GroupVar="study_id", NewVar="inf_rv_lag", slideBy = 1))
cluster_34$rv_yn[cluster_34$inf_rv_yn=="Y" | cluster_34$inf_rv_lag=="Y"] <- "Yes"
cluster_34$rv_yn[cluster_34$inf_rv_yn=="N" & (cluster_34$inf_rv_lag=="N" | is.na(cluster_34$inf_rv_lag))] <- "No"
cluster_34$study_id <- paste0(cluster_34$study_id, "_34")
cluster_34 <- cluster_34[,c("study_id", "timing", "cluster", "rv_yn")]

cluster_45 <- subset(cluster_rv, (month=="4" | month=="5") & !is.na(inf_rv_yn))
cluster_45 <- cluster_45[order(cluster_45$study_id, cluster_45$month),]
cluster_45 <- subset(cluster_45, month!="4" | inf_rv_yn=="N")
cluster_45 <- cluster_45 %>% group_by(study_id) %>% filter(n() == 2)
cluster_45$timing <- ""
cluster_45$timing[cluster_45$month=="4"] <- "Pre"
cluster_45$timing[cluster_45$month=="5"] <- "Post"
cluster_45 <- suppressWarnings(slide(cluster_45, "inf_rv_yn", TimeVar="month", GroupVar="study_id", NewVar="inf_rv_lag", slideBy = 1))
cluster_45$rv_yn[cluster_45$inf_rv_yn=="Y" | cluster_45$inf_rv_lag=="Y"] <- "Yes"
cluster_45$rv_yn[cluster_45$inf_rv_yn=="N" & (cluster_45$inf_rv_lag=="N" | is.na(cluster_45$inf_rv_lag))] <- "No"
cluster_45$study_id <- paste0(cluster_45$study_id, "_45")
cluster_45 <- cluster_45[,c("study_id", "timing", "cluster", "rv_yn")]

cluster_56 <- subset(cluster_rv, (month=="5" | month=="6") & !is.na(inf_rv_yn))
cluster_56 <- cluster_56[order(cluster_56$study_id, cluster_56$month),]
cluster_56 <- subset(cluster_56, month!="5" | inf_rv_yn=="N")
cluster_56 <- cluster_56 %>% group_by(study_id) %>% filter(n() == 2)
cluster_56$timing <- ""
cluster_56$timing[cluster_56$month=="5"] <- "Pre"
cluster_56$timing[cluster_56$month=="6"] <- "Post"
cluster_56 <- suppressWarnings(slide(cluster_56, "inf_rv_yn", TimeVar="month", GroupVar="study_id", NewVar="inf_rv_lag", slideBy = 1))
cluster_56$rv_yn[cluster_56$inf_rv_yn=="Y" | cluster_56$inf_rv_lag=="Y"] <- "Yes"
cluster_56$rv_yn[cluster_56$inf_rv_yn=="N" & (cluster_56$inf_rv_lag=="N" | is.na(cluster_56$inf_rv_lag))] <- "No"
cluster_56$study_id <- paste0(cluster_56$study_id, "_56")
cluster_56 <- cluster_56[,c("study_id", "timing", "cluster", "rv_yn")]

cluster_68 <- subset(cluster_rv, (month=="6" | month=="8") & !is.na(inf_rv_yn))
cluster_68 <- cluster_68[order(cluster_68$study_id, cluster_68$month),]
cluster_68 <- subset(cluster_68, month!="6" | inf_rv_yn=="N")
cluster_68 <- cluster_68 %>% group_by(study_id) %>% filter(n() == 2)
cluster_68$timing <- ""
cluster_68$timing[cluster_68$month=="6"] <- "Pre"
cluster_68$timing[cluster_68$month=="8"] <- "Post"
cluster_68 <- suppressWarnings(slide(cluster_68, "inf_rv_yn", TimeVar="month", GroupVar="study_id", NewVar="inf_rv_lag", slideBy = 1))
cluster_68$rv_yn[cluster_68$inf_rv_yn=="Y" | cluster_68$inf_rv_lag=="Y"] <- "Yes"
cluster_68$rv_yn[cluster_68$inf_rv_yn=="N" & (cluster_68$inf_rv_lag=="N" | is.na(cluster_68$inf_rv_lag))] <- "No"
cluster_68$study_id <- paste0(cluster_68$study_id, "_68")
cluster_68 <- cluster_68[,c("study_id", "timing", "cluster", "rv_yn")]

cluster_810 <- subset(cluster_rv, (month=="8" | month=="10") & !is.na(inf_rv_yn))
cluster_810 <- cluster_810[order(cluster_810$study_id, cluster_810$month),]
cluster_810 <- subset(cluster_810, month!="8" | inf_rv_yn=="N")
cluster_810 <- cluster_810 %>% group_by(study_id) %>% filter(n() == 2)
cluster_810$timing <- ""
cluster_810$timing[cluster_810$month=="8"] <- "Pre"
cluster_810$timing[cluster_810$month=="10"] <- "Post"
cluster_810 <- suppressWarnings(slide(cluster_810, "inf_rv_yn", TimeVar="month", GroupVar="study_id", NewVar="inf_rv_lag", slideBy = 1))
cluster_810$rv_yn[cluster_810$inf_rv_yn=="Y" | cluster_810$inf_rv_lag=="Y"] <- "Yes"
cluster_810$rv_yn[cluster_810$inf_rv_yn=="N" & (cluster_810$inf_rv_lag=="N" | is.na(cluster_810$inf_rv_lag))] <- "No"
cluster_810$study_id <- paste0(cluster_810$study_id, "_810")
cluster_810 <- cluster_810[,c("study_id", "timing", "cluster", "rv_yn")]

cluster_1012 <- subset(cluster_rv, (month=="10" | month=="12") & !is.na(inf_rv_yn))
cluster_1012 <- cluster_1012[order(cluster_1012$study_id, cluster_1012$month),]
cluster_1012 <- subset(cluster_1012, month!="10" | inf_rv_yn=="N")
cluster_1012 <- cluster_1012 %>% group_by(study_id) %>% filter(n() == 2)
cluster_1012$timing <- ""
cluster_1012$timing[cluster_1012$month=="10"] <- "Pre"
cluster_1012$timing[cluster_1012$month=="12"] <- "Post"
cluster_1012 <- suppressWarnings(slide(cluster_1012, "inf_rv_yn", TimeVar="month", GroupVar="study_id", NewVar="inf_rv_lag", slideBy = 1))
cluster_1012$rv_yn[cluster_1012$inf_rv_yn=="Y" | cluster_1012$inf_rv_lag=="Y"] <- "Yes"
cluster_1012$rv_yn[cluster_1012$inf_rv_yn=="N" & (cluster_1012$inf_rv_lag=="N" | is.na(cluster_1012$inf_rv_lag))] <- "No"
cluster_1012$study_id <- paste0(cluster_1012$study_id, "_1012")
cluster_1012 <- cluster_1012[,c("study_id", "timing", "cluster", "rv_yn")]

cluster_all <- rbind(cluster_01, cluster_12, cluster_23, cluster_34, cluster_45, cluster_56, cluster_68, cluster_810, cluster_1012)

counts <- cluster_all %>% dplyr::count(rv_yn) 
cluster_all <- merge(cluster_all, counts, by="rv_yn")
cluster_all$pct <- (1/(cluster_all$n))*2
cluster_all$timing <- factor(cluster_all$timing, levels=c("Pre","Post"))

theme_barplot3 <-   theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                          axis.text.y = element_text(size=8, color="black"), 
                          axis.title.y = element_text(angle=90, size=9, color="black"), 
                          axis.text.x = element_text(size=8.5, color="black"),
                          axis.title.x = element_text(size=9, color="black"),
                          legend.position = "right", legend.text=element_text(size=8, color="black"), legend.title=element_blank(),
                          strip.background =element_rect(fill="grey20"), 
                          strip.text = element_text(size=8, colour = 'white', margin = ggplot2::margin(-0.01,0,-0.01,0, "cm")), 
                          plot.title = element_text(size = 9, hjust = 0.5), plot.margin = unit(c(0.3, 0, 0, 0), "cm"),
                          legend.key.size = unit(0.4, 'cm'), legend.box.spacing = unit(0, "pt"), panel.spacing = unit(0,'lines'))

rv_alluvium <- ggplot(cluster_all, aes(x = timing, stratum = cluster, alluvium = study_id, y = pct, fill = cluster)) +
  scale_fill_manual(values=cluster_cols_8) +
  geom_flow() +
  geom_stratum(width=0.38, color="gray20") + ylab("Proportion of samples") + xlab("") +
  theme_barplot3 + facet_wrap(~rv_yn) + ggtitle("Respiratory virus infection")

png(file="R_Plots/Figure_3/Figure_3d.png", width = 3.375, height = 3.2, units = 'in', res = 1200)
plot(rv_alluvium)
dev.off()

maaslin_results_rv <- read_csv('Pixu_Analyses/results_v2/maaslin_impact/inf_rv_yn/res_w_taxa.csv', show_col_types = FALSE)
maaslin_results_rv <- subset(maaslin_results_rv, metadata=="inf_rv_yn")
maaslin_results_rv <- maaslin_results_rv[,c("ASV", "Genus", "Species_eHOMD", "metadata","coef","coef_dir","stderr","pval","qval")]
maaslin_results_rv <- arrange(maaslin_results_rv, coef)
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV3"] <- "Staphylococcus sp. (3)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV6"] <- "C. accolens/macginleyi (6)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV4"] <- "C. propinquum/pseudodiphtheriticum (4)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV18"] <- "C. propinquum/pseudodiphtheriticum (18)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV14"] <- "Neisseriaceae [G-1] HMT-174 (14)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV15"] <- "S. thermophilus/vestibularis/salivarius (15)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV26"] <- "Micrococcus luteus (26)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV86"] <- "C.  tuberculostearicum (86)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV31"] <- "Neisseria sp. (31)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV1"] <- "M. catarrhalis/nonliquefaciens (1)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV30"] <- "Gemella sp. (30)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV8"] <- "Haemophilus sp. (8)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV7"] <- "Haemophilus sp. (7)"
maaslin_results_rv$ASV[maaslin_results_rv$ASV=="ASV5"] <- "S. pneumoniae/pseudopneumoniae (5)"

title_rv <- expression(paste("Respiratory virus infection"))
maaslin_rv <- maaslin_results_rv %>% ggplot(aes(x = coef, y = reorder(ASV, -coef), label = ASV, fill = coef_dir)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("darkslateblue","seagreen")) +
  theme(axis.text.y = element_text(size=8, color="black", face="italic", angle=5), plot.title = element_text(size=9, hjust = 0.5),
        axis.text.x = element_text(size=8, color="black"), axis.title.x = element_text(size=8.5, color="black"),
        legend.position = "none", axis.line.x = element_line(linewidth = 0.1, linetype = "solid", colour = "black"),
        panel.background = element_blank(), panel.grid.major.y = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(0.3, 0, 0, 0), "cm")) +
  labs(x = "MaAsLin2 coefficient", y = "", fill = "") + scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5), limits=c(-1.15, 1.77)) + ggtitle(title_rv)

png(file="R_Plots/Figure_3/Figure_3e.png", width = 5.5, height = 3, units = 'in', res = 1200)
plot(maaslin_rv)
dev.off()

fig_3ac <- plot_grid(cluster_nmds, NULL, cluster_alluvium, rel_widths=c(0.42, 0.03, 0.55), nrow=1, labels=c("a","", "c"), label_size=10)
fig_3bde <- plot_grid(cluster_barplot, NULL, rv_alluvium, maaslin_rv, rel_widths=c(0.35, 0.01, 0.28, 0.36), nrow=1, labels=c("b","","d","e"), label_size=10) 

png(file="R_Plots/Figure_3.png", width = 11.5, height = 7, units = 'in', res = 1200)
plot_grid(fig_3ac, fig_3bde, labels=NULL, nrow=2, rel_heights=c(1, 1)) 
dev.off()

# Save files as a Source Data file
source_data <- list('Fig3b'=relative_inf, 'Fig3c'=cluster_month, 'Fig3d'=cluster_all, 'Fig3e'=maaslin_results_rv) 
openxlsx::write.xlsx(source_data, file="Source_Data/Figure_3.xlsx")