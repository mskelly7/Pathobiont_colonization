# Botswana Infant Microbiome Study - RV-Bacterial Analyses
# Matthew Kelly, MD, MPH 
# Supplementary Figure 3
# Last updated: October 18, 2024

remove(list=ls())
setwd("G:/My Drive/Research/SAS_Bots_Microbiome/RV_Pathogen") 
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

countdata <- read.csv("G:/My Drive/Research/SAS_Bots_Microbiome/RV_Pathogen/Pixu_Analyses/countdata_proc.csv", header=T, row.names=1)
metadata <- read.csv("G:/My Drive/Research/SAS_Bots_Microbiome/RV_Pathogen/Pixu_Analyses/metadata_proc.csv", header=T, row.names=1)
var_infect <- c("inf_rv_yn", "inf_multi_hi", "inf_multi_mc", "inf_multi_sa", "inf_multi_sp")
pathogen_name <- c("Respiratory virus", "H. influenzae", "M. catarrhalis", "S. aureus", "S. pneumoniae")

propdata <- countdata/rowSums(countdata)
tab_interest <- cbind(propdata[,paste0("ASV",c(1,3,5,7,8,9,15))], metadata[,c("inf_sp", var_infect)])
interest_pairs <- rbind(c("ASV7", var_infect[2], pathogen_name[2]),
                        c("ASV8", var_infect[2], pathogen_name[2]),
                        c("ASV1", var_infect[3], pathogen_name[3]),               
                        c("ASV3", var_infect[4], pathogen_name[4]),
                        c("ASV5", var_infect[5], pathogen_name[5]),
                        c("ASV15", var_infect[5], pathogen_name[5]))

wilcox.test(ASV7 ~ inf_multi_hi, data=tab_interest)
hi <- expression(paste(italic("H. influenzae"), " colonization"))
asv7_title <- expression(paste(italic("Haemophilus"), " sp. (ASV7)"))
asv7 <- ggplot(tab_interest, aes(x=inf_multi_hi, y=ASV7)) + geom_violin(scale = "width", width = 0.7) + xlab(hi) + ylab("Relative abundance") +
  geom_boxplot(width=0.05, fill="indianred1", outlier.size=0.7) + ggtitle(asv7_title) + scale_y_sqrt() +
  scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + 
  theme_bw() + theme(axis.title.x = element_text(size=6), axis.title.y = element_text(size=6),
                     plot.title = element_text(size=6, hjust = 0.5, vjust = -2.0), ggplot2::margin(0,0,0,0), plot.margin = unit(c(0.2, 0, 0, 0), "cm"), legend.position="none",
                     axis.text.y = element_text(size=5.5, color="black"), axis.text.x = element_text(size=5.5, color="black")) +
  annotate("text", x=1.5, y=1.00, size=2, label= "p<0.0001")

wilcox.test(ASV8 ~ inf_multi_hi, data=tab_interest)
hi <- expression(paste(italic("H. influenzae"), " colonization"))
asv8_title <- expression(paste(italic("Haemophilus"), " sp. (ASV8)"))
asv8 <- ggplot(tab_interest, aes(x=inf_multi_hi, y=ASV8)) + geom_violin(scale = "width", width = 0.7) + xlab(hi) + ylab("") +
  geom_boxplot(width=0.05, fill="indianred1", outlier.size=0.7) + ggtitle(asv8_title) + scale_y_sqrt() +
  scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + 
  theme_bw() + theme(axis.title.x = element_text(size=6), axis.title.y = element_text(size=6),
                     plot.title = element_text(size=6, hjust = 0.5, vjust = -2.0), ggplot2::margin(0,0,0,0), plot.margin = unit(c(0.2, 0, 0, 0), "cm"), legend.position="none",
                     axis.text.y = element_text(size=5.5, color="black"), axis.text.x = element_text(size=5.5, color="black")) +
  annotate("text", x=1.5, y=1.00, size=2, label= "p<0.0001")

wilcox.test(ASV1 ~ inf_multi_mc, data=tab_interest)
mc <- expression(paste(italic("M. catarrhalis"), " colonization"))
asv1_title <- expression(paste(italic("M. catarrhalis/nonliquefaciens"), " (ASV1)"))
asv1 <- ggplot(tab_interest, aes(x=inf_multi_mc, y=ASV1)) + geom_violin(scale = "width", width = 0.7) + xlab(mc) + ylab("") +
  geom_boxplot(width=0.05, fill="indianred1", outlier.size=0.7) + ggtitle(asv1_title) + scale_y_sqrt() +
  scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + 
  theme_bw() + theme(axis.title.x = element_text(size=6), axis.title.y = element_text(size=6),
                     plot.title = element_text(size=6, hjust = 0.5, vjust = -2.0), ggplot2::margin(0,0,0,0), plot.margin = unit(c(0.2, 0, 0, 0), "cm"), legend.position="none",
                     axis.text.y = element_text(size=5.5, color="black"), axis.text.x = element_text(size=5.5, color="black")) +
  annotate("text", x=1.5, y=1.00, size=2, label= "p<0.0001")

wilcox.test(ASV3 ~ inf_multi_sa, data=tab_interest)
sa <- expression(paste(italic("S. aureus"), " colonization"))
asv3_title <- expression(paste(italic("Staphylococcus"), " sp. (ASV3)"))
asv3 <- ggplot(tab_interest, aes(x=inf_multi_sa, y=ASV3)) + geom_violin(scale = "width", width = 0.7) + xlab(sa) + ylab("Relative abundance") +
  geom_boxplot(width=0.05, fill="indianred1", outlier.size=0.7) + ggtitle(asv3_title) + scale_y_sqrt() +
  scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + 
  theme_bw() + theme(axis.title.x = element_text(size=6), axis.title.y = element_text(size=6),
                     plot.title = element_text(size=6, hjust = 0.5, vjust = -2.0), ggplot2::margin(0,0,0,0), plot.margin = unit(c(0.2, 0, 0, 0), "cm"), legend.position="none",
                     axis.text.y = element_text(size=5.5, color="black"), axis.text.x = element_text(size=5.5, color="black")) +
  annotate("text", x=1.5, y=1.00, size=2, label= "p<0.0001") 

wilcox.test(ASV5 ~ inf_multi_sp, data=tab_interest)
sp <- expression(paste(italic("S. pneumoniae"), " colonization"))
asv5_title <- expression(paste(italic("S. pneumoniae/pseudopneumoniae"), " (ASV5)"))
asv5 <- ggplot(tab_interest, aes(x=inf_multi_sp, y=ASV5)) + geom_violin(scale = "width", width = 0.7) + xlab(sp) + ylab("") +
  geom_boxplot(width=0.05, fill="indianred1", outlier.size=0.7) + ggtitle(asv5_title) + scale_y_sqrt() +
  scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + 
  theme_bw() + theme(axis.title.x = element_text(size=6), axis.title.y = element_text(size=6),
                     plot.title = element_text(size=6, hjust = 0.5, vjust = -2.0), ggplot2::margin(0,0,0,0), plot.margin = unit(c(0.2, 0, 0, 0), "cm"), legend.position="none",
                     axis.text.y = element_text(size=5.5, color="black"), axis.text.x = element_text(size=5.5, color="black")) +
  annotate("text", x=1.5, y=1.00, size=2, label= "p<0.0001")

cor(tab_interest$ASV5, tab_interest$inf_sp, method="spearman")
cor.test(tab_interest$ASV5, tab_interest$inf_sp, method="spearman", exact=FALSE)
sp_qpcr <- expression(paste(italic("S. pneumoniae"), " qPCR (log"[10], " copies/mL)"))
tab_interest$inf_sp[tab_interest$inf_sp>10] <- 10
asv5_scatter <- ggplot(tab_interest, aes(x=inf_sp, y=ASV5)) +
    geom_point(alpha=0.3, stroke=NA, size=1.5) +
    scale_y_sqrt() + scale_x_continuous(limits = c(0, 10), breaks = 0:10) + 
    labs(x=sp_qpcr, y="") + ggtitle(asv5_title) + 
  theme_bw() + theme(axis.title.x = element_text(size=6), axis.title.y = element_text(size=6),
                     plot.title = element_text(size=6, hjust = 0.5, vjust = -2.0), ggplot2::margin(0,0,0,0), plot.margin = unit(c(0.2, 0, 0, 0), "cm"), legend.position="none",
                     axis.text.y = element_text(size=5.5, color="black"), axis.text.x = element_text(size=5.5, color="black")) +
  annotate("text", x=2.95, y=0.98, size=2, label= expression(rho)) +
  annotate("text", x=5, y=1.00, size=2, label= "= 0.83, p<0.0001")

wilcox.test(ASV15 ~ inf_multi_sp, data=tab_interest)
sp <- expression(paste(italic("S. pneumoniae"), " colonization"))
asv15_title <- expression(paste(italic("S. thermophilus/vestibularis/salivarius"), " (ASV15)"))
asv15 <- ggplot(tab_interest, aes(x=inf_multi_sp, y=ASV15)) + geom_violin(scale = "width", width = 0.7) + xlab(sp) + ylab("") +
  geom_boxplot(width=0.05, fill="indianred1", outlier.size=0.7) + ggtitle(asv15_title) + scale_y_sqrt() +
  scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + 
  theme_bw() + theme(axis.title.x = element_text(size=6), axis.title.y = element_text(size=6),
                     plot.title = element_text(size=6, hjust = 0.5, vjust = -2.0), ggplot2::margin(0,0,0,0), plot.margin = unit(c(0.2, 0, 0, 0), "cm"), legend.position="none",
                     axis.text.y = element_text(size=5.5, color="black"), axis.text.x = element_text(size=5.5, color="black")) +
  annotate("text", x=1.5, y=1.00, size=2, label= "p<0.0001")

cor(tab_interest$ASV15, tab_interest$inf_sp, method="spearman")
cor.test(tab_interest$ASV15, tab_interest$inf_sp, method="spearman", exact=FALSE)
tab_interest$inf_sp[tab_interest$inf_sp>10] <- 10
asv15_scatter <- ggplot(tab_interest, aes(x=inf_sp, y=ASV15)) +
  geom_point(alpha=0.3, stroke=NA, size=1.5) +
  scale_y_sqrt() + scale_x_continuous(limits = c(0, 10), breaks = 0:10) + 
  labs(x=sp_qpcr, y="") + ggtitle(asv15_title) + 
  theme_bw() + theme(axis.title.x = element_text(size=6), axis.title.y = element_text(size=6),
                     plot.title = element_text(size=6, hjust = 0.5, vjust = -2.0), ggplot2::margin(0,0,0,0), plot.margin = unit(c(0.2, 0, 0, 0), "cm"), legend.position="none",
                     axis.text.y = element_text(size=5.5, color="black"), axis.text.x = element_text(size=5.5, color="black")) +
  annotate("text", x=2.85, y=0.98, size=2, label= expression(rho)) +
  annotate("text", x=5, y=1.00, size=2, label= "= -0.22, p<0.0001")

fig_S3 <- plot_grid(asv7, asv8, asv5, asv5_scatter, asv1, asv3, asv15, asv15_scatter, labels=c("a","","","b","","","",""), label_size=7.5, nrow=2) 

png(file="G:/My Drive/Research/SAS_Bots_Microbiome/RV_Pathogen/R_Plots/Figure_S3.png", width = 8.5, height = 4, units = 'in', res = 1200)
plot(fig_S3)
dev.off()