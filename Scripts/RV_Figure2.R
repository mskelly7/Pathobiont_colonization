# Botswana Infant Microbiome Study - RV-Bacterial Analyses
# Matthew Kelly, MD, MPH 
# Figure 2
# Last updated: October 18, 2024

remove(list=ls())
setwd("G:/My Drive/Research/SAS_Bots_Microbiome/RV_Pathogen") 
set.seed(1234)

version
library(phyloseq)
packageVersion("phyloseq")
library(tidyverse)
library(dplyr)
library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(cowplot)
library(magick)
library(ggtext)
library(EBImage)
theme_set(theme_light()) 

# Create color palettes
library(RColorBrewer)
virus_cols <- c("navy", "white", "green4")
virus_cols <- colorRampPalette(virus_cols)
virus_cols <- virus_cols(7)
barplot(1:7, col=virus_cols)
virus_cols

bacteria_cols <- c("darkslateblue", "firebrick", "forestgreen", "mediumorchid4")
bacteria_cols_rev <- c("mediumorchid4", "forestgreen", "firebrick", "darkslateblue")

# Load required datasets
metadata_inf_np_RV <- read.csv("metadata_inf_np_RV.csv")
metadata_virus <- metadata_inf_np_RV %>% filter(!is.na(inf_rv))
metadata_bact <- metadata_inf_np_RV %>% filter(!is.na(inf_multi_hi))

# Results from survival analyses for pathobiont acquisition

# Respiratory virus infection

df <- tibble(pathogen = c("H. influenzae", "M. catarrhalis", "S. aureus", "S. pneumoniae"),
             HR = c(1.29, 1.39, 0.45, 1.55), 
             lower = c(0.98, 1.10, 0.28, 1.23),
             upper = c(1.72, 1.77, 0.73, 1.96))

df$pathogen <- factor(df$pathogen, levels=c( "S. pneumoniae", "S. aureus", "M. catarrhalis", "H. influenzae"))
levels(df$pathogen)

forest_rv <- ggplot(df, aes(x=pathogen, y=HR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=bacteria_cols_rev, guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=bacteria_cols_rev, guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="Hazard ratio", limits = c(0, 2.5)) +
  coord_flip() + xlab("") + ggtitle("Coincident respiratory virus infection") +
  theme(legend.text=element_text(size=6, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.1, 'cm'), 
        legend.box.spacing = unit(0.18, "pt"), axis.title.x = element_text(size=6.5),  
        plot.title = element_text(size=6.5, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size=5.5, color="black")) 

# H. influenzae colonization

df <- tibble(pathogen = c("M. catarrhalis", "S. aureus", "S. pneumoniae"),
             HR = c(1.28, 1.25, 1.29), 
             lower = c(0.97, 0.76, 0.98),
             upper = c(1.70, 2.08, 1.71))

df$pathogen <- factor(df$pathogen, levels=c( "S. pneumoniae", "S. aureus", "M. catarrhalis"))
levels(df$pathogen)

title_hi <- expression(paste("Preceding ", italic("H. influenzae"), " colonization"))
forest_hi <- ggplot(df, aes(x=pathogen, y=HR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=c("mediumorchid4", "forestgreen", "firebrick"), guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=c("mediumorchid4", "forestgreen", "firebrick"), guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="  ", limits = c(0, 2.5)) +
  coord_flip() + xlab("") + ggtitle(title_hi) +
  theme(legend.text=element_text(size=6, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.1, 'cm'), 
        legend.box.spacing = unit(0.18, "pt"), axis.title.x = element_text(size=6.5),  
        plot.title = element_text(size=6.5, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size=5.5, color="black")) 

# M. catarrhalis colonization

df <- tibble(pathogen = c("H. influenzae", "S. aureus", "S. pneumoniae"),
             HR = c(1.55, 0.80, 1.37), 
             lower = c(1.11, 0.52, 1.07),
             upper = c(2.16, 1.23, 1.76))

df$pathogen <- factor(df$pathogen, levels=c("S. pneumoniae", "S. aureus", "H. influenzae"))
levels(df$pathogen)

title_mc <- expression(paste("Preceding ", italic("M. catarrhalis"), " colonization"))
forest_mc <- ggplot(df, aes(x=pathogen, y=HR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=c("mediumorchid4", "forestgreen", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=c("mediumorchid4", "forestgreen", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="  ", limits = c(0, 2.5)) +
  coord_flip() + xlab("") + ggtitle(title_mc) +
  theme(legend.text=element_text(size=6, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.1, 'cm'), 
        legend.box.spacing = unit(0.18, "pt"), axis.title.x = element_text(size=6.5),  
        plot.title = element_text(size=6.5, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size=5.5, color="black")) 

# S. aureus colonization

df <- tibble(pathogen = c("H. influenzae", "M. catarrhalis", "S. pneumoniae"),
             HR = c(0.93, 0.84, 0.93), 
             lower = c(0.64, 0.61, 0.68),
             upper = c(1.36, 1.15, 1.27))

df$pathogen <- factor(df$pathogen, levels=c("S. pneumoniae", "M. catarrhalis", "H. influenzae"))
levels(df$pathogen)

title_sa <- expression(paste("Preceding ", italic("S. aureus"), " colonization"))
forest_sa <- ggplot(df, aes(x=pathogen, y=HR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=c("mediumorchid4", "firebrick", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=c("mediumorchid4", "firebrick", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="Hazard ratio", limits = c(0, 2.5)) +
  coord_flip() + xlab("") + ggtitle(title_sa) +
  theme(legend.text=element_text(size=6, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.1, 'cm'), 
        legend.box.spacing = unit(0.18, "pt"), axis.title.x = element_text(size=6.5),  
        plot.title = element_text(size=6.5, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size=5.5, color="black")) 

# S. pneumoniae colonization

df <- tibble(pathogen = c("H. influenzae", "M. catarrhalis", "S. aureus"),
             HR = c(1.29, 1.32, 1.06), 
             lower = c(0.96, 1.01, 0.64),
             upper = c(1.74, 1.73, 1.77))

df$pathogen <- factor(df$pathogen, levels=c( "S. aureus", "M. catarrhalis", "H. influenzae"))
levels(df$pathogen)

title_sp <- expression(paste("Preceding ", italic("S. pneumoniae"), " colonization"))
forest_sp <- ggplot(df, aes(x=pathogen, y=HR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=c("forestgreen", "firebrick", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=c("forestgreen", "firebrick", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="Hazard ratio", limits = c(0, 2.5)) +
  coord_flip() + xlab("") + ggtitle(title_sp) +
  theme(legend.text=element_text(size=6, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.1, 'cm'), 
        legend.box.spacing = unit(0.18, "pt"), axis.title.x = element_text(size=6.5),  
        plot.title = element_text(size=6.5, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size=5.5, color="black")) 

fig_2a <- plot_grid(forest_rv, NULL, forest_hi, NULL, forest_mc, NULL, forest_sa, NULL, forest_sp, 
                    rel_heights=c(0.63, -0.07, 0.5, -0.07, 0.5, -0.07, 0.5, -0.07, 0.5), labels=c("a","","","","","","","",""), label_size=7,  nrow=9)

png(file="R_Plots/Figure_2/Figure_2a.png", width = 2.25, height = 2.75, units = 'in', res = 600)
print(fig_2a)
dev.off()

# S. pneumoniae colonization density by respiratory virus infection

library(ggpubr)
library(rstatix)

ylab_sp <- expression(paste(italic("S. pneumoniae"), " log"[10], " copies/mL"))
metadata_inf_spy <- subset(metadata_virus, inf_sp_yn=="Y" & inf_multi_sp=="Y")
tapply(metadata_inf_spy$inf_sp, metadata_inf_spy$inf_rv_yn, summary)
# For graphical purposes, set samples with Sp PCR > 10 log copies/mL to maximum of 10 log copies/mL
metadata_inf_spy$inf_sp[metadata_inf_spy$inf_sp>10] <- 10
metadata_inf_spy$inf_rv_yn <- revalue(metadata_inf_spy$inf_rv_yn, c("N"="No", "Y"="Yes"))
compare_means(inf_sp ~ inf_rv_yn, data = metadata_inf_spy)
fig_2b_rv <- ggplot(metadata_inf_spy, aes(x=inf_rv_yn, y=inf_sp, fill=inf_rv_yn)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.15, width=0.3) +
  scale_fill_manual(values=c("mediumorchid4", "mediumorchid4")) + scale_y_continuous(breaks=seq(2, 10.2, 2), limits=c(2, 10.2)) +
  ylab(ylab_sp) + xlab("Respiratory virus infection") + ggtitle("  ") +
  theme(axis.title.y = element_text(size=6.75), axis.title.x = element_text(size=6.5), 
        plot.title = element_text(size=7, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0.3), "cm"), legend.position="none",
        axis.text.y = element_text(size=6, color="black"), axis.text.x = element_text(size=6, color="black")) + 
  annotate("text", x=1.5, y=10.2, size=1.8, label= "p=0.006") 

title_hi <- expression(paste(italic("H. influenzae"), " colonization"))
metadata_inf_spy$inf_multi_hi <- revalue(metadata_inf_spy$inf_multi_hi, c("N"="No", "Y"="Yes"))
compare_means(inf_sp ~ inf_multi_hi, data = metadata_inf_spy)
fig_2b_hi <- ggplot(metadata_inf_spy, aes(x=inf_multi_hi, y=inf_sp, fill=inf_multi_hi)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.15, width=0.3) +
  scale_fill_manual(values=c("mediumorchid4", "mediumorchid4")) + scale_y_continuous(breaks=seq(2, 10.2, 2), limits=c(2, 10.2)) +
  ylab(ylab_sp) + xlab(title_hi) + ggtitle("  ") +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=6.5), 
        plot.title = element_text(size=7, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0.3), "cm"), legend.position="none",
        axis.text.y = element_text(size=6, color="black"), axis.text.x = element_text(size=6, color="black")) + 
  annotate("text", x=1.5, y=10.2, size=1.8, label= "p<0.0001") 

fig_2b <- plot_grid(fig_2b_rv, fig_2b_hi, rel_widths=c(2,1.81), nrow=1)

png(file="R_Plots/Figure_2/Figure_2b.png", width = 3, height = 2.25, units = 'in', res = 600)
print(fig_2b)
dev.off()

# Venn diagram of bacterial pathobionts

library(ggVennDiagram)
venn <- metadata_virus[,c("SampleID", "inf_multi_sa","inf_multi_sp","inf_multi_mc","inf_multi_hi")]
sa <- subset(venn, inf_multi_sa=="Y")
sa <- sa$SampleID
sp <- subset(venn, inf_multi_sp=="Y")
sp <- sp$SampleID
mc <- subset(venn, inf_multi_mc=="Y")
mc <- mc$SampleID
hi <- subset(venn, inf_multi_hi=="Y")
hi <- hi$SampleID
venn <- list(hi, mc, sp, sa)
remove(hi, mc, sa, sp)

venn_plot <- ggVennDiagram(venn, label_alpha = 0, edge_size=0.4, set_size=2.3, label_size=2.1, 
                           category.names = c("", "", "", "")) +
  scale_fill_gradient(low="white", high="indianred3") + scale_color_grey(start = 0.0, end = 0.0) + theme(legend.position = "none") +
  scale_x_continuous(expand = expansion(mult = .1)) 
venn_plot <- venn_plot + annotate("text", x=0.12, y=0.77, size=2.2, label="H. influenzae", fontface="italic") 
venn_plot <- venn_plot + annotate("text", x=0.30, y=0.86, size=2.2, label="M. catarrhalis", fontface="italic") 
venn_plot <- venn_plot + annotate("text", x=0.62, y=0.85, size=2.2, label="S. pneumoniae", fontface="italic") 
venn_plot <- venn_plot + annotate("text", x=0.84, y=0.78, size=2.2, label="S. aureus", fontface="italic") 

png(file="R_Plots/Figure_2/Figure_2c.png", width = 4, height = 2.75, units = 'in', res = 300)
plot(venn_plot)
dev.off()

venn_plot <- venn_plot + theme(plot.margin = unit(c(0.2, -0.5, -0.25, 0.3), "cm"))
venn_plot <- plot_grid(NULL, venn_plot, NULL, rel_widths = c(-0.10, 0.6, -0.06), nrow=1)
fig_2b <- plot_grid(NULL, fig_2b, NULL, rel_widths = c(0.05, 0.85,0.10), labels=NULL, nrow=1) 
fig_2bc <- plot_grid(fig_2b, venn_plot, rel_heights = c(0.45,0.55), labels=c("b","c"), label_size=7, ncol=1) 

png(file="R_Plots/Figure_2.png", width = 7, height = 4.75, units = 'in', res = 1200)
plot_grid(fig_2a, fig_2bc, labels=NULL, ncol=2, rel_widths=c(1.8,2)) 
dev.off()