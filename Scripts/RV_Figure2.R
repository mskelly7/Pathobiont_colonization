# Botswana Infant Microbiome Study - RV-Bacterial Analyses
# Matthew Kelly, MD, MPH 
# Figure 2
# Last updated: October 21, 2024

remove(list=ls())
setwd("___________________") 
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
             OR = c(1.44, 1.34, 0.45, 1.83), 
             lower = c(1.06, 0.96, 0.29, 1.33),
             upper = c(1.96, 1.87, 0.70, 2.52))

df$pathogen <- factor(df$pathogen, levels=c( "S. pneumoniae", "S. aureus", "M. catarrhalis", "H. influenzae"))
levels(df$pathogen)

forest_rv <- ggplot(df, aes(x=pathogen, y=OR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=bacteria_cols_rev, guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=bacteria_cols_rev, guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="Odds ratio", limits = c(0, 2.75)) +
  coord_flip() + xlab("") + ggtitle("Coincident respiratory virus infection") +
  theme(legend.text=element_text(size=6, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.1, 'cm'), 
        legend.box.spacing = unit(0.18, "pt"), axis.title.x = element_text(size=6.5),  
        plot.title = element_text(size=6.5, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size=5.5, color="black")) 

# H. influenzae colonization

df <- tibble(pathogen = c("M. catarrhalis", "S. aureus", "S. pneumoniae"),
             OR = c(1.29, 0.92, 1.42), 
             lower = c(0.80, 0.55, 0.90),
             upper = c(2.08, 1.51, 2.24))

df$pathogen <- factor(df$pathogen, levels=c( "S. pneumoniae", "S. aureus", "M. catarrhalis"))
levels(df$pathogen)

title_hi <- expression(paste("Preceding ", italic("H. influenzae"), " colonization"))
forest_hi <- ggplot(df, aes(x=pathogen, y=OR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=c("mediumorchid4", "forestgreen", "firebrick"), guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=c("mediumorchid4", "forestgreen", "firebrick"), guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="  ", limits = c(0, 2.75)) +
  coord_flip() + xlab("") + ggtitle(title_hi) +
  theme(legend.text=element_text(size=6, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.1, 'cm'), 
        legend.box.spacing = unit(0.18, "pt"), axis.title.x = element_text(size=6.5),  
        plot.title = element_text(size=6.5, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size=5.5, color="black")) 

# M. catarrhalis colonization

df <- tibble(pathogen = c("H. influenzae", "S. aureus", "S. pneumoniae"),
             OR = c(1.82, 0.75, 1.58), 
             lower = c(1.28, 0.48, 1.12),
             upper = c(2.60, 1.18, 2.24))

df$pathogen <- factor(df$pathogen, levels=c("S. pneumoniae", "S. aureus", "H. influenzae"))
levels(df$pathogen)

title_mc <- expression(paste("Preceding ", italic("M. catarrhalis"), " colonization"))
forest_mc <- ggplot(df, aes(x=pathogen, y=OR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=c("mediumorchid4", "forestgreen", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=c("mediumorchid4", "forestgreen", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="  ", limits = c(0, 2.75)) +
  coord_flip() + xlab("") + ggtitle(title_mc) +
  theme(legend.text=element_text(size=6, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.1, 'cm'), 
        legend.box.spacing = unit(0.18, "pt"), axis.title.x = element_text(size=6.5),  
        plot.title = element_text(size=6.5, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size=5.5, color="black")) 

# S. aureus colonization

df <- tibble(pathogen = c("H. influenzae", "M. catarrhalis", "S. pneumoniae"),
             OR = c(0.94, 1.27, 1.02), 
             lower = c(0.64, 0.89, 0.71),
             upper = c(1.40, 1.80, 1.48))

df$pathogen <- factor(df$pathogen, levels=c("S. pneumoniae", "M. catarrhalis", "H. influenzae"))
levels(df$pathogen)

title_sa <- expression(paste("Preceding ", italic("S. aureus"), " colonization"))
forest_sa <- ggplot(df, aes(x=pathogen, y=OR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=c("mediumorchid4", "firebrick", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=c("mediumorchid4", "firebrick", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="Odds ratio", limits = c(0, 2.75)) +
  coord_flip() + xlab("") + ggtitle(title_sa) +
  theme(legend.text=element_text(size=6, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.1, 'cm'), 
        legend.box.spacing = unit(0.18, "pt"), axis.title.x = element_text(size=6.5),  
        plot.title = element_text(size=6.5, hjust = 0.5, margin = ggplot2::margin(0,0,0,0)), plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size=5.5, color="black")) 

# S. pneumoniae colonization

df <- tibble(pathogen = c("H. influenzae", "M. catarrhalis", "S. aureus"),
             OR = c(1.42, 1.07, 0.66), 
             lower = c(1.01, 0.71, 0.41),
             upper = c(1.98, 1.60, 1.06))

df$pathogen <- factor(df$pathogen, levels=c( "S. aureus", "M. catarrhalis", "H. influenzae"))
levels(df$pathogen)

title_sp <- expression(paste("Preceding ", italic("S. pneumoniae"), " colonization"))
forest_sp <- ggplot(df, aes(x=pathogen, y=OR, ymin=lower, ymax=upper, col=pathogen, fill=pathogen)) + 
  geom_linerange(linewidth=3,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) + geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=c("forestgreen", "firebrick", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=c("forestgreen", "firebrick", "darkslateblue"), guide = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(name="Odds ratio", limits = c(0, 2.75)) +
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