# Botswana Infant Microbiome Study - RV-Bacterial Analyses
# Matthew Kelly, MD, MPH 
# Supplementary Figure 1
# Last updated: March 9, 2025

remove(list=ls())
setwd("_________________") 
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
metadata_inf_np_16s$kmeans_cluster <- factor(metadata_inf_np_16s$kmeans_cluster, levels=c("OTH","STA","COR","STR","CD","HAE","CDM","MOR"))
sample_data(phy.inf.np.16s) <- metadata_inf_np_16s
phy.agglom <- tax_glom(phy.inf.np.16s, taxrank = 'Genus')
ntaxa(phy.agglom)
phy.relative <- transform_sample_counts(phy.agglom, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.relative))

# Compare classification of samples using the two clustering algorithms

table(metadata_inf_np_16s$cluster, metadata_inf_np_16s$kmeans_cluster)

# Fig. S1a and S1d - create NMDS plots on Bray-Curtis distances

ord <- ordinate(phy.relative, method="NMDS", distance="bray")

ord_theme   <-  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                      axis.text.y = element_text(size=8, color="black"), 
                      axis.title.y = element_text(size=9, color="black"),
                      axis.text.x  = element_text(size=8, color="black"),
                      axis.title.x = element_text(size=9, color="black"), plot.title = element_text(size=9.5, hjust=0.5),
                      legend.position = "right", legend.text = element_text(size=8, color="black"), 
                      legend.title = element_blank(), legend.key = element_blank(), plot.margin = unit(c(0.5, -0.2, 0, 0), "cm"),
                      legend.box.spacing = unit(0, "pt"), legend.key.size = unit(0.4, 'cm'))

fig_S1a <- plot_ordination(phy.relative, ord, color="cluster") +
  geom_point(shape=16, size=1.5) + ord_theme + scale_color_manual(values=cluster_cols_8) +
  stat_ellipse(geom="polygon", alpha=0, type="t", level=0.9, linewidth=0.5) + ggtitle("K-medoids clustering")

png(file="R_Plots/Figure_S1/Figure_S1a.png", width = 3.375, height = 3, units = 'in', res = 1200)
plot(fig_S1a)
dev.off()

fig_S1d <- plot_ordination(phy.relative, ord, color="kmeans_cluster") +
  geom_point(shape=16, size=1.5) + ord_theme + scale_color_manual(values=cluster_cols_8) +
  stat_ellipse(geom="polygon", alpha=0, type="t", level=0.9, linewidth=0.5) + ggtitle("K-means clustering") 

png(file="R_Plots/Figure_S1/Figure_S1d.png", width = 3.375, height = 3, units = 'in', res = 1200)
plot(fig_S1d)
dev.off()

# Fig. S1b & S1e - cluster barplots

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
relative_inf$Genus <- as.factor(relative_inf$Genus)
table(relative_inf$Genus)
relative_inf$Genus <- factor(relative_inf$Genus, levels=c("Acinetobacter", "Corynebacterium", "Dolosigranulum", "Gemella", 
                                                          "Haemophilus", "Lactobacillus", "Micrococcus", "Moraxella", "Neisseria", "Prevotella", 
                                                          "Pseudomonas", "Staphylococcus", "Streptococcus", "Veillonella", "Other"))
relative_inf <- subset(relative_inf, Abundance!=0)
relative_inf <- relative_inf[,c("Sample","Genus","cluster","kmeans_cluster","Abundance")]

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

fig_S1b <- ggplot(arrange(relative_inf, Genus), aes(x=cluster, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="fill") + theme_barplot + 
  scale_fill_manual(values=earthy_cols_15, 
                    labels=c(expression(italic("Acinetobacter")), expression(italic("Corynebacterium")), expression(italic("Dolosigranulum")), expression(italic("Gemella")), 
                             expression(italic("Haemophilus")), expression(italic("Lactobacillus")), expression(italic("Micrococcus")), expression(italic("Moraxella")), 
                             expression(italic("Neisseria")), expression(italic("Prevotella")), expression(italic("Pseudomonas")), expression(italic("Staphylococcus")), 
                             expression(italic("Streptococcus")), expression(italic("Veillonella")), "Other")) + 
  xlab("K-medoids microbiota type") + ylab("Relative abundance") 

png(file="R_Plots/Figure_S1/Figure_S1b.png", width = 4.125, height = 2.8, units = 'in', res = 1200)
plot(fig_S1b)
dev.off()

fig_S1e <- ggplot(arrange(relative_inf, Genus), aes(x=kmeans_cluster, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="fill") + theme_barplot + 
  scale_fill_manual(values=earthy_cols_15, 
                    labels=c(expression(italic("Acinetobacter")), expression(italic("Corynebacterium")), expression(italic("Dolosigranulum")), expression(italic("Gemella")), 
                             expression(italic("Haemophilus")), expression(italic("Lactobacillus")), expression(italic("Micrococcus")), expression(italic("Moraxella")), 
                             expression(italic("Neisseria")), expression(italic("Prevotella")), expression(italic("Pseudomonas")), expression(italic("Staphylococcus")), 
                             expression(italic("Streptococcus")), expression(italic("Veillonella")), "Other")) + 
  xlab("K-means microbiota type") + ylab("Relative abundance") 

png(file="R_Plots/Figure_S1/Figure_S1e.png", width = 4.125, height = 2.8, units = 'in', res = 1200)
plot(fig_S1e)
dev.off()

# Fig. S1c and S1f - overlap of k-medoids and k-means clustering algorithms

library(patchwork)

theme_barplot_clusters <-   theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                         axis.title.x = element_blank(), axis.text.x = element_blank(),
                         axis.text.y = element_text(size=8, color="black"), 
                         axis.title.y = element_text(angle=90, size=9, color="black"), 
                         legend.position = "left", legend.text=element_text(size=8, color="black"), legend.text.align = 0,
                         legend.box.spacing = unit(0, "pt"),
                         legend.title=element_blank(), legend.key.size = unit(0.35, 'cm')) 

fig_S1c <- ggplot(relative_inf, aes(x=cluster, y=1, fill=kmeans_cluster)) +
  geom_bar(stat="identity", position="fill") + theme_barplot_clusters +
  scale_fill_manual(values=cluster_cols_8) + ylab("K-medoids microbiota type") 

fig_S1c_x <- ggplot(relative_inf, aes(x=cluster, y=1, fill=cluster)) +
  geom_col(width = 1, position = "identity", color = "#00000020") +
  scale_fill_manual(values = cluster_cols_8, guide = "none") +
  coord_cartesian(expand = FALSE) +
  theme_minimal() + xlab("K-means microbiota type") +
  theme(axis.text.x = element_text(hjust = 0.5, size = 8, color="black"), axis.title.x = element_text(size=9, color="black"),
        axis.text.y = element_blank(), axis.title.y = element_blank())

fig_S1c <- fig_S1c + fig_S1c_x + plot_layout(ncol = 1, nrow = 2, heights = c(20, 1))

png(file="R_Plots/Figure_S1/Figure_S1c.png", width = 4.0, height = 3.2, units = 'in', res = 1200)
plot(fig_S1c)
dev.off()

fig_S1f <- ggplot(relative_inf, aes(x=kmeans_cluster, y=1, fill=cluster)) +
  geom_bar(stat="identity", position="fill") + theme_barplot_clusters +
  scale_fill_manual(values=cluster_cols_8) + ylab("K-means microbiota type") 

fig_S1f_x <- ggplot(relative_inf, aes(x=kmeans_cluster, y=1, fill=cluster)) +
  geom_col(width = 1, position = "identity", color = "#00000020") +
  scale_fill_manual(values = cluster_cols_8, guide = "none") +
  coord_cartesian(expand = FALSE) +
  theme_minimal() + xlab("K-medoids microbiota type") +
  theme(axis.text.x = element_text(hjust = 0.5, size = 8, color="black"), axis.title.x = element_text(size=9, color="black"),
        axis.text.y = element_blank(), axis.title.y = element_blank())

fig_S1f <- fig_S1f + fig_S1f_x + plot_layout(ncol = 1, nrow = 2, heights = c(20, 1))

png(file="R_Plots/Figure_S1/Figure_S1f.png", width = 4.0, height = 3.2, units = 'in', res = 1200)
plot(fig_S1f)
dev.off()

fig_S1abc <- plot_grid(fig_S1a, NULL, fig_S1b, fig_S1c, nrow=1, rel_widths=c(0.45, 0.03, 0.55, 0.45), 
                        labels=c("a","","b","c"), label_size=10)
fig_S1def <- plot_grid(fig_S1d, NULL, fig_S1e, fig_S1f, nrow=1, rel_widths=c(0.45, 0.03, 0.55, 0.45), 
                        labels=c("d","","e","f"), label_size=10)

png(file="R_Plots/Figure_S1.png", width = 13, height = 7, units = 'in', res = 1200)
plot_grid(fig_S1abc, fig_S1def, nrow=2)
dev.off()

# Save files as a Source Data file
source_data <- list('FigS1bcef'=relative_inf) 
openxlsx::write.xlsx(source_data, file="Source_Data/Figure_S1.xlsx")