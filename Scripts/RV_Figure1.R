# Botswana Infant Microbiome Study - RV-Bacterial Analyses
# Matthew Kelly, MD, MPH 
# Figure 1
# Last updated: March 9, 2025

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
library(wesanderson)
budapest1 <- wes_palette("GrandBudapest1")
budapest2 <- wes_palette("GrandBudapest2")

library(RColorBrewer)
virus_cols <- c("midnightblue", "white", "darkgreen")
virus_cols <- colorRampPalette(virus_cols)
virus_cols <- virus_cols(8)
barplot(1:8, col=virus_cols)
virus_cols

bacteria_cols <- c("darkslateblue", "firebrick", "forestgreen", "mediumorchid4")
bacteria_cols_rev <- c("mediumorchid4", "forestgreen", "firebrick", "darkslateblue")

# Load required datasets
metadata_inf_np_RV <- read.csv("metadata_inf_np_RV.csv")
phy.inf.np.16s <- read_rds("phy.inf.np.16s.rds")
metadata_inf_np_16s <- data.frame(sample_data(phy.inf.np.16s))

metadata_virus <- metadata_inf_np_RV %>% filter(!is.na(inf_rv))
metadata_bact <- metadata_inf_np_RV %>% filter(!is.na(inf_multi_hi))

# Overview of sample microbiological testing

infant_virus <- length(unique(metadata_virus$study_id))
infant_bact <- length(unique(metadata_bact$study_id))
infant_16s <- length(unique(metadata_inf_np_16s$study_id))
infant_n <- data.frame(variable = c("Respiratory viruses", "Bacterial pathobionts", "URT microbiota"), value = c(infant_virus, infant_bact, infant_16s), 
                       x = c("Infants", "Infants", "Infants"))
infant_n$variable <- factor(infant_n$variable, levels=c("Respiratory viruses", "Bacterial pathobionts", "URT microbiota"))

dot_infant <- ggplot(infant_n, aes(x = x, y = variable)) + 
  geom_point(aes(size = value, fill = variable), shape = 21) + 
  geom_text(aes(label = value), parse = TRUE, size=3.5, color="black") +
  scale_size_continuous(limits = c(0, 300), range = c(1,11)) + 
  labs( x= "", y = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 8), 
        axis.text.y = element_text(colour = "black", size = 8),
        panel.background = element_blank(), panel.border = element_rect(colour = "gray70", fill = NA, linewidth=0.8), 
        legend.position = "none", plot.margin=unit(c(0,0,0,-0.3), "cm")) +  
  scale_fill_manual(values = budapest2, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(infant_n$variable))) +
  scale_x_discrete(position="top")

total_virus <- nrow(metadata_virus)
total_bact <- nrow(metadata_bact)
total_16s <- nrow(metadata_inf_np_16s)
total_n <- data.frame(variable = c("Respiratory viruses", "Bacterial pathobionts", "URT microbiota"), value = c(total_virus, total_bact, total_16s), 
                      x = c("Samples", "Samples", "Samples"))
total_n$variable <- factor(total_n$variable, levels=c("Respiratory viruses", "Bacterial pathobionts", "URT microbiota"))

dot_total <- ggplot(total_n, aes(x = x, y = variable)) + 
  geom_point(aes(size = value, fill = variable), shape = 21) + 
  geom_text(aes(label = value), parse = TRUE, size=3.5, color="black") +
  scale_size_continuous(limits = c(0, 2409), range = c(1,13)) + 
  labs( x= "", y = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 8), 
        axis.text.y = element_text(colour = "black", size = 8),
        panel.background = element_blank(), panel.border = element_rect(colour = "gray70", fill = NA, size=0.8), 
        legend.position = "none", plot.margin=unit(c(0,0,0,-0.5), "cm")) +  
  scale_fill_manual(values = budapest2, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(total_n$variable))) +
  scale_x_discrete(position="top")

month_virus <- data.frame(unclass(table(metadata_virus$month)))
month_bact <- data.frame(unclass(table(metadata_bact$month)))
month_16s <- data.frame(unclass(table(metadata_inf_np_16s$month)))
month_n <- data.frame(t(cbind(month_virus, month_bact, month_16s)))
colnames(month_n) <- gsub("X","m",colnames(month_n))
month_n <- data.frame(t(month_n))
month_n$month <- row.names(month_n)
colnames(month_n) <- c("Respiratory viruses", "Bacterial pathobionts", "URT microbiota", "month")
month_n <- melt(month_n, id = c("month"))
month_n$variable <- factor(month_n$variable, levels=c("Respiratory viruses", "Bacterial pathobionts", "URT microbiota"))
month_n$month <- factor(month_n$month, levels=c("m0","m1","m2","m3","m4","m5","m6","m8","m10","m12"))

dot_month <- ggplot(month_n, aes(x = month, y = variable)) + 
  geom_point(aes(size = value, fill = variable), shape = 21) + 
  geom_text(aes(label = value), parse = TRUE, size=3, color="black") +
  scale_size_continuous(limits = c(0, 300), range = c(1,10)) + 
  labs( x= "", y = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 8), 
        axis.text.y = element_text(colour = "black", size = 8), 
        panel.background = element_blank(), panel.border = element_rect(colour = "gray70", fill = NA, linewidth=0.8), 
        legend.position = "none", plot.margin=unit(c(0,0,0,-0.5), "cm")) +  
  scale_fill_manual(values = budapest2, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(month_n$variable))) +
  scale_x_discrete(position="top")

dot_total <- dot_total + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
dot_month <- dot_month + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

dot_final <- cowplot::plot_grid(dot_infant, dot_total, dot_month, nrow=1, rel_widths = c(0.36, 0.14, 0.80))

infant <- readImage("infant.png")
infant <- ggdraw() + draw_image(infant, scale=1)

fig_1a <- plot_grid(infant, dot_final, rel_widths = c(0.2,0.8), labels="a", label_size=9, nrow=1) 

png(file="R_Plots/Figure_1/Figure_1a.png", width = 7.5, height = 2.0, units = 'in', res = 300)
plot(fig_1a)
dev.off()

# Barplot of detection of respiratory viruses by month

virus_n <- data.frame(t(unclass(table(metadata_virus$month, metadata_virus$inf_rv_cat))))
colnames(virus_n) <- gsub("X","m",colnames(virus_n))
# Convert counts of respiratory viral detections to % detection
virus_num <- as.vector(colSums(virus_n))
virus_n <- as.matrix(virus_n)
library(sweep)
virus_n <- sweep(virus_n, 2, virus_num, '/')
virus_n <- virus_n[!(row.names(virus_n)=="No viruses"),]
virus_n <- data.frame(t(virus_n))
virus_n$month <- row.names(virus_n)
virus_n$month <- factor(virus_n$month, levels=c("m0","m1","m2","m3","m4","m5","m6","m8","m10","m12"))
virus_n <- melt(virus_n, id = c("month"))
virus_n$variable <- as.character(virus_n$variable)
virus_n$variable[virus_n$variable=="X.1.virus"] <- ">1 virus"
virus_n$variable[virus_n$variable=="Rhinovirus.enterovirus"] <- "Rhinovirus/enterovirus"
virus_n$variable[virus_n$variable=="SARS2"] <- "SARS-CoV-2"
virus_n$variable <- factor(virus_n$variable, levels=c("Adenovirus", "HMPV", "Influenza", "Parainfluenza", "Rhinovirus/enterovirus", "RSV", 
                                                      "SARS-CoV-2", ">1 virus"))

virus_plot <- ggplot(virus_n, aes(x=month, y=value, fill=variable)) + 
  geom_col(color = "black", linewidth = 0.25) + xlab("") + ylab("Proportion of samples") + 
  scale_x_discrete(limits=c("m0","m1","m2","m3","m4","m5","m6","m8","m10","m12")) + scale_y_continuous(breaks = seq(0, 0.4, by = 0.1)) + 
  scale_fill_manual(values=virus_cols, labels = c("Adenovirus", "Human metapneumovirus", "Influenza viruses A/B", 
                                                  "Parainfluenza viruses", "Rhinovirus/enterovirus", "Respiratory syncytial virus", "SARS-CoV-2", ">1 virus")) + 
  theme(legend.text=element_text(size=7), legend.title=element_blank(), legend.key.size = unit(0.3, 'cm'), legend.box.spacing = unit(0, "pt"),
        axis.title.y = element_text(angle=90, size=8, margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.y = element_text(size=7.5, color = "black"), axis.text.x = element_text(size=7.5, color = "black"), 
        plot.margin = unit(c(0.5, 0, 0, 0), "cm")) 

png(file="R_Plots/Figure_1/Figure_1b.png", width = 4, height = 2.75, units = 'in', res = 300)
plot(virus_plot)
dev.off()

# Prevalence of bacterial pathogen colonization by age

metadata_bact$inf_multi_sp_pos[metadata_bact$inf_multi_sp=="N"] <- 0
metadata_bact$inf_multi_sp_pos[metadata_bact$inf_multi_sp=="Y"] <- 1
metadata_bact$inf_multi_sa_pos[metadata_bact$inf_multi_sa=="N"] <- 0
metadata_bact$inf_multi_sa_pos[metadata_bact$inf_multi_sa=="Y"] <- 1
metadata_bact$inf_multi_hi_pos[metadata_bact$inf_multi_hi=="N"] <- 0
metadata_bact$inf_multi_hi_pos[metadata_bact$inf_multi_hi=="Y"] <- 1
metadata_bact$inf_multi_mc_pos[metadata_bact$inf_multi_mc=="N"] <- 0
metadata_bact$inf_multi_mc_pos[metadata_bact$inf_multi_mc=="Y"] <- 1
pathogens <- metadata_bact %>% group_by(month) %>% 
  summarise_at(vars(inf_multi_sp_pos, inf_multi_hi_pos, inf_multi_mc_pos, inf_multi_sa_pos), list(name = mean))
pathogens <- melt(pathogens, id.vars = "month")
pathogens$month <- as.factor(pathogens$month)
pathogens$pathogen[pathogens$variable=="inf_multi_sp_pos_name"] <- "S. pneumoniae"
pathogens$pathogen[pathogens$variable=="inf_multi_hi_pos_name"] <- "H. influenzae"
pathogens$pathogen[pathogens$variable=="inf_multi_mc_pos_name"] <- "M. catarrhalis"
pathogens$pathogen[pathogens$variable=="inf_multi_sa_pos_name"] <- "S. aureus"
pathogens$month2[pathogens$month=="0"] <- "m0"
pathogens$month2[pathogens$month=="1"] <- "m1"
pathogens$month2[pathogens$month=="2"] <- "m2"
pathogens$month2[pathogens$month=="3"] <- "m3"
pathogens$month2[pathogens$month=="4"] <- "m4"
pathogens$month2[pathogens$month=="5"] <- "m5"
pathogens$month2[pathogens$month=="6"] <- "m6"
pathogens$month2[pathogens$month=="8"] <- "m8"
pathogens$month2[pathogens$month=="10"] <- "m10"
pathogens$month2[pathogens$month=="12"] <- "m12"

pathogen_plot <- ggplot(pathogens, aes(x = month2, y = value, fill = pathogen, group=pathogen)) +  
  geom_line(aes(color=pathogen), linewidth=1.5) + geom_point(size = 2, shape = 21) + 
  scale_x_discrete(limits=c("m0","m1","m2","m3","m4","m5","m6","m8","m10","m12")) +
  scale_fill_manual(values=bacteria_cols) + ylim(0,1) + 
  scale_color_manual(values=bacteria_cols) + ylab("Colonization prevalence") + xlab("") +
  theme(legend.text=element_text(size=7, face="italic"), legend.title=element_blank(), legend.key.size = unit(0.27, 'cm'), 
        legend.box.spacing = unit(0, "pt"),
        axis.title.y = element_text(size=8, margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.y = element_text(size=7.5, color = "black"), axis.text.x = element_text(size=7.5, color = "black"), 
        plot.margin = unit(c(0.5, 0, 0, 0), "cm")) 

png(file="R_Plots/Figure_1/Figure_1c.png", width = 3.75, height = 2.75, units = 'in', res = 300)
plot(pathogen_plot)
dev.off()

fig_1bc <- plot_grid(virus_plot, pathogen_plot, rel_widths = c(0.53,0.47), labels=c("b","c"), label_size=9, nrow=1)  

png(file="R_Plots/Figure_1.png", width = 7.5, height = 4.75, units = 'in', res = 1200)
plot_grid(fig_1a, NULL, fig_1bc, labels=NULL, nrow=3, rel_heights=c(2, 0.08, 2.75)) 
dev.off()

# Save files as a Source Data file
source_data <- list('Fig1a_infant'=infant_n, 'Fig1a_samples'=total_n, 'Fig1a_months'=month_n, 
                    'Fig1b'=virus_n, 'Fig1c'=pathogens)
openxlsx::write.xlsx(source_data, file="Source_Data/Figure_1.xlsx")