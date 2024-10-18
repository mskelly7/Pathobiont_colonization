# Botswana Infant Microbiome Study - RV-Bacterial Pathogen Analyses
# Matthew Kelly, MD, MPH 
# Dataset Preprocessing
# Last update: October 18, 2024

remove(list=ls())
setwd("__________________") 
set.seed(1234)

version
library(tidyverse)
library(DataCombine)
library(data.table)
library(phyloseq)
packageVersion("phyloseq")
library(vegan)
library(cluster)
library(clusterSim)
library(ade4)
library(decontam)
packageVersion("decontam")

##################################
# EXAMINE NEGATIVE CONTROL SAMPLES
##################################

phy.bots.nps <- readRDS("phy.bots.nps.rds")
sample_data <- data.frame(sample_data(phy.bots.nps))
total_reads <- sum(sample_sums(phy.bots.nps))
phy.controls <- subset_samples(phy.bots.nps, control=="Y")
phy.controls <- filter_taxa(phy.controls, function(x) sum(x) > 0, TRUE)
sample_sums(phy.controls)
ntaxa(phy.controls)
control_asvs <- data.frame(t(otu_table(phy.controls)))
control_asvs$ASV <- rownames((control_asvs))
taxtable <- data.frame(tax_table(phy.bots.nps))
taxtable$ASV <- rownames((taxtable))
control_asvs <- merge(taxtable, control_asvs, by="ASV", all.x=FALSE)
rownames(control_asvs) <- control_asvs[,1]
control_asvs <- control_asvs[,-1]
phy.controls.agglom <- tax_glom(phy.controls, taxrank = 'Genus')
ntaxa(phy.controls.agglom)
phy.controls.relative <- transform_sample_counts(phy.controls.agglom, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.controls.relative))  
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling
relative_controls <- psmelt(phy.controls.relative)
relative_controls$Genus <- as.character(relative_controls$Genus)
control_genera <- aggregate(relative_controls$Abundance, by=list(Phylum=relative_controls$Phylum, Genus=relative_controls$Genus,
                                                                 OTU=relative_controls$OTU), FUN=mean)
control_genera <- aggregate(control_genera$x, by=list(Phylum=control_genera$Phylum,
                                                      Genus=control_genera$Genus), FUN=sum)
colnames(control_genera)[3] <- "genus_Ab"
sum(control_genera$genus_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(control_genera)            # Corresponds to # of unique genera
control_genera <- control_genera[order(-control_genera$genus_Ab),]
head(control_genera)
remove(phy.controls, phy.controls.agglom, phy.controls.relative, relative_controls)

#########################################
# REMOVAL OF CONTAMINANT ASVS IN DECONTAM
#########################################

# Identification of contaminant ASVs using the FREQUENCY METHOD
# Exclude samples from which library prep DNA concentration is unavailable
phy.dna <- subset_samples(phy.bots.nps, !is.na(dna))
nsamples(phy.dna)
# Exclude 0 month samples from decontam as these have vastly different microbiota and, as expected, lower DNA concentrations 
tapply(sample_data$dna, sample_data$month, summary)
phy.dna <- subset_samples(phy.dna, month!="0")
nsamples(phy.dna)
sample_data_dna <- data.frame(sample_data(phy.dna))
table(sample_data_dna$control)
sample_data_dna$counts <- sample_sums(phy.dna)
sample_data_dna <- sample_data_dna[order(sample_data_dna$counts),]
sample_data_dna$index <- seq(nrow(sample_data_dna))
ggplot(data=sample_data_dna, aes(x=index, y=counts, color=control)) + geom_point()

# Identify contaminants based on frequency using library preparation concentrations
suppressWarnings(freq <- isContaminant(phy.dna, method="frequency", conc="dna", threshold=0.1))
table(freq$contaminant)
freq_prev50 <- freq
freq_prev50$contaminant[freq_prev50$prev<50] <- FALSE
plot_frequency(phy.dna, taxa_names(phy.dna)[sample(which(freq_prev50$contaminant))], conc="dna") + xlab("DNA Concentration")

# Review contaminant plot for highly prevalent ASVs to confirm that these represent contaminants
contam_freq <- subset(freq, contaminant==TRUE)
contam_freq <- merge(taxtable, contam_freq, by="row.names")
names(contam_freq)[names(contam_freq)=="Row.names"] <- "Taxon"
# Remove the following ASVs from the list of contaminants: 
# ASV14 (Neisseriaceae), ASV26 (M. luteus), ASV38 (Veillonella atypica), ASV107 (Kocuria sp.), 
# ASV67 (Moraxella osloensis), ASV113 (Haemophilus parainfluenzae), ASV69 (Streptococcus sp.)
table(freq$contaminant)
freq$contaminant[row.names(freq)=="ASV14" | row.names(freq)=="ASV26" | row.names(freq)=="ASV38" |
                          row.names(freq)=="ASV107" | row.names(freq)=="ASV67" | row.names(freq)=="ASV113" |
                          row.names(freq)=="ASV69"] <- FALSE
table(freq$contaminant)
contaminants <- subset(freq, contaminant==TRUE)
contaminants$ASV <- rownames((contaminants))
decontam_asvs_freq <- merge(taxtable, contaminants, by="ASV", all.x=FALSE)
rownames(decontam_asvs_freq) <- decontam_asvs_freq[,1]
ntaxa(phy.dna)
phy.decontam.freq <- prune_taxa(!freq$contaminant, phy.bots.nps)
phy.decontam.freq <- prune_samples(sample_sums(phy.decontam.freq)>0, phy.decontam.freq)
nsamples(phy.decontam.freq)
ntaxa(phy.decontam.freq)
remove(sample_data_dna, phy.dna, contam_freq, contaminants, freq, freq_prev50)

# Identification of contaminant ASVs using the PREVALENCE METHOD
sample_data(phy.decontam.freq)$is.neg <- sample_data(phy.decontam.freq)$Sample_or_Control == "Control Sample"
contam_prev <- isContaminant(phy.decontam.freq, method="prevalence", neg="is.neg", threshold=0.10)
table(contam_prev$contaminant)
ps.pa <- transform_sample_counts(phy.decontam.freq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contam_prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
contaminants <- subset(contam_prev, contaminant==TRUE)
contaminants$ASV <- rownames((contaminants))
decontam_asvs_prev <- merge(taxtable, contaminants, by="ASV", all.x=FALSE)
rownames(decontam_asvs_prev) <- decontam_asvs_prev[,1]
nsamples(phy.decontam.freq)
ntaxa(phy.decontam.freq)
phy.decontam.prev <- prune_taxa(!contam_prev$contaminant, phy.decontam.freq)
phy.decontam.prev <- prune_samples(sample_sums(phy.decontam.prev)>0, phy.decontam.prev)
nsamples(phy.decontam.prev)
ntaxa(phy.decontam.prev)
decontam_reads <- sum(sample_sums(phy.decontam.prev))
decontam_reads/total_reads
saveRDS(phy.decontam.prev, "phy.bots.nps.decontam.rds")
remove(contaminants, contam_prev, df.pa, ps.pa, ps.pa.neg, ps.pa.pos)

##################################################################################
# BATCH CONTAMINANT REMOVAL (see Moosavia et al. Microbiome. 2021 Feb 10;9(1):41.)
##################################################################################

library(RColorBrewer) 
library(vegan) 
library(zCompositions)
library(irr)
library(plyr)
library(data.table)

phy.decontam.prev <- readRDS("phy.bots.nps.decontam.rds")
decontam_reads <- sum(sample_sums(phy.decontam.prev))
ntaxa(phy.decontam.prev)

# Remove control samples
nsamples(phy.decontam.prev)
phy.decontam.prev <- subset_samples(phy.decontam.prev, control=="N")
nsamples(phy.decontam.prev)

# Calculate proportion of microbiome variation explained by batch prior to batch filtering
#phy.ordination <- prune_taxa(names(sort(taxa_sums(phy.decontam.prev), TRUE)[1:50]), phy.decontam.prev)
#phy.ordination <- prune_samples(sample_sums(phy.ordination)>0, phy.ordination)
#bray_pcoa <- ordinate(phy.ordination, "PCoA", "bray")
#plot_ordination(phy.ordination, bray_pcoa, "samples", color="batch")

#D_BC <- phyloseq::distance(phy.ordination, "bray")
#D_JC <- phyloseq::distance(phy.ordination, "jaccard")
#dist_list = list('Bray Curtis' = D_BC, 'Jaccard' = D_JC)
#for (i in 1:length(dist_list)) {res_batch = vegan::adonis2(dist_list[[i]] ~ phyloseq::sample_data(phy.ordination)$batch)
#  cat(names(dist_list[i]))
#  print(res_batch)
#  cat("\n", "--------------------------------------------------------------------------", "\n")}

# Let x = the expected ('reference') prevalence/proportion of an ASV. If we think the prevalence of the ASV is this, given the sample sizes of 
# both batches (N1 and N2), what is the lowest proportion "expected" in the other batch?
Structure_comp <- function(x) {(x - sqrt(x*(1-x)/N1) - sqrt(x*(1-x)/N2))*k}
# Let k = a constant stringency factor provided by the user (value between 0 and 1, the lower the value, the more relaxed the test is) 
# How similar is the prevalence expected to be between batches?
k = 0.05
# If set to 1, expect prevalences to be exactly the same; Moosavia et al. used k=1/15 (0.067)

phy.decontam.batch <- transform_sample_counts(phy.decontam.prev, function(abund) 1*(abund>0))
phy.set1 <- prune_samples(sample_data(phy.decontam.batch)$batch =="nps_set1", phy.decontam.batch) 
phy.set1
N1 <- 329
phy.set2 <- prune_samples(sample_data(phy.decontam.batch)$batch =="nps_set2", phy.decontam.batch) 
phy.set2
N2 <- 327
phy.set3 <- prune_samples(sample_data(phy.decontam.batch)$batch =="nps_set3", phy.decontam.batch) 
phy.set3
N3 <- 383
phy.set4 <- prune_samples(sample_data(phy.decontam.batch)$batch =="nps_set4", phy.decontam.batch) 
phy.set4
N4 <- 357
phy.set5 <- prune_samples(sample_data(phy.decontam.batch)$batch =="nps_set5", phy.decontam.batch) 
phy.set5
N5 <- 222
phy.set6 <- prune_samples(sample_data(phy.decontam.batch)$batch =="nps_set6", phy.decontam.batch) 
phy.set6
N6 <- 381
phy.set7 <- prune_samples(sample_data(phy.decontam.batch)$batch =="nps_set7", phy.decontam.batch) 
phy.set7
N7 <- 376
phy.set8 <- prune_samples(sample_data(phy.decontam.batch)$batch =="nps_set8", phy.decontam.batch) 
phy.set8
N8 <- 278

# Set 1 batch comparisons

phy.set1_vs_set2 <- data.frame(phy.set1=taxa_sums(phy.set1)/N1, phy.set2=taxa_sums(phy.set2)/N2, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set1_vs_set2$phy.set1>0 | phy.set1_vs_set2$phy.set2>0)
phy.set1_vs_set2$formula_set1 <- Structure_comp(phy.set1_vs_set2$phy.set1) # Given prevalence in set1
phy.set1_vs_set2$formula_set2 <- Structure_comp(phy.set1_vs_set2$phy.set2) # Given prevalence in set2
batch_asvs_set1C_vs_set2 <- subset(phy.set1_vs_set2, phy.set2 < formula_set1) # Potential contaminants in set2 (shaded area on x-axis)
batch_asvs_set1C_vs_set2 <- batch_asvs_set1C_vs_set2[,-c(1,2,10,11)]
batch_asvs_set1C_vs_set2 <- setDT(batch_asvs_set1C_vs_set2, keep.rownames = "ASV")
batch_asvs_set1_vs_set2C <- subset(phy.set1_vs_set2, phy.set1 < formula_set2) # Potential contaminants in set2 (shaded area on x-axis)
batch_asvs_set1_vs_set2C <- batch_asvs_set1_vs_set2C[,-c(1,2,10,11)]
batch_asvs_set1_vs_set2C <- setDT(batch_asvs_set1_vs_set2C, keep.rownames = "ASV")

phy.set1_vs_set3 <- data.frame(phy.set1=taxa_sums(phy.set1)/N1, phy.set3=taxa_sums(phy.set3)/N3, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set1_vs_set3$phy.set1>0 | phy.set1_vs_set3$phy.set3>0)
phy.set1_vs_set3$formula_set1 <- Structure_comp(phy.set1_vs_set3$phy.set1) # Given prevalence in set1
phy.set1_vs_set3$formula_set3 <- Structure_comp(phy.set1_vs_set3$phy.set3) # Given prevalence in set3
batch_asvs_set1C_vs_set3 <- subset(phy.set1_vs_set3, phy.set3 < formula_set1) # Potential contaminants in set3 (shaded area on x-axis)
batch_asvs_set1C_vs_set3 <- batch_asvs_set1C_vs_set3[,-c(1,2,10,11)]
batch_asvs_set1C_vs_set3 <- setDT(batch_asvs_set1C_vs_set3, keep.rownames = "ASV")
batch_asvs_set1_vs_set3C <- subset(phy.set1_vs_set3, phy.set1 < formula_set3) # Potential contaminants in set3 (shaded area on x-axis)
batch_asvs_set1_vs_set3C <- batch_asvs_set1_vs_set3C[,-c(1,2,10,11)]
batch_asvs_set1_vs_set3C <- setDT(batch_asvs_set1_vs_set3C, keep.rownames = "ASV")

phy.set1_vs_set3 <- data.frame(phy.set1=taxa_sums(phy.set1)/N1, phy.set3=taxa_sums(phy.set3)/N3, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set1_vs_set3$phy.set1>0 | phy.set1_vs_set3$phy.set3>0)
phy.set1_vs_set3$formula_set1 <- Structure_comp(phy.set1_vs_set3$phy.set1) # Given prevalence in set1
phy.set1_vs_set3$formula_set3 <- Structure_comp(phy.set1_vs_set3$phy.set3) # Given prevalence in set3
batch_asvs_set1C_vs_set3 <- subset(phy.set1_vs_set3, phy.set3 < formula_set1) # Potential contaminants in set3 (shaded area on x-axis)
batch_asvs_set1C_vs_set3 <- batch_asvs_set1C_vs_set3[,-c(1,2,10,11)]
batch_asvs_set1C_vs_set3 <- setDT(batch_asvs_set1C_vs_set3, keep.rownames = "ASV")
batch_asvs_set1_vs_set3C <- subset(phy.set1_vs_set3, phy.set1 < formula_set3) # Potential contaminants in set3 (shaded area on x-axis)
batch_asvs_set1_vs_set3C <- batch_asvs_set1_vs_set3C[,-c(1,2,10,11)]
batch_asvs_set1_vs_set3C <- setDT(batch_asvs_set1_vs_set3C, keep.rownames = "ASV")

phy.set1_vs_set4 <- data.frame(phy.set1=taxa_sums(phy.set1)/N1, phy.set4=taxa_sums(phy.set4)/N4, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set1_vs_set4$phy.set1>0 | phy.set1_vs_set4$phy.set4>0)
phy.set1_vs_set4$formula_set1 <- Structure_comp(phy.set1_vs_set4$phy.set1) # Given prevalence in set1
phy.set1_vs_set4$formula_set4 <- Structure_comp(phy.set1_vs_set4$phy.set4) # Given prevalence in set4
batch_asvs_set1C_vs_set4 <- subset(phy.set1_vs_set4, phy.set4 < formula_set1) # Potential contaminants in set4 (shaded area on x-axis)
batch_asvs_set1C_vs_set4 <- batch_asvs_set1C_vs_set4[,-c(1,2,10,11)]
batch_asvs_set1C_vs_set4 <- setDT(batch_asvs_set1C_vs_set4, keep.rownames = "ASV")
batch_asvs_set1_vs_set4C <- subset(phy.set1_vs_set4, phy.set1 < formula_set4) # Potential contaminants in set4 (shaded area on x-axis)
batch_asvs_set1_vs_set4C <- batch_asvs_set1_vs_set4C[,-c(1,2,10,11)]
batch_asvs_set1_vs_set4C <- setDT(batch_asvs_set1_vs_set4C, keep.rownames = "ASV")

phy.set1_vs_set5 <- data.frame(phy.set1=taxa_sums(phy.set1)/N1, phy.set5=taxa_sums(phy.set5)/N5, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set1_vs_set5$phy.set1>0 | phy.set1_vs_set5$phy.set5>0)
phy.set1_vs_set5$formula_set1 <- Structure_comp(phy.set1_vs_set5$phy.set1) # Given prevalence in set1
phy.set1_vs_set5$formula_set5 <- Structure_comp(phy.set1_vs_set5$phy.set5) # Given prevalence in set5
batch_asvs_set1C_vs_set5 <- subset(phy.set1_vs_set5, phy.set5 < formula_set1) # Potential contaminants in set5 (shaded area on x-axis)
batch_asvs_set1C_vs_set5 <- batch_asvs_set1C_vs_set5[,-c(1,2,10,11)]
batch_asvs_set1C_vs_set5 <- setDT(batch_asvs_set1C_vs_set5, keep.rownames = "ASV")
batch_asvs_set1_vs_set5C <- subset(phy.set1_vs_set5, phy.set1 < formula_set5) # Potential contaminants in set5 (shaded area on x-axis)
batch_asvs_set1_vs_set5C <- batch_asvs_set1_vs_set5C[,-c(1,2,10,11)]
batch_asvs_set1_vs_set5C <- setDT(batch_asvs_set1_vs_set5C, keep.rownames = "ASV")

phy.set1_vs_set6 <- data.frame(phy.set1=taxa_sums(phy.set1)/N1, phy.set6=taxa_sums(phy.set6)/N6, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set1_vs_set6$phy.set1>0 | phy.set1_vs_set6$phy.set6>0)
phy.set1_vs_set6$formula_set1 <- Structure_comp(phy.set1_vs_set6$phy.set1) # Given prevalence in set1
phy.set1_vs_set6$formula_set6 <- Structure_comp(phy.set1_vs_set6$phy.set6) # Given prevalence in set6
batch_asvs_set1C_vs_set6 <- subset(phy.set1_vs_set6, phy.set6 < formula_set1) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set1C_vs_set6 <- batch_asvs_set1C_vs_set6[,-c(1,2,10,11)]
batch_asvs_set1C_vs_set6 <- setDT(batch_asvs_set1C_vs_set6, keep.rownames = "ASV")
batch_asvs_set1_vs_set6C <- subset(phy.set1_vs_set6, phy.set1 < formula_set6) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set1_vs_set6C <- batch_asvs_set1_vs_set6C[,-c(1,2,10,11)]
batch_asvs_set1_vs_set6C <- setDT(batch_asvs_set1_vs_set6C, keep.rownames = "ASV")

phy.set1_vs_set7 <- data.frame(phy.set1=taxa_sums(phy.set1)/N1, phy.set7=taxa_sums(phy.set7)/N7, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set1_vs_set7$phy.set1>0 | phy.set1_vs_set7$phy.set7>0)
phy.set1_vs_set7$formula_set1 <- Structure_comp(phy.set1_vs_set7$phy.set1) # Given prevalence in set1
phy.set1_vs_set7$formula_set7 <- Structure_comp(phy.set1_vs_set7$phy.set7) # Given prevalence in set7
batch_asvs_set1C_vs_set7 <- subset(phy.set1_vs_set7, phy.set7 < formula_set1) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set1C_vs_set7 <- batch_asvs_set1C_vs_set7[,-c(1,2,10,11)]
batch_asvs_set1C_vs_set7 <- setDT(batch_asvs_set1C_vs_set7, keep.rownames = "ASV")
batch_asvs_set1_vs_set7C <- subset(phy.set1_vs_set7, phy.set1 < formula_set7) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set1_vs_set7C <- batch_asvs_set1_vs_set7C[,-c(1,2,10,11)]
batch_asvs_set1_vs_set7C <- setDT(batch_asvs_set1_vs_set7C, keep.rownames = "ASV")

phy.set1_vs_set8 <- data.frame(phy.set1=taxa_sums(phy.set1)/N1, phy.set8=taxa_sums(phy.set8)/N8, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set1_vs_set8$phy.set1>0 | phy.set1_vs_set8$phy.set8>0)
phy.set1_vs_set8$formula_set1 <- Structure_comp(phy.set1_vs_set8$phy.set1) # Given prevalence in set1
phy.set1_vs_set8$formula_set8 <- Structure_comp(phy.set1_vs_set8$phy.set8) # Given prevalence in set8
batch_asvs_set1C_vs_set8 <- subset(phy.set1_vs_set8, phy.set8 < formula_set1) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set1C_vs_set8 <- batch_asvs_set1C_vs_set8[,-c(1,2,10,11)]
batch_asvs_set1C_vs_set8 <- setDT(batch_asvs_set1C_vs_set8, keep.rownames = "ASV")
batch_asvs_set1_vs_set8C <- subset(phy.set1_vs_set8, phy.set1 < formula_set8) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set1_vs_set8C <- batch_asvs_set1_vs_set8C[,-c(1,2,10,11)]
batch_asvs_set1_vs_set8C <- setDT(batch_asvs_set1_vs_set8C, keep.rownames = "ASV")

set1_asvs <- unique(rbind(batch_asvs_set1C_vs_set2, batch_asvs_set1_vs_set2C, batch_asvs_set1C_vs_set3, batch_asvs_set1_vs_set3C,
                          batch_asvs_set1C_vs_set4, batch_asvs_set1_vs_set4C, batch_asvs_set1C_vs_set5, batch_asvs_set1_vs_set5C,
                          batch_asvs_set1C_vs_set6, batch_asvs_set1_vs_set6C, batch_asvs_set1C_vs_set7, batch_asvs_set1_vs_set7C,
                          batch_asvs_set1C_vs_set8, batch_asvs_set1_vs_set8C))

# Set 2 batch contaminants

phy.set2_vs_set3 <- data.frame(phy.set2=taxa_sums(phy.set2)/N2, phy.set3=taxa_sums(phy.set3)/N3, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set2_vs_set3$phy.set2>0 | phy.set2_vs_set3$phy.set3>0)
phy.set2_vs_set3$formula_set2 <- Structure_comp(phy.set2_vs_set3$phy.set2) # Given prevalence in set2
phy.set2_vs_set3$formula_set3 <- Structure_comp(phy.set2_vs_set3$phy.set3) # Given prevalence in set3
batch_asvs_set2C_vs_set3 <- subset(phy.set2_vs_set3, phy.set3 < formula_set2) # Potential contaminants in set3 (shaded area on x-axis)
batch_asvs_set2C_vs_set3 <- batch_asvs_set2C_vs_set3[,-c(1,2,10,11)]
batch_asvs_set2C_vs_set3 <- setDT(batch_asvs_set2C_vs_set3, keep.rownames = "ASV")
batch_asvs_set2_vs_set3C <- subset(phy.set2_vs_set3, phy.set2 < formula_set3) # Potential contaminants in set3 (shaded area on x-axis)
batch_asvs_set2_vs_set3C <- batch_asvs_set2_vs_set3C[,-c(1,2,10,11)]
batch_asvs_set2_vs_set3C <- setDT(batch_asvs_set2_vs_set3C, keep.rownames = "ASV")

phy.set2_vs_set4 <- data.frame(phy.set2=taxa_sums(phy.set2)/N2, phy.set4=taxa_sums(phy.set4)/N4, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set2_vs_set4$phy.set2>0 | phy.set2_vs_set4$phy.set4>0)
phy.set2_vs_set4$formula_set2 <- Structure_comp(phy.set2_vs_set4$phy.set2) # Given prevalence in set2
phy.set2_vs_set4$formula_set4 <- Structure_comp(phy.set2_vs_set4$phy.set4) # Given prevalence in set4
batch_asvs_set2C_vs_set4 <- subset(phy.set2_vs_set4, phy.set4 < formula_set2) # Potential contaminants in set4 (shaded area on x-axis)
batch_asvs_set2C_vs_set4 <- batch_asvs_set2C_vs_set4[,-c(1,2,10,11)]
batch_asvs_set2C_vs_set4 <- setDT(batch_asvs_set2C_vs_set4, keep.rownames = "ASV")
batch_asvs_set2_vs_set4C <- subset(phy.set2_vs_set4, phy.set2 < formula_set4) # Potential contaminants in set4 (shaded area on x-axis)
batch_asvs_set2_vs_set4C <- batch_asvs_set2_vs_set4C[,-c(1,2,10,11)]
batch_asvs_set2_vs_set4C <- setDT(batch_asvs_set2_vs_set4C, keep.rownames = "ASV")

phy.set2_vs_set5 <- data.frame(phy.set2=taxa_sums(phy.set2)/N2, phy.set5=taxa_sums(phy.set5)/N5, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set2_vs_set5$phy.set2>0 | phy.set2_vs_set5$phy.set5>0)
phy.set2_vs_set5$formula_set2 <- Structure_comp(phy.set2_vs_set5$phy.set2) # Given prevalence in set2
phy.set2_vs_set5$formula_set5 <- Structure_comp(phy.set2_vs_set5$phy.set5) # Given prevalence in set5
batch_asvs_set2C_vs_set5 <- subset(phy.set2_vs_set5, phy.set5 < formula_set2) # Potential contaminants in set5 (shaded area on x-axis)
batch_asvs_set2C_vs_set5 <- batch_asvs_set2C_vs_set5[,-c(1,2,10,11)]
batch_asvs_set2C_vs_set5 <- setDT(batch_asvs_set2C_vs_set5, keep.rownames = "ASV")
batch_asvs_set2_vs_set5C <- subset(phy.set2_vs_set5, phy.set2 < formula_set5) # Potential contaminants in set5 (shaded area on x-axis)
batch_asvs_set2_vs_set5C <- batch_asvs_set2_vs_set5C[,-c(1,2,10,11)]
batch_asvs_set2_vs_set5C <- setDT(batch_asvs_set2_vs_set5C, keep.rownames = "ASV")

phy.set2_vs_set6 <- data.frame(phy.set2=taxa_sums(phy.set2)/N2, phy.set6=taxa_sums(phy.set6)/N6, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set2_vs_set6$phy.set2>0 | phy.set2_vs_set6$phy.set6>0)
phy.set2_vs_set6$formula_set2 <- Structure_comp(phy.set2_vs_set6$phy.set2) # Given prevalence in set2
phy.set2_vs_set6$formula_set6 <- Structure_comp(phy.set2_vs_set6$phy.set6) # Given prevalence in set6
batch_asvs_set2C_vs_set6 <- subset(phy.set2_vs_set6, phy.set6 < formula_set2) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set2C_vs_set6 <- batch_asvs_set2C_vs_set6[,-c(1,2,10,11)]
batch_asvs_set2C_vs_set6 <- setDT(batch_asvs_set2C_vs_set6, keep.rownames = "ASV")
batch_asvs_set2_vs_set6C <- subset(phy.set2_vs_set6, phy.set2 < formula_set6) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set2_vs_set6C <- batch_asvs_set2_vs_set6C[,-c(1,2,10,11)]
batch_asvs_set2_vs_set6C <- setDT(batch_asvs_set2_vs_set6C, keep.rownames = "ASV")

phy.set2_vs_set7 <- data.frame(phy.set2=taxa_sums(phy.set2)/N2, phy.set7=taxa_sums(phy.set7)/N7, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set2_vs_set7$phy.set2>0 | phy.set2_vs_set7$phy.set7>0)
phy.set2_vs_set7$formula_set2 <- Structure_comp(phy.set2_vs_set7$phy.set2) # Given prevalence in set2
phy.set2_vs_set7$formula_set7 <- Structure_comp(phy.set2_vs_set7$phy.set7) # Given prevalence in set7
batch_asvs_set2C_vs_set7 <- subset(phy.set2_vs_set7, phy.set7 < formula_set2) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set2C_vs_set7 <- batch_asvs_set2C_vs_set7[,-c(1,2,10,11)]
batch_asvs_set2C_vs_set7 <- setDT(batch_asvs_set2C_vs_set7, keep.rownames = "ASV")
batch_asvs_set2_vs_set7C <- subset(phy.set2_vs_set7, phy.set2 < formula_set7) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set2_vs_set7C <- batch_asvs_set2_vs_set7C[,-c(1,2,10,11)]
batch_asvs_set2_vs_set7C <- setDT(batch_asvs_set2_vs_set7C, keep.rownames = "ASV")

phy.set2_vs_set8 <- data.frame(phy.set2=taxa_sums(phy.set2)/N2, phy.set8=taxa_sums(phy.set8)/N8, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set2_vs_set8$phy.set2>0 | phy.set2_vs_set8$phy.set8>0)
phy.set2_vs_set8$formula_set2 <- Structure_comp(phy.set2_vs_set8$phy.set2) # Given prevalence in set2
phy.set2_vs_set8$formula_set8 <- Structure_comp(phy.set2_vs_set8$phy.set8) # Given prevalence in set8
batch_asvs_set2C_vs_set8 <- subset(phy.set2_vs_set8, phy.set8 < formula_set2) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set2C_vs_set8 <- batch_asvs_set2C_vs_set8[,-c(1,2,10,11)]
batch_asvs_set2C_vs_set8 <- setDT(batch_asvs_set2C_vs_set8, keep.rownames = "ASV")
batch_asvs_set2_vs_set8C <- subset(phy.set2_vs_set8, phy.set2 < formula_set8) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set2_vs_set8C <- batch_asvs_set2_vs_set8C[,-c(1,2,10,11)]
batch_asvs_set2_vs_set8C <- setDT(batch_asvs_set2_vs_set8C, keep.rownames = "ASV")

set2_asvs <- unique(rbind(batch_asvs_set2C_vs_set3, batch_asvs_set2_vs_set3C, batch_asvs_set2C_vs_set4, batch_asvs_set2_vs_set4C, 
                          batch_asvs_set2C_vs_set5, batch_asvs_set2_vs_set5C, batch_asvs_set2C_vs_set6, batch_asvs_set2_vs_set6C, 
                          batch_asvs_set2C_vs_set7, batch_asvs_set2_vs_set7C, batch_asvs_set2C_vs_set8, batch_asvs_set2_vs_set8C))

# Set 3 batch contaminants

phy.set3_vs_set4 <- data.frame(phy.set3=taxa_sums(phy.set3)/N3, phy.set4=taxa_sums(phy.set4)/N4, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set3_vs_set4$phy.set3>0 | phy.set3_vs_set4$phy.set4>0)
phy.set3_vs_set4$formula_set3 <- Structure_comp(phy.set3_vs_set4$phy.set3) # Given prevalence in set3
phy.set3_vs_set4$formula_set4 <- Structure_comp(phy.set3_vs_set4$phy.set4) # Given prevalence in set4
batch_asvs_set3C_vs_set4 <- subset(phy.set3_vs_set4, phy.set4 < formula_set3) # Potential contaminants in set4 (shaded area on x-axis)
batch_asvs_set3C_vs_set4 <- batch_asvs_set3C_vs_set4[,-c(1,2,10,11)]
batch_asvs_set3C_vs_set4 <- setDT(batch_asvs_set3C_vs_set4, keep.rownames = "ASV")
batch_asvs_set3_vs_set4C <- subset(phy.set3_vs_set4, phy.set3 < formula_set4) # Potential contaminants in set4 (shaded area on x-axis)
batch_asvs_set3_vs_set4C <- batch_asvs_set3_vs_set4C[,-c(1,2,10,11)]
batch_asvs_set3_vs_set4C <- setDT(batch_asvs_set3_vs_set4C, keep.rownames = "ASV")

phy.set3_vs_set5 <- data.frame(phy.set3=taxa_sums(phy.set3)/N3, phy.set5=taxa_sums(phy.set5)/N5, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set3_vs_set5$phy.set3>0 | phy.set3_vs_set5$phy.set5>0)
phy.set3_vs_set5$formula_set3 <- Structure_comp(phy.set3_vs_set5$phy.set3) # Given prevalence in set3
phy.set3_vs_set5$formula_set5 <- Structure_comp(phy.set3_vs_set5$phy.set5) # Given prevalence in set5
batch_asvs_set3C_vs_set5 <- subset(phy.set3_vs_set5, phy.set5 < formula_set3) # Potential contaminants in set5 (shaded area on x-axis)
batch_asvs_set3C_vs_set5 <- batch_asvs_set3C_vs_set5[,-c(1,2,10,11)]
batch_asvs_set3C_vs_set5 <- setDT(batch_asvs_set3C_vs_set5, keep.rownames = "ASV")
batch_asvs_set3_vs_set5C <- subset(phy.set3_vs_set5, phy.set3 < formula_set5) # Potential contaminants in set5 (shaded area on x-axis)
batch_asvs_set3_vs_set5C <- batch_asvs_set3_vs_set5C[,-c(1,2,10,11)]
batch_asvs_set3_vs_set5C <- setDT(batch_asvs_set3_vs_set5C, keep.rownames = "ASV")

phy.set3_vs_set6 <- data.frame(phy.set3=taxa_sums(phy.set3)/N3, phy.set6=taxa_sums(phy.set6)/N6, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set3_vs_set6$phy.set3>0 | phy.set3_vs_set6$phy.set6>0)
phy.set3_vs_set6$formula_set3 <- Structure_comp(phy.set3_vs_set6$phy.set3) # Given prevalence in set3
phy.set3_vs_set6$formula_set6 <- Structure_comp(phy.set3_vs_set6$phy.set6) # Given prevalence in set6
batch_asvs_set3C_vs_set6 <- subset(phy.set3_vs_set6, phy.set6 < formula_set3) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set3C_vs_set6 <- batch_asvs_set3C_vs_set6[,-c(1,2,10,11)]
batch_asvs_set3C_vs_set6 <- setDT(batch_asvs_set3C_vs_set6, keep.rownames = "ASV")
batch_asvs_set3_vs_set6C <- subset(phy.set3_vs_set6, phy.set3 < formula_set6) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set3_vs_set6C <- batch_asvs_set3_vs_set6C[,-c(1,2,10,11)]
batch_asvs_set3_vs_set6C <- setDT(batch_asvs_set3_vs_set6C, keep.rownames = "ASV")

phy.set3_vs_set7 <- data.frame(phy.set3=taxa_sums(phy.set3)/N3, phy.set7=taxa_sums(phy.set7)/N7, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set3_vs_set7$phy.set3>0 | phy.set3_vs_set7$phy.set7>0)
phy.set3_vs_set7$formula_set3 <- Structure_comp(phy.set3_vs_set7$phy.set3) # Given prevalence in set3
phy.set3_vs_set7$formula_set7 <- Structure_comp(phy.set3_vs_set7$phy.set7) # Given prevalence in set7
batch_asvs_set3C_vs_set7 <- subset(phy.set3_vs_set7, phy.set7 < formula_set3) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set3C_vs_set7 <- batch_asvs_set3C_vs_set7[,-c(1,2,10,11)]
batch_asvs_set3C_vs_set7 <- setDT(batch_asvs_set3C_vs_set7, keep.rownames = "ASV")
batch_asvs_set3_vs_set7C <- subset(phy.set3_vs_set7, phy.set3 < formula_set7) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set3_vs_set7C <- batch_asvs_set3_vs_set7C[,-c(1,2,10,11)]
batch_asvs_set3_vs_set7C <- setDT(batch_asvs_set3_vs_set7C, keep.rownames = "ASV")

phy.set3_vs_set8 <- data.frame(phy.set3=taxa_sums(phy.set3)/N3, phy.set8=taxa_sums(phy.set8)/N8, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set3_vs_set8$phy.set3>0 | phy.set3_vs_set8$phy.set8>0)
phy.set3_vs_set8$formula_set3 <- Structure_comp(phy.set3_vs_set8$phy.set3) # Given prevalence in set3
phy.set3_vs_set8$formula_set8 <- Structure_comp(phy.set3_vs_set8$phy.set8) # Given prevalence in set8
batch_asvs_set3C_vs_set8 <- subset(phy.set3_vs_set8, phy.set8 < formula_set3) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set3C_vs_set8 <- batch_asvs_set3C_vs_set8[,-c(1,2,10,11)]
batch_asvs_set3C_vs_set8 <- setDT(batch_asvs_set3C_vs_set8, keep.rownames = "ASV")
batch_asvs_set3_vs_set8C <- subset(phy.set3_vs_set8, phy.set3 < formula_set8) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set3_vs_set8C <- batch_asvs_set3_vs_set8C[,-c(1,2,10,11)]
batch_asvs_set3_vs_set8C <- setDT(batch_asvs_set3_vs_set8C, keep.rownames = "ASV")

set3_asvs <- unique(rbind(batch_asvs_set3C_vs_set4, batch_asvs_set3_vs_set4C, batch_asvs_set3C_vs_set5, batch_asvs_set3_vs_set5C, 
                          batch_asvs_set3C_vs_set6, batch_asvs_set3_vs_set6C, batch_asvs_set3C_vs_set7, batch_asvs_set3_vs_set7C, 
                          batch_asvs_set3C_vs_set8, batch_asvs_set3_vs_set8C))

# Set 4 batch contaminants

phy.set4_vs_set5 <- data.frame(phy.set4=taxa_sums(phy.set4)/N4, phy.set5=taxa_sums(phy.set5)/N5, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set4_vs_set5$phy.set4>0 | phy.set4_vs_set5$phy.set5>0)
phy.set4_vs_set5$formula_set4 <- Structure_comp(phy.set4_vs_set5$phy.set4) # Given prevalence in set4
phy.set4_vs_set5$formula_set5 <- Structure_comp(phy.set4_vs_set5$phy.set5) # Given prevalence in set5
batch_asvs_set4C_vs_set5 <- subset(phy.set4_vs_set5, phy.set5 < formula_set4) # Potential contaminants in set5 (shaded area on x-axis)
batch_asvs_set4C_vs_set5 <- batch_asvs_set4C_vs_set5[,-c(1,2,10,11)]
batch_asvs_set4C_vs_set5 <- setDT(batch_asvs_set4C_vs_set5, keep.rownames = "ASV")
batch_asvs_set4_vs_set5C <- subset(phy.set4_vs_set5, phy.set4 < formula_set5) # Potential contaminants in set5 (shaded area on x-axis)
batch_asvs_set4_vs_set5C <- batch_asvs_set4_vs_set5C[,-c(1,2,10,11)]
batch_asvs_set4_vs_set5C <- setDT(batch_asvs_set4_vs_set5C, keep.rownames = "ASV")

phy.set4_vs_set6 <- data.frame(phy.set4=taxa_sums(phy.set4)/N4, phy.set6=taxa_sums(phy.set6)/N6, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set4_vs_set6$phy.set4>0 | phy.set4_vs_set6$phy.set6>0)
phy.set4_vs_set6$formula_set4 <- Structure_comp(phy.set4_vs_set6$phy.set4) # Given prevalence in set4
phy.set4_vs_set6$formula_set6 <- Structure_comp(phy.set4_vs_set6$phy.set6) # Given prevalence in set6
batch_asvs_set4C_vs_set6 <- subset(phy.set4_vs_set6, phy.set6 < formula_set4) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set4C_vs_set6 <- batch_asvs_set4C_vs_set6[,-c(1,2,10,11)]
batch_asvs_set4C_vs_set6 <- setDT(batch_asvs_set4C_vs_set6, keep.rownames = "ASV")
batch_asvs_set4_vs_set6C <- subset(phy.set4_vs_set6, phy.set4 < formula_set6) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set4_vs_set6C <- batch_asvs_set4_vs_set6C[,-c(1,2,10,11)]
batch_asvs_set4_vs_set6C <- setDT(batch_asvs_set4_vs_set6C, keep.rownames = "ASV")

phy.set4_vs_set7 <- data.frame(phy.set4=taxa_sums(phy.set4)/N4, phy.set7=taxa_sums(phy.set7)/N7, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set4_vs_set7$phy.set4>0 | phy.set4_vs_set7$phy.set7>0)
phy.set4_vs_set7$formula_set4 <- Structure_comp(phy.set4_vs_set7$phy.set4) # Given prevalence in set4
phy.set4_vs_set7$formula_set7 <- Structure_comp(phy.set4_vs_set7$phy.set7) # Given prevalence in set7
batch_asvs_set4C_vs_set7 <- subset(phy.set4_vs_set7, phy.set7 < formula_set4) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set4C_vs_set7 <- batch_asvs_set4C_vs_set7[,-c(1,2,10,11)]
batch_asvs_set4C_vs_set7 <- setDT(batch_asvs_set4C_vs_set7, keep.rownames = "ASV")
batch_asvs_set4_vs_set7C <- subset(phy.set4_vs_set7, phy.set4 < formula_set7) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set4_vs_set7C <- batch_asvs_set4_vs_set7C[,-c(1,2,10,11)]
batch_asvs_set4_vs_set7C <- setDT(batch_asvs_set4_vs_set7C, keep.rownames = "ASV")

phy.set4_vs_set8 <- data.frame(phy.set4=taxa_sums(phy.set4)/N4, phy.set8=taxa_sums(phy.set8)/N8, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set4_vs_set8$phy.set4>0 | phy.set4_vs_set8$phy.set8>0)
phy.set4_vs_set8$formula_set4 <- Structure_comp(phy.set4_vs_set8$phy.set4) # Given prevalence in set4
phy.set4_vs_set8$formula_set8 <- Structure_comp(phy.set4_vs_set8$phy.set8) # Given prevalence in set8
batch_asvs_set4C_vs_set8 <- subset(phy.set4_vs_set8, phy.set8 < formula_set4) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set4C_vs_set8 <- batch_asvs_set4C_vs_set8[,-c(1,2,10,11)]
batch_asvs_set4C_vs_set8 <- setDT(batch_asvs_set4C_vs_set8, keep.rownames = "ASV")
batch_asvs_set4_vs_set8C <- subset(phy.set4_vs_set8, phy.set4 < formula_set8) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set4_vs_set8C <- batch_asvs_set4_vs_set8C[,-c(1,2,10,11)]
batch_asvs_set4_vs_set8C <- setDT(batch_asvs_set4_vs_set8C, keep.rownames = "ASV")

set4_asvs <- unique(rbind(batch_asvs_set4C_vs_set5, batch_asvs_set4_vs_set5C, batch_asvs_set4C_vs_set6, batch_asvs_set4_vs_set6C, 
                          batch_asvs_set4C_vs_set7, batch_asvs_set4_vs_set7C, batch_asvs_set4C_vs_set8, batch_asvs_set4_vs_set8C))

# Set 5 batch contaminants

phy.set5_vs_set6 <- data.frame(phy.set5=taxa_sums(phy.set5)/N5, phy.set6=taxa_sums(phy.set6)/N6, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set5_vs_set6$phy.set5>0 | phy.set5_vs_set6$phy.set6>0)
phy.set5_vs_set6$formula_set5 <- Structure_comp(phy.set5_vs_set6$phy.set5) # Given prevalence in set5
phy.set5_vs_set6$formula_set6 <- Structure_comp(phy.set5_vs_set6$phy.set6) # Given prevalence in set6
batch_asvs_set5C_vs_set6 <- subset(phy.set5_vs_set6, phy.set6 < formula_set5) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set5C_vs_set6 <- batch_asvs_set5C_vs_set6[,-c(1,2,10,11)]
batch_asvs_set5C_vs_set6 <- setDT(batch_asvs_set5C_vs_set6, keep.rownames = "ASV")
batch_asvs_set5_vs_set6C <- subset(phy.set5_vs_set6, phy.set5 < formula_set6) # Potential contaminants in set6 (shaded area on x-axis)
batch_asvs_set5_vs_set6C <- batch_asvs_set5_vs_set6C[,-c(1,2,10,11)]
batch_asvs_set5_vs_set6C <- setDT(batch_asvs_set5_vs_set6C, keep.rownames = "ASV")

phy.set5_vs_set7 <- data.frame(phy.set5=taxa_sums(phy.set5)/N5, phy.set7=taxa_sums(phy.set7)/N7, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set5_vs_set7$phy.set5>0 | phy.set5_vs_set7$phy.set7>0)
phy.set5_vs_set7$formula_set5 <- Structure_comp(phy.set5_vs_set7$phy.set5) # Given prevalence in set5
phy.set5_vs_set7$formula_set7 <- Structure_comp(phy.set5_vs_set7$phy.set7) # Given prevalence in set7
batch_asvs_set5C_vs_set7 <- subset(phy.set5_vs_set7, phy.set7 < formula_set5) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set5C_vs_set7 <- batch_asvs_set5C_vs_set7[,-c(1,2,10,11)]
batch_asvs_set5C_vs_set7 <- setDT(batch_asvs_set5C_vs_set7, keep.rownames = "ASV")
batch_asvs_set5_vs_set7C <- subset(phy.set5_vs_set7, phy.set5 < formula_set7) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set5_vs_set7C <- batch_asvs_set5_vs_set7C[,-c(1,2,10,11)]
batch_asvs_set5_vs_set7C <- setDT(batch_asvs_set5_vs_set7C, keep.rownames = "ASV")

phy.set5_vs_set8 <- data.frame(phy.set5=taxa_sums(phy.set5)/N5, phy.set8=taxa_sums(phy.set8)/N8, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set5_vs_set8$phy.set5>0 | phy.set5_vs_set8$phy.set8>0)
phy.set5_vs_set8$formula_set5 <- Structure_comp(phy.set5_vs_set8$phy.set5) # Given prevalence in set5
phy.set5_vs_set8$formula_set8 <- Structure_comp(phy.set5_vs_set8$phy.set8) # Given prevalence in set8
batch_asvs_set5C_vs_set8 <- subset(phy.set5_vs_set8, phy.set8 < formula_set5) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set5C_vs_set8 <- batch_asvs_set5C_vs_set8[,-c(1,2,10,11)]
batch_asvs_set5C_vs_set8 <- setDT(batch_asvs_set5C_vs_set8, keep.rownames = "ASV")
batch_asvs_set5_vs_set8C <- subset(phy.set5_vs_set8, phy.set5 < formula_set8) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set5_vs_set8C <- batch_asvs_set5_vs_set8C[,-c(1,2,10,11)]
batch_asvs_set5_vs_set8C <- setDT(batch_asvs_set5_vs_set8C, keep.rownames = "ASV")

set5_asvs <- unique(rbind(batch_asvs_set5C_vs_set6, batch_asvs_set5_vs_set6C, batch_asvs_set5C_vs_set7, batch_asvs_set5_vs_set7C, 
                          batch_asvs_set5C_vs_set8, batch_asvs_set5_vs_set8C))

# Set 6 batch contaminants

phy.set6_vs_set7 <- data.frame(phy.set6=taxa_sums(phy.set6)/N6, phy.set7=taxa_sums(phy.set7)/N7, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set6_vs_set7$phy.set6>0 | phy.set6_vs_set7$phy.set7>0)
phy.set6_vs_set7$formula_set6 <- Structure_comp(phy.set6_vs_set7$phy.set6) # Given prevalence in set6
phy.set6_vs_set7$formula_set7 <- Structure_comp(phy.set6_vs_set7$phy.set7) # Given prevalence in set7
batch_asvs_set6C_vs_set7 <- subset(phy.set6_vs_set7, phy.set7 < formula_set6) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set6C_vs_set7 <- batch_asvs_set6C_vs_set7[,-c(1,2,10,11)]
batch_asvs_set6C_vs_set7 <- setDT(batch_asvs_set6C_vs_set7, keep.rownames = "ASV")
batch_asvs_set6_vs_set7C <- subset(phy.set6_vs_set7, phy.set6 < formula_set7) # Potential contaminants in set7 (shaded area on x-axis)
batch_asvs_set6_vs_set7C <- batch_asvs_set6_vs_set7C[,-c(1,2,10,11)]
batch_asvs_set6_vs_set7C <- setDT(batch_asvs_set6_vs_set7C, keep.rownames = "ASV")

phy.set6_vs_set8 <- data.frame(phy.set6=taxa_sums(phy.set6)/N6, phy.set8=taxa_sums(phy.set8)/N8, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set6_vs_set8$phy.set6>0 | phy.set6_vs_set8$phy.set8>0)
phy.set6_vs_set8$formula_set6 <- Structure_comp(phy.set6_vs_set8$phy.set6) # Given prevalence in set6
phy.set6_vs_set8$formula_set8 <- Structure_comp(phy.set6_vs_set8$phy.set8) # Given prevalence in set8
batch_asvs_set6C_vs_set8 <- subset(phy.set6_vs_set8, phy.set8 < formula_set6) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set6C_vs_set8 <- batch_asvs_set6C_vs_set8[,-c(1,2,10,11)]
batch_asvs_set6C_vs_set8 <- setDT(batch_asvs_set6C_vs_set8, keep.rownames = "ASV")
batch_asvs_set6_vs_set8C <- subset(phy.set6_vs_set8, phy.set6 < formula_set8) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set6_vs_set8C <- batch_asvs_set6_vs_set8C[,-c(1,2,10,11)]
batch_asvs_set6_vs_set8C <- setDT(batch_asvs_set6_vs_set8C, keep.rownames = "ASV")

set6_asvs <- unique(rbind(batch_asvs_set6C_vs_set7, batch_asvs_set6_vs_set7C, batch_asvs_set6C_vs_set8, batch_asvs_set6_vs_set8C))

# Set 7 batch contaminants

phy.set7_vs_set8 <- data.frame(phy.set7=taxa_sums(phy.set7)/N7, phy.set8=taxa_sums(phy.set8)/N8, 
                               taxonomy=as(tax_table(phy.decontam.batch), "matrix"))
sum(phy.set7_vs_set8$phy.set7>0 | phy.set7_vs_set8$phy.set8>0)
phy.set7_vs_set8$formula_set7 <- Structure_comp(phy.set7_vs_set8$phy.set7) # Given prevalence in set7
phy.set7_vs_set8$formula_set8 <- Structure_comp(phy.set7_vs_set8$phy.set8) # Given prevalence in set8
batch_asvs_set7C_vs_set8 <- subset(phy.set7_vs_set8, phy.set8 < formula_set7) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set7C_vs_set8 <- batch_asvs_set7C_vs_set8[,-c(1,2,10,11)]
batch_asvs_set7C_vs_set8 <- setDT(batch_asvs_set7C_vs_set8, keep.rownames = "ASV")
batch_asvs_set7_vs_set8C <- subset(phy.set7_vs_set8, phy.set7 < formula_set8) # Potential contaminants in set8 (shaded area on x-axis)
batch_asvs_set7_vs_set8C <- batch_asvs_set7_vs_set8C[,-c(1,2,10,11)]
batch_asvs_set7_vs_set8C <- setDT(batch_asvs_set7_vs_set8C, keep.rownames = "ASV")

set7_asvs <- unique(rbind(batch_asvs_set7C_vs_set8, batch_asvs_set7_vs_set8C))

batch_asvs <- unique(rbind(set1_asvs, set2_asvs, set3_asvs, set4_asvs, set5_asvs, set6_asvs, set7_asvs))
nrow(batch_asvs)

batch_asv_ids <- batch_asvs$ASV
all_asvs <- taxa_names(phy.decontam.prev) 
keep_asvs <- all_asvs[!(all_asvs %in% batch_asv_ids)] 
phy.decontam.batch <- prune_taxa(keep_asvs, phy.decontam.prev)
ntaxa(phy.decontam.batch)
batch_reads <- sum(sample_sums(phy.decontam.batch))
batch_reads/decontam_reads
batch_reads/total_reads
saveRDS(phy.decontam.batch, "phy.bots.nps.batch.rds")

##########################
# PRUNING LOW-READ SAMPLES
##########################

# Check for any duplicate SampleIDs
phy.decontam.batch <- readRDS("phy.bots.nps.batch.rds")
duplicates <- data.frame(sample_data(phy.decontam.batch))
n_occur <- data.frame(table(duplicates$SampleID))
n_occur[n_occur$Freq > 1,]
# This should return <0 rows> if there are no duplicates
remove(duplicates, n_occur)

# Examine final phyloseq object
nsamples(phy.decontam.batch)
ntaxa(phy.decontam.batch)
metadata <- data.frame(sample_data(phy.decontam.batch))
table(metadata$batch, useNA="always")
metadata$Source <- as.factor(metadata$Source)
summary(metadata$Source, useNA="always")
metadata$SampleType <- as.factor(metadata$SampleType)
summary(metadata$SampleType, useNA="always")
summary(metadata$month, useNA="always")
table(metadata$batch, useNA="always")

# Creates dataframes with # of reads for NP samples
np_counts <- data.frame(samples=sample_names(phy.decontam.batch), 
                        counts=sample_sums(phy.decontam.batch))
np_richness <- estimate_richness(phy.decontam.batch, measure=c("Shannon", "Observed"))
np_richness$SampleID <- rownames(np_richness)
np_richness <- cbind(np_richness, np_counts)

# Plots of alpha diversity by number of reads for NP SWABS
nps_shannon <- ggplot(np_richness, aes(counts, Shannon))+
  geom_point() +
  scale_x_log10()
png(file="bots_nps_shannon.png", 
    width = 10, height = 6, units = 'in', res = 600)
print(nps_shannon)
dev.off()
nps_observed <- ggplot(np_richness, aes(counts, Observed))+
  geom_point() +
  scale_x_log10()
png(file="bots_nps_observed.png", 
    width = 10, height = 6, units = 'in', res = 600)
print(nps_observed)
dev.off()
# BASED ON REVIEW OF THESE DIVERSITY PLOTS, USE A CUTOFF OF 1000 READS FOR FILTERING OF NP SWABS

# Prune NP samples that don't have at least 1000 reads
phy.np.less1000 <- prune_samples(sample_sums(phy.decontam.batch)<1000, phy.decontam.batch)
nsamples(phy.np.less1000)
metadata_np_less1000 <- data.frame(sample_data(phy.np.less1000))
table(metadata_np_less1000$SampleType, metadata_np_less1000$month)
# Samples from 0 months of age are markedly overrepresented 
phy.np.1000 <- prune_samples(sample_sums(phy.decontam.batch)>=1000, phy.decontam.batch)
nsamples(phy.np.1000)

phy.bots.nps.pruned <- subset_samples(phy.np.1000, Source=="INF" | month=="0") # removes maternal samples from visits other than birth
# Remove samples and taxa with 0 reads 
nsamples(phy.bots.nps.pruned)
ntaxa(phy.bots.nps)
phy.bots.nps.pruned <- prune_samples(sample_sums(phy.bots.nps.pruned)>0, phy.bots.nps.pruned)
phy.bots.nps.pruned <- filter_taxa(phy.bots.nps.pruned, function(x) sum(x) > 0, TRUE)
nsamples(phy.bots.nps.pruned)
ntaxa(phy.bots.nps.pruned)
metadata_np <- data.frame(sample_data(phy.bots.nps.pruned))
table(metadata_np$Source, metadata_np$month)
taxtable_final <- data.frame(tax_table(phy.bots.nps.pruned))
remove(phy.np.1000, phy.np.less1000, np_counts, np_richness, metadata_np_less1000, metadata_np)

saveRDS(phy.bots.nps.pruned, "phy.bots.nps.pruned.rds")

###############################################
# PREPROCESSING OF VIRAL AND BACTERIAL METADATA
###############################################

metadata_inf_np_RV <- read.csv("metadata_inf_np.csv")
metadata_inf_np_RV$inf_sp_yn[metadata_inf_np_RV$inf_sp>0] <- "Y"
metadata_inf_np_RV$inf_sp_yn[metadata_inf_np_RV$inf_sp==0] <- "N"
# Create variable for number of PCV-13 doses received (do not include doses administered on the day of the study visit)
metadata_inf_np_RV$pcv[is.na(metadata_inf_np_RV$pcv1) & is.na(metadata_inf_np_RV$pcv2) & is.na(metadata_inf_np_RV$pcv3)] <- "0"
metadata_inf_np_RV$pcv[(metadata_inf_np_RV$age_days)<=(metadata_inf_np_RV$pcv1+1)] <- "0"  
metadata_inf_np_RV$pcv[(metadata_inf_np_RV$age_days)>(metadata_inf_np_RV$pcv1+1) & is.na(metadata_inf_np_RV$pcv2)] <- "1"  
metadata_inf_np_RV$pcv[metadata_inf_np_RV$age_days>(metadata_inf_np_RV$pcv1+1) & metadata_inf_np_RV$age_days<=(metadata_inf_np_RV$pcv2+1)] <- "1"
metadata_inf_np_RV$pcv[(metadata_inf_np_RV$age_days)>(metadata_inf_np_RV$pcv2+1) & is.na(metadata_inf_np_RV$pcv3)] <- "2"  
metadata_inf_np_RV$pcv[metadata_inf_np_RV$age_days>(metadata_inf_np_RV$pcv2+1) & metadata_inf_np_RV$age_days<=(metadata_inf_np_RV$pcv3+1)] <- "2"
metadata_inf_np_RV$pcv[metadata_inf_np_RV$age_days>(metadata_inf_np_RV$pcv3+1)] <- "3"
metadata_inf_np_RV$pcv[metadata_inf_np_RV$age_days<60 & is.na(metadata_inf_np_RV$pcv)] <- "0"
metadata_inf_np_RV$pcv <- as.numeric(metadata_inf_np_RV$pcv)
table(metadata_inf_np_RV$month, metadata_inf_np_RV$pcv, useNA="always")
# Create variable for number of Hib doses received (do not include doses administered on the day of the study visit)
metadata_inf_np_RV$hib[is.na(metadata_inf_np_RV$dpt1) & is.na(metadata_inf_np_RV$dpt2) & is.na(metadata_inf_np_RV$dpt3)] <- "0"
metadata_inf_np_RV$hib[(metadata_inf_np_RV$age_days)<=(metadata_inf_np_RV$dpt1+1)] <- "0"  
metadata_inf_np_RV$hib[(metadata_inf_np_RV$age_days)>(metadata_inf_np_RV$dpt1+1) & is.na(metadata_inf_np_RV$dpt2)] <- "1"  
metadata_inf_np_RV$hib[metadata_inf_np_RV$age_days>(metadata_inf_np_RV$dpt1+1) & metadata_inf_np_RV$age_days<=(metadata_inf_np_RV$dpt2+1)] <- "1"
metadata_inf_np_RV$hib[(metadata_inf_np_RV$age_days)>(metadata_inf_np_RV$dpt2+1) & is.na(metadata_inf_np_RV$dpt3)] <- "2"  
metadata_inf_np_RV$hib[metadata_inf_np_RV$age_days>(metadata_inf_np_RV$dpt2+1) & metadata_inf_np_RV$age_days<=(metadata_inf_np_RV$dpt3+1)] <- "2"
metadata_inf_np_RV$hib[metadata_inf_np_RV$age_days>(metadata_inf_np_RV$dpt3+1)] <- "3"
metadata_inf_np_RV$hib[metadata_inf_np_RV$age_days<60 & is.na(metadata_inf_np_RV$dpt)] <- "0"
metadata_inf_np_RV$hib <- as.numeric(metadata_inf_np_RV$hib)
table(metadata_inf_np_RV$month, metadata_inf_np_RV$hib, useNA="always")
metadata_inf_np_RV$recent_uri[metadata_inf_np_RV$current_uri=="Y"] <- "Y"
table(metadata_inf_np_RV$recent_uri, metadata_inf_np_RV$current_uri)

table(metadata_inf_np_RV$inf_rv)
metadata_inf_np_RV$inf_rv_yn[metadata_inf_np_RV$inf_rv=="NEG"] <- "N"
metadata_inf_np_RV$inf_rv_yn[metadata_inf_np_RV$inf_rv!="NEG" & !is.na(metadata_inf_np_RV$inf_rv)] <- "Y"
metadata_inf_np_RV$inf_rv_cat[metadata_inf_np_RV$inf_rv=="NEG"] <- "No viruses"
metadata_inf_np_RV$inf_rv_cat[metadata_inf_np_RV$inf_rv=="ADENO"] <- "Adenovirus"
metadata_inf_np_RV$inf_rv_cat[metadata_inf_np_RV$inf_rv=="ER"] <- "Rhinovirus/enterovirus"
metadata_inf_np_RV$inf_rv_cat[metadata_inf_np_RV$inf_rv=="FLU A" | metadata_inf_np_RV$inf_rv=="FLU B"] <- "Influenza"
metadata_inf_np_RV$inf_rv_cat[metadata_inf_np_RV$inf_rv=="MPV"] <- "HMPV"
metadata_inf_np_RV$inf_rv_cat[metadata_inf_np_RV$inf_rv=="P1" | metadata_inf_np_RV$inf_rv=="P2" | metadata_inf_np_RV$inf_rv=="P3"] <- "Parainfluenza"
metadata_inf_np_RV$inf_rv_cat[metadata_inf_np_RV$inf_rv=="RSV"] <- "RSV"
metadata_inf_np_RV$inf_rv_cat[metadata_inf_np_RV$inf_rv=="SARS2"] <- "SARS2"
metadata_inf_np_RV$inf_rv_cat[metadata_inf_np_RV$inf_rv=="ADENO/MPV" | metadata_inf_np_RV$inf_rv=="ADENO/P3" | 
                                metadata_inf_np_RV$inf_rv=="ER/ADENO" | metadata_inf_np_RV$inf_rv=="ER/P1" | metadata_inf_np_RV$inf_rv=="ER/P2" | 
                                metadata_inf_np_RV$inf_rv=="ER/P3" | metadata_inf_np_RV$inf_rv=="ER/RSV" | metadata_inf_np_RV$inf_rv=="ER/RSV/MPV" | 
                                metadata_inf_np_RV$inf_rv=="FLU A/MPV" | metadata_inf_np_RV$inf_rv=="MPV/P3" | metadata_inf_np_RV$inf_rv=="RSV/ADENO"] <- ">1 virus"
metadata_inf_np_RV$inf_rv_cat <- factor(metadata_inf_np_RV$inf_rv_cat, levels = c("No viruses", "Adenovirus", "HMPV", "Influenza", 
                                                                                  "Parainfluenza", "Rhinovirus/enterovirus", "RSV", "SARS2", ">1 virus"))
table(metadata_inf_np_RV$inf_rv_cat, metadata_inf_np_RV$inf_rv)

metadata_inf_np_RV$adeno[metadata_inf_np_RV$inf_rv=="ADENO" | metadata_inf_np_RV$inf_rv=="ADENO/MPV" | metadata_inf_np_RV$inf_rv=="ADENO/P3" |
                           metadata_inf_np_RV$inf_rv=="ER/ADENO" | metadata_inf_np_RV$inf_rv=="RSV/ADENO"] <- "Y"
metadata_inf_np_RV$adeno[metadata_inf_np_RV$inf_rv!="ADENO" & metadata_inf_np_RV$inf_rv!="ADENO/MPV" & metadata_inf_np_RV$inf_rv!="ADENO/P3" &
                           metadata_inf_np_RV$inf_rv!="ER/ADENO" & metadata_inf_np_RV$inf_rv!="RSV/ADENO"] <- "N"
metadata_inf_np_RV$rhino[metadata_inf_np_RV$inf_rv=="ER" | metadata_inf_np_RV$inf_rv=="ER/ADENO" | metadata_inf_np_RV$inf_rv=="ER/P2" |
                           metadata_inf_np_RV$inf_rv=="ER/P3" | metadata_inf_np_RV$inf_rv=="ER/RSV" | metadata_inf_np_RV$inf_rv=="ER/RSV/MPV"] <- "Y"
metadata_inf_np_RV$rhino[metadata_inf_np_RV$inf_rv!="ER" & metadata_inf_np_RV$inf_rv!="ER/ADENO" & metadata_inf_np_RV$inf_rv!="ER/P2" &
                           metadata_inf_np_RV$inf_rv!="ER/P3" & metadata_inf_np_RV$inf_rv!="ER/RSV" & metadata_inf_np_RV$inf_rv!="ER/RSV/MPV"] <- "N"
metadata_inf_np_RV$flu[metadata_inf_np_RV$inf_rv=="FLU A" | metadata_inf_np_RV$inf_rv=="FLU A/MPV" | metadata_inf_np_RV$inf_rv=="FLU B"] <- "Y"
metadata_inf_np_RV$flu[metadata_inf_np_RV$inf_rv!="FLU A" & metadata_inf_np_RV$inf_rv!="FLU A/MPV" & metadata_inf_np_RV$inf_rv!="FLU B"] <- "N"
metadata_inf_np_RV$hmpv[metadata_inf_np_RV$inf_rv=="ADENO/MPV" | metadata_inf_np_RV$inf_rv=="ER/RSV/MPV" | metadata_inf_np_RV$inf_rv=="FLU A/MPV" |
                          metadata_inf_np_RV$inf_rv=="MPV" | metadata_inf_np_RV$inf_rv=="MPV/P3"] <- "Y"
metadata_inf_np_RV$hmpv[metadata_inf_np_RV$inf_rv!="ADENO/MPV" & metadata_inf_np_RV$inf_rv!="ER/RSV/MPV" & metadata_inf_np_RV$inf_rv!="FLU A/MPV" &
                          metadata_inf_np_RV$inf_rv!="MPV" & metadata_inf_np_RV$inf_rv!="MPV/P3"] <- "N"
metadata_inf_np_RV$paraflu[metadata_inf_np_RV$inf_rv=="ADENO/P3" | metadata_inf_np_RV$inf_rv=="ER/P1" | metadata_inf_np_RV$inf_rv=="ER/P2" | metadata_inf_np_RV$inf_rv=="ER/P3" |
                             metadata_inf_np_RV$inf_rv=="MPV/P3" | metadata_inf_np_RV$inf_rv=="P1" | metadata_inf_np_RV$inf_rv=="P2" | metadata_inf_np_RV$inf_rv=="P3"] <- "Y"
metadata_inf_np_RV$paraflu[metadata_inf_np_RV$inf_rv!="ADENO/P3" & metadata_inf_np_RV$inf_rv!="ER/P1" & metadata_inf_np_RV$inf_rv!="ER/P2" & metadata_inf_np_RV$inf_rv!="ER/P3" &
                             metadata_inf_np_RV$inf_rv!="MPV/P3" & metadata_inf_np_RV$inf_rv!="P1" & metadata_inf_np_RV$inf_rv!="P2" & metadata_inf_np_RV$inf_rv!="P3"] <- "N"
metadata_inf_np_RV$rsv[metadata_inf_np_RV$inf_rv=="ER/RSV" | metadata_inf_np_RV$inf_rv=="ER/RSV/MPV" | metadata_inf_np_RV$inf_rv=="RSV" |
                         metadata_inf_np_RV$inf_rv=="RSV/ADENO"] <- "Y"
metadata_inf_np_RV$rsv[metadata_inf_np_RV$inf_rv!="ER/RSV" & metadata_inf_np_RV$inf_rv!="ER/RSV/MPV" & metadata_inf_np_RV$inf_rv!="RSV" &
                         metadata_inf_np_RV$inf_rv!="RSV/ADENO"] <- "N"
metadata_inf_np_RV$sars2[metadata_inf_np_RV$inf_rv=="SARS2"] <- "Y"
metadata_inf_np_RV$sars2[metadata_inf_np_RV$inf_rv!="SARS2"] <- "N"

metadata_inf_np_RV$month <- as.numeric(metadata_inf_np_RV$month)
metadata_inf_np_RV <- metadata_inf_np_RV[order(metadata_inf_np_RV$study_id, metadata_inf_np_RV$month),]
# Create variables for acquisition of bacterial respiratory pathogens
metadata_inf_np_RV <- suppressWarnings(slide(metadata_inf_np_RV, "inf_multi_hi", TimeVar="month", GroupVar="study_id", NewVar="inf_multi_hi_lag", slideBy = -1))
metadata_inf_np_RV$inf_hi_new[metadata_inf_np_RV$inf_multi_hi=="Y" & (metadata_inf_np_RV$inf_multi_hi_lag=="N" | is.na(metadata_inf_np_RV$inf_multi_hi_lag))] <- "Y"
metadata_inf_np_RV$inf_hi_new[metadata_inf_np_RV$inf_multi_hi=="Y" & metadata_inf_np_RV$inf_multi_hi_lag=="Y"] <- "N"
metadata_inf_np_RV$inf_hi_new[metadata_inf_np_RV$inf_multi_hi=="Y" & metadata_inf_np_RV$month==0] <- "Y"
metadata_inf_np_RV$inf_hi_new[metadata_inf_np_RV$inf_multi_hi=="N"] <- "N"
table(metadata_inf_np_RV$inf_hi_new, metadata_inf_np_RV$inf_multi_hi, useNA="always")
metadata_inf_np_RV <- suppressWarnings(slide(metadata_inf_np_RV, "inf_multi_mc", TimeVar="month", GroupVar="study_id", NewVar="inf_multi_mc_lag", slideBy = -1))
metadata_inf_np_RV$inf_mc_new[metadata_inf_np_RV$inf_multi_mc=="Y" & (metadata_inf_np_RV$inf_multi_mc_lag=="N" | is.na(metadata_inf_np_RV$inf_multi_mc_lag))] <- "Y"
metadata_inf_np_RV$inf_mc_new[metadata_inf_np_RV$inf_multi_mc=="Y" & metadata_inf_np_RV$inf_multi_mc_lag=="Y"] <- "N"
metadata_inf_np_RV$inf_mc_new[metadata_inf_np_RV$inf_multi_mc=="Y" & metadata_inf_np_RV$month==0] <- "Y"
metadata_inf_np_RV$inf_mc_new[metadata_inf_np_RV$inf_multi_mc=="N"] <- "N"
table(metadata_inf_np_RV$inf_mc_new, metadata_inf_np_RV$inf_multi_mc, useNA="always")
metadata_inf_np_RV <- suppressWarnings(slide(metadata_inf_np_RV, "inf_multi_sa", TimeVar="month", GroupVar="study_id", NewVar="inf_multi_sa_lag", slideBy = -1))
metadata_inf_np_RV$inf_sa_new[metadata_inf_np_RV$inf_multi_sa=="Y" & (metadata_inf_np_RV$inf_multi_sa_lag=="N" | is.na(metadata_inf_np_RV$inf_multi_sa_lag))] <- "Y"
metadata_inf_np_RV$inf_sa_new[metadata_inf_np_RV$inf_multi_sa=="Y" & metadata_inf_np_RV$inf_multi_sa_lag=="Y"] <- "N"
metadata_inf_np_RV$inf_sa_new[metadata_inf_np_RV$inf_multi_sa=="Y" & metadata_inf_np_RV$month==0] <- "Y"
metadata_inf_np_RV$inf_sa_new[metadata_inf_np_RV$inf_multi_sa=="N"] <- "N"
table(metadata_inf_np_RV$inf_sa_new, metadata_inf_np_RV$inf_multi_sa, useNA="always")
metadata_inf_np_RV <- suppressWarnings(slide(metadata_inf_np_RV, "inf_multi_sp", TimeVar="month", GroupVar="study_id", NewVar="inf_multi_sp_lag", slideBy = -1))
metadata_inf_np_RV$inf_sp_new[metadata_inf_np_RV$inf_multi_sp=="Y" & (metadata_inf_np_RV$inf_multi_sp_lag=="N" | is.na(metadata_inf_np_RV$inf_multi_sp_lag))] <- "Y"
metadata_inf_np_RV$inf_sp_new[metadata_inf_np_RV$inf_multi_sp=="Y" & metadata_inf_np_RV$inf_multi_sp_lag=="Y"] <- "N"
metadata_inf_np_RV$inf_sp_new[metadata_inf_np_RV$inf_multi_sp=="Y" & metadata_inf_np_RV$month==0] <- "Y"
metadata_inf_np_RV$inf_sp_new[metadata_inf_np_RV$inf_multi_sp=="N"] <- "N"
table(metadata_inf_np_RV$inf_sp_new, metadata_inf_np_RV$inf_multi_sp, useNA="always")

# Create variable for infections caused by rhino/entero and other respiratory viruses
metadata_inf_np_RV$inf_rv_cat2[metadata_inf_np_RV$inf_rv_cat=="Rhinovirus/enterovirus"] <- "Rhino"
metadata_inf_np_RV$inf_rv_cat2[metadata_inf_np_RV$inf_rv_cat=="No viruses"] <- "None"
metadata_inf_np_RV$inf_rv_cat2[metadata_inf_np_RV$inf_rv_cat!="Rhinovirus/enterovirus" & metadata_inf_np_RV$inf_rv_cat!="No viruses"] <- "Other"
table(metadata_inf_np_RV$inf_rv_cat, metadata_inf_np_RV$inf_rv_cat2)

metadata_inf_np_RV <- metadata_inf_np_RV[,c("SampleID", "study_id", "month", "age_days", "calendar_mo", "season", "sex", "bw", "lbw", "inf_hiv", "mat_hiv",
                                            "residence", "mat_educ", "water", "electric", "fridge", "wood", "num_total", "num_kids", "num_adol", "num_adults", 
                                            "breastmilk", "current_uri", "recent_uri", "pneumonia", "inf_abx_any", "inf_abx_amox", "inf_abx_metro", "inf_abx_cotrim", 
                                            "pcv", "hib", "inf_rv_yn", "inf_rv_cat", "inf_rv_cat2", "inf_rv", "adeno", "flu", "hmpv", "paraflu", "rhino", "rsv", "sars2",
                                            "inf_multi_hi", "inf_hi_new", "inf_multi_mc", "inf_mc_new", "inf_multi_sa", "inf_sa_new", "inf_multi_sp", "inf_sp_new",
                                            "inf_sp_yn", "inf_sp", "batch")]

study_ids <- unique(metadata_inf_np_RV$study_id)
nrow(metadata_inf_np_RV)
write.csv(metadata_inf_np_RV, "metadata_inf_np_RV.csv")

###############################
# PREPROCESSING OF 16S METADATA
###############################

phy.bots.nps.pruned <- readRDS("phy.bots.nps.pruned.rds")
phy.bots.nps.pruned.inf <- subset_samples(phy.bots.nps.pruned, Source=="INF")
nsamples(phy.bots.nps.pruned.inf)
ntaxa(phy.bots.nps.pruned.inf)
summary(sample_sums(phy.bots.nps.pruned.inf))

# Create dataframe with Shannon index and observed ASVs for each specimen
diversity <- estimate_richness(phy.bots.nps.pruned.inf, measures = c("Shannon", "Observed"))
setDT(diversity, keep.rownames = TRUE)[]
diversity$SampleID <- diversity$rn
diversity$rn <- NULL
diversity$SampleID <- gsub("[.]", "-", diversity$SampleID)
metadata_16s <- merge(metadata_inf_np_RV, diversity, by="SampleID")
summary(metadata_16s$Observed)
summary(metadata_16s$Shannon)
row.names(metadata_16s) <- metadata_16s$SampleID
sample_data(phy.bots.nps.pruned.inf) <- metadata_16s

# Create file for BLAST searches for specific ASVs in order to provide further information on species/species groups
refseq <- data.frame(refseq(phy.bots.nps.pruned.inf))
refseq <- tibble::rownames_to_column(refseq, "ASV")
write.csv(refseq, "ASV_refseqs.csv")

# Do not rely upon species assignments from eHOMD but instead base these on BLASTn searches
taxtable <- data.frame(tax_table(phy.bots.nps.pruned.inf))
names(taxtable)[names(taxtable) == "Species"] <- "Species_eHOMD"
taxtable$ASV <- as.factor(row.names(taxtable))
tax_table(phy.bots.nps.pruned.inf) <- as.matrix(taxtable)

# Use k-medoids clustering algorithm based on Bray-Curtis distances

# Agglomerate taxa into genera 
phy.agglom <- tax_glom(phy.bots.nps.pruned.inf, taxrank = 'Genus')
ntaxa(phy.agglom)

# Transform to relative abundances
phy.relative <- transform_sample_counts(phy.agglom, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.relative))

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)}
  return(as(OTU, "matrix"))}

pam.clustering <- function(x,k) { # x is a distance matrix and k is the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)}

noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe -> Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)}

# Use Bray-Curtis as the distance measure
data <- t(as.data.frame(vegan_otu(phy.relative)))
data.dist <- vegdist(vegan_otu(phy.relative), method = "bray")

# Determine the optimal number of clusters
data.cluster_temp <- pam.clustering(data.dist, k=5)

# Determine the optimal number of clusters
nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data), data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")}}
par(mar = c(1, 1, 1, 1))
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
axis(1,at=1:20)

# Set the number of clusters = to k based on above plot
data.cluster <- pam.clustering(data.dist, k=8)
table(data.cluster)
sample_data(phy.bots.nps.pruned.inf)$cluster_num <- as.factor(data.cluster)
metadata_cluster <- data.frame(sample_data(phy.bots.nps.pruned.inf))
metadata_cluster$cluster <- ""

# Rename clusters based on composition
sample_data(phy.bots.nps.pruned.inf) <- metadata_cluster
# Agglomerate taxa into genera for compositional comparisons
phy.agglom <- tax_glom(phy.bots.nps.pruned.inf, taxrank = 'Genus')
ntaxa(phy.agglom)
# Transform to relative abundances
phy.relative <- transform_sample_counts(phy.agglom, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.relative))

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
TOPGenera <- unique(genera_abundances$Genus[1:16])
genus_df <- genera_abundances[genera_abundances$Genus %in% TOPGenera,]
genus_df$Genus <- factor(genus_df$Genus, levels = genus_df$Genus[order(-genus_df$genus_Ab)])
head(genera_abundances, 25)
# Rename genera other than top genera as "Other" in creating dataframe relative_inf 
relative_inf$Genus[!(relative_inf$Genus %in% TOPGenera)] <- "Other"
TOPGenera
relative_inf$Genus[relative_inf$Genus=="Bacteria;p;c;o;f;g"] <- "Other"
relative_inf$Genus[relative_inf$Genus=="Enterobacteriaceae;g"] <- "Other"
unique(relative_inf$Genus)
sum(relative_inf$Abundance)
relative_inf$Genus <- as.factor(relative_inf$Genus)
table(relative_inf$Genus)

theme_barplot <-   theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                         axis.title.x = element_text(size=9, color="black"), 
                         axis.text.y = element_text(size=8, color="black"), 
                         axis.title.y = element_text(angle=90, size=9, color="black"), 
                         axis.text.x = element_text(size=8, color="black"),
                         legend.position = "right", legend.text=element_text(size=8, color="black"), legend.text.align = 0,
                         legend.box.spacing = unit(0, "pt"),
                         legend.title=element_blank(), legend.key.size = unit(0.35, 'cm'), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) 

earthy_cols_15 <- c("rosybrown2", "gray50", "darkslateblue", "indianred4", "forestgreen", "cornflowerblue", "peru", "indianred1", "darkolivegreen4", 
                    "dodgerblue3", "plum3", "navajowhite3",  "mediumorchid4",  "coral3", "gray10") 

ggplot(arrange(relative_inf, Genus), aes(x=cluster_num, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="fill") + theme_barplot + 
  scale_fill_manual(values=earthy_cols_15) + xlab("Microbiota profile") + ylab("Relative abundance") 

metadata_cluster$cluster[metadata_cluster$cluster_num=="1"] <- "COR"
metadata_cluster$cluster[metadata_cluster$cluster_num=="2"] <- "STR"
metadata_cluster$cluster[metadata_cluster$cluster_num=="3"] <- "STA"
metadata_cluster$cluster[metadata_cluster$cluster_num=="4"] <- "CD"
metadata_cluster$cluster[metadata_cluster$cluster_num=="5"] <- "MOR"
metadata_cluster$cluster[metadata_cluster$cluster_num=="6"] <- "CDM"
metadata_cluster$cluster[metadata_cluster$cluster_num=="7"] <- "HAE"
metadata_cluster$cluster[metadata_cluster$cluster_num=="8"] <- "OTH"
# Order microbiome clusters based on mean age
tapply(metadata_cluster$age_days, metadata_cluster$cluster, mean)
metadata_cluster$cluster <- factor(metadata_cluster$cluster, levels=c("OTH","STA","COR","STR","CD","HAE","CDM","MOR"))
metadata_cluster$cluster_num <- NULL

# Use of k-means clustering based on the default Euclidean distances as an alternative algorithm for identification of microbiota profiles
otutable <- data.frame(otu_table(phy.agglom))
otutable_norm <- decostand(otutable, "normalize") 
otutable.kmeans <- kmeans(otutable_norm, centers = 8, nstart = 100)
SampleID <- row.names(otutable_norm)
kmeans_cluster_num <- otutable.kmeans$cluster
kmeans <- data.frame(SampleID, kmeans_cluster_num)

relative_inf <- merge(relative_inf, kmeans, by="SampleID")

ggplot(arrange(relative_inf, Genus), aes(x=kmeans_cluster_num, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="fill") + theme_barplot + 
  scale_fill_manual(values=earthy_cols_15) + xlab("Microbiota profile") + ylab("Relative abundance") 

kmeans$kmeans_cluster[kmeans$kmeans_cluster_num=="1"] <- "HAE"
kmeans$kmeans_cluster[kmeans$kmeans_cluster_num=="2"] <- "STA"
kmeans$kmeans_cluster[kmeans$kmeans_cluster_num=="3"] <- "STR"
kmeans$kmeans_cluster[kmeans$kmeans_cluster_num=="4"] <- "COR"
kmeans$kmeans_cluster[kmeans$kmeans_cluster_num=="5"] <- "OTH"
kmeans$kmeans_cluster[kmeans$kmeans_cluster_num=="6"] <- "MOR"
kmeans$kmeans_cluster[kmeans$kmeans_cluster_num=="7"] <- "CDM"
kmeans$kmeans_cluster[kmeans$kmeans_cluster_num=="8"] <- "CD"

metadata_cluster <- merge(metadata_cluster, kmeans, by="SampleID")
metadata_cluster$kmeans_cluster_num <- NULL
table(metadata_cluster$cluster, metadata_cluster$kmeans_cluster)
row.names(metadata_cluster) <- metadata_cluster[,1]
sample_data(phy.bots.nps.pruned.inf) <- metadata_cluster
saveRDS(phy.bots.nps.pruned.inf, "phy.inf.np.16s.rds")

# Write fasta file for PICRUSt2 analyses
phy.inf.np.16s <- readRDS("phy.inf.np.16s.rds")
library(seqateurs)
seqateurs::ps_to_fasta(phy.inf.np.16s, "asv_refseqs_picrust2.fasta", "ASV")

# Write ASV table for PICRUSt2 analyses
otu_table <- data.frame(t(otu_table(phy.inf.np.16s)))
write.table(otu_table, file='asv_table_picrust2.tsv', quote=FALSE, sep='\t', col.names=NA)

# Transfer these files to the DCC and run the following command (first need to install picrust2 into conda environment)
# See https://github.com/picrust/picrust2/wiki/Installation for installation instructions
# Run the following code to generate picrust2 output
# picrust2_pipeline.py -s asv_refseqs_picrust2.fasta -i asv_table_picrust2.tsv -o picrust2 -p 1
# add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv.gz

# Read in PICRUSt2 output file
picrust <- read.table("PICRUSt2/pathways_out/path_abun_unstrat_descrip.tsv", sep = '\t', header = TRUE)
names(picrust) <- gsub(x = names(picrust), pattern = "\\.", replacement = "-") 
row.names(picrust) <- picrust[,1]
picrust_tax <- as.matrix(picrust[,c(1:2)])
picrust_otu <- picrust[,c(3:2237)]
phy.picrust <- phyloseq(otu_table(picrust_otu, taxa_are_rows=TRUE), sample_data(metadata_16s), tax_table(picrust_tax))
nsamples(phy.picrust)
ntaxa(phy.picrust)
phy.picrust <- prune_taxa(taxa_sums(phy.picrust) > 0, phy.picrust)
ntaxa(phy.picrust)
# All pathways are present in at least 1 sample
# Transform to relative pathway abundances
phy.picrust <- transform_sample_counts(phy.picrust, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.picrust))
saveRDS(phy.picrust, "phy.inf.np.picrust.rds")