# Botswana Infant Microbiome Study - RV-Bacterial Analyses
# Matthew Kelly, MD, MPH 
# Figure 4
# Last updated: October 18, 2024

remove(list=ls())
setwd("_____________________________") 
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
library(xlsx)
library(pROC)
packageVersion("pROC")
library(pheatmap)
library(EBImage)
library(grid)
library(gtable)
library(viridis)

# Summary of results of random forest analyses
roc_auc_month <- read.xlsx(file="randomForest_month_only.xlsx", sheetN="AUC")
roc_auc_month
# Age on its own can be used to predict colonization with accuracy better than chance for all four bacterial pathobionts
roc_auc_clinical <- read.xlsx(file="randomForest_clinical.xlsx", sheetN="AUC")
roc_auc_clinical
# The addition of other clinical variables improves the accuracy of all models
# The most important variables are month (particularly for Sa) and number of children in the household (particularly for Hi, Mc, and Sp)
# Surprisingly, respiratory virus detection was only the 3rd most important variable for prediction in all models (4th after PCV for Sp)
roc_auc_asv <- read.xlsx(file="randomForest_asv.xlsx", sheetN="AUC")
roc_auc_asv
# Models built only with microbiome data from the preceding sample perform as well as or better than models with clinical variables for Hi, Sa, and Sp (not Mc)
roc_auc_asv_pathway <- read.xlsx(file="randomForest_asv_pathway.xlsx", sheetN="AUC")
roc_auc_asv_pathway
# The addition of functional pathway data does not substantively improve the performance of any models based on microbiome data
roc_auc_asv_clinical <- read.xlsx(file="randomForest_asv_clinical.xlsx", sheetN="AUC")
roc_auc_asv_clinical
# The addition of clinical variables data does not substantively improve the performance of any models based on microbiome data
roc_auc_all <- read.xlsx(file="randomForest_all.xlsx", sheetN="AUC")
roc_auc_all
# Models based on only microbiome data have similar accuracy to a model incorporating clinical variables and functional pathways
# Variability in ability to predict colonization based on specific bacterial pathobiont

rf_palette <- c('#377eb8', '#e41a1c', '#984ea3', '#4caf4a', '#a65628', 'gray30')

pathogen_name <- c("Respiratory virus", "H. influenzae", "M. catarrhalis", "S. aureus", "S. pneumoniae")
tab_all_auc <- as.data.frame(rbind(roc_auc_month,
                                   roc_auc_clinical,
                                   roc_auc_asv,
                                   roc_auc_asv_pathway,
                                   roc_auc_asv_clinical,
                                   roc_auc_all))
vec_methods <- c("Age", "Clinical variables", "URT microbiota taxonomy", "Taxonomy + microbial pathways", "Taxonomy + clinical variables", "Full model")
tab_all_auc$method <- factor(rep(vec_methods, each=nrow(roc_auc_asv)), level=vec_methods)
tab_all_auc$pathogen <- rep(pathogen_name[-1], length(vec_methods))

fig_4a <- tab_all_auc %>%
  ggplot(aes(x=pathogen, color=method, y=AUC)) +
  geom_point(position=position_dodge(.75)) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.2,
                position=position_dodge(.75)) + ylim(0.5,0.8) +
  theme_bw() +  theme(legend.title=element_blank(), axis.text.x = element_text(size=7.5, face="italic", color="black"),
                      axis.text.y = element_text(size=7.5, color="black"), axis.title.y = element_text(size=8, color="black"),
                      legend.text = element_text(size=8, color="black"), legend.key.size = unit(0.3, "cm")) +
  labs(y="AUC-ROC", x="") +
  geom_hline(yintercept=0.5, linetype="dotted") + 
  scale_color_manual(values = rf_palette)

png(file="R_Plots/Figure_4/Figure_4a.png", width = 6.5, height = 3, units = 'in', res = 1200)
plot(fig_4a)
dev.off()

# Figure 4b - heatmap of importance of clinical variables and ASVs

var_infect <- c("inf_multi_hi", "inf_multi_mc", "inf_multi_sa", "inf_multi_sp")
var_clinical <- c("num_kids_trunc", "inf_rv_yn", "mat_hiv", "breastmilk", "inf_abx_any", "season", "lbw", "wood", "residence", "sex", "pcv")
count_filter <- read.csv("Pixu_Analyses/count_filter.csv", header=T, row.names=1)
var_ASV <- colnames(count_filter)
pathogen_name <- c("Respiratory virus", "H. influenzae", "M. catarrhalis", "S. aureus", "S. pneumoniae")

res_rf_asv_clinical <- readRDS("res_rf_asv_clinical.RDS")
roc_auc_asv_clinical <- matrix(NA, length(var_infect), 3)
colnames(roc_auc_asv_clinical) <- c("CI_lower", "AUC", "CI_upper")
rownames(roc_auc_asv_clinical) <- var_infect
importance_asv_clinical <- matrix(NA, 1+length(var_clinical)+length(var_ASV), length(var_infect))
rownames(importance_asv_clinical) <- c("month", var_clinical, var_ASV)
colnames(importance_asv_clinical) <- var_infect
for (ii in 1:length(var_infect)){
  res <- res_rf_asv_clinical[[ii]]$predprob
  predprob <- res$predprob
  vec_id <- unique(res$study_id)
  ind_Y <- res$observed=="Y"
  set.seed(1234)
  roc_auc_asv_clinical[ii,] <- ci.auc(roc(cases=predprob[ind_Y],
                                          controls=predprob[!ind_Y]), 
                                      method="bootstrap")
  if (var_infect[ii]=="inf_multi_sp"){
    importance_asv_clinical[,ii] <- apply(res_rf_asv_clinical[[ii]]$importance, 2, mean)
  }else{
    importance_asv_clinical[-length(var_clinical)-1,ii] <- apply(res_rf_asv_clinical[[ii]]$importance, 2, mean)
  }
}

tab_importance <- importance_asv_clinical
rownames(tab_importance)[1:12] <- c("Age", "Number of children", "Respiratory virus infection", "Maternal HIV infection", "Breastfeeding", "Receipt of antibiotics", 
                                    "Season", "Low birth weight", "Use of solid fuels", "Location of residence", "Sex", "PCV-13")
rownames(tab_importance)[13:44] <- c("M. catarrhalis/nonliquefaciens (1)","Dolosigranulum pigrum (2)","Staphylococcus sp. (3)","C. propinquum/pseudodiphtheriticum (4)",
                                     "S. pneumoniae/pseudopneumoniae (5)","C. accolens/macginleyi (6)","Haemophilus sp. (7)","Haemophilus sp. (8)",
                                     "Streptococcus mitis/oralis (9)","C. tuberculostearicum (10)","Neisseriaceae [G-1] HMT-174 (14)",
                                     "S. thermophilus/vestibularis/salivarius (15)", "Unclassified Enterobacteriaceae (16)", "C. propinquum/pseudodiphtheriticum (18)",
                                     "C. accolens/macginleyi (23)", "Dolosigranulum pigrum (24)", "Micrococcus luteus (26)", "Dolosigranulum pigrum (28)",
                                     "Gemella sp. (30)", "Neisseria sp. (31)", "Veillonella dispar (38)", "Streptococcus lactarius/perioris (40)",
                                     "Alloprevotella sp. (45)", "Veillonella massiliensis (50)", "C. amycolatum (53)", "Moraxella osloensis (67)",
                                     "Streptococcus (69)","C. tuberculostearicum (86)","Acinetobacter sp. (93)","Kocuria sp. (107)",
                                     "Granulicatella elegans (115)","Other ASVs")
colnames(tab_importance) <- pathogen_name[-1]
it_patho_names <- lapply(colnames(tab_importance), function(x) bquote(italic(.(x))))
it_asv_names <- c(as.list(rownames(tab_importance)[(1:12)]), lapply(rownames(tab_importance)[-c(1:12,44)], function(x){bquote(italic(.(x)))}),rownames(tab_importance)[44])
legend_breaks <- seq(0,7, by=1)
fig_4b <- pheatmap(tab_importance, cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 labels_col="", labels_row=as.expression(it_asv_names), breaks = legend_breaks, 
                 angle_row = 45, color=viridis(n = 7, alpha = 1, begin = 0, end = 1, option = "viridis"))

save_pheatmap_png <- function(x, filename, width=2200, height=2000, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(fig_4b, "R_Plots/Figure_4/Figure_4b.png")

fig_4b <- readImage("R_Plots/Figure_4/Figure_4b.png")
fig_4b <- ggdraw() + draw_image(fig_4b, scale=1)

# Figure 4c: barplots of ASVs associated with pathobiont acquisition

theme_lower <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size=7.5, hjust = 0.5),
                     axis.text.x = element_text(size=7, color="black"), axis.title.x = element_text(size=7.5, color="black"),
                     legend.position = "none", axis.line.x = element_line(linewidth = 0.1, linetype = "solid", colour = "black"),
                     panel.background = element_blank(), panel.grid.major.y = element_blank(), panel.border = element_blank(),
                     plot.margin = unit(c(0.3, 0, 0, 0), "cm"))

theme_upper <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size=7.5, hjust = 0.5),
                     axis.text.x = element_text(size=7, color="black"), axis.title.x = element_blank(),
                     legend.position = "none", axis.line.x = element_line(linewidth = 0.1, linetype = "solid", colour = "black"),
                     panel.background = element_blank(), panel.grid.major.y = element_blank(), panel.border = element_blank(),
                     plot.margin = unit(c(0.3, 0, 0.2, 0), "cm"))

log_reg_hi <- read.xlsx('glmer_ASV_acquisition_adj_month.xlsx', sheetN="inf_multi_hi")
log_reg_hi <- subset(log_reg_hi, qvalue<0.20)
names(log_reg_hi)[names(log_reg_hi) == "NA."] <- "ASV"
names(log_reg_hi)[names(log_reg_hi) == "Estimate"] <- "coef"
log_reg_hi <- arrange(log_reg_hi, coef)
log_reg_hi$coef_dir[log_reg_hi$coef<0] <- "Neg"
log_reg_hi$coef_dir[log_reg_hi$coef>0] <- "Pos"
log_reg_hi <- subset(log_reg_hi, ASV!="ASV7" & ASV!="ASV8") # exclude these ASVs as they potentially contain the pathobiont of interest
log_reg_hi$ASV[log_reg_hi$ASV=="ASV9"] <- "Streptococcus mitis/oralis (9)"
log_reg_hi$ASV[log_reg_hi$ASV=="ASV15"] <- "S. thermophilus/vestibularis/salivarius (15)"
log_reg_hi$ASV[log_reg_hi$ASV=="ASV1"] <- "M. catarrhalis/nonliquefaciens (1)"
log_reg_hi$ASV[log_reg_hi$ASV=="ASV5"] <- "S. pneumoniae/pseudopneumoniae (5)"

title_hi <- expression(paste(italic("H. influenzae"), " acquisition"))
hi_plot <- log_reg_hi %>% ggplot(aes(x = coef, y = reorder(ASV, -coef), label = ASV, fill = coef_dir)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("darkslateblue","seagreen")) +
  geom_text(data = log_reg_hi[log_reg_hi$coef<0,], aes(x = 0, hjust = 0), position=position_nudge(x=0.005),  
            size = 2.4, color = "black", lineheight=1, fontface="italic") +
  geom_text(data = log_reg_hi[log_reg_hi$coef>0,], aes(x = 0, hjust = 1), position=position_nudge(x=-0.005), 
            size = 2.4, color = "black", lineheight=1, fontface="italic") +
  theme_upper + labs(x = "Logistic regression coefficient", y = "", fill = "") + scale_x_continuous(breaks=c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15), limits=c(-0.15, 0.15)) +
  ggtitle(title_hi)
hi_plot

log_reg_mc <- read.xlsx('glmer_ASV_acquisition_adj_month.xlsx', sheetN="inf_multi_mc")
log_reg_mc <- subset(log_reg_mc, qvalue<0.20)
names(log_reg_mc)[names(log_reg_mc) == "NA."] <- "ASV"
names(log_reg_mc)[names(log_reg_mc) == "Estimate"] <- "coef"
log_reg_mc <- arrange(log_reg_mc, coef)
log_reg_mc$coef_dir[log_reg_mc$coef<0] <- "Neg"
log_reg_mc$coef_dir[log_reg_mc$coef>0] <- "Pos"
log_reg_mc$ASV[log_reg_mc$ASV=="ASV6"] <- "C. accolens/macginleyi (6)"
log_reg_mc$ASV[log_reg_mc$ASV=="ASV9"] <- "Streptococcus mitis/oralis (9)"
log_reg_mc$ASV[log_reg_mc$ASV=="ASV2"] <- "Dolosigranulum pigrum (2)"
log_reg_mc$ASV[log_reg_mc$ASV=="ASV4"] <- "C. propinquum/pseudodiphtheriticum (4)"
log_reg_mc$ASV[log_reg_mc$ASV=="ASV8"] <- "Haemophilus sp. (8)"

title_mc <- expression(paste(italic("M. catarrhalis"), " acquisition"))
mc_plot <- log_reg_mc %>% ggplot(aes(x = coef, y = reorder(ASV, -coef), label = ASV, fill = coef_dir)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("darkslateblue","seagreen")) +
  geom_text(data = log_reg_mc[log_reg_mc$coef<0,], aes(x = 0, hjust = 0), position=position_nudge(x=0.005),  
            size = 2.4, color = "black", lineheight=1, fontface="italic") +
  geom_text(data = log_reg_mc[log_reg_mc$coef>0,], aes(x = 0, hjust = 1), position=position_nudge(x=-0.005), 
            size = 2.4, color = "black", lineheight=1, fontface="italic") +
  theme_lower + labs(x = "Logistic regression coefficient", y = "", fill = "") + scale_x_continuous(breaks=c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15), limits=c(-0.15, 0.15)) +
  ggtitle(title_mc)
mc_plot

log_reg_sa <- read.xlsx('glmer_ASV_acquisition_adj_month.xlsx', sheetN="inf_multi_sa")
log_reg_sa <- subset(log_reg_sa, qvalue<0.20)
names(log_reg_sa)[names(log_reg_sa) == "NA."] <- "ASV"
names(log_reg_sa)[names(log_reg_sa) == "Estimate"] <- "coef"
log_reg_sa <- arrange(log_reg_sa, coef)
log_reg_sa$coef_dir[log_reg_sa$coef<0] <- "Neg"
log_reg_sa$coef_dir[log_reg_sa$coef>0] <- "Pos"
log_reg_sa <- subset(log_reg_sa, ASV!="ASV3") # exclude these ASVs as they potentially contain the pathobiont of interest
log_reg_sa$ASV[log_reg_sa$ASV=="ASV2"] <- "Dolosigranulum pigrum (2)"
log_reg_sa$ASV[log_reg_sa$ASV=="ASV4"] <- "C. propinquum/pseudodiphtheriticum (4)"
log_reg_sa$ASV[log_reg_sa$ASV=="ASV16"] <- "Unclassified Enterobacteriaceae (16)"

title_sa <- expression(paste(italic("S. aureus"), " acquisition"))
sa_plot <- log_reg_sa %>% ggplot(aes(x = coef, y = reorder(ASV, -coef), label = ASV, fill = coef_dir)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("darkslateblue","seagreen")) +
  geom_text(data = log_reg_sa[log_reg_sa$coef<0,], aes(x = 0, hjust = 0), position=position_nudge(x=0.005),  
            size = 2.4, color = "black", lineheight=1, fontface="italic") +
  geom_text(data = log_reg_sa[log_reg_sa$coef>0,], aes(x = 0, hjust = 1), position=position_nudge(x=-0.005), 
            size = 2.4, color = "black", lineheight=1, fontface="italic") +
  theme_upper + labs(x = "Logistic regression coefficient", y = "", fill = "") + scale_x_continuous(breaks=c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15), limits=c(-0.15, 0.15)) +
  ggtitle(title_sa)
sa_plot

log_reg_sp <- read.xlsx('glmer_ASV_acquisition_adj_month.xlsx', sheetN="inf_multi_sp")
log_reg_sp <- subset(log_reg_sp, qvalue<0.20)
names(log_reg_sp)[names(log_reg_sp) == "NA."] <- "ASV"
names(log_reg_sp)[names(log_reg_sp) == "Estimate"] <- "coef"
log_reg_sp <- arrange(log_reg_sp, coef)
log_reg_sp$coef_dir[log_reg_sp$coef<0] <- "Neg"
log_reg_sp$coef_dir[log_reg_sp$coef>0] <- "Pos"
log_reg_sp <- subset(log_reg_sp, ASV!="ASV5") # exclude these ASVs as they potentially contain the pathobiont of interest
log_reg_sp$ASV[log_reg_sp$ASV=="ASV50"] <- "Veillonella massiliensis (50)"
log_reg_sp$ASV[log_reg_sp$ASV=="ASV53"] <- "C. amycolatum (53)"
log_reg_sp$ASV[log_reg_sp$ASV=="ASV9"] <- "Streptococcus mitis/oralis (9)"
log_reg_sp$ASV[log_reg_sp$ASV=="ASV3"] <- "Staphylococcus sp. (3)"
log_reg_sp$ASV[log_reg_sp$ASV=="ASV10"] <- "C. tuberculostearicum (10)"
log_reg_sp$ASV[log_reg_sp$ASV=="ASV4"] <- "C. propinquum/pseudodiphtheriticum (4)"
log_reg_sp$ASV[log_reg_sp$ASV=="ASV2"] <- "Dolosigranulum pigrum (2)"
log_reg_sp$ASV[log_reg_sp$ASV=="ASV1"] <- "M. catarrhalis/nonliquefaciens (1)"

title_sp <- expression(paste(italic("S. pneumoniae"), " acquisition"))
sp_plot <- log_reg_sp %>% ggplot(aes(x = coef, y = reorder(ASV, -coef), label = ASV, fill = coef_dir)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("darkslateblue","seagreen")) +
  geom_text(data = log_reg_sp[log_reg_sp$coef<0,], aes(x = 0, hjust = 0), position=position_nudge(x=0.005),  
            size = 2.4, color = "black", lineheight=1, fontface="italic") +
  geom_text(data = log_reg_sp[log_reg_sp$coef>0,], aes(x = 0, hjust = 1), position=position_nudge(x=-0.005), 
            size = 2.4, color = "black", lineheight=1, fontface="italic") +
  theme_lower + labs(x = "Logistic regression coefficient", y = "", fill = "") + scale_x_continuous(breaks=c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15), limits=c(-0.15, 0.15)) +
  ggtitle(title_sp)
sp_plot

fig_4c <- plot_grid(hi_plot, sa_plot, NULL, NULL, mc_plot, sp_plot, labels=NULL, nrow=3, rel_heights = c(0.7,-0.03,1)) 

png(file="R_Plots/Figure_4/Figure_4c.png", width = 6.5, height = 3, units = 'in', res = 1200)
plot(fig_4c)
dev.off()

fig_4ac <- plot_grid(fig_4a, NULL, fig_4c, labels=c("a","","c"), label_size=9, nrow=3, rel_heights = c(0.85,-0.03,1)) 

fig_4 <- plot_grid(fig_4ac, NULL, fig_4b, labels=c("","b",""), label_size=9, nrow=1, ncol=3, rel_widths=c(1,0.02,0.68))
fig_4 <- fig_4 + annotate("text", x=0.632, y=0.982, size=2.6, label=expression(italic("H. influenzae")))
fig_4 <- fig_4 + annotate("text", x=0.688, y=0.982, size=2.6, label=expression(italic("M. catarrhalis")))
fig_4 <- fig_4 + annotate("text", x=0.744, y=0.982, size=2.6, label=expression(italic("S. aureus")))
fig_4 <- fig_4 + annotate("text", x=0.801, y=0.980, size=2.6, label=expression(italic("S. pneumoniae")))

png(file="R_Plots/Figure_4.png", width = 13.2, height = 5, units = 'in', res = 1200)
suppressWarnings(plot(fig_4))
dev.off()