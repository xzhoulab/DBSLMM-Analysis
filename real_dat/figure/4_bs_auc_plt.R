rm(list = ls())
library(plyr)
library(data.table)
load("/net/mulan/home/yasheng/Biobank/code/summ_res/bs_list_b.RData")
load("/net/mulan/home/yasheng/Biobank/code/summ_res/auc_list_b.RData")

col_dat <- read.table("/net/mulan/home/yasheng/Biobank/code/figures/col.txt")
method_ord <- c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")
col_dat <- col_dat[match(method_ord, col_dat[, 2]), ]

pheno <- factor(rep(c("TA", "BI", "HT", "EZ", "MP", "CO"), each = 30),
                levels = c("BI", "CO", "TA", "MP", "HT", "EZ"))
pheno_uni <- c("BI", "CO", "TA", "MP", "HT", "EZ")

## 
library(ggplot2)
bs_dat <- ldply(bs_list, function (a){
  data.frame(Method = factor(rep(c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T"), each = 5), 
                             level = c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")),
             cross = c(1:5),
             Type = factor(c(rep("SP", 5), rep("P", 5), rep("P", 5), rep("SP", 5), rep("S", 5), rep("S", 5)), 
                           levels = c("SP", "P", "S")),
             bs = matrix(a, ncol = 1)[, 1])
})
bs_dat <- data.frame(pheno = pheno, bs_dat)
max(bs_dat[bs_dat[, 1] == "EZ", 5], na.rm = T)
min(bs_dat[bs_dat[, 1] == "EZ", 5], na.rm = T)
source("/net/mulan/home/yasheng/Biobank/code/figures/r2_est/override.R")
est_bs_plt <- ggplot(bs_dat, aes(Method, bs)) + 
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2)) +  
  facet_wrap_custom(~pheno, scales = "free", ncol = 3, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0.2, 0.47, 0.09), limits = c(0.2, 0.47))), 
    scale_override(2, scale_y_continuous(breaks = seq(0.2, 0.47, 0.09), limits = c(0.2, 0.47))), 
    scale_override(3, scale_y_continuous(breaks = seq(0.2, 0.8, 0.2), limits = c(0.2, 0.8))), 
    
    scale_override(4, scale_y_continuous(breaks = seq(0.2, 0.47, 0.09), limits = c(0.2, 0.47))), 
    scale_override(5, scale_y_continuous(breaks = seq(0.2, 0.5, 0.1), limits = c(0.2, 0.5))), 
    scale_override(6, scale_y_continuous(breaks = seq(0.2, 0.47, 0.09), limits = c(0.2, 0.47)))
  ))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -2, aes(label=round(..y.., digits=3)), size = 8) +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 40, face = "bold"),
        legend.position = "bottom",
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        strip.text = element_text(face = 'bold', size = rel(3)), 
        axis.line.y.left = element_line(size = 1, color = "#E0E0E0"), 
        axis.line.x = element_line(size = 1, color = "#E0E0E0"))+
  scale_color_manual(values=as.character(col_dat[, 1])) + 
  ylab("Prediction Brier Score")
tiff("/net/mulan/disk2/yasheng/figures/real_dat/est_bs.tiff", width = 1800, height = 1000)
est_bs_plt
dev.off()

# mse data
auc_dat <- ldply(auc_list, function (a){
  data.frame(Method = factor(rep(c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T"), each = 5), 
                             level = c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")),
             cross = c(1:5),
             Type = factor(c(rep("SP", 5), rep("P", 5), rep("P", 5), rep("SP", 5), rep("S", 5), rep("S", 5)), 
                           levels = c("SP", "P", "S")),
             auc = matrix(a, ncol = 1)[, 1])
})
auc_dat <- data.frame(pheno = pheno, auc_dat)

for (p in 1: 6){
  cat(min(auc_dat[auc_dat[, 1] == pheno_uni[p], 5], na.rm = T), "\n")
}
est_auc_plt <- ggplot(auc_dat, aes(Method, auc, fill = Method)) + 
  geom_boxplot(width=0.7, position=position_dodge(0.8)) +  
  facet_wrap_custom(~pheno, scales = "free", ncol = 3, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0.5, 0.74, 0.08), limits = c(0.5, 0.74))), 
    scale_override(2, scale_y_continuous(breaks = seq(0.5, 0.74, 0.08), limits = c(0.5, 0.74))), 
    scale_override(3, scale_y_continuous(breaks = seq(0.5, 0.74, 0.08), limits = c(0.5, 0.74))), 
    
    scale_override(4, scale_y_continuous(breaks = seq(0.5, 0.62, 0.04), limits = c(0.5, 0.62))), 
    scale_override(5, scale_y_continuous(breaks = seq(0.5, 0.62, 0.04), limits = c(0.5, 0.62))), 
    scale_override(6, scale_y_continuous(breaks = seq(0.5, 0.62, 0.04), limits = c(0.5, 0.62)))
  ))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -2, aes(label=round(..y.., digits=3)), size = 8) +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 40, face = "bold"),
        legend.position = "bottom",
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        strip.text = element_text(face = 'bold', size = rel(3)), 
        axis.line.y.left = element_line(size = 1, color = "#E0E0E0"), 
        axis.line.x = element_line(size = 1, color = "#E0E0E0"))+
  scale_fill_manual(values=as.character(col_dat[, 1])) + 
  ylab("Prediction AUC" )

tiff("/net/mulan/disk2/yasheng/figures/real_dat/est_auc.tiff", width = 1800, height = 1000)
est_auc_plt
dev.off()