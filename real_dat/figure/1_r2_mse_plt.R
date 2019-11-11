rm(list = ls())
library(plyr)
library(data.table)
load("/net/mulan/home/yasheng/Biobank/code/summ_res/r2_list_c.RData")
load("/net/mulan/home/yasheng/Biobank/code/summ_res/mse_list_c.RData")

col_dat <- read.table("/net/mulan/home/yasheng/Biobank/code/figures/col.txt")
method_ord <- c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")
col_dat <- col_dat[match(method_ord, col_dat[, 2]), ]

pheno <- factor(rep(c("SH", "BMD", "BMR", "BMI", "WHR", "AM", "FVC", "FFR", 
                      "SBP", "PLT", "RBC", "EOS", "WBC", "RDW", "NS", "YE"), each = 30),
                levels = c("SH", "PLT", "BMD", "BMR", "BMI", "RBC", "AM", "RDW",
                           "EOS", "WBC", "FVC", "FFR", "WHR", "NS", "SBP", "YE"))
pheno_uni <- c("SH", "BMD", "BMR", "BMI", "WHR", "AM", "FVC", "FFR", 
               "SBP", "PLT", "RBC", "EOS", "WBC", "RDW", "NS", "YE")
## ldsc
ldsc_tot <- vector()
for (chr in 1: 22){
  ldsc_chr <- data.frame(fread(paste0("/net/mulan/disk2/yasheng/sample500/ldsc/", chr, ".l2.ldscore.gz")))[, 4]
  ldsc_tot <- c(ldsc_tot, ldsc_chr)
}
ldsc_mean <- mean(ldsc_tot)

## herit
E_r2 <- vector()
min_ref <- vector("numeric", 16)
max_ref <- vector("numeric", 16)
for (p in 1:16){
  E_r2_p <- vector("numeric", 5)
  for (cross in 1:5){
    
    herit_str <- paste0("/net/mulan/disk2/yasheng/pheno", p, "/heritability/h2_cross", cross, "_ukb.log")
    herit_file <- read.delim(herit_str)[, 1]
    herit_dat <- as.character(herit_file)[24]
    herit <- as.numeric(strsplit(strsplit(herit_dat, ": ")[[1]][2], " \\(")[[1]][1])
    M_dat <- as.character(herit_file)[22]
    M <- as.numeric(strsplit(strsplit(M_dat, ", ")[[1]][2], " ")[[1]][1]) / ldsc_mean
    N_dat <- read.table(paste0("/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_", p, "_n_train.txt"))
    N <- N_dat[cross, 1]
    E_r2_p[cross] <- herit/(1 + M/(N*herit))
  }
  E_r2_pp <- rep(E_r2_p, length(method_ord))
  min_ref[p] <- min(E_r2_p)
  max_ref[p] <- max(E_r2_p)
  E_r2 <- c(E_r2, E_r2_pp)
}

## 
library(ggplot2)
r2_dat <- ldply(r2_list, function (a){
  data.frame(Method = factor(rep(c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T"), each = 5), 
                             level = c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")),
             cross = c(1:5),
             Type = factor(c(rep("SP", 5), rep("P", 5), rep("P", 5), rep("SP", 5), rep("S", 5), rep("S", 5)), 
                           levels = c("SP", "P", "S")),
             r2 = matrix(a, ncol = 1)[, 1])
})
r2_dat <- data.frame(pheno = pheno, r2_dat, E_r2 = E_r2)
min_dat <- data.frame(pheno = pheno_uni, min_line = min_ref)
max_dat <- data.frame(pheno = pheno_uni, max_line = max_ref)


source("/net/mulan/home/yasheng/Biobank/code/figures/r2_est/override.R")
est_r2_plt <- ggplot(r2_dat, aes(Method, r2)) + 
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2)) +  
  facet_wrap_custom(~pheno, scales = "free", ncol = 4, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0.1, 0.4, 0.1), limits = c(0.1, 0.4))), 
    scale_override(2, scale_y_continuous(breaks = seq(0.05, 0.2, 0.05), limits = c(0.05, 0.22))), 
    scale_override(3, scale_y_continuous(breaks = seq(0.05, 0.2, 0.05), limits = c(0.05, 0.2))), 
    scale_override(4, scale_y_continuous(breaks = seq(0.05, 0.2, 0.05), limits = c(0.05, 0.2))), 
    
    scale_override(5, scale_y_continuous(breaks = seq(0.06, 0.15, 0.03), limits = c(0.06, 0.15))), 
    scale_override(6, scale_y_continuous(breaks = seq(0.04, 0.16, 0.04), limits = c(0.04, 0.16))), 
    scale_override(7, scale_y_continuous(breaks = seq(0.03, 0.075, 0.015), limits = c(0.03, 0.075))), 
    scale_override(8, scale_y_continuous(breaks = seq(0.03, 0.15, 0.04), limits = c(0.03, 0.15))), 
    
    scale_override(9, scale_y_continuous(breaks = seq(0.03, 0.12, 0.03), limits = c(0.03, 0.12))), 
    scale_override(10, scale_y_continuous(breaks = seq(0.025, 0.1, 0.025), limits = c(0.025, 0.1))), 
    scale_override(11, scale_y_continuous(breaks = seq(0.03, 0.12, 0.03), limits = c(0.03, 0.12))), 
    scale_override(12, scale_y_continuous(breaks = seq(0.025, 0.1, 0.025), limits = c(0.025, 0.1))), 
    
    scale_override(13, scale_y_continuous(breaks = seq(0.02, 0.08, 0.02), limits = c(0.02, 0.08))), 
    scale_override(14, scale_y_continuous(breaks = seq(0.005, 0.035, 0.01), limits = c(0.005, 0.035))), 
    scale_override(15, scale_y_continuous(breaks = seq(0.015, 0.051, 0.012), limits = c(0.015, 0.051))), 
    scale_override(16, scale_y_continuous(breaks = seq(0.015, 0.033, 0.006), limits = c(0.015, 0.033)))
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
  ylab(expression("Prediction " ~ italic(R^2)))+
  geom_hline(aes(yintercept = min_line), min_dat, color="#EFC000FF", 
             linetype="dashed", size=1.5, na.rm = T) +
  geom_hline(aes(yintercept = max_line), max_dat, color="#EFC000FF", 
             linetype="dashed", size=1.5, na.rm = T)

tiff("/net/mulan/disk2/yasheng/figures/real_dat/est_r2.tiff", width = 2000, height = 1600)
est_r2_plt
dev.off()

# mse data
mse_dat <- ldply(mse_list, function (a){
  data.frame(Method = factor(rep(c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T"), each = 5), 
                             level = c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")),
             cross = c(1:5),
             Type = factor(c(rep("SP", 5), rep("P", 5), rep("P", 5), rep("SP", 5), rep("S", 5), rep("S", 5)), 
                           levels = c("SP", "P", "S")),
             mse = matrix(a, ncol = 1)[, 1])
})
mse_dat <- data.frame(pheno = pheno, mse_dat)

est_mse_plt <- ggplot(mse_dat, aes(Method, mse, fill = Method)) + 
  geom_boxplot(width=0.7, position=position_dodge(0.8)) +  
  facet_wrap_custom(~pheno, scales = "free", ncol = 4, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0, 480, 120), limits = c(0, 480))), 
    scale_override(2, scale_y_continuous(breaks = seq(0, 52, 13), limits = c(0, 53))), 
    scale_override(3, scale_y_continuous(breaks = seq(0, 8, 2), limits = c(0, 8))), 
    scale_override(4, scale_y_continuous(breaks = seq(0, 60, 15), limits = c(0, 60))), 
    
    scale_override(5, scale_y_continuous(breaks = seq(0, 16, 4), limits = c(0, 16))), 
    scale_override(6, scale_y_continuous(breaks = seq(0, 60, 15), limits = c(0, 60))), 
    scale_override(7, scale_y_continuous(breaks = seq(0, 6, 1.5), limits = c(0, 6))), 
    scale_override(8, scale_y_continuous(breaks = seq(0, 600, 150), limits = c(0, 600))), 
    
    scale_override(9, scale_y_continuous(breaks = seq(0, 120, 30), limits = c(0, 120))), 
    scale_override(10, scale_y_continuous(breaks = seq(0, 320, 80), limits = c(0, 325))), 
    scale_override(11, scale_y_continuous(breaks = seq(0, 28, 7), limits = c(0, 28))), 
    scale_override(12, scale_y_continuous(breaks = seq(0, 60, 15), limits = c(0, 60))), 
    
    scale_override(13, scale_y_continuous(breaks = seq(0, 10, 2.5), limits = c(0, 10))), 
    scale_override(14, scale_y_continuous(breaks = seq(0, 28, 7), limits = c(0, 28))), 
    scale_override(15, scale_y_continuous(breaks = seq(0, 12, 3), limits = c(0, 13))), 
    scale_override(16, scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 80)))
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
  ylab("Prediction MSE" )

tiff("/net/mulan/disk2/yasheng/figures/real_dat/est_mse.tiff", width = 2000, height = 1600)
est_mse_plt
dev.off()