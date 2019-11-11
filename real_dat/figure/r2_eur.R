rm(list = ls())
library(ggplot2)
library(Metrics)
load("/net/mulan/disk2/yasheng/out_sample/summ.RData")

### color
col_dat <- read.table("/net/mulan/home/yasheng/Biobank/code/figures/col.txt")
method_ord <- c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")
col_dat <- col_dat[match(method_ord, col_dat[, 2]), ]

pheno <- c("SH", "BMI", "PLT", "RBC", "EOS", "WBC")
pheno2 <- factor(rep(c("SH", "BMI", "PLT", "RBC", "EOS", "WBC"), each = 30),
                levels = c("SH", "BMI", "PLT", "RBC", "EOS", "WBC"))
cross <- c(1:5)
method_label <- c("SBSLMM", "SBLUP", "ldpredinf", "ldpredpbest", "lassosum", "PT")

## ldsc
ldsc_tot <- vector()
for (chr in 1: 22){
  ldsc_chr <- data.frame(fread(paste0("/net/mulan/disk2/yasheng/sub_snp/ldsc/", chr, ".l2.ldscore.gz")))[, 4]
  ldsc_tot <- c(ldsc_tot, ldsc_chr)
}
ldsc_mean <- mean(ldsc_tot)

## herit
phenoc <- c(1, 4, 10, 11, 12, 13)
E_r2 <- vector()
min_ref <- vector("numeric", 6)
max_ref <- vector("numeric", 6)
for (p in 1: length(phenoc)){
  E_r2_p <- vector("numeric", 5)
  for (cross in 1:5){
    
    herit_str <- paste0("/net/mulan/disk2/yasheng/out_sample/pheno", phenoc[p], "/herit/h2_cross", cross, ".log")
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


library(ggplot2)
r2_dat <- ldply(eff_summ, function (a){
  data.frame(Method = factor(rep(c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T"), each = 5), 
                             level = c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")),
             cross = c(1:5),
             Type = factor(c(rep("SP", 5), rep("P", 5), rep("P", 5), rep("SP", 5), rep("S", 5), rep("S", 5)), 
                           levels = c("SP", "P", "S")),
             r2 = matrix(a, ncol = 1)[, 1])
})
r2_dat <- data.frame(pheno = pheno2, r2_dat, E_r2 = E_r2)
min_dat <- data.frame(pheno = pheno, min_line = min_ref)
max_dat <- data.frame(pheno = pheno, max_line = max_ref)

min(r2_dat[r2_dat[, 1] == "WBC", 5])
max(r2_dat[r2_dat[, 1] == "WBC", 5])

source("/net/mulan/home/yasheng/Biobank/code/figures/r2_est/override.R")
est_r2_plt <- ggplot(r2_dat, aes(Method, r2)) + 
  geom_jitter(aes(color = Method), size = 5, shape = 17, position=position_jitter(0.2)) +  
  facet_wrap_custom(~pheno, scales = "free", ncol = 3, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0.0, 0.3, 0.1), limits = c(0.0, 0.3))), 
    scale_override(2, scale_y_continuous(breaks = seq(0.06, 0.12, 0.02), limits = c(0.06, 0.12))), 
    scale_override(3, scale_y_continuous(breaks = seq(0.0, 0.21, 0.07), limits = c(0.0, 0.21))), 
    scale_override(4, scale_y_continuous(breaks = seq(0.06, 0.12, 0.02), limits = c(0.06, 0.12))), 
    scale_override(5, scale_y_continuous(breaks = seq(0.015, 0.075, 0.02), limits = c(0.015, 0.075))), 
    scale_override(6, scale_y_continuous(breaks = seq(0.03, 0.12, 0.03), limits = c(0.03, 0.12)))
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

tiff("/net/mulan/disk2/yasheng/figures/external/r2_eur.tiff", width = 1800, height = 1000)
est_r2_plt
dev.off()