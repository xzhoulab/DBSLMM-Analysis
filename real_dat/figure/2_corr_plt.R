rm(list=ls())
library(plyr)
load("/net/mulan/home/yasheng/Biobank/code/summ_res/r2_list_c.RData")

col_dat <- read.table("/net/mulan/home/yasheng/Biobank/code/figures/col.txt")
method_ord <- c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_pbest", "lassosum", "P+T")
col_dat <- col_dat[match(method_ord, col_dat[, 2]), ]

## herit
herit_summ <- matrix(NA, 5, 16)
for (cross in 1:5){
  for (pheno in 1:16){
    herit_str <- paste0("/net/mulan/disk2/yasheng/pheno", pheno, "/heritability/h2_cross", cross, "_ukb.log")
    herit_dat <- as.character(read.delim(herit_str)[, 1])[24]
    herit_summ[cross, pheno] <- as.numeric(strsplit(strsplit(herit_dat, ": ")[[1]][2], " \\(")[[1]][1])
  }
}

# herit <- vector()
# for (i in 1: 16){
#   herit_tmp <- rep(herit_summ[, i], length(method_ord))
#   herit <- c(herit, herit_tmp)
# }
herit <- rep(colMeans(herit_summ), each = length(method_ord))

## r2
# r2_dat <- ldply(r2_list, function (a){
#   data.frame(Method = rep(c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_pbest", "lassosum", "P+T"), each = 5),
#              cross = c(1:5),
#              r2 = matrix(a, ncol = 1)[, 1])
# })
# r2_dat[, 1] <- factor(r2_dat[, 1],
#                       level = c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_pbest", "lassosum", "P+T"))
r2_dat <- ldply(r2_list, function (a){
  data.frame(Method = c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_pbest", "lassosum", "P+T"),
             r2 = matrix(colMeans(a, na.rm = T), ncol = 1)[, 1])
})
r2_dat <- data.frame(r2_dat, 
                     Pheno  = factor(rep(c("SH", "BMD", "BMR", "BMI", "WHR", "AM", "FVC", "FFR", 
                                           "SBP", "PLT", "RBC", "EOS", "WBC", "RDW", "NS", "YE"), each = 6),
                                     levels = c("SH", "PLT", "BMD", "BMR", "BMI", "RBC", "AM", "RDW",
                                                "EOS", "WBC", "FVC", "FFR", "WHR", "NS", "SBP", "YE")), 
                     Class = rep(c("Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", 
                                   "Anthropometric", "Other", "Other", "Other",  
                                   "Anthropometric", "Blood cell", "Blood cell", "Blood cell", 
                                   "Blood cell", "Blood cell", "Other", "Other"), each = 6))
r2_dat[, 1] <- factor(r2_dat[, 1],
                      level = c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_pbest", "lassosum", "P+T"))
herit_r2_dat <- cbind(r2_dat, herit)

cor_res <- vector("numeric", length(method_ord))
for (i in 1: length(method_ord)){
  tmp_dat <- herit_r2_dat[herit_r2_dat$Method == method_ord[i], ]
  if(sum(is.na(tmp_dat$r2)) != 0){
    tmp_dat <- tmp_dat[-which(is.na(tmp_dat[, 4])), ]
  }
  cor_res[i] <- cor(tmp_dat$r2, tmp_dat$herit, method = "pearson")
  # cor_res[i] <- summary(lm(tmp_dat$herit ~ tmp_dat$r2))$coef[2, 1]
}

dat_text <- data.frame(
  Method = c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_pbest", "lassosum", "P+T")
)
dat_text$label <- sprintf(
  "Cor=%s",
  matrix(round(cor_res, 3), ncol = 1)
)

lm_plt <- ggplot(herit_r2_dat, aes(herit, r2)) + 
  geom_point(aes(colour = Class, shape = Pheno), size = 10) +  
  facet_wrap(~Method, ncol = 3) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 30, face = "bold"), 
        axis.title.x = element_text(size = 60, face = "bold"),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        strip.text = element_text(face = 'bold', size = rel(3)), 
        axis.line.y.left = element_line(size = 1, color = "#E0E0E0"), 
        axis.line.x = element_line(size = 1, color = "#E0E0E0"),
        legend.text = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 40, face = "bold"),
        legend.position = "bottom")+
  scale_shape_manual(values = c(0: 15)) +
  # scale_size_manual(values = c(10, 10, 10)) +
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF", "#868686FF"))+
  geom_text(size = 12, data = dat_text,
            mapping = aes(x = Inf, y = Inf, label = label), hjust = 1.4, vjust = 1.5)+
  xlab(expression("Prediction " ~ italic(R^2)))+
  ylab("Estimation heritability" ) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.1), limits = c(0, 0.4)) +
  scale_x_continuous(breaks = seq(0, 0.4, 0.1), limits = c(0, 0.4)) +
  geom_abline(intercept = 0, slope = 1, color = "#EFC000FF", 
              linetype = "dashed", size = 1.5, na.rm = T)

tiff("/net/mulan/disk2/yasheng/figures/real_dat/corrlationx.tiff", width = 2000, height = 1600)
lm_plt
dev.off()
