rm(list = ls())
library(ggplot2)
load("/net/mulan/home/yasheng/Biobank/code/summ_res/r2_group_dat.RData")
load("/net/mulan/home/yasheng/Biobank/code/summ_res/mse_group_dat.RData")


r2_group_dat[, 1] <- as.character(r2_group_dat[, 1])

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

col_dat <- read.table("/net/mulan/home/yasheng/Biobank/code/figures/col.txt")
method_ord <- c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")
col_dat <- col_dat[match(method_ord, col_dat[, 2]), ]
pheno_label <- c("SH", "BMD", "BMR", "BMI", "WHR", "AM", "FVC", "FFR", 
                 "SBP", "PLT", "RBC", "EOS", "WBC", "RDW", "NS", "YE")
pheno_ord <- c("SH", "PLT", "BMD", "BMR", "BMI", "RBC", "AM", "RDW",
               "EOS", "WBC", "FVC", "FFR", "WHR", "NS", "SBP", "YE")
pheno <- c(1:16)
min_r2_vec <- vector("numeric", length(pheno))
max_r2_vec <- vector("numeric", length(pheno))
min_mse_vec <- vector("numeric", length(pheno))
max_mse_vec <- vector("numeric", length(pheno))
for (i in 1: length(pheno)){
  min_r2_vec[i] <- min(r2_group_dat[r2_group_dat$pheno == pheno_ord[i], 5])
  max_r2_vec[i] <- max(r2_group_dat[r2_group_dat$pheno == pheno_ord[i], 5])
  min_mse_vec[i] <- min(mse_group_dat[mse_group_dat$pheno == pheno_ord[i], 5])
  max_mse_vec[i] <- max(mse_group_dat[mse_group_dat$pheno == pheno_ord[i], 5])
}

r2_group_SE_dat <- summarySE(r2_group_dat, measurevar = "r2", 
                             groupvars = c("Method", "group", "pheno"))
r2_group_SE_dat[, 3] <- factor(r2_group_SE_dat[, 3], levels = pheno_ord)
for(p in 1:16){
  temp <- r2_group_SE_dat[r2_group_SE_dat[, 3] == pheno_ord[p], ]
  minx <- min(temp$r2 - 1.96*temp$se)
  maxx <- max(temp$r2 + 1.96*temp$se)
  cat( pheno_ord[p], "min: ", minx, "max: ", maxx, "\n")
}

source("/net/mulan/home/yasheng/Biobank/code/figures/real_dat/override.R")
r2_curve <- ggplot(r2_group_SE_dat, aes(x=group, y=r2, group=Method, color=Method)) + 
  geom_errorbar(aes(ymin=r2-1.96*se, ymax=r2+1.96*se), width=.1, size = 1.5) +
  geom_line(size = 1.5) +
  geom_point(size = 1.5) +
  facet_wrap_custom(~pheno, scales = "free", ncol = 4, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0, 0.09, 0.03), limits = c(0, 0.09))), 
    scale_override(2, scale_y_continuous(breaks = seq(0, 0.0081, 0.0027), limits = c(0, 0.0081))), 
    scale_override(3, scale_y_continuous(breaks = seq(0, 0.006, 0.002), limits = c(0, 0.0063))), 
    scale_override(4, scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(0, 0.03))), 
    
    scale_override(5, scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(0, 0.03))), 
    scale_override(6, scale_y_continuous(breaks = seq(0, 0.024, 0.008), limits = c(0, 0.024))), 
    scale_override(7, scale_y_continuous(breaks = seq(0, 0.27, 0.09), limits = c(0, 0.27))), 
    scale_override(8, scale_y_continuous(breaks = seq(0, 0.024, 0.008), limits = c(0, 0.025))), 
    
    scale_override(9, scale_y_continuous(breaks = seq(0, 0.027, 0.009), limits = c(0, 0.027))), 
    scale_override(10, scale_y_continuous(breaks = seq(0, 0.012, 0.004), limits = c(0, 0.012))), 
    scale_override(11, scale_y_continuous(breaks = seq(0, 0.009, 0.003), limits = c(0, 0.01))), 
    scale_override(12, scale_y_continuous(breaks = seq(0, 0.006, 0.002), limits = c(0, 0.006))), 
    
    scale_override(13, scale_y_continuous(breaks = seq(0, 0.009, 0.003), limits = c(0, 0.009))), 
    scale_override(14, scale_y_continuous(breaks = seq(0, 0.27, 0.09), limits = c(0, 0.27))), 
    scale_override(15, scale_y_continuous(breaks = seq(0, 0.009, 0.003), limits = c(0, 0.009))), 
    scale_override(16, scale_y_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.32)))
  ))+
  scale_color_manual(values=as.character(col_dat[, 1]))+
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 40, face = "bold"),
        legend.position = "bottom",
        title = element_text(size = 20, face = "bold"), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+
  ylab(expression("Prediction " ~ italic(R^2)))

tiff("/net/mulan/disk2/yasheng/figures/real_dat/r2_curvex.tiff", height = 1600, width = 2000)
r2_curve
dev.off()
