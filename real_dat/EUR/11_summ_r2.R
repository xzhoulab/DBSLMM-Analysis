#! /usr/bin/env Rscript

library(data.table)
library(plyr)
###
### paramters
method_ord <- c("DBSLMM", "SBLUP", "ldpredinf", "ldpredpbest", "lassosum", "PT" )
pheno <- c(1, 4, 10, 11, 12, 13)
eff_summ <- rep(list(matrix(NA, ncol = length(method_ord), nrow = 5)), length(pheno))
for (p in 1: length(pheno)){
  for (m in 1: length(method_ord)){
    for (cross in 1:5){
      eff_list <- alply(c(1: 22), 1, function (a) {
                        read.table(paste0("/net/mulan/disk2/yasheng/out_sample/pheno", pheno[p], 
                                "/r2/", method_ord[m], "_cross", cross, "_chr", a, ".txt"))
                        })
      eff_dat <- do.call(rbind, eff_list)
      eff_summ[[p]][cross, m] <- (sum(eff_dat[, 1], na.rm = T) / sqrt(sum(eff_dat[, 2], na.rm = T)))^2
    }
  }
}
save(eff_summ, file = "/net/mulan/disk2/yasheng/out_sample/summ.RData")

load("/net/mulan/disk2/yasheng/out_sample/summ.RData")
library(plyr)
r2_mean <- laply(eff_summ, function(a) colMeans(a, na.rm = T))
r2_best_num_DBSLMM <- apply(r2_mean, 1, function(a) which.max(a) == 1)
cat("DBSLMM best traits: ", sum(r2_best_num_DBSLMM), "\n")

## imporve 
r2_mean_best <- r2_mean[r2_best_num_DBSLMM, ]
r2_mean_best_imp <- aaply(r2_mean_best, 1, function(a) sort(a, decreasing = T)[c(1, 2)])
r2_imp <- r2_mean_best_imp[, 1] / r2_mean_best_imp[, 2] - 1 
max(r2_imp)
min(r2_imp)
mean(r2_imp)