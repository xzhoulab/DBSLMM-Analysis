rm(list=ls())
library(data.table)
library(plyr)
snp_inter <- data.frame(fread("/net/mulan/disk2/yasheng/out_sample/raw_data/snp_UKB_inter.txt", header = T))[, 1]
pheno=10
for (cross in 1: 5){
  for (chr in 1: 22){
    summ_chr <- data.frame(fread(paste0("/net/mulan/disk2/yasheng/pheno", pheno, 
                                        "/output/summary_cross", cross, "_chr", chr, ".assoc.txt")))
    summ_chr_sub <- summ_chr[summ_chr[, 2] %in% snp_inter, ]
    write.table(summ_chr_sub, file = paste0("/net/mulan/disk2/yasheng/out_sample/pheno", 
                                            pheno, "/summ/summary_cross", cross, "_chr", chr, ".assoc.txt"),
                col.names = F, row.names = F, quote = F, sep = "\t")   
  }
}

