#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--pred", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--pheno", type="character", default=NULL,
              help="INPUT: phenotype number", metavar="character"), 
  make_option("--r2", type="character", default=NULL,
              help="OUTPUT: r2 file", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)
# if (opt$pheno == "SH") {
#   pheno_num <- 1
# }
# if (opt$pheno == "BMI") {
#   pheno_num <- 4
# }
# if (opt$pheno == "AM") {
#   pheno_num <- 6
# }
# if (opt$pheno == "PC") {
#   pheno_num <- 10
# }
# if (opt$pheno == "RBC") {
#   pheno_num <- 11
# }
# if (opt$pheno == "NS") {
#   pheno_num <- 15
# }
# if (opt$pheno == "YE") {
#   pheno_num <- 16
# }
pred_pheno <- data.frame(fread(opt$pred, header = T))[, 6]
# pheno <- data.frame(fread("/net/mulan/disk2/yasheng/plink_file/subsample/pheno_sub_b_v2.txt"))[, as.numeric(opt$pheno)]
pheno <- data.frame(fread("/net/mulan/disk2/yasheng/plink_file/subsample/pheno_sub_v2.txt"))[, as.numeric(opt$pheno)]
r2 <- cor(pred_pheno, pheno)^2
write.table(r2, file = opt$r2, col.names = F, row.names = F, quote = F)
