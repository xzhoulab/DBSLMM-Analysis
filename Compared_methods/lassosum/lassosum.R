#! /usr/bin/env Rscript
library(plyr)
library(data.table)
library(optparse)
library(lassosum, lib.loc = "/net/mulan/home/yasheng/R/x86_64-pc-linux-gnu-library/3.3")
library(Metrics, lib.loc = "/net/mulan/home/yasheng/R/x86_64-pc-linux-gnu-library/3.3")

## Parameter setting
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: summary data prefix", metavar="character"), 
  make_option("--ref", type="character", default=NULL,
              help="INPUT: reference data prefix", metavar="character"), 
  make_option("--valid", type="character", default=NULL,
              help="INPUT: validation data prefix", metavar="character"),
  make_option("--test", type="character", default=NULL,
              help="INPUT: test data prefix", metavar="character"),
  make_option("--res", type="character", default=NULL,
              help="OUTPUT: output mse and r2", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(summ = "/exports/disk2/yasheng/simulation2/SLM/output/summary_block5_her0.5_cross8_dist1_ps0.001.assoc.txt",
#             ref = "/exports/disk2/yasheng/sample500/frq/chr1",
#             valid = "/exports/disk2/yasheng/simulation2/SLM/data/val_block5_her0.5_cross8_dist1_ps0.001",
#             test = "/exports/disk2/yasheng/simulation2/SLM/data/test_block5_her0.5_cross8_dist1_ps0.001",
#             res = "/exports/disk2/yasheng/simulation2/SLM/lassosum/res_block5_her0.5_cross8_dist1_ps0.001.profile")

summ <- data.frame(fread(opt$summ))

## sample size
n <- as.numeric(summ[1, 4]) + as.numeric(summ[1, 5])

## p to cor
pval <- ifelse(summ[, 11] == 0, 1e-100, summ[, 11]) 
cor <- p2cor(p = pval, n = n, sign = summ[, 9])

## block information
LDblocks <- "EUR.hg19"

## estimation
out <- lassosum.pipeline(cor = cor, chr = summ[, 1], pos = summ[, 3], 
                         A1 = summ[, 6], A2 = summ[, 7], 
                         ref.bfile = opt$ref, test.bfile = opt$valid, 
                         LDblocks = LDblocks)
## validation
valid <- validate(out, plot = F)

## testing
out_best <- subset(out, s = valid$best.s, lambda = valid$best.lambda)
valid_test <- validate(out_best, test.bfile = opt$test, plot = F)

## results
mse_test <- mse(valid_test$results.table$pheno, valid_test$results.table$best.pgs)
r2_test <- cor(valid_test$results.table$pheno, valid_test$results.table$best.pgs)^2
save(out_best, file = paste0(opt$res, ".RData"))
write.table(c(r2_test, mse_test), file = paste0(opt$res, ".txt"), row.names = F,
            col.names = F, quote = F)
