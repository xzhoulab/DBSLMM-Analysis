#!/bin/bash
rm(list=ls())
library(data.table)
library(doParallel)
library(foreach)
library(optparse)

## Parameter setting
args_list <- list(
  make_option("--chr", type = "character", default = NULL,
              help = "INPUT: chromosome", metavar = "character"))
EAS_ref_str <- "/net/mulan/disk2/yasheng/EAS/genotype/chr"
emeraLD_str <- "/net/mulan/home/yasheng/Biobank/program/emeraLD/bin/emeraLD"
chr_block_str <- "/net/mulan/home/yasheng/Biobank/data/LDblock_EAS/chr"
path <- "/net/mulan/disk2/yasheng/EAS/LDmat/"

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)
chr <- opt$chr

# for (chr in 1: 22){
#   plink_sel_cmd <- paste0("plink-1.9 --bfile ", EAS_ref_str, chr, " --recode vcf --out ", EAS_ref_str, chr)
#   system(plink_sel_cmd)
#   system(paste0("bgzip ", EAS_ref_str, chr, ".vcf"))
#   system(paste0("rm ", EAS_ref_str, chr, ".log"))
#   system(paste0("tabix -p vcf ", EAS_ref_str, chr, ".vcf.gz"))
#   chr_block <- read.table(paste0(chr_block_str, chr, ".bed"), header = T)
#   chr_block[, 1] <- sub("chr", "", chr_block[, 1])
#   block_string <- paste0(chr_block[, 1], ":", chr_block[, 2], "-", chr_block[, 3])
#   cl <- makeCluster(4)
#   registerDoParallel(cl)
#   LD_par <-  foreach(b = c(1: length(block_string)), .combine = c) %dopar% {
#     emeraLD_str <- "/net/mulan/home/yasheng/Biobank/program/emeraLD/bin/emeraLD"
#     block_string <- paste0(chr_block[, 1], ":", chr_block[, 2], "-", chr_block[, 3])
#     emeraLD_cmd <- paste0(emeraLD_str, " -i ", EAS_ref_str, chr, ".vcf.gz --phased --region ", block_string[b],
#                           " --stdout | bgzip -c > ", path, "/LD_chr", chr, "_block", b, ".txt.gz")
#     system(emeraLD_cmd)
#   }
# }

snp_info <- data.frame(fread(paste0(EAS_ref_str, chr, ".bim")))
chr_block <- read.table(paste0(chr_block_str, chr, ".bed"))
cl <- makeCluster(6)
registerDoParallel(cl)
LD_par <-  foreach(b = c(1: dim(chr_block)[1]), .combine = c) %dopar% {
  require(data.table)
  ld_mat_str <- paste0(path, "LD_chr", chr, "_block", b, ".txt.gz")
  ld_mat_long <- data.frame(fread(ld_mat_str))
  pos_uni <- unique(c(ld_mat_long[, 2], ld_mat_long[, 3]))
  ld_mat <- matrix(0, length(pos_uni), length(pos_uni))

  for(i in 1: c(length(pos_uni)-1)){
    r_tmp <- ld_mat_long[ld_mat_long[, 2] == pos_uni[i], ]
    pos_res <- data.frame(POS2 = pos_uni[c((i+1):length(pos_uni))])
    r_tmp <- merge(pos_res, r_tmp, by = "POS2", all.x = T)
    ld_mat[(i), c((i+1): dim(ld_mat)[2])] <- ifelse(is.na(r_tmp[, 4]), 0, r_tmp[, 4])
  }
  ld_mat <- ld_mat + t(ld_mat)
  diag(ld_mat) <- 1
  snp_sub <- snp_info[match(pos_uni, snp_info[, 4]), 2]
  save(ld_mat, file = paste0(path, "LD_chr", chr, "_block", b, ".RData"))
  save(snp_sub, file = paste0(path, "snplist_chr", chr, "_block", b, ".RData"))
}