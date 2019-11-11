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
sel_snp_str <- "/net/mulan/disk2/yasheng/out_sample/raw_data/snp_UKB_inter.txt"
UKB_ref_str <- "/net/mulan/disk2/yasheng/sample500/frq/chr"
sel_ref_str <- "/net/mulan/disk2/yasheng/sample500/out_sample/sel_chr"
imp_ref_str <- "/net/mulan/disk2/yasheng/sample500/out_sample/imp_chr"
emeraLD_str <- "/net/mulan/home/yasheng/Biobank/program/emeraLD/bin/emeraLD"
chr_block_str <- "/net/mulan/home/yasheng/Biobank/data/LDblock/chr"
path <- "/net/mulan/disk2/yasheng/sample500/out_sample/"

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)
chr <- opt$chr
# ## select snp
snp1 <- data.frame(fread("/net/mulan/disk2/yasheng/sample500/frq/merge.bim"))
load("/net/mulan/disk2/yasheng/out_sample/raw_data/snp_inter.RData")

write.table(intersect(snp_inter, snp1[, 2]), file = "/net/mulan/disk2/yasheng/out_sample/raw_data/snp_UKB_inter.txt",
            col.names = F, row.names = F, quote = F)

# for (chr in 11: 21){
#   plink_sel_cmd <- paste0("plink-1.9 --bfile ", UKB_ref_str, chr, ".imp --extract ", sel_snp_str,
#                         " --recode vcf --out ", sel_ref_str, chr)
#   system(plink_sel_cmd)
#   system(paste0("bgzip ", sel_ref_str, chr, ".vcf"))
#   system(paste0("rm ", sel_ref_str, chr, ".log"))
#   system(paste0("tabix -p vcf ", sel_ref_str, chr, ".vcf.gz"))
#   chr_block <- read.table(paste0(chr_block_str, chr, ".bed"))
#   chr_block[, 1] <- sub("chr", "", chr_block[, 1])
#   block_string <- paste0(chr_block[, 1], ":", chr_block[, 2], "-", chr_block[, 3])
#   cl <- makeCluster(4)
#   registerDoParallel(cl)
#   LD_par <-  foreach(b = c(1: length(block_string)), .combine = c) %dopar% {
#     emeraLD_str <- "/net/mulan/home/yasheng/Biobank/program/emeraLD/bin/emeraLD"
#     imp_ref_str <- "/net/mulan/disk2/yasheng/sample500/out_sample/imp_chr"
#     block_string <- paste0(chr_block[, 1], ":", chr_block[, 2], "-", chr_block[, 3])
#     emeraLD_cmd <- paste0(emeraLD_str, " -i ", sel_ref_str, chr, ".vcf.gz --phased --region ", block_string[b],
#                           " --stdout | bgzip -c > ", path, "/LDmat/LD_chr", chr, "_block", b, ".txt.gz")
#     system(emeraLD_cmd)
#   }
# }

snp_info <- data.frame(fread(paste0(UKB_ref_str, chr, ".bim")))
chr_block <- read.table(paste0(chr_block_str, chr, ".bed"))
cl <- makeCluster(6)
registerDoParallel(cl)
LD_par <-  foreach(b = c(1: dim(chr_block)[1]), .combine = c) %dopar% {
  require(data.table)
  ld_mat_str <- paste0(path, "/LDmat/LD_chr", chr, "_block", b, ".txt.gz")
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
  save(ld_mat, file = paste0(path, "/LDmat/LD_chr", chr, "_block", b, ".RData"))
  save(snp_sub, file = paste0(path, "/LDmat/snplist_chr", chr, "_block", b, ".RData"))
}