#! /usr/bin/env Rscript

library(data.table)
library(optparse)

## Parameter setting
args_list = list(
  make_option("--method", type="character", default=NULL,
              help="INPUT: method", metavar="character"),
  make_option("--pheno", type="character", default=NULL,
              help="INPUT: method", metavar="character"),
  make_option("--cross", type="character", default=NULL,
              help="INPUT: cross number", metavar="character"), 
  make_option("--chr", type="character", default=NULL,
              help="INPUT: chr", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt = list (method = "ldpredinf", pheno = 1, cross = 1, chr = "22")

chr <- as.numeric(opt$chr)
eff_PT_str <- paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/PT/summary_cross", opt$cross, "_best_chr")
eff_SBSLMM_str <- paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/slmm/esteff_summary_cross", opt$cross, "_chr")
eff_SBLUP_str <- paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/sblup/esteff_cross", opt$cross, "_chr")
eff_ldpredinf_str <- paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/ldpred/esteff_cross", opt$cross, "_ldr200_inf_chr")
eff_ldpredpbest_str <- paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/ldpred/esteff_cross", opt$cross, "_ldr200_pbest_chr")
# eff_SBLUP_str <- paste0("/net/mulan/disk2/yasheng/pheno", opt$pheno, "/sblup_ukb/esteff_cross", opt$cross, "_chr")
# eff_ldpredinf_str <- paste0("/net/mulan/disk2/yasheng/pheno", opt$pheno, "/ldpred_ukb/esteff_cross", opt$cross, "_ldr200_inf_chr")
# eff_ldpredpbest_str <- paste0("/net/mulan/disk2/yasheng/pheno", opt$pheno, "/ldpred_ukb/esteff_cross", opt$cross, "_ldr200_pbest_chr")
LD_path <- "/net/mulan/disk2/yasheng/sample500/out_sample/LDmat/"
chr_block_str <- "/net/mulan/home/yasheng/Biobank/data/LDblock/chr"
### out summary statistics
summ_out <- data.frame(fread(paste0("/net/mulan/disk2/yasheng/out_sample/raw_data/pheno", opt$pheno, "_scale.txt.gz")))
# colnames(summ_out) <- c("MarkerName", "A1", "A2", "beta")
# summ_out$beta_s <- summ_out$beta * sqrt(2*summ_out$af*(1-summ_out$af))
# summ_out <- summ_out[, c(1, 2, 3, 6)]

### in summary statistics
summ_in <- data.frame(fread(paste0("/net/mulan/disk2/yasheng/pheno", opt$pheno, "/summ/summary_cross", opt$cross, "_chr", chr, ".assoc.txt")))
summ_in <- summ_in[, c(2, 8)]

### P+T
if (opt$method == "PT"){
  eff_PT <- data.frame(fread(paste0(eff_PT_str, chr, ".txt")))
  chr_block <- read.table(paste0(chr_block_str, chr, ".bed"))
  beta_inter <- merge(eff_PT, summ_out, by.x = "V1", by.y = "MarkerName")
  beta_inter <- merge(beta_inter, summ_in, by.x = "V1", by.y = "V2")
  beta_inter$V3 <- beta_inter$V3 * sqrt(2*beta_inter$V8*(1-beta_inter$V8))
  beta_inter[, 6] <- ifelse(beta_inter[,2] == beta_inter[, 4], beta_inter[, 6], -beta_inter[, 6])
  nume <- vector("numeric", dim(chr_block)[1])
  deno <- vector("numeric", dim(chr_block)[1])
  for(b in 1: dim(chr_block)[1]){
    snp_load <- try(load(paste0(LD_path, "snplist_chr", chr, "_block", b, ".RData")), silent = T)
    if(inherits(snp_load, "try-error") == T){
      nume[b] <- NA
      deno[b] <- NA
    }else{
      snp_inter <- Reduce(intersect, list(snp_sub, eff_PT[, 1], summ_out[, 1]))
      if (length(snp_inter) == 0){
        nume[b] <- NA
        deno[b] <- NA
      }else{
        ld_load <- try(load(paste0(LD_path, "LD_chr", chr, "_block", b, ".RData")), silent = T)
        beta_inter_sub <- beta_inter[match(snp_inter, beta_inter[, 1]), ]
        beta_inter_sub <- beta_inter_sub[!is.na(beta_inter_sub$beta_s), ]
        ld_mat_inter <- ld_mat[which(snp_sub %in% beta_inter_sub[, 1]), which(snp_sub %in% beta_inter_sub[, 1])]
        nume[b] <- (t(beta_inter_sub$V3) %*% beta_inter_sub$beta_s)[1, 1]
        deno[b] <- (t(beta_inter_sub$V3) %*% ld_mat_inter %*% beta_inter_sub$V3)[1, 1]
      }
    }
  }
  cat("chr: ", chr, "\n")
  write.table(cbind(nume, deno),
              file = paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/r2/PT_cross", opt$cross, "_chr", chr, ".txt"),
              row.names = F, col.names = F, quote = F)
}


### SBSLMM
if (opt$method == "DBSLMM"){
  eff_DBSLMM <- data.frame(fread(paste0(eff_SBSLMM_str, chr, ".assoc.sbslmm.txt")))
  chr_block <- read.table(paste0(chr_block_str, chr, ".bed"))
  beta_inter <- merge(eff_DBSLMM, summ_out, by.x = "V1", by.y = "MarkerName")
  beta_inter[, 8] <- ifelse(beta_inter[,2] == beta_inter[, 7], -beta_inter[, 8], beta_inter[, 8])
  nume <- vector("numeric", dim(chr_block)[1])
  deno <- vector("numeric", dim(chr_block)[1])
  for(b in 1: dim(chr_block)[1]){
    snp_load <- try(load(paste0(LD_path, "snplist_chr", chr, "_block", b, ".RData")), silent = T)
    if(inherits(snp_load, "try-error") == T){
      nume[b] <- NA
      deno[b] <- NA
    }else{
      snp_inter <- Reduce(intersect, list(snp_sub, eff_DBSLMM[, 1], summ_out[, 1]))
      if (length(snp_inter) == 0){
        nume[b] <- NA
        deno[b] <- NA
      }else{
        ld_load <- try(load(paste0(LD_path, "LD_chr", chr, "_block", b, ".RData")), silent = T)
        snp_inter <- Reduce(intersect, list(snp_sub, eff_DBSLMM[, 1], summ_out[, 1]))
        beta_inter_sub <- beta_inter[match(snp_inter, beta_inter[, 1]), ]
        beta_inter_sub <- beta_inter_sub[!is.na(beta_inter_sub$beta_s), ]
        ld_mat_inter <- ld_mat[which(snp_sub %in% beta_inter_sub[, 1]), which(snp_sub %in% beta_inter_sub[, 1])]
        
        nume[b] <- (t(beta_inter_sub$V3) %*% beta_inter_sub$beta_s)[1, 1]
        deno[b] <- (t(beta_inter_sub$V3) %*% ld_mat_inter %*% beta_inter_sub$V3)[1, 1]
      }
    }
  }

  write.table(cbind(nume, deno),
              file = paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/r2/DBSLMM_cross", opt$cross, "_chr", chr, ".txt"),
              row.names = F, col.names = F, quote = F)
}

### SBLUP
if (opt$method == "SBLUP"){

  eff_SBLUP <- data.frame(fread(paste0(eff_SBLUP_str, chr, ".sblup.cojo")))
  chr_block <- read.table(paste0(chr_block_str, chr, ".bed"))
  beta_inter <- merge(eff_SBLUP, summ_out, by.x = "V1", by.y = "MarkerName")
  beta_inter <- merge(beta_inter, summ_in, by.x = "V1", by.y = "V2")
  beta_inter$V4 <- beta_inter$V4 * sqrt(2*beta_inter$V8*(1-beta_inter$V8))
  beta_inter[, 7] <- ifelse(beta_inter[,2] == beta_inter[, 6], -beta_inter[, 7], beta_inter[, 7])
  nume <- vector("numeric", dim(chr_block)[1])
  deno <- vector("numeric", dim(chr_block)[1])
  for(b in 1: dim(chr_block)[1]){
    snp_load <- try(load(paste0(LD_path, "snplist_chr", chr, "_block", b, ".RData")), silent = T)
    if(inherits(snp_load, "try-error") == T){
      nume[b] <- NA
      deno[b] <- NA
    }else{
      snp_inter <- Reduce(intersect, list(snp_sub, eff_SBLUP[, 1], summ_out[, 1]))
      if (length(snp_inter) == 0){
        nume[b] <- NA
        deno[b] <- NA
      }else{
        ld_load <- try(load(paste0(LD_path, "LD_chr", chr, "_block", b, ".RData")), silent = T)
        beta_inter_sub <- beta_inter[match(snp_inter, beta_inter[, 1]), ]
        beta_inter_sub <- beta_inter_sub[!is.na(beta_inter_sub$beta_s), ]
        ld_mat_inter <- ld_mat[which(snp_sub %in% beta_inter_sub[, 1]), which(snp_sub %in% beta_inter_sub[, 1])]
        nume[b] <- (t(beta_inter_sub$V4) %*% beta_inter_sub$beta_s)[1, 1]
        deno[b] <- (t(beta_inter_sub$V4) %*% ld_mat_inter %*% beta_inter_sub$V4)[1, 1]
      }
    }
  }
  
  write.table(cbind(nume, deno),
              file = paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/r2/SBLUP_cross", opt$cross, "_chr", chr, ".txt"),
              row.names = F, col.names = F, quote = F)
}

### lassosum
if (opt$method == "lassosum"){
  load(paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/lassosum/esteff_cross", opt$cross, ".RData"))
  eff_tot_lassosum <- data.frame(chr = out_best$sumstats$chr, pos = out_best$sumstats$pos,
                                 beta = out_best$beta[[1]])
  eff_tot_lassosum <- eff_tot_lassosum[-which(eff_tot_lassosum[, 3] == 0), ]
  eff_lassosum <- eff_tot_lassosum[eff_tot_lassosum$chr == chr, ]
  bim_file <- data.frame(fread(paste0("/net/mulan/disk2/yasheng/plink_file/subsample/frq/chr", chr, ".bim")))
  eff_lassosum <- merge(eff_lassosum, bim_file, by.x = "pos", by.y = "V4")
  eff_lassosum <- data.frame(V1 = eff_lassosum$V2, V2 = eff_lassosum$V5, V3 = eff_lassosum$beta)
  chr_block <- read.table(paste0(chr_block_str, chr, ".bed"))
  beta_inter <- merge(eff_lassosum, summ_out, by.x = "V1", by.y = "MarkerName")
  beta_inter <- merge(beta_inter, summ_in, by.x = "V1", by.y = "V2")
  beta_inter$V3 <- beta_inter$V3 * sqrt(2*beta_inter$V8*(1-beta_inter$V8))
  beta_inter[, 6] <- ifelse(beta_inter[,2] == beta_inter[, 4], beta_inter[, 6], -beta_inter[, 6])
  nume <- vector("numeric", dim(chr_block)[1])
  deno <- vector("numeric", dim(chr_block)[1])
  for(b in 1: dim(chr_block)[1]){
    snp_load <- try(load(paste0(LD_path, "snplist_chr", chr, "_block", b, ".RData")), silent = T)
    if(inherits(snp_load, "try-error") == T){
      nume[b] <- NA
      deno[b] <- NA
    }else{
      snp_inter <- Reduce(intersect, list(snp_sub, eff_lassosum[, 1], summ_out[, 1]))
      if (length(snp_inter) == 0){
        nume[b] <- NA
        deno[b] <- NA
      }else{
        ld_load <- try(load(paste0(LD_path, "LD_chr", chr, "_block", b, ".RData")), silent = T)
        beta_inter_sub <- beta_inter[match(snp_inter, beta_inter[, 1]), ]
        beta_inter_sub <- beta_inter_sub[!is.na(beta_inter_sub$beta_s), ]
        ld_mat_inter <- ld_mat[which(snp_sub %in% beta_inter_sub[, 1]), which(snp_sub %in% beta_inter_sub[, 1])]
        nume[b] <- (t(beta_inter_sub$V3) %*% beta_inter_sub$beta_s)[1, 1]
        deno[b] <- (t(beta_inter_sub$V3) %*% ld_mat_inter %*% beta_inter_sub$V3)[1, 1]
      }
    }
  }
  cat("chr: ", chr, "\n")
  write.table(cbind(nume, deno),
              file = paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/r2/lassosum_cross", opt$cross, "_chr", chr, ".txt"),
              row.names = F, col.names = F, quote = F)
}

### LDpred inf
if (opt$method == "ldpredinf"){

  eff_ldpredinf <- data.frame(fread(paste0(eff_ldpredinf_str, chr, ".txt")))
  chr_block <- read.table(paste0(chr_block_str, chr, ".bed"))
  beta_inter <- merge(eff_ldpredinf, summ_out, by.x = "V1", by.y = "MarkerName")
  beta_inter <- merge(beta_inter, summ_in, by.x = "V1", by.y = "V2")
  beta_inter$V3 <- beta_inter$V3 * sqrt(2*beta_inter$V8*(1-beta_inter$V8))
  beta_inter[, 6] <- ifelse(beta_inter[,2] == beta_inter[, 5], -beta_inter[, 6], beta_inter[, 6])
  nume <- vector("numeric", dim(chr_block)[1])
  deno <- vector("numeric", dim(chr_block)[1])
  for(b in 1: dim(chr_block)[1]){
    snp_load <- try(load(paste0(LD_path, "snplist_chr", chr, "_block", b, ".RData")), silent = T)
    if(inherits(snp_load, "try-error") == T){
      nume[b] <- NA
      deno[b] <- NA
    }else{
      snp_inter <- Reduce(intersect, list(snp_sub, eff_ldpredinf[, 1], summ_out[, 1]))
      if (length(snp_inter) == 0){
        nume[b] <- NA
        deno[b] <- NA
      }else{
          ld_load <- try(load(paste0(LD_path, "LD_chr", chr, "_block", b, ".RData")), silent = T)
          beta_inter_sub <- beta_inter[match(snp_inter, beta_inter[, 1]), ]
          beta_inter_sub <- beta_inter_sub[!is.na(beta_inter_sub$beta_s), ]
          ld_mat_inter <- ld_mat[which(snp_sub %in% beta_inter_sub[, 1]), which(snp_sub %in% beta_inter_sub[, 1])]
          nume[b] <- (t(beta_inter_sub$V3) %*% beta_inter_sub$beta_s)[1, 1]
          deno[b] <- (t(beta_inter_sub$V3) %*% ld_mat_inter %*% beta_inter_sub$V3)[1, 1]
      }
    }
  }
  cat("chr: ", chr, "\n")
  write.table(cbind(nume, deno),
              file = paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/r2/ldpredinf_cross", opt$cross, "_chr", chr, ".txt"),
              row.names = F, col.names = F, quote = F)
}

### LDpred pbest
if (opt$method == "ldpredpbest"){

  eff_ldpredinf <- data.frame(fread(paste0(eff_ldpredpbest_str, chr, ".txt")))
  chr_block <- read.table(paste0(chr_block_str, chr, ".bed"))
  beta_inter <- merge(eff_ldpredinf, summ_out, by.x = "V1", by.y = "MarkerName")
  beta_inter <- merge(beta_inter, summ_in, by.x = "V1", by.y = "V2")
  beta_inter$V3 <- beta_inter$V3 * sqrt(2*beta_inter$V8*(1-beta_inter$V8))
  beta_inter[, 6] <- ifelse(beta_inter[,2] == beta_inter[, 5], -beta_inter[, 6], beta_inter[, 6])
  nume <- vector("numeric", dim(chr_block)[1])
  deno <- vector("numeric", dim(chr_block)[1])
  for(b in 1: dim(chr_block)[1]){
    snp_load <- try(load(paste0(LD_path, "snplist_chr", chr, "_block", b, ".RData")), silent = T)
    if(inherits(snp_load, "try-error") == T){
      nume[b] <- NA
      deno[b] <- NA
    }else{
      snp_inter <- Reduce(intersect, list(snp_sub, eff_ldpredinf[, 1], summ_out[, 1]))
      if (length(snp_inter) == 0){
        nume[b] <- NA
        deno[b] <- NA
      }else{
        ld_load <- try(load(paste0(LD_path, "LD_chr", chr, "_block", b, ".RData")), silent = T)
        beta_inter_sub <- beta_inter[match(snp_inter, beta_inter[, 1]), ]
        beta_inter_sub <- beta_inter_sub[!is.na(beta_inter_sub$beta_s), ]
        ld_mat_inter <- ld_mat[which(snp_sub %in% beta_inter_sub[, 1]), which(snp_sub %in% beta_inter_sub[, 1])]
        nume[b] <- (t(beta_inter_sub$V3) %*% beta_inter_sub$beta_s)[1, 1]
        deno[b] <- (t(beta_inter_sub$V3) %*% ld_mat_inter %*% beta_inter_sub$V3)[1, 1]
      }
    }
  }
  cat("chr: ", chr, "\n")

  write.table(cbind(nume, deno),
              file = paste0("/net/mulan/disk2/yasheng/out_sample/pheno", opt$pheno, "/r2/ldpredpbest_cross", opt$cross, "_chr", chr, ".txt"),
              row.names = F, col.names = F, quote = F)
}
