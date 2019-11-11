rm(list=ls())
library(data.table)
library(plyr)
setwd("/net/mulan/disk2/yasheng/out_bbj/raw_data/")

### load datasets of 6 traits
SH_dat <- data.frame(fread("SH_raw.txt.gz"))
BMI_dat <- data.frame(fread("BMI_raw.txt.gz"))
PLT_dat <- data.frame(fread("PLT_raw.txt.gz"))
RBC_dat <- data.frame(fread("RBC_raw.txt.gz"))
EOS_dat <- data.frame(fread("EOS_raw.txt.gz"))
WBC_dat <- data.frame(fread("WBC_raw.txt.gz"))
snp_rs <- PLT_dat[, c(1:3)]

###
snp_pos <- vector()
for (chr in 1: 22){
  
  snp_rs_chr <- snp_rs[snp_rs[, 2] == chr, ]
  SH_dat_chr <- SH_dat[SH_dat$CHR == chr, ]
  snp_ex_chr <- merge(SH_dat_chr, snp_rs_chr, by.x = "POS", by.y = "POS")
  snp_pos_chr <- cbind(snp_ex_chr$Variants, snp_ex_chr$SNP)
  snp_pos <- rbind(snp_pos, snp_pos_chr)
}
colnames(snp_pos) <- c("Variants", "rs")
SH_dat <- merge(snp_pos, SH_dat, by.x = "Variants", by.y = "Variants")
save(SH_dat, file = "SH_dat.RData")

snp_inter <- intersect(SH_dat$rs, BMI_dat$SNP)
write.table(snp_inter, file = "/net/mulan/disk2/yasheng/EAS/snp.txt", 
            col.names = F, row.names = F, quote = F)
save(snp_inter, file = "snp_inter.RData")


###
setwd("/net/mulan/disk2/yasheng/out_bbj/raw_data/")
# SH_dat <- data.frame(fread("SH_raw.txt.gz"))
SH_dat <- SH_dat[SH_dat$rs %in% snp_inter, ]
SH_dat$beta_s <- SH_dat$BETA * sqrt(2*SH_dat$MAF*(1 - SH_dat$MAF))
SH_dat_x <- data.frame(MarkerName = SH_dat$rs, 
                       A1 = SH_dat$ALT, 
                       A2 = SH_dat$REF, 
                       beta_s = SH_dat$beta_s)
SH_dat_ldsc <- data.frame(SNP = SH_dat$rs,
                          N = 159095,
                          Z = as.numeric(SH_dat$BETA) / as.numeric(SH_dat$SE),
                          A1 = SH_dat$ALT,
                          A2 = SH_dat$REF)
write.table(SH_dat_x, file = "pheno1_scale.txt", row.names = F, quote = F)
write.table(SH_dat_ldsc, file = "pheno1.ldsc", row.names = F, quote = F)
system("gzip pheno1_scale.txt")
system("gzip pheno1.ldsc")

###
# load("snp_inter.RData")
# BMI_dat <- data.frame(fread("BMI_raw.txt.gz"))
BMI_dat <- BMI_dat[BMI_dat$SNP %in% snp_inter, ]
BMI_dat$beta_s <- BMI_dat$BETA * sqrt(2*BMI_dat$Frq*(1-BMI_dat$Frq))
BMI_dat_x <- data.frame(MarkerName = BMI_dat$SNP, 
                        A1 = BMI_dat$ALT, 
                        A2 = BMI_dat$REF, 
                        beta_s = BMI_dat$beta_s)
BMI_dat_ldsc <- data.frame(SNP = BMI_dat$SNP,
                           N = 158284,
                           Z = as.numeric(BMI_dat$BETA) / as.numeric(BMI_dat$SE),
                           A1 = BMI_dat$ALT,
                           A2 = BMI_dat$REF)
write.table(BMI_dat_x, file = "pheno4_scale.txt", row.names = F, quote = F)
write.table(BMI_dat_ldsc, file = "pheno4.ldsc", row.names = F, quote = F)
system("gzip pheno4_scale.txt")
system("gzip pheno4.ldsc")

###
PLT_dat <- PLT_dat[PLT_dat$SNP %in% snp_inter, ]
PLT_dat_x <- data.frame(MarkerName = PLT_dat$SNP,
                        A1 = PLT_dat$ALT,
                        A2 = PLT_dat$REF,
                        beta_s = PLT_dat$BETA * sqrt(2*PLT_dat$Frq*(1-PLT_dat$Frq)))
write.table(PLT_dat_x, file = "pheno10_scale.txt", row.names = F, quote = F)
system("gzip pheno10_scale.txt")
PLT_dat_ldsc <- data.frame(SNP = PLT_dat$SNP,
                           N = 108208,
                           Z = PLT_dat$BETA / PLT_dat$SE,
                           A1 = PLT_dat$ALT,
                           A2 = PLT_dat$REF)
write.table(PLT_dat_ldsc, file = "pheno10.ldsc", row.names = F, quote = F)
system("gzip pheno10.ldsc")

###
RBC_dat <- RBC_dat[RBC_dat$SNP %in% snp_inter, ]
RBC_dat_x <- data.frame(MarkerName = RBC_dat$SNP,
                        A1 = RBC_dat$ALT,
                        A2 = RBC_dat$REF,
                        beta_s = RBC_dat$BETA * sqrt(2*RBC_dat$Frq*(1-RBC_dat$Frq)))
write.table(RBC_dat_x, file = "pheno11_scale.txt", row.names = F, quote = F)
system("gzip pheno11_scale.txt")
RBC_dat_ldsc <- data.frame(SNP = RBC_dat$SNP,
                           N = 108794,
                           Z = RBC_dat$BETA / RBC_dat$SE,
                           A1 = RBC_dat$ALT,
                           A2 = RBC_dat$REF)
write.table(RBC_dat_ldsc, file = "pheno11.ldsc", row.names = F, quote = F)
system("gzip pheno11.ldsc")

# ###
EOS_dat <- EOS_dat[EOS_dat$SNP %in% snp_inter, ]
EOS_dat_x <- data.frame(MarkerName = EOS_dat$SNP,
                        A1 = EOS_dat$ALT,
                        A2 = EOS_dat$REF,
                        beta_s = EOS_dat$BETA * sqrt(2*EOS_dat$Frq*(1-EOS_dat$Frq)))
write.table(EOS_dat_x, file = "pheno12_scale.txt", row.names = F, quote = F)
system("gzip pheno12_scale.txt")
EOS_dat_ldsc <- data.frame(SNP = EOS_dat$SNP,
                           N = 62076,
                           Z = EOS_dat$BETA / EOS_dat$SE,
                           A1 = EOS_dat$ALT,
                           A2 = EOS_dat$REF)
write.table(EOS_dat_ldsc, file = "pheno12.ldsc", row.names = F, quote = F)
system("gzip pheno12.ldsc")

###
WBC_dat <- WBC_dat[WBC_dat$SNP %in% snp_inter, ]
WBC_dat_x <- data.frame(MarkerName = WBC_dat$SNP,
                        A1 = WBC_dat$ALT,
                        A2 = WBC_dat$REF,
                        beta_s = WBC_dat$BETA * sqrt(2*WBC_dat$Frq*(1-WBC_dat$Frq)))
write.table(WBC_dat_x, file = "pheno13_scale.txt", row.names = F, quote = F)
system("gzip pheno13_scale.txt")
WBC_dat_ldsc <- data.frame(SNP = WBC_dat$SNP,
                           N = 107964,
                           Z = WBC_dat$BETA / WBC_dat$SE,
                           A1 = WBC_dat$ALT,
                           A2 = WBC_dat$REF)
write.table(WBC_dat_ldsc, file = "pheno13.ldsc", row.names = F, quote = F)
system("gzip pheno13.ldsc")
