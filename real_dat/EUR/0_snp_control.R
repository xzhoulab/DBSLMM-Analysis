rm(list=ls())
library(data.table)
library(plyr)
setwd("/net/mulan/disk2/yasheng/out_sample/raw_data/")

### load datasets of 9 traits
SH_dat <- data.frame(fread("SH_raw.txt.gz"))
BMI_dat <- data.frame(fread("BMI_raw.txt.gz"))
# AM_dat <- data.frame(fread("AM_raw.txt.gz"))
PLT_dat <- data.frame(fread("PLT_raw.txt.gz"))
RBC_dat <- data.frame(fread("RBC_raw.txt.gz"))
EOS_dat <- data.frame(fread("EOS_raw.txt.gz"))
WBC_dat <- data.frame(fread("WBC_raw.txt.gz"))
NS_dat <- data.frame(fread("NS_raw.txt.gz"))
# SBP_dat <- data.frame(fread("SBP_raw.txt.gz"))
YE_dat <- data.frame(fread("YE_raw.txt.gz"))

snp_inter <- Reduce(intersect, list(SH_dat[, 1], BMI_dat[, 1], PLT_dat[, 3], 
                                    NS_dat[, 1], YE_dat[, 1]))
save(snp_inter, file = "snp_inter.RData")


###
setwd("/net/mulan/disk2/yasheng/out_sample/raw_data/")
load("snp_inter.RData")
# SH_dat <- data.frame(fread("SH_raw.txt.gz"))
SH_dat <- SH_dat[SH_dat[, 1] %in% snp_inter, ]
SH_dat$beta_s <- SH_dat$b * sqrt(2*SH_dat$Freq.Allele1.HapMapCEU*(1 - SH_dat$Freq.Allele1.HapMapCEU))
SH_dat_x <- data.frame(MarkerName = SH_dat$MarkerName, 
                       A1 = SH_dat$Allele1, 
                       A2 = SH_dat$Allele2, 
                       beta_s = SH_dat$beta_s)
SH_dat_ldsc <- data.frame(SNP = SH_dat$MarkerName,
                          N = SH_dat$N,
                          Z = as.numeric(SH_dat$b) / as.numeric(SH_dat$SE),
                          A1 = SH_dat$Allele1,
                          A2 = SH_dat$Allele2)[SH_dat[, 1] %in% snp_inter, ]
write.table(SH_dat_x, file = "pheno1_scale.txt", row.names = F, quote = F)
write.table(SH_dat_ldsc, file = "pheno1.ldsc", row.names = F, quote = F)
system("gzip pheno1_scale.txt")
system("gzip pheno1.ldsc")

###
# load("snp_inter.RData")
# BMI_dat <- data.frame(fread("BMI_raw.txt.gz"))
BMI_dat <- BMI_dat[BMI_dat[, 1] %in% snp_inter, ]
BMI_dat$beta_s <- BMI_dat$b * sqrt(2*BMI_dat$Freq1.Hapmap*(1-BMI_dat$Freq1.Hapmap))
BMI_dat_x <- data.frame(MarkerName = BMI_dat$SNP, 
                        A1 = BMI_dat$A1, 
                        A2 = BMI_dat$A2, 
                        beta_s = BMI_dat$beta_s)
BMI_dat_ldsc <- data.frame(SNP = BMI_dat$SNP,
                           N = BMI_dat$N,
                           Z = as.numeric(BMI_dat$b) / as.numeric(BMI_dat$se),
                           A1 = BMI_dat$A1,
                           A2 = BMI_dat$A2)[BMI_dat[, 1] %in% snp_inter, ]
write.table(BMI_dat_x, file = "BMI_scale.txt", row.names = F, quote = F)
write.table(BMI_dat_ldsc, file = "BMI.ldsc", row.names = F, quote = F)
system("gzip pheno4_scale.txt")
system("gzip pheno4.ldsc")

# ###
# AM_dat <- AM_rs_dat[AM_rs_dat[, 1] %in% snp_inter, ]
# AM_dat <- merge(AM_dat, SH_dat, by.x = "Markername", by.y = "MarkerName")
# AM_dat$beta_s <- AM_dat$Effect* sqrt(2*AM_dat$Freq.Allele1.HapMapCEU*(1-AM_dat$Freq.Allele1.HapMapCEU))
# AM_dat_x <- data.frame(MarkerName = AM_dat$Markername,
#                        A1 = toupper(AM_dat$Allele1.x),
#                        A2 = toupper(AM_dat$Allele2.x),
#                        beta_s = AM_dat$beta_s)
# AM_dat_ldsc <- data.frame(SNP = AM_dat$Markername,
#                           N =
#                           A1 = toupper(AM_dat$Allele1.x),
#                           A2 = toupper(AM_dat$Allele2.x),
#                           beta_s = AM_dat$beta_s)
# write.table(AM_dat_x, file = "AM_scale_summ.txt", row.names = F, quote = F)
# write.table(BMI_dat_ldsc, file = "AM.ldsc", row.names = F, quote = F)
# system("gzip AM_scale_summ.txt")
# system("gzip AM.ldsc")

###
PLT_dat <- data.frame(fread("PLT_raw.txt.gz"))
PLT_dat <- PLT_dat[PLT_dat[, 3] %in% snp_inter, 1:6]
PLT_dat <- merge(PLT_dat, SH_dat, by.x = "SNP", by.y = "MarkerName")
PLT_dat$A2 <- ifelse(PLT_dat$A1 == PLT_dat$Allele1, PLT_dat$Allele2, PLT_dat$Allele1)
PLT_dat_x <- data.frame(MarkerName = PLT_dat$SNP,
                        A1 = PLT_dat$A1,
                        A2 = PLT_dat$A2,
                        beta_s = PLT_dat$BETA * sqrt(2*PLT_dat$Freq.Allele1.HapMapCEU*(1-PLT_dat$Freq.Allele1.HapMapCEU)))
write.table(PLT_dat_x, file = "pheno10_scale.txt", row.names = F, quote = F)
system("gzip pheno10_scale.txt")
PLT_dat_ldsc <- data.frame(SNP = PLT_dat$SNP,
                           N = 4250,
                           Z = PLT_dat$BETA / PLT_dat$SE.x,
                           A1 = PLT_dat$A1,
                           A2 = PLT_dat$A2)
write.table(PLT_dat_ldsc, file = "pheno10.ldsc", row.names = F, quote = F)
system("gzip pheno10.ldsc")

###
RBC_dat <- data.frame(fread("RBC_raw.txt.gz"))
RBC_dat <- RBC_dat[RBC_dat[, 3] %in% snp_inter, 1:6]
RBC_dat <- merge(RBC_dat, SH_dat, by.x = "SNP", by.y = "MarkerName")
RBC_dat$A2 <- ifelse(RBC_dat$A1 == RBC_dat$Allele1, RBC_dat$Allele2, RBC_dat$Allele1)
RBC_dat_x <- data.frame(MarkerName = RBC_dat$SNP,
                        A1 = RBC_dat$A1,
                        A2 = RBC_dat$A2,
                        beta_s = RBC_dat$BETA * sqrt(2*RBC_dat$Freq.Allele1.HapMapCEU*(1-RBC_dat$Freq.Allele1.HapMapCEU)))
write.table(RBC_dat_x, file = "pheno11_scale.txt", row.names = F, quote = F)
system("gzip pheno11_scale.txt")
RBC_dat_ldsc <- data.frame(SNP = RBC_dat$SNP,
                           N = 4250,
                           Z = RBC_dat$BETA / RBC_dat$SE.x,
                           A1 = RBC_dat$A1,
                           A2 = RBC_dat$A2)
write.table(RBC_dat_ldsc, file = "pheno11.ldsc", row.names = F, quote = F)
system("gzip pheno11.ldsc")

# ###
EOS_dat <- data.frame(fread("EOS_raw.txt.gz"))
EOS_dat <- EOS_dat[EOS_dat[, 3] %in% snp_inter, 1:6]
EOS_dat <- merge(EOS_dat, SH_dat, by.x = "SNP", by.y = "MarkerName")
EOS_dat$A2 <- ifelse(EOS_dat$A1 == EOS_dat$Allele1, EOS_dat$Allele2, EOS_dat$Allele1)
EOS_dat_x <- data.frame(MarkerName = EOS_dat$SNP,
                        A1 = EOS_dat$A1,
                        A2 = EOS_dat$A2,
                        beta_s = EOS_dat$BETA * sqrt(2*EOS_dat$Freq.Allele1.HapMapCEU*(1-EOS_dat$Freq.Allele1.HapMapCEU)))
write.table(EOS_dat_x, file = "pheno12_scale.txt", row.names = F, quote = F)
system("gzip pheno12_scale.txt")
EOS_dat_ldsc <- data.frame(SNP = EOS_dat$SNP,
                           N = 4250,
                           Z = EOS_dat$BETA / EOS_dat$SE.x,
                           A1 = EOS_dat$A1,
                           A2 = EOS_dat$A2)
write.table(EOS_dat_ldsc, file = "pheno12.ldsc", row.names = F, quote = F)
system("gzip pheno12.ldsc")
#
# ###
WBC_dat <- data.frame(fread("WBC_raw.txt.gz"))
WBC_dat <- WBC_dat[WBC_dat[, 3] %in% snp_inter, 1:6]
WBC_dat <- merge(WBC_dat, SH_dat, by.x = "SNP", by.y = "MarkerName")
WBC_dat$A2 <- ifelse(WBC_dat$A1 == WBC_dat$Allele1, WBC_dat$Allele2, WBC_dat$Allele1)
WBC_dat_x <- data.frame(MarkerName = WBC_dat$SNP,
                        A1 = WBC_dat$A1,
                        A2 = WBC_dat$A2,
                        beta_s = WBC_dat$BETA * sqrt(2*WBC_dat$Freq.Allele1.HapMapCEU*(1-WBC_dat$Freq.Allele1.HapMapCEU)))
write.table(WBC_dat_x, file = "pheno13_scale.txt", row.names = F, quote = F)
system("gzip pheno13_scale.txt")
WBC_dat_ldsc <- data.frame(SNP = WBC_dat$SNP,
                           N = 4250,
                           Z = WBC_dat$BETA / WBC_dat$SE.x,
                           A1 = WBC_dat$A1,
                           A2 = WBC_dat$A2)
write.table(WBC_dat_ldsc, file = "pheno13.ldsc", row.names = F, quote = F)
system("gzip pheno13.ldsc")

###
NS_dat <- data.frame(fread("NS_raw.txt.gz"))
NS_dat <- NS_dat[NS_dat[, 1] %in% snp_inter, ]
NS_dat_x <- data.frame(MarkerName = NS_dat[, 1], 
                       A1 = NS_dat$A1,
                       A2 = NS_dat$A2, 
                       beta_s = NS_dat$Beta * sqrt(2*NS_dat$EAF*(1-NS_dat$EAF)))
write.table(NS_dat_x, file = "pheno15_scale.txt", row.names = F, quote = F)
system("gzip pheno15_scale.txt")
NS_dat_ldsc <- data.frame(SNP = NS_dat[, 1], 
                          N = 170911,
                          Z = NS_dat$Beta / NS_dat$SE, 
                          A1 = NS_dat$A1,
                          A2 = NS_dat$A2)
write.table(NS_dat_ldsc, file = "pheno15.ldsc", row.names = F, quote = F)
system("gzip pheno15.ldsc")

# ###
# SBP_dat <- data.frame(fread("SBP_raw.txt.gz"))
# SBP_dat <- SBP_dat[SBP_dat$SNP.ID %in% snp_inter, ]
# SBP_dat <- merge(SBP_dat, SH_dat, by.x = "SNP.ID", by.y = "MarkerName")
# YE_dat_x <- data.frame(MarkerName = YE_dat[, 1], 
#                        A1 = YE_dat$A1,
#                        A2 = YE_dat$A2, 
#                        beta_s = YE_dat$Beta * sqrt(2*YE_dat$EAF*(1-YE_dat$EAF)))

###
YE_dat <- data.frame(fread("YE_raw.txt.gz"))
YE_dat <- YE_dat[YE_dat[, 1] %in% snp_inter, ]
YE_dat_x <- data.frame(MarkerName = YE_dat[, 1], 
                       A1 = YE_dat$A1,
                       A2 = YE_dat$A2, 
                       beta_s = YE_dat$Beta * sqrt(2*YE_dat$EAF*(1-YE_dat$EAF)))
YE_dat_ldsc <- data.frame(SNP = YE_dat[, 1], 
                          N = 170911,
                          Z = YE_dat$Beta / YE_dat$SE, 
                          A1 = YE_dat$A1,
                          A2 = YE_dat$A2)
write.table(YE_dat_x, file = "pheno16_scale.txt", row.names = F, quote = F)
system("gzip pheno16_scale.txt")
write.table(YE_dat_ldsc, file = "pheno16.ldsc", row.names = F, quote = F)
system("gzip pheno16.ldsc")