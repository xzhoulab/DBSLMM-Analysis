rm(list = ls())
library(Metrics)

### paramters
herit <- c(0.1, 0.2, 0.5)
method_label <- c("DBSLMM", "P+T", "LDpred_inf", "LDpred_sparse", "SBLUP")
method <- c("slmm", "PT", "ldpred", "ldpred", "sblup")
method_ord <- c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")
status <- c("", "", "_inf", "_pbest", "_200")

### estmation expectation of R2
ldsc_chr1 <- data.frame(fread(paste0("/net/mulan/disk2/yasheng/sample500/ldsc/1.l2.ldscore.gz")))[, 4]
ldsc_mean <- mean(ldsc_chr1[1:100000])
M <- 100000 / ldsc_mean
N <- 10000
E_r2 <- vector("numeric", 3)
for (i in 1:3){
  E_r2[i] <- herit[i]/(1 + M /(N*herit[i]))
}


### SLM datasets R2 and mse
r2_SLM <- vector()
mse_SLM <- vector()
for (h in 1: length(herit)){
  for (m in 1: length(method)){
    for (c in 1: 10){
      pred_str <- paste0("/net/mulan/disk2/yasheng/simulation2/SLM/", method[m], "/pheno_her", herit[h], 
                         "_cross", c, "_dist1_ps0.001", status[m], ".profile")
      pred_tmp <- try(read.table(pred_str, header = T), silent = T)
      if (inherits(pred_tmp, "try-error")){
        r2 <- NA
        mse <- NA
      }else{
        r2 <- cor(pred_tmp[, 3], pred_tmp[, 6])^2
        mse <- mse(pred_tmp[, 3], pred_tmp[, 6])
      }
      r2_v <- c(paste0("heritability=", herit[h]), method_label[m], c, r2)
      mse_v <- c(paste0("heritability=", herit[h]), method_label[m], c, mse)
      r2_SLM <- rbind(r2_SLM, r2_v)
      mse_SLM <- rbind(mse_SLM, mse_v)
    }
  }
}

for (h in 1:length(herit)){
  for (c in 1: 10){
    pred_str <- paste0("/net/mulan/disk2/yasheng/simulation2/SLM/lassosum/res_her", herit[h], 
                       "_cross", c, "_dist1_ps0.001.txt.txt")
    pred_tmp <- try(read.table(pred_str, header = F), silent = T)
    if (inherits(pred_tmp, "try-error")){
      r2 <- NA
      mse <- NA
    }else{
      r2 <- pred_tmp[1, 1]
      mse <- pred_tmp[2, 1]
    }
    r2_v <- c(paste0("heritability=", herit[h]), "lassosum", c, r2)
    mse_v <- c(paste0("heritability=", herit[h]), "lassosum", c, mse)
    r2_SLM <- rbind(r2_SLM, r2_v)
    mse_SLM <- rbind(mse_SLM, mse_v)
  }
}

r2_SLM_dat <- data.frame(heritability = factor(r2_SLM[, 1], 
                                               levels = c("heritability=0.1", "heritability=0.2", "heritability=0.5")), 
                         Method = factor(r2_SLM[, 2], levels = method_ord), 
                         cross = r2_SLM[, 3], 
                         R2 = as.numeric(r2_SLM[, 4]))
mse_SLM_dat <- data.frame(heritability = factor(mse_SLM[, 1], 
                                                levels = c("heritability=0.1", "heritability=0.2", "heritability=0.5")), 
                          Method = factor(mse_SLM[, 2], levels = method_ord), 
                          cross = mse_SLM[, 3], 
                          MSE = as.numeric(mse_SLM[, 4]))
r2_SLM_dat1 <- r2_SLM_dat[r2_SLM_dat$heritability == "heritability=0.1", ]
r2_SLM_dat2 <- r2_SLM_dat[r2_SLM_dat$heritability == "heritability=0.2", ]
r2_SLM_dat3 <- r2_SLM_dat[r2_SLM_dat$heritability == "heritability=0.5", ]


### LMM assumption
r2_LMM <- vector()
mse_LMM <- vector()
for (h in 1: length(herit)){
  for (m in 1: length(method)){
    for (c in 1: 10){
      pred_str <- paste0("/net/mulan/disk2/yasheng/simulation2/LMM/", method[m], "/pheno_her", herit[h], 
                         "_cross", c, "_dist1", status[m], ".profile")
      pred_tmp <- try(read.table(pred_str, header = T), silent = T)
      if (inherits(pred_tmp, "try-error")){
        r2 <- NA
        mse <- NA
      }else{
        r2 <- cor(pred_tmp[, 3], pred_tmp[, 6])^2
        mse <- mse(pred_tmp[, 3], pred_tmp[, 6])
      }
      r2_v <- c(paste0("heritability=", herit[h]), method_label[m], c, r2)
      mse_v <- c(paste0("heritability=", herit[h]), method_label[m], c, mse)
      r2_LMM <- rbind(r2_LMM, r2_v)
      mse_LMM <- rbind(mse_LMM, mse_v)
    }
  }
}

#lassosum
for (h in 1: length(herit)){
  for (c in 1: 10){
    pred_str <- paste0("/net/mulan/disk2/yasheng/simulation2/LMM/lassosum/res_her", herit[h], 
                       "_cross", c, "_dist1.txt")
    pred_tmp <- try(read.table(pred_str, header = F), silent = T)
    if (inherits(pred_tmp, "try-error")){
      r2 <- NA
      mse <- NA
    }else{
      r2 <- pred_tmp[1, 1]
      mse <- pred_tmp[2, 1]
    }
    r2_v <- c(paste0("heritability=", herit[h]), "lassosum", c, r2)
    mse_v <- c(paste0("heritability=", herit[h]), "lassosum", c, mse)
    r2_LMM <- rbind(r2_LMM, r2_v)
    mse_LMM <- rbind(mse_LMM, mse_v)
  }
}


r2_LMM_dat <- data.frame(heritability = factor(r2_LMM[, 1], 
                                               levels = c("heritability=0.1", "heritability=0.2", "heritability=0.5")), 
                         Method = factor(r2_LMM[, 2], levels = method_ord), 
                         cross = r2_LMM[, 3], 
                         R2 = as.numeric(r2_LMM[, 4]))
mse_LMM_dat <- data.frame(heritability = factor(mse_LMM[, 1], levels = c("heritability=0.1", "heritability=0.2", "heritability=0.5")), 
                          Method = factor(mse_LMM[, 2], levels = method_ord), 
                          cross = mse_LMM[, 3], 
                          MSE = as.numeric(mse_LMM[, 4]))
r2_LMM_dat1 <- r2_LMM_dat[r2_LMM_dat$heritability == "heritability=0.1", ]
r2_LMM_dat2 <- r2_LMM_dat[r2_LMM_dat$heritability == "heritability=0.2", ]
r2_LMM_dat3 <- r2_LMM_dat[r2_LMM_dat$heritability == "heritability=0.5", ]


### SLMM
prop <- c(0.2, 0.5)
r2_SLMM <- vector()
mse_SLMM <- vector()
for (pp in 1: length(prop)){
  for (h in 1: length(herit)){
    for (m in 1: length(method)){
      for (c in 1: 10){
        pred_str <- paste0("/net/mulan/disk2/yasheng/simulation2/SLMM/", method[m], "/pheno_her", herit[h],
                           "_cross", c, "_dist1", "_ps0.001", "_prop", prop[pp], status[m], ".profile")
        pred_tmp <- try(read.table(pred_str, header = T), silent = T)
        if (inherits(pred_tmp, "try-error")){
          r2 <- NA
          mse <- NA
        }else{
          r2 <- cor(pred_tmp[, 3], pred_tmp[, 6])^2
          mse <- mse(pred_tmp[, 3], pred_tmp[, 6])
        }
        r2_v <- c(paste0("heritability=", herit[h]), paste0("prop=", prop[pp]),
                  method_label[m], c, r2)
        mse_v <- c(paste0("heritability=", herit[h]), paste0("prop=", prop[pp]),
                   method_label[m], c, mse)
        r2_SLMM <- rbind(r2_SLMM, r2_v)
        mse_SLMM <- rbind(mse_SLMM, mse_v)
      }
    }
  }
}
for (pp in 1: length(prop)){
  for (h in 1: length(herit)){
    for (c in 1: 10){
      pred_str <- paste0("/net/mulan/disk2/yasheng/simulation2/SLMM/lassosum/res_her", herit[h], 
                         "_cross", c, "_dist1", "_ps0.001", "_prop", prop[pp], ".txt.txt")
      pred_tmp <- try(read.table(pred_str, header = F), silent = T)
      if (inherits(pred_tmp, "try-error")){
        r2 <- NA
        mse <- NA
      }else{
        r2 <- pred_tmp[1, 1]
        mse <- pred_tmp[2, 1]
      }
      r2_v <- c(paste0("heritability=", herit[h]), paste0("prop=", prop[pp]),
                "lassosum", c, r2)
      mse_v <- c(paste0("heritability=", herit[h]), paste0("prop=", prop[pp]),
                 "lassosum", c, mse)
      r2_SLMM <- rbind(r2_SLMM, r2_v)
      mse_SLMM <- rbind(mse_SLMM, mse_v)
    }
  }
}

r2_SLMM_dat <- data.frame(heritability = factor(r2_SLMM[, 1], 
                                                levels = c("heritability=0.1", "heritability=0.2", "heritability=0.5")), 
                          PS = r2_SLMM[, 2],
                          Method = factor(r2_SLMM[, 3], levels = method_ord), 
                          cross = r2_SLMM[, 4], 
                          R2 = as.numeric(r2_SLMM[, 5]))
mse_SLMM_dat <- data.frame(heritability = factor(mse_SLMM[, 1], levels = c("heritability=0.1", "heritability=0.2", "heritability=0.5")), 
                           PS = mse_SLMM[, 2],
                           Method = factor(mse_SLMM[, 3], levels = method_ord), 
                           cross = mse_SLMM[, 4], 
                           MSE = as.numeric(mse_SLMM[, 5]))
r2_SLMM_dat11 <- r2_SLMM_dat[r2_SLMM_dat$PS == "prop=0.2"&r2_SLMM_dat$heritability == "heritability=0.1", ]
r2_SLMM_dat12 <- r2_SLMM_dat[r2_SLMM_dat$PS == "prop=0.2"&r2_SLMM_dat$heritability == "heritability=0.2", ]
r2_SLMM_dat13 <- r2_SLMM_dat[r2_SLMM_dat$PS == "prop=0.2"&r2_SLMM_dat$heritability == "heritability=0.5", ]
r2_SLMM_dat21 <- r2_SLMM_dat[r2_SLMM_dat$PS == "prop=0.5"&r2_SLMM_dat$heritability == "heritability=0.1", ]
r2_SLMM_dat22 <- r2_SLMM_dat[r2_SLMM_dat$PS == "prop=0.5"&r2_SLMM_dat$heritability == "heritability=0.2", ]
r2_SLMM_dat23 <- r2_SLMM_dat[r2_SLMM_dat$PS == "prop=0.5"&r2_SLMM_dat$heritability == "heritability=0.5", ]
mse_SLMM_dat1 <- mse_SLMM_dat[mse_SLMM_dat$PS == "prop=0.2", ]
mse_SLMM_dat2 <- mse_SLMM_dat[mse_SLMM_dat$PS == "prop=0.5", ]

save.image("/net/mulan/home/yasheng/Biobank/code/figures/simulation/normal_summary_result.RData")