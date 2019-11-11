rm(list=ls())
library(ggplot2)
library(plyr)

mem_time <- read.table("/net/mulan/disk2/yasheng/simulation_bslmm/Comp/Memory.txt", header = T)
DBSLMM_dat <- mem_time[mem_time[, 2] == "DBSLMM", ]
BSLMM_dat <- mem_time[mem_time[, 2] == "BSLMM", ]

DBSLMM_dat$Time_s / BSLMM_dat$Time_s
DBSLMM_dat$Memory / BSLMM_dat$Memory

DBSLMM_dat$Time_s/60
BSLMM_dat$Time_s/60

DBSLMM_dat$Memory/1024/1024
BSLMM_dat$Memory/1024/1024