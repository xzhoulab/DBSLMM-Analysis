rm(list=ls())
library(ggplot2)
library(plyr)
library(data.table)
load("/net/mulan/home/yasheng/Biobank/code/figures/simulation_bslmm/summary_result.RData")

## corelation and accuracy difference
SLM_sub <- r2_SLM_dat[r2_SLM_dat$Method %in% c("DBSLMM", "BSLMM"), ]
LMM_sub <- r2_LMM_dat[r2_LMM_dat$Method %in% c("DBSLMM", "BSLMM"), ]
SLMM_sub <- r2_SLMM_dat[r2_SLMM_dat$Method %in% c("DBSLMM", "BSLMM"), -2]

tot_sub <- rbind(rbind(SLM_sub, LMM_sub), SLMM_sub)
tot_sub <- data.frame(Settings = factor(c(rep("Sparse", nrow(SLM_sub)), rep("Polygenic", nrow(LMM_sub)), 
                                          rep("Hybrid I", nrow(SLM_sub)), 
                                          rep("Hybrid II", nrow(SLM_sub))), 
                                        levels = c("Sparse", "Polygenic", "Hybrid I", "Hybrid II")), 
                      tot_sub)     
DBSLMM_sub <- tot_sub[tot_sub$Method == "DBSLMM", -3]
BSLMM_sub <- tot_sub[tot_sub$Method == "BSLMM", -3]
lm_sub <- data.frame(Settings = DBSLMM_sub$Settings, 
                     PVE = DBSLMM_sub$PVE,
                     DBSLMM = DBSLMM_sub[, 4], 
                     BSLMM = BSLMM_sub[, 4])
cor(lm_sub$DBSLMM, lm_sub$BSLMM)

lm_plt <- ggplot(lm_sub, aes(x=DBSLMM, y=BSLMM, group = Settings))+
  geom_point(aes(color=Settings), size = 5)+
  scale_color_manual(values=c("#8A4198FF", "#008EA0FF", "#FF6F00FF", "#C71000FF"))+ 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 30),  
        axis.title.x = element_text(size = 60),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5))+
  scale_x_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5))+
  geom_abline(intercept = 0, slope = 1, color="#EFC000FF", 
              linetype="dashed", size=1.5)+
  ylab(expression(italic(R[BSLMM]^2)))+
  xlab(expression(italic(R[DBSLMM]^2)))+
  annotate(geom="text", x=0.1, y=0.4, label="italic(R) ^ 2 == 0.969",
           color="black", parse = TRUE, size = 20)+
  labs(tag = "A", title = "Correlation")

PVE_s <- as.numeric(laply(strsplit(as.character(lm_sub$PVE), "="), function (a) a[2]))
lm_sub$diff <- (lm_sub$DBSLMM - lm_sub$BSLMM) / PVE_s

loss_per_mean <- mean((lm_sub$BSLMM - lm_sub$DBSLMM))

diff_plt <- ggplot(lm_sub, aes(x=Settings, y=diff, fill=Settings))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  scale_fill_manual(values=c("#8A4198FF", "#008EA0FF", "#FF6F00FF", "#C71000FF"))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_text(size = 60),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.6, 0.9),
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-1, 0.5, 0.5), limits = c(-1, 0.5))+
  geom_abline(intercept = 0, slope = 0, color="#EFC000FF", 
              linetype="dashed", size=1.5)+
  ylab("Difference")+
  xlab("Scenario")+
  labs(tag = "B", title = "Accuracy Difference")

## color
col_dat <- read.table("/net/mulan/home/yasheng/Biobank/code/figures/col.txt")
method_ord <- c("DBSLMM", "BSLMM")
col_dat <- col_dat[match(method_ord, col_dat[, 2]), ]

## memory and time
mem_time <- read.table("/net/mulan/disk2/yasheng/simulation_bslmm/Comp/Memory.txt", header = T)
mem_time[, 3] <- mem_time[, 3]/60
mem_time[, 4] <- mem_time[, 4]/1024/1024
mem_time[, 1] <- factor(as.character(mem_time[, 1]), levels = c("200", "500", "1000", "2000", "5000", "10000"))
mem_time[, 2] <- factor(mem_time[, 2], levels = method_ord)

time_plt <- ggplot(mem_time, aes(x=SampleSize, y=Time_s, group=Method, color=Method))+
  geom_point(aes(shape = Method), size = 10)+
  geom_line(linetype="dashed", size=1.2) +
  theme_bw() + 
  scale_color_manual(values=as.character(col_dat[, 1]))+
  theme(axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 60),
        axis.text.y = element_text(size = 28),
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1, color = "#E0E0E0"),
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+
  scale_y_continuous(breaks = seq(0, 80, 20), limits = c(0, 81))+
  ylab("CPU time (min)")+
  xlab(expression(italic(N) ~ " = # of Samples"))+ 
  labs(tag = "C", title = "Time Usage")

mem_plt <- ggplot(mem_time, aes(x=SampleSize, y=Memory, group=Method, color=Method))+
  geom_point(aes(shape = Method), size = 10)+
  geom_line(linetype="dashed", size=1.2) +
  theme_bw() + 
  scale_color_manual(values=as.character(col_dat[, 1]))+
  theme(axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 60),
        axis.text.y = element_text(size = 28),
        axis.title.y = element_text(size = 60),
        legend.text = element_text(size = 25, face = "bold"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.3, 0.8),
        title = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1, color = "#E0E0E0"),
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+
  scale_y_continuous(breaks = seq(0, 12, 3), limits = c(0, 12.2))+
  ylab("Memory (GB)")+
  xlab(expression(italic(N) ~ " = # of Samples"))+ 
  labs(tag = "D", title = "CPU Usage")

library(patchwork)
tiff("/net/mulan/disk2/yasheng/figures/simulation_bslmm/Comp_DBSLMM_BSLMM.tiff", height = 1600, width = 1600)
(lm_plt | diff_plt) /
  (time_plt | mem_plt)
dev.off()