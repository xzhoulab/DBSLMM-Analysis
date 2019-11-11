rm(list=ls())
library(ggplot2)
library(data.table)
load("/net/mulan/home/yasheng/Biobank/code/figures/simulation_bslmm/summary_result.RData")

## color
col_dat <- read.table("/net/mulan/home/yasheng/Biobank/code/figures/col.txt")
method_ord <- c("DBSLMM", "BSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "lasso", "P+T")
col_dat <- col_dat[match(method_ord, col_dat[, 2]), ]

## Sparse scenario difference
r2_SLM_diff_dat1 <- r2_SLM_dat1
r2_SLM_diff_dat1[, 4] <- r2_SLM_dat1[, 4] - as.numeric(r2_SLM_dat1[r2_SLM_dat1[, 2] == "DBSLMM", 4])
r2_SLM_diff_dat1 <- r2_SLM_diff_dat1[r2_SLM_diff_dat1[, 2] != "DBSLMM", ]
r2_SLM_diff_dat2 <- r2_SLM_dat2
r2_SLM_diff_dat2[, 4] <- r2_SLM_dat2[, 4] - as.numeric(r2_SLM_dat2[r2_SLM_dat2[, 2] == "DBSLMM", 4])
r2_SLM_diff_dat2 <- r2_SLM_diff_dat2[r2_SLM_diff_dat2[, 2] != "DBSLMM", ]
r2_SLM_diff_dat3 <- r2_SLM_dat3
r2_SLM_diff_dat3[, 4] <- r2_SLM_dat3[, 4] - as.numeric(r2_SLM_dat3[r2_SLM_dat3[, 2] == "DBSLMM", 4])
r2_SLM_diff_dat3 <- r2_SLM_diff_dat3[r2_SLM_diff_dat3[, 2] != "DBSLMM", ]

r2_SLM_diff_plt1 <- ggplot(r2_SLM_diff_dat1, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.04, 0.04, 0.02), limits = c(-0.05, 0.05))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "A", title = "Sparse")
r2_SLM_diff_plt2 <- ggplot(r2_SLM_diff_dat2, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.1, 0.1, 0.05), limits = c(-0.1, 0.1))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "A", title = "Sparse")

r2_SLM_diff_plt3 <- ggplot(r2_SLM_diff_dat3, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.3, 0.3, 0.15), limits = c(-0.3, 0.3))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "A", title = "Sparse")

## Polygenic scenario difference
r2_LMM_diff_dat1 <- r2_LMM_dat1
r2_LMM_diff_dat1[, 4] <- r2_LMM_dat1[, 4] - as.numeric(r2_LMM_dat1[r2_LMM_dat1[, 2] == "DBSLMM", 4])
r2_LMM_diff_dat1 <- r2_LMM_diff_dat1[r2_LMM_diff_dat1[, 2] != "DBSLMM", ]
r2_LMM_diff_dat2 <- r2_LMM_dat2
r2_LMM_diff_dat2[, 4] <- r2_LMM_dat2[, 4] - as.numeric(r2_LMM_dat2[r2_LMM_dat2[, 2] == "DBSLMM", 4])
r2_LMM_diff_dat2 <- r2_LMM_diff_dat2[r2_LMM_diff_dat2[, 2] != "DBSLMM", ]
r2_LMM_diff_dat3 <- r2_LMM_dat3
r2_LMM_diff_dat3[, 4] <- r2_LMM_dat3[, 4] - as.numeric(r2_LMM_dat3[r2_LMM_dat3[, 2] == "DBSLMM", 4])
r2_LMM_diff_dat3 <- r2_LMM_diff_dat3[r2_LMM_diff_dat3[, 2] != "DBSLMM", ]

r2_LMM_diff_plt1 <- ggplot(r2_LMM_diff_dat1, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.9),
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.03, 0.03, 0.015), limits = c(-0.03, 0.03))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "B", title = "Polygenic")

r2_LMM_diff_plt2 <- ggplot(r2_LMM_diff_dat2, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.9),
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.1, 0.1, 0.05), limits = c(-0.1, 0.1))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "B", title = "Polygenic")

r2_LMM_diff_plt3 <- ggplot(r2_LMM_diff_dat3, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.9),
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1), limits = c(-0.2, 0.2))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "B", title = "Polygenic")

## Hybrid scenrio difference
r2_SLMM_diff_dat11 <- r2_SLMM_dat11
r2_SLMM_diff_dat11[, 5] <- r2_SLMM_dat11[, 5] - as.numeric(r2_SLMM_dat11[r2_SLMM_dat11[, 3] == "DBSLMM", 5])
r2_SLMM_diff_dat11 <- r2_SLMM_diff_dat11[r2_SLMM_diff_dat11[, 3] != "DBSLMM", ]
r2_SLMM_diff_dat12 <- r2_SLMM_dat12
r2_SLMM_diff_dat12[, 5] <- r2_SLMM_diff_dat12[, 5] - as.numeric(r2_SLMM_dat12[r2_SLMM_dat12[, 3] == "DBSLMM", 5])
r2_SLMM_diff_dat12 <- r2_SLMM_diff_dat12[r2_SLMM_diff_dat12[, 3] != "DBSLMM", ]
r2_SLMM_diff_dat13 <- r2_SLMM_dat13
r2_SLMM_diff_dat13[, 5] <- r2_SLMM_dat13[, 5] - as.numeric(r2_SLMM_dat13[r2_SLMM_dat13[, 3] == "DBSLMM", 5])
r2_SLMM_diff_dat13 <- r2_SLMM_diff_dat13[r2_SLMM_diff_dat13[, 3] != "DBSLMM", ]

r2_SLMM_diff_dat21 <- r2_SLMM_dat21
r2_SLMM_diff_dat21[, 5] <- r2_SLMM_dat21[, 5] - as.numeric(r2_SLMM_dat21[r2_SLMM_dat21[, 3] == "DBSLMM", 5])
r2_SLMM_diff_dat21 <- r2_SLMM_diff_dat21[r2_SLMM_diff_dat21[, 3] != "DBSLMM", ]
r2_SLMM_diff_dat22 <- r2_SLMM_dat22
r2_SLMM_diff_dat22[, 5] <- r2_SLMM_diff_dat22[, 5] - as.numeric(r2_SLMM_diff_dat22[r2_SLMM_dat22[, 3] == "DBSLMM", 5])
r2_SLMM_diff_dat22 <- r2_SLMM_diff_dat22[r2_SLMM_diff_dat22[, 3] != "DBSLMM", ]
r2_SLMM_diff_dat23 <- r2_SLMM_dat23
r2_SLMM_diff_dat23[, 5] <- r2_SLMM_dat23[, 5] - as.numeric(r2_SLMM_dat23[r2_SLMM_dat23[, 3] == "DBSLMM", 5])
r2_SLMM_diff_dat23 <- r2_SLMM_diff_dat23[r2_SLMM_diff_dat23[, 3] != "DBSLMM", ]

# max(r2_SLMM_diff_dat11$R2)
# min(r2_SLMM_diff_dat11$R2)
r2_SLMM_diff_plt11 <- ggplot(r2_SLMM_diff_dat11, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.1, 0.1, 0.05), limits = c(-0.1, 0.1))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "C", title = "Hybrid I")

# max(r2_SLMM_diff_dat12$R2)
# min(r2_SLMM_diff_dat12$R2)
r2_SLMM_diff_plt12 <- ggplot(r2_SLMM_diff_dat12, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.1, 0.1, 0.05), limits = c(-0.1, 0.1))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "C", title = "Hybrid I")

# max(r2_SLMM_diff_dat13$R2)
# min(r2_SLMM_diff_dat13$R2)
r2_SLMM_diff_plt13 <- ggplot(r2_SLMM_diff_dat13, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.15, 0.15, 0.075), limits = c(-0.15, 0.15))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "C", title = "Hybrid I")

# max(r2_SLMM_diff_dat21$R2)
# min(r2_SLMM_diff_dat21$R2)
r2_SLMM_diff_plt21 <- ggplot(r2_SLMM_diff_dat21, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.05, 0.05, 0.025), limits = c(-0.05, 0.05))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "D", title = "Hybrid II")

# max(r2_SLMM_diff_dat22$R2)
# min(r2_SLMM_diff_dat22$R2)
r2_SLMM_diff_plt22 <- ggplot(r2_SLMM_diff_dat22, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1), limits = c(-0.2, 0.2))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "D", title = "Hybrid II")

# max(r2_SLMM_diff_dat23$R2)
# min(r2_SLMM_diff_dat23$R2)
r2_SLMM_diff_plt23 <- ggplot(r2_SLMM_diff_dat23, aes(x=Method, y=R2, fill=Method))+
  geom_boxplot(size = 1)+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -6, aes(label=round(..y.., digits=4)), size = 8) +
  scale_fill_manual(values=as.character(col_dat[-1, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 60, face = "bold"),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(-0.15, 0.15, 0.075), limits = c(-0.15, 0.15))+
  geom_hline(aes(yintercept = 0), color="#EFC000FF",
             linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression(italic(R^2) ~ "Difference"))+
  labs(tag = "D", title = "Hybrid II")

## Output
library(patchwork)
tiff("/net/mulan/disk2/yasheng/figures/simulation_bslmm/normal_diff_herit0.1.tiff", height = 1200, width = 1600)
(r2_SLM_diff_plt1 | r2_LMM_diff_plt1) /
  (r2_SLMM_diff_plt11 | r2_SLMM_diff_plt21)
dev.off()
tiff("/net/mulan/disk2/yasheng/figures/simulation_bslmm/normal_diff_herit0.2.tiff", height = 1200, width = 1600)
(r2_SLM_diff_plt2 | r2_LMM_diff_plt2) /
  (r2_SLMM_diff_plt12 | r2_SLMM_diff_plt22)
dev.off()
tiff("/net/mulan/disk2/yasheng/figures/simulation_bslmm/normal_diff_herit0.5.tiff", height = 1200, width = 1600)
(r2_SLM_diff_plt3 | r2_LMM_diff_plt3) /
  (r2_SLMM_diff_plt13 | r2_SLMM_diff_plt23)
dev.off()
