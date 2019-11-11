rm(list=ls())
library(ggplot2)
library(data.table)
load("/net/mulan/home/yasheng/Biobank/code/figures/simulation/laplace_summary_result.RData")

### color
col_dat <- read.table("/net/mulan/home/yasheng/Biobank/code/figures/col.txt")
method_ord <- c("DBSLMM", "SBLUP", "LDpred_inf", "LDpred_sparse", "lassosum", "P+T")
col_dat <- col_dat[match(method_ord, col_dat[, 2]), ]

## Sparse setting R2
r2_SLM_plt1 <- ggplot(r2_SLM_dat1, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE,
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.1, 0.025), limits = c(0, 0.1))+
  ylab(expression("Prediction " ~ italic(R^2)))+ 
  # geom_hline(aes(yintercept = E_r2[1]), color="#EFC000FF",
  #            linetype="dashed", size=1.5, na.rm = T) +
  labs(tag = "A", title = 'Sparse') 
r2_SLM_plt2 <- ggplot(r2_SLM_dat2, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.2, 0.05), limits = c(0, 0.2))+
  ylab(expression("Prediction " ~ italic(R^2)))+ 
  # geom_hline(aes(yintercept = E_r2[2]), color="#EFC000FF",
  #            linetype="dashed", size=1.5, na.rm = T) +
  labs(tag = "A", title = 'Sparse') 

r2_SLM_plt3 <- ggplot(r2_SLM_dat3, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5))+
  ylab(expression("Prediction " ~ italic(R^2)))+ 
  # geom_hline(aes(yintercept = E_r2[3]), color="#EFC000FF",
  #            linetype="dashed", size=1.5, na.rm = T) +
  labs(tag = "A", title = 'Sparse') 

## Sparse setting MSE
source("/net/mulan/home/yasheng/Biobank/code/figures/r2_est/override.R")
mse_SLM_plt <- ggplot(mse_SLM_dat, aes(x=Method, y=MSE, fill=Method))+
  geom_boxplot(width=0.7, position=position_dodge(0.8))+
  facet_wrap_custom(~heritability, scales = "free", ncol = 1, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0, 8, 4), limits = c(0, 8))), 
    scale_override(2, scale_y_continuous(breaks = seq(0, 12, 3), limits = c(0, 12))), 
    scale_override(3, scale_y_continuous(breaks = seq(0, 300, 75), limits = c(0, 300)))
  ))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -2, aes(label=round(..y.., digits=3)), size = 10) +
  scale_fill_manual(values=as.character(col_dat[, 1]))+
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 20))+
  ylab("Prediction MSE")+
  labs(tag = "A", title = 'Sparse') 

r2_LMM_plt1 <- ggplot(r2_LMM_dat1, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE,
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF",
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
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
  scale_y_continuous(breaks = seq(0, 0.1, 0.025), limits = c(0, 0.1))+
  # geom_hline(aes(yintercept = E_r2[1]), color="#EFC000FF", 
  #            linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression("Prediction " ~ italic(R^2)))+
  labs(tag = "B", title = "Polygenic")

r2_LMM_plt2 <- ggplot(r2_LMM_dat2, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
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
  scale_y_continuous(breaks = seq(0, 0.2, 0.05), limits = c(0, 0.2))+
  # geom_hline(aes(yintercept = E_r2[2]), color="#EFC000FF", 
  #            linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression("Prediction " ~ italic(R^2)))+
  labs(tag = "B", title = "Polygenic")

r2_LMM_plt3 <- ggplot(r2_LMM_dat3, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
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
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5))+
  # geom_hline(aes(yintercept = E_r2[3]), color="#EFC000FF", 
  #            linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression("Prediction " ~ italic(R^2)))+
  labs(tag = "B", title = "Polygenic")

mse_LMM_plt <- ggplot(mse_LMM_dat, aes(x=Method, y=MSE, fill=Method))+
  geom_boxplot(width=0.7, position=position_dodge(0.8))+
  geom_boxplot(width=0.7, position=position_dodge(0.8))+
  facet_wrap_custom(~heritability, scales = "free", ncol = 1, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0, 600, 150), limits = c(0, 600))), 
    scale_override(2, scale_y_continuous(breaks = seq(0, 12, 3), limits = c(0, 12))), 
    scale_override(3, scale_y_continuous(breaks = seq(0, 160, 40), limits = c(0, 160)))
  ))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -2, aes(label=round(..y.., digits=3)), size = 10) +
  scale_fill_manual(values=as.character(col_dat[, 1]))+
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.4, 0.9),
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  ylab("Prediction MSE")+
  labs(tag = "B", title = "Polygenic")

r2_SLMM_plt11 <- ggplot(r2_SLMM_dat11, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.1, 0.025), limits = c(0, 0.1))+
  # geom_hline(aes(yintercept = E_r2[1]), color="#EFC000FF",
  #            linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression("Prediction " ~ italic(R^2)))+
  labs(tag = "C", title = "Hybrid I")

r2_SLMM_plt12 <- ggplot(r2_SLMM_dat12, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.2, 0.05), limits = c(0, 0.2))+
  # geom_hline(aes(yintercept = E_r2[2]), color="#EFC000FF", 
  #            linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression("Prediction " ~ italic(R^2)))+
  labs(tag = "C", title = "Hybrid I")

r2_SLMM_plt13 <- ggplot(r2_SLMM_dat13, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5))+
  # geom_hline(aes(yintercept = E_r2[3]), color="#EFC000FF", 
  #            linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression("Prediction " ~ italic(R^2)))+
  labs(tag = "C", title = "Hybrid I")

r2_SLMM_plt21 <- ggplot(r2_SLMM_dat21, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.1, 0.025), limits = c(0, 0.1))+
  # geom_hline(aes(yintercept = E_r2[1]), color="#EFC000FF", 
  #            linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression("Prediction " ~ italic(R^2)))+
  labs(tag = "D", title = "Hybrid II")

r2_SLMM_plt22 <- ggplot(r2_SLMM_dat22, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.2, 0.05), limits = c(0, 0.2))+
  # geom_hline(aes(yintercept = E_r2[2]), color="#EFC000FF", 
  #            linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression("Prediction " ~ italic(R^2)))+
  labs(tag = "D", title = "Hybrid II")

r2_SLMM_plt23 <- ggplot(r2_SLMM_dat23, aes(x=Method, y=R2, color=Method))+
  geom_jitter(aes(color = Method), size = 5, position=position_jitter(0.2))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -4, aes(label=round(..y.., digits=3)), size = 10) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, colour = "#EFC000FF", 
               geom = "crossbar", width = 0.5) +
  scale_color_manual(values=as.character(col_dat[, 1]))+ 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5))+
  # geom_hline(aes(yintercept = E_r2[3]), color="#EFC000FF", 
  #            linetype="dashed", size=1.5, na.rm = T) +
  ylab(expression("Prediction " ~ italic(R^2)))+
  labs(tag = "D", title = "Hybrid II")
mse_SLMM_plt1 <- ggplot(mse_SLMM_dat1, aes(x=Method, y=MSE, fill=Method))+
  geom_boxplot(width=0.7, position=position_dodge(0.8))+
  facet_wrap_custom(~heritability, scales = "free", ncol = 1, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0, 800, 200), limits = c(0, 800))), 
    scale_override(2, scale_y_continuous(breaks = seq(0, 10, 2.5), limits = c(0, 10))), 
    scale_override(3, scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100)))
  ))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -2, aes(label=round(..y.., digits=3)), size = 10) +
  scale_fill_manual(values=as.character(col_dat[, 1]))+
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  ylab("Prediction MSE")+
  labs(tag = "C", title = "Hybrid I")

mse_SLMM_plt2 <- ggplot(mse_SLMM_dat2, aes(x=Method, y=MSE, fill=Method))+
  geom_boxplot(width=0.7, position=position_dodge(0.8))+
  facet_wrap_custom(~heritability, scales = "free", ncol = 1, scale_overrides = list(
    scale_override(1, scale_y_continuous(breaks = seq(0, 400, 100), limits = c(0, 400))), 
    scale_override(2, scale_y_continuous(breaks = seq(0, 8, 4), limits = c(0, 8))), 
    scale_override(3, scale_y_continuous(breaks = seq(0, 240, 60), limits = c(0, 240)))
  ))+
  stat_summary(fun.y = mean, colour = "black", geom="text", show_guide = FALSE, 
               vjust = -2, aes(label=round(..y.., digits=3)), size = 10) +
  scale_fill_manual(values=as.character(col_dat[, 1]))+
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 60),
        legend.position = "none",
        title = element_text(size = 20), 
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, color = "#E0E0E0"), 
        plot.tag = element_text(size = 50, face = "bold"),
        strip.text = element_text(face = 'bold', size = rel(3)))+ 
  ylab("Prediction MSE")+
  labs(tag = "D", title = "Hybrid II")

library(patchwork)
tiff("/net/mulan/disk2/yasheng/figures/simulation/r2_laplace_herit0.1.tiff", height = 1200, width = 1600)
(r2_SLM_plt1 | r2_LMM_plt1) /
  (r2_SLMM_plt11 | r2_SLMM_plt21)
dev.off()
tiff("/net/mulan/disk2/yasheng/figures/simulation/r2_laplace_herit0.2.tiff", height = 1200, width = 1600)
(r2_SLM_plt2 | r2_LMM_plt2) /
  (r2_SLMM_plt12 | r2_SLMM_plt22)
dev.off()
tiff("/net/mulan/disk2/yasheng/figures/simulation/r2_laplace_herit0.5.tiff", height = 1200, width = 1600)
(r2_SLM_plt3 | r2_LMM_plt3) /
  (r2_SLMM_plt13 | r2_SLMM_plt23)
dev.off()

tiff("/net/mulan/disk2/yasheng/figures/simulation/mse_laplace.tiff", height = 1600, width = 2000)
(mse_SLM_plt | mse_LMM_plt) /
  (mse_SLMM_plt1 | mse_SLMM_plt2)
dev.off()