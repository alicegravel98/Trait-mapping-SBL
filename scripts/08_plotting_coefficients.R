### Plot PLSR coefficients

## Load libraries----
library(tidyverse)
library(patchwork)
library(ggpubr) #annotate_figure()
library(egg)
library(grid) #textGrob()
library(RColorBrewer)

## Set working directory----
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Load functions----
source("scripts/00_useful_functions.R")


## Import data----
all.jack.coefs.list <- readRDS("output_plsr/all_jack_coefs_list.rds")
wv <- readRDS("output_plsr/wavelength.rds")

## Set up data to plot----
all.jack.coefs.list.summ <-lapply(all.jack.coefs.list, function(coef.list) {
  coef.mat <- matrix(unlist(coef.list), nrow = length(coef.list), byrow = T)
  return(coef.mat)
})

coef.plot.df.list<-lapply(all.jack.coefs.list.summ,function(coef.mat){
  coef.means <- colMeans(coef.mat)
  coef.025 <- apply(coef.mat, 2, quantile, probs = 0.025)
  coef.975 <- apply(coef.mat, 2, quantile, probs = 0.975)
  coef.summ.df <- data.frame(mean = coef.means,
                           per.025 = coef.025,
                           per.975 = coef.975)
  return(coef.summ.df)
})

## Plot coefficients----

# SLA----
SLA_coefs <- ggplot(coef.plot.df.list$SLA[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("SLA")

SLA_coefs

# LMA----
LMA_coefs <- ggplot(coef.plot.df.list$LMA[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("LMA")

LMA_coefs

# LDMC----
LDMC_coefs <- ggplot(coef.plot.df.list$LDMC[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("LDMC")

LDMC_coefs

# EWT----
EWT_coefs <- ggplot(coef.plot.df.list$EWT[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("EWT")

EWT_coefs

# Hemicellulose----
hemi_coefs <- ggplot(coef.plot.df.list$Hemi[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("Hemicellulose")

hemi_coefs

# Cellulose----
cell_coefs <- ggplot(coef.plot.df.list$Cell[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("Cellulose")

cell_coefs

# Lignin----
lignin_coefs <- ggplot(coef.plot.df.list$Lignin[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("Lignin")

lignin_coefs

# Chlorophyll a----
chla_coefs <- ggplot(coef.plot.df.list$"Chl a"[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle(expression(Chl ~ italic(a)))

chla_coefs

# Chlorophyll b----
chlb_coefs <- ggplot(coef.plot.df.list$"Chl b"[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle(expression(Chl ~ italic(b)))

chlb_coefs

# Carotenoids----
carot_coefs <- ggplot(coef.plot.df.list$Car[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("Carotenoids")

carot_coefs

# Carbon----
carbon_coefs <- ggplot(coef.plot.df.list$Carbon[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("Carbon")

carbon_coefs

# Nitrogen----
nitrogen_coefs <- ggplot(coef.plot.df.list$Nitrogen[-1, ])+
  geom_ribbon(aes(x = wv,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "#1f78b4") +
  geom_line(aes(x = wv,
                y = mean),
            linewidth = 1,color = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey33", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(390,2410)) +
  ggtitle("Nitrogen")

nitrogen_coefs

## Save plots----
coefs <- egg::ggarrange(plots = list(LMA_coefs, EWT_coefs, SLA_coefs,
                               carot_coefs, chla_coefs, nitrogen_coefs,
                               chlb_coefs, cell_coefs, lignin_coefs,
                               carbon_coefs, LDMC_coefs, hemi_coefs),
               ncol = 2, nrow = 6)
coefs_axis <- annotate_figure(coefs, left = textGrob("PLSR coefficient", rot = 90, vjust = 1, gp = gpar(cex = 2)),
                bottom = textGrob("Wavelength (nm)", gp = gpar(cex = 2)))

pdf("output_plsr/PLSR_coef_plot.pdf", width = 18, height = 20, onefile = F)
coefs_axis
dev.off()

png("output_plsr/PLSR_coef_plot.png", width = 18, height = 20, units = "in", res = 300)
coefs_axis
dev.off()
