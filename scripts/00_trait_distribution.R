### Trait distribution analysis

## Load libraries----
library(tidyverse)
library(terra)
library(raster)
library(egg)
library(ggpubr) #function get_legend()

## Set working directory----
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Import data----
mosaic_pred <- rast("mapping/flg01/pred_traits01.tif") #trait predictions
band_names <- c("SLA", "LMA", "LDMC", "EWT", "hemicellulose", "cellulose",
                "lignin", "chla", "chlb", "carotenoid", "carbon", "nitrogen") #define band names
names(mosaic_pred) <- band_names
names(mosaic_pred) #check results
twi_clip <- rast("C:/Users/alice/Documents/ArcGIS/TWI/Resample_TWI_merged_Clip.tif") #TWI values
twi <- as.data.frame(twi_clip, xy = T, na.rm = NA) #change raster as dataframe
pred_traits <- as.data.frame(mosaic_pred, xy = T, na.rm = NA)

rm(mosaic_pred, twi_clip)
gc()

join <- right_join(pred_traits, twi, by = c("x","y"))
join <- join[rowSums(is.na(join[, 3:14])) < 12, ] #remove rows with traits as complete NAs

join$Milieux <- ifelse(join$Resample_TWI_merged_Clip > 8, "hydriques",
                             ifelse(join$Resample_TWI_merged_Clip <= 1.5, "xériques", "mésiques"))

## Trait distribution----

dist_SLA <- ggplot(join, aes(x = SLA)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(SLA, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = expression("SLA (m"^2*" kg"^-1*")"), y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none",
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26))
dist_SLA

dist_LMA <- ggplot(join, aes(x = LMA)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(LMA, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = expression("LMA (g m"^-2*")"), y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_LMA

dist_LDMC <- ggplot(join, aes(x = LDMC)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(LDMC, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = expression("LDMC (mg g"^-1*")"), y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_LDMC

dist_EWT <- ggplot(join, aes(x = EWT)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(EWT, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = "EWT (mm)", y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_EWT

dist_hemi <- ggplot(join, aes(x = hemicellulose)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(hemicellulose, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = "Hémicellulose (%)", y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_hemi

dist_cell <- ggplot(join, aes(x = cellulose)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(cellulose, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = "Cellulose (%)", y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_cell

dist_lignin <- ggplot(join, aes(x = lignin)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(lignin, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = "Lignine (%)", y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_lignin

dist_chla <- ggplot(join, aes(x = chla)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(chla, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = expression("Chl" ~ italic(a) ~ "(mg g"^-1*")"), y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_chla

dist_chlb <- ggplot(join, aes(x = chlb)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(chlb, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = expression("Chl" ~ italic(b) ~ "(mg g"^-1*")"), y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_chlb

dist_carot <- ggplot(join, aes(x = carotenoid)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(carotenoid, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = expression("Car (mg g"^-1*")"), y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_carot

dist_carbon <- ggplot(join, aes(x = carbon)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(carbon, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = "Carbone (%)", y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_carbon

dist_nitrogen <- ggplot(join, aes(x = nitrogen)) +
  geom_density(aes(color = Milieux), linewidth = 1.5) +
  geom_vline(aes(xintercept = mean(nitrogen, na.rm = T)), color = "#414141", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(values = c("#1f77b4", "#2ca02c", "#d62728")) +
  labs(x = "Azote (%)", y = "Densité") +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        legend.position = "none")
dist_nitrogen

## Save plots----
common_legend <- ggpubr::get_legend(dist_SLA, position = "bottom")

all_plots <- egg::ggarrange(plots = list(dist_nitrogen, dist_carot, dist_carbon,
                                         dist_cell, dist_chla, dist_chlb,
                                         dist_EWT, dist_hemi, dist_LDMC,
                                         dist_lignin, dist_LMA, dist_SLA),
                            ncol = 3, nrow = 4, bottom = common_legend)

pdf("output_topography/trait_distribution.pdf", width = 18, height = 18, onefile = F)
all_plots
dev.off()

png("output_topography/trait_distribution.png", width = 18, height = 18, units = "in", res = 300)
all_plots
dev.off()