### Build PLSR models

## Install all packages and load libraries----
install.packages("reshape2")
library(tidyverse)
library(spectrolab)
library(pls)
library(reshape2) #melt function
library(egg)

## Set working directory----
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Load functions----
# Download VIP.R from Bjorn-Helge Mevik's website: https://mevik.net/work/software/pls.html
source("scripts/00_VIP.R")
source("scripts/00_useful_functions.R")

## Import train/test data----
ref.train.nonbn <- readRDS("output_plsr/ref_train.rds")
ref.test.nonbn <- readRDS("output_plsr/ref_test.rds")
wv <- readRDS("output_plsr/wavelength.rds")

## Vector-normalization of spectra----
ref.train <- t(apply(as.matrix(ref.train.nonbn), 1, function(x) x/sqrt(sum(x^2))))
ref.train <- spectra(ref.train, bands = bands(ref.train.nonbn), names = names(ref.train.nonbn))
meta(ref.train) <- meta(ref.train.nonbn)

ref.test <- t(apply(as.matrix(ref.test.nonbn), 1, function(x) x/sqrt(sum(x^2))))
ref.test <- spectra(ref.test, bands = bands(ref.test.nonbn), names = names(ref.test.nonbn))
meta(ref.test) <- meta(ref.test.nonbn)

plot_interactive(ref.train)
plot_interactive(ref.test)

## Building first-pass models using K-fold cross validation to get the optimal number of components----
# Specific leaf area (m2/kg)
SLA_CVmodel<-plsr(meta(ref.train)$SLA~as.matrix(ref.train),
                  ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_SLA_CVmodel <- selectNcomp(SLA_CVmodel, method = "onesigma", plot = T)
SLA_valid <- which(!is.na(meta(ref.train)$SLA))
SLA_pred <- data.frame(ID = meta(ref.train)$sample_id[SLA_valid],
                     Species = meta(ref.train)$specie[SLA_valid],
                     measured = meta(ref.train)$SLA[SLA_valid],
                     val_pred = SLA_CVmodel$validation$pred[, , ncomp_SLA_CVmodel])
SLA_pred$measured <- as.numeric(SLA_pred$measured)
ggplot(SLA_pred, aes(y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(0, 20),ylim = c(0, 20)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting SLA from canopy spectra")

# Leaf mass per area (g/m2)
LMA_CVmodel <- plsr(meta(ref.train)$LMA~as.matrix(ref.train),
                  ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_LMA_CVmodel <- selectNcomp(LMA_CVmodel, method = "onesigma", plot = T)
LMA_valid <- which(!is.na(meta(ref.train)$LMA))
LMA_pred <- data.frame(ID = meta(ref.train)$sample_id[LMA_valid],
                     Species = meta(ref.train)$specie[LMA_valid],
                     measured = meta(ref.train)$LMA[LMA_valid],
                     val_pred = LMA_CVmodel$validation$pred[, , ncomp_LMA_CVmodel])
LMA_pred$measured <- as.numeric(LMA_pred$measured)
ggplot(LMA_pred, aes(y = measured, x = val_pred))+
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1,intercept = 0,linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(0, 400), ylim = c(0, 400)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting LMA from canopy spectra")

# Leaf dry matter content (mg/g)
LDMC_CVmodel <- plsr(meta(ref.train)$LDMC~as.matrix(ref.train),
                   ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_LDMC_CVmodel <- selectNcomp(LDMC_CVmodel, method = "onesigma", plot = T) #returns 0 (we have to force one)
ncomp_LDMC_CVmodel <- 8 #chose the second minimum RMSEP
LDMC_valid <- which(!is.na(meta(ref.train)$LDMC))
LDMC_pred <- data.frame(ID = meta(ref.train)$sample_id[LDMC_valid],
                      Species = meta(ref.train)$specie[LDMC_valid],
                      measured = meta(ref.train)$LDMC[LDMC_valid],
                      val_pred = LDMC_CVmodel$validation$pred[, , ncomp_LDMC_CVmodel])
LDMC_pred$measured <- as.numeric(LDMC_pred$measured)
ggplot(LDMC_pred, aes(y = measured, x = val_pred))+
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept=0, linetype="dashed", linewidth=2) +
  coord_cartesian(xlim = c(375, 600), ylim = c(300, 600)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting LDMC from canopy spectra")

# Equivalent water thickness (mm)
EWT_CVmodel <- plsr(meta(ref.train)$EWT~as.matrix(ref.train),
                  ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_EWT_CVmodel <- selectNcomp(EWT_CVmodel, method = "onesigma", plot = T)
EWT_valid <- which(!is.na(meta(ref.train)$EWT))
EWT_pred <- data.frame(ID = meta(ref.train)$sample_id[EWT_valid],
                     Species = meta(ref.train)$specie[EWT_valid],
                     measured = meta(ref.train)$EWT[EWT_valid],
                     val_pred = EWT_CVmodel$validation$pred[, , ncomp_EWT_CVmodel])
EWT_pred$measured <- as.numeric(EWT_pred$measured)
ggplot(EWT_pred, aes(y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting EWT from canopy spectra")

# Hemicellulose (%)
hemi_CVmodel <- plsr(meta(ref.train)$hemicellulose~as.matrix(ref.train),
                   ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_hemi_CVmodel <- selectNcomp(hemi_CVmodel, method = "onesigma", plot = T) #returns 0 (we have to force one)
ncomp_hemi_CVmodel <- 5 #chose near the second minimum RMSEP
hemi_valid <- which(!is.na(meta(ref.train)$hemicellulose))
hemi_pred <- data.frame(ID = meta(ref.train)$sample_id[hemi_valid],
                      Species = meta(ref.train)$specie[hemi_valid],
                      measured = meta(ref.train)$hemicellulose[hemi_valid],
                      val_pred = hemi_CVmodel$validation$pred[, , ncomp_hemi_CVmodel])
hemi_pred$measured <- as.numeric(hemi_pred$measured)
ggplot(hemi_pred, aes( y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(0, 25), ylim = c(0, 25)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting hemicellulose from canopy spectra")

# Cellulose (%)
cell_CVmodel <- plsr(meta(ref.train)$cellulose~as.matrix(ref.train),
                   ncomp = 30, method = "oscorespls", validation = "CV", segments = 10) #number too small
ncomp_cell_CVmodel <- selectNcomp(cell_CVmodel, method = "onesigma", plot = T)
ncomp_cell_CVmodel <- 5 #chose the second minimum RMSEP
cell_valid <- which(!is.na(meta(ref.train)$cellulose))
cell_pred <- data.frame(ID = meta(ref.train)$sample_id[cell_valid],
                      Species = meta(ref.train)$specie[cell_valid],
                      measured = meta(ref.train)$cellulose[cell_valid],
                      val_pred = cell_CVmodel$validation$pred[, , ncomp_cell_CVmodel])
cell_pred$measured <- as.numeric(cell_pred$measured)
ggplot(cell_pred, aes(y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F)+
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2)+
  coord_cartesian(xlim = c(6, 16), ylim = c(6, 16)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting cellulose from canopy spectra")

# Lignin (%)
lignin_CVmodel <- plsr(meta(ref.train)$lignin~as.matrix(ref.train),
                     ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_lignin_CVmodel <- selectNcomp(lignin_CVmodel, method = "onesigma", plot = T)
lignin_valid <- which(!is.na(meta(ref.train)$lignin))
lignin_pred <- data.frame(ID = meta(ref.train)$sample_id[lignin_valid],
                        Species = meta(ref.train)$specie[lignin_valid],
                        measured = meta(ref.train)$lignin[lignin_valid],
                        val_pred = lignin_CVmodel$validation$pred[, , ncomp_lignin_CVmodel])
lignin_pred$measured <- as.numeric(lignin_pred$measured)
ggplot(lignin_pred, aes(y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(0, 18), ylim = c(0, 18)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2))+
  labs(y = "Measured", x = "Predicted")+
  ggtitle("Predicting lignin from canopy spectra")

# Chlorophyll a (mg/g)
chla_CVmodel <- plsr(meta(ref.train)$chla~as.matrix(ref.train),
                   ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_chla_CVmodel <- selectNcomp(chla_CVmodel, method = "onesigma", plot = T)
chla_valid <- which(!is.na(meta(ref.train)$chla))
chla_pred <- data.frame(ID = meta(ref.train)$sample_id[chla_valid],
                      Species = meta(ref.train)$specie[chla_valid],
                      measured = meta(ref.train)$chla[chla_valid],
                      val_pred = chla_CVmodel$validation$pred[, , ncomp_chla_CVmodel])
chla_pred$measured <- as.numeric(chla_pred$measured)
ggplot(chla_pred, aes(y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting chlorophyll a from canopy spectra")

# Chlorophyll b (mg/g)
chlb_CVmodel <- plsr(meta(ref.train)$chlb~as.matrix(ref.train),
                   ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_chlb_CVmodel <- selectNcomp(chlb_CVmodel, method = "onesigma", plot = T) #number too small
ncomp_chlb_CVmodel <- 6 #chose the absolute minimum RMSEP
chlb_valid <- which(!is.na(meta(ref.train)$chlb))
chlb_pred <- data.frame(ID = meta(ref.train)$sample_id[chlb_valid],
                      Species = meta(ref.train)$specie[chlb_valid],
                      measured = meta(ref.train)$chlb[chlb_valid],
                      val_pred = chlb_CVmodel$validation$pred[, , ncomp_chlb_CVmodel])
chlb_pred$measured <- as.numeric(chlb_pred$measured)
ggplot(chlb_pred,aes(y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(0, 0.9), ylim = c(0, 0.9)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting chlorophyll b from canopy spectra")

# Carotenoids (mg/g)
carot_CVmodel <- plsr(meta(ref.train)$carotenoids~as.matrix(ref.train),
                    ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_carot_CVmodel <- selectNcomp(carot_CVmodel, method = "onesigma", plot = T)
carot_valid <- which(!is.na(meta(ref.train)$carotenoids))
carot_pred <- data.frame(ID = meta(ref.train)$sample_id[carot_valid],
                       Species = meta(ref.train)$specie[carot_valid],
                       measured = meta(ref.train)$carotenoids[carot_valid],
                       val_pred = carot_CVmodel$validation$pred[, , ncomp_carot_CVmodel])
carot_pred$measured <- as.numeric(carot_pred$measured)
ggplot(carot_pred, aes(y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm",se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(0, 0.7), ylim = c(0, 0.7)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting carotenoids a from canopy spectra")

# Carbon (%)
carbon_CVmodel <- plsr(meta(ref.train)$carbon~as.matrix(ref.train),
                     ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_carbon_CVmodel <- selectNcomp(carbon_CVmodel, method = "onesigma", plot = T) #number too small
ncomp_carbon_CVmodel <- 7 #chose near the second minimum RMSEP
carbon_valid <- which(!is.na(meta(ref.train)$carbon))
carbon_pred <- data.frame(ID = meta(ref.train)$sample_id[carbon_valid],
                        Species = meta(ref.train)$specie[carbon_valid],
                        measured = meta(ref.train)$carbon[carbon_valid],
                        val_pred = carbon_CVmodel$validation$pred[, , ncomp_carbon_CVmodel])
carbon_pred$measured <- as.numeric(carbon_pred$measured)
ggplot(carbon_pred, aes(y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(47, 51), ylim = c(47, 51)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting carbon a from canopy spectra")

# Nitrogen (%)
nitrogen_CVmodel <- plsr(meta(ref.train)$nitrogen~as.matrix(ref.train),
                       ncomp = 30, method = "oscorespls", validation = "CV", segments = 10)
ncomp_nitrogen_CVmodel <- selectNcomp(nitrogen_CVmodel, method = "onesigma", plot = T)
nitrogen_valid <- which(!is.na(meta(ref.train)$nitrogen))
nitrogen_pred<-data.frame(ID = meta(ref.train)$sample_id[nitrogen_valid],
                          Species = meta(ref.train)$specie[nitrogen_valid],
                          measured = meta(ref.train)$nitrogen[nitrogen_valid],
                          val_pred = nitrogen_CVmodel$validation$pred[, , ncomp_nitrogen_CVmodel])
nitrogen_pred$measured <- as.numeric(nitrogen_pred$measured)
ggplot(nitrogen_pred, aes(y = measured, x = val_pred)) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8, 0.2)) +
  labs(y = "Measured", x = "Predicted") +
  ggtitle("Predicting nitrogen a from canopy spectra")

## VIP from calibration models----
# Extract VIP values
VIP1.df <- data.frame(EWT = VIP(EWT_CVmodel)[ncomp_EWT_CVmodel, ],
                      LDMC = VIP(LDMC_CVmodel)[ncomp_LDMC_CVmodel, ],
                      LMA = VIP(LMA_CVmodel)[ncomp_LMA_CVmodel, ],
                      SLA = VIP(SLA_CVmodel)[ncomp_SLA_CVmodel, ],
                    wavelength = wv)

rownames(VIP1.df) <- 1:220 #change rownames

VIP2.df <- data.frame(cell = VIP(cell_CVmodel)[ncomp_cell_CVmodel, ],
                    hemi = VIP(hemi_CVmodel)[ncomp_hemi_CVmodel, ],
                    lignin = VIP(lignin_CVmodel)[ncomp_lignin_CVmodel, ],
                    wavelength = wv)


rownames(VIP2.df) <- 1:220

VIP3.df <- data.frame(car = VIP(carot_CVmodel)[ncomp_carot_CVmodel, ],
                      chla = VIP(chla_CVmodel)[ncomp_chla_CVmodel, ],
                    chlb = VIP(chlb_CVmodel)[ncomp_chlb_CVmodel, ],
                    wavelength = wv)

rownames(VIP3.df) <- 1:220

VIP4.df <- data.frame(carbon = VIP(carbon_CVmodel)[ncomp_carbon_CVmodel, ],
                      nitrogen = VIP(nitrogen_CVmodel)[ncomp_nitrogen_CVmodel, ],
                      wavelength = wv)

rownames(VIP4.df) <- 1:220

saveRDS(VIP1.df, "output_plsr/VIP1_df.rds")
saveRDS(VIP2.df, "output_plsr/VIP2_df.rds")
saveRDS(VIP3.df, "output_plsr/VIP3_df.rds")
saveRDS(VIP4.df, "output_plsr/VIP4_df.rds")

VIP1.long <- melt(VIP1.df, id.vars = "wavelength")
VIP2.long <- melt(VIP2.df, id.vars = "wavelength")
VIP3.long <- melt(VIP3.df, id.vars = "wavelength")
VIP4.long <- melt(VIP4.df, id.vars = "wavelength")

# Choose a color palette
focal_palette <- c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd" )

# Plots
VIP1.plot <- ggplot(VIP1.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 0.8) + theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4)) +
  scale_color_manual(values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")
VIP1.plot

VIP2.plot <- ggplot(VIP2.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 0.8) + theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0,4)) +
  scale_color_manual(labels = c("Cell", "Hemi", "Lignin"),
                     values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")
VIP2.plot

VIP3.plot <- ggplot(VIP3.long,aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 0.8) + theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4)) +
  scale_color_manual(labels = expression("Car", Chl ~ italic(a), Chl ~ italic(b)),
                     values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")
VIP3.plot

VIP4.plot <- ggplot(VIP4.long,aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 0.8) + theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4)) +
  scale_color_manual(labels = c("Carbon", "Nitrogen"),
                     values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")
VIP4.plot

# Arrange all plots together
VIP_all <- egg::ggarrange(plots = list(VIP1.plot, VIP2.plot, VIP3.plot, VIP4.plot),
                          ncol = 1, nrow = 4, labels = c("(A)", "(B)", "(C)", "(D)"))

# Save plots
pdf("output_plsr/VIPs_all.pdf", height = 15, width = 10, onefile = F)
VIP_all
dev.off()

png("output_plsr/VIPs_all.png", width = 12, height = 10, units = "in", res = 300)
VIP_all
dev.off()

## Jackknife analyses----
# Change traits of interest to numeric
meta(ref.train)$SLA <- as.numeric(meta(ref.train)$SLA)
meta(ref.test)$SLA <- as.numeric(meta(ref.test)$SLA)
meta(ref.train)$LMA <- as.numeric(meta(ref.train)$LMA)
meta(ref.test)$LMA <- as.numeric(meta(ref.test)$LMA)
meta(ref.train)$LDMC <- as.numeric(meta(ref.train)$LDMC)
meta(ref.test)$LDMC <- as.numeric(meta(ref.test)$LDMC)
meta(ref.train)$EWT <- as.numeric(meta(ref.train)$EWT)
meta(ref.test)$EWT <- as.numeric(meta(ref.test)$EWT)
meta(ref.train)$hemi <- as.numeric(meta(ref.train)$hemi)
meta(ref.test)$hemi <- as.numeric(meta(ref.test)$hemi)
meta(ref.train)$cellulose <- as.numeric(meta(ref.train)$cellulose)
meta(ref.test)$cellulose <- as.numeric(meta(ref.test)$cellulose)
meta(ref.train)$lignin <- as.numeric(meta(ref.train)$lignin)
meta(ref.test)$lignin <- as.numeric(meta(ref.test)$lignin)
meta(ref.train)$chla <- as.numeric(meta(ref.train)$chla)
meta(ref.test)$chla <- as.numeric(meta(ref.test)$chla)
meta(ref.train)$chlb <- as.numeric(meta(ref.train)$chlb)
meta(ref.test)$chlb <- as.numeric(meta(ref.test)$chlb)
meta(ref.train)$carotenoids <- as.numeric(meta(ref.train)$carotenoids)
meta(ref.test)$carotenoids <- as.numeric(meta(ref.test)$carotenoids)
meta(ref.train)$carbon <- as.numeric(meta(ref.train)$carbon)
meta(ref.test)$carbon <- as.numeric(meta(ref.test)$carbon)
meta(ref.train)$nitrogen <- as.numeric(meta(ref.train)$nitrogen)
meta(ref.test)$nitrogen <- as.numeric(meta(ref.test)$nitrogen)

SLA.jack.coefs <- list()
LMA.jack.coefs <- list()
LDMC.jack.coefs <- list()
EWT.jack.coefs <- list()
hemi.jack.coefs <- list()
cell.jack.coefs <- list()
lignin.jack.coefs <- list()
chla.jack.coefs <- list()
chlb.jack.coefs <- list()
carot.jack.coefs <- list()
carbon.jack.coefs <- list()
nitrogen.jack.coefs <- list()

SLA.jack.stats <- list()
LMA.jack.stats <- list()
LDMC.jack.stats <- list()
EWT.jack.stats <- list()
hemi.jack.stats <- list()
cell.jack.stats <- list()
lignin.jack.stats <- list()
chla.jack.stats <- list()
chlb.jack.stats <- list()
carot.jack.stats <- list()
carbon.jack.stats <- list()
nitrogen.jack.stats <- list()

nreps <- 100

for(i in 1:nreps){
  print(i)
  
  n.cal.spec <- nrow(ref.train)
  train.jack <- sample(1:n.cal.spec, floor(0.7*n.cal.spec))
  test.jack <- setdiff(1:n.cal.spec, train.jack)
  
  calib.jack <- ref.train[train.jack]
  val.jack <- ref.train[test.jack]
  
  SLA.jack <- plsr(meta(calib.jack)$SLA~as.matrix(calib.jack),
                 ncomp = 30, method = "oscorespls", validation = "none")
  LMA.jack <- plsr(meta(calib.jack)$LMA~as.matrix(calib.jack),
                 ncomp = 30, method = "oscorespls", validation = "none")
  LDMC.jack <- plsr(meta(calib.jack)$LDMC~as.matrix(calib.jack),
                  ncomp = 30, method = "oscorespls", validation = "none")
  EWT.jack <- plsr(meta(calib.jack)$EWT~as.matrix(calib.jack),
                 ncomp = 30, method = "oscorespls", validation = "none")
  hemi.jack <- plsr(meta(calib.jack)$hemi~as.matrix(calib.jack),
                    ncomp = 30, method = "oscorespls", validation = "none")
  cell.jack <- plsr(meta(calib.jack)$cellulose~as.matrix(calib.jack),
                  ncomp = 30, method = "oscorespls", validation = "none")
  lignin.jack <- plsr(meta(calib.jack)$lignin~as.matrix(calib.jack),
                    ncomp = 30, method = "oscorespls", validation = "none")
  chla.jack <- plsr(meta(calib.jack)$chla~as.matrix(calib.jack),
                  ncomp = 30,method = "oscorespls", validation = "none")
  chlb.jack <- plsr(meta(calib.jack)$chlb~as.matrix(calib.jack),
                  ncomp = 30, method = "oscorespls", validation = "none")
  carot.jack <- plsr(meta(calib.jack)$carotenoids~as.matrix(calib.jack),
                   ncomp = 30, method = "oscorespls", validation = "none")
  carbon.jack <- plsr(meta(calib.jack)$carbon~as.matrix(calib.jack),
                    ncomp = 30, method = "oscorespls", validation = "none")
  nitrogen.jack <- plsr(meta(calib.jack)$nitrogen~as.matrix(calib.jack),
                      ncomp = 30, method = "oscorespls", validation = "none")
  
  SLA.jack.val.pred <- as.vector(predict(SLA.jack, newdata = as.matrix(val.jack), ncomp = ncomp_SLA_CVmodel)[, , 1])
  SLA.jack.val.fit <- lm(SLA.jack.val.pred~meta(val.jack)$SLA)
  SLA.jack.stats[[i]] <- c(R2 = summary(SLA.jack.val.fit)$r.squared,
                         RMSE = RMSD(meta(val.jack)$SLA, SLA.jack.val.pred),
                         perRMSE = percentRMSD(meta(val.jack)$SLA, SLA.jack.val.pred, 0.025, 0.975),
                         bias = mean(SLA.jack.val.pred, na.rm = T) - mean(meta(val.jack)$SLA, na.rm = T))
  
  LMA.jack.val.pred <- as.vector(predict(LMA.jack,newdata = as.matrix(val.jack), ncomp = ncomp_LMA_CVmodel)[, , 1])
  LMA.jack.val.fit <- lm(LMA.jack.val.pred~meta(val.jack)$LMA)
  LMA.jack.stats[[i]] <- c(R2 = summary(LMA.jack.val.fit)$r.squared,
                         RMSE = RMSD(meta(val.jack)$LMA, LMA.jack.val.pred),
                         perRMSE = percentRMSD(meta(val.jack)$LMA, LMA.jack.val.pred, 0.025, 0.975),
                         bias = mean(LMA.jack.val.pred, na.rm = T) - mean(meta(val.jack)$LMA, na.rm = T))
  
  LDMC.jack.val.pred <- as.vector(predict(LDMC.jack, newdata = as.matrix(val.jack), ncomp = ncomp_LDMC_CVmodel)[, , 1])
  LDMC.jack.val.fit <- lm(LDMC.jack.val.pred~meta(val.jack)$LDMC)
  LDMC.jack.stats[[i]] <- c(R2 = summary(LDMC.jack.val.fit)$r.squared,
                          RMSE = RMSD(meta(val.jack)$LDMC, LDMC.jack.val.pred),
                          perRMSE = percentRMSD(meta(val.jack)$LDMC, LDMC.jack.val.pred, 0.025, 0.975),
                          bias = mean(LDMC.jack.val.pred, na.rm = T) - mean(meta(val.jack)$LDMC, na.rm = T))
  
  EWT.jack.val.pred <- as.vector(predict(EWT.jack, newdata = as.matrix(val.jack), ncomp = ncomp_EWT_CVmodel)[, , 1])
  EWT.jack.val.fit <- lm(EWT.jack.val.pred~meta(val.jack)$EWT)
  EWT.jack.stats[[i]] <- c(R2 = summary(EWT.jack.val.fit)$r.squared,
                         RMSE = RMSD(meta(val.jack)$EWT, EWT.jack.val.pred),
                         perRMSE = percentRMSD(meta(val.jack)$EWT, EWT.jack.val.pred, 0.025, 0.975),
                         bias = mean(EWT.jack.val.pred,na.rm = T) - mean(meta(val.jack)$EWT, na.rm = T))
  
  hemi.jack.val.pred <- as.vector(predict(hemi.jack, newdata = as.matrix(val.jack), ncomp = ncomp_hemi_CVmodel)[, , 1])
  hemi.jack.val.fit <- lm(hemi.jack.val.pred~meta(val.jack)$hemi)
  hemi.jack.stats[[i]] <- c(R2 = summary(hemi.jack.val.fit)$r.squared,
                            RMSE = RMSD(meta(val.jack)$hemi, hemi.jack.val.pred),
                            perRMSE = percentRMSD(meta(val.jack)$hemi, hemi.jack.val.pred, 0.025, 0.975),
                            bias = mean(hemi.jack.val.pred, na.rm = T) - mean(meta(val.jack)$hemi, na.rm = T))
  
  cell.jack.val.pred <- as.vector(predict(cell.jack, newdata = as.matrix(val.jack), ncomp = ncomp_cell_CVmodel)[, , 1])
  cell.jack.val.fit <- lm(cell.jack.val.pred~meta(val.jack)$cellulose)
  cell.jack.stats[[i]] <- c(R2 = summary(cell.jack.val.fit)$r.squared,
                          RMSE = RMSD(meta(val.jack)$cellulose, cell.jack.val.pred),
                          perRMSE = percentRMSD(meta(val.jack)$cellulose, cell.jack.val.pred, 0.025, 0.975),
                          bias = mean(cell.jack.val.pred, na.rm = T) - mean(meta(val.jack)$cellulose, na.rm = T))
  
  lignin.jack.val.pred <- as.vector(predict(lignin.jack, newdata = as.matrix(val.jack), ncomp = ncomp_lignin_CVmodel)[, , 1])
  lignin.jack.val.fit <- lm(lignin.jack.val.pred~meta(val.jack)$lignin)
  lignin.jack.stats[[i]] <- c(R2 = summary(lignin.jack.val.fit)$r.squared,
                            RMSE = RMSD(meta(val.jack)$lignin, lignin.jack.val.pred),
                            perRMSE = percentRMSD(meta(val.jack)$lignin,lignin.jack.val.pred, 0.025, 0.975),
                            bias = mean(lignin.jack.val.pred, na.rm = T) - mean(meta(val.jack)$lignin, na.rm = T))
  
  chla.jack.val.pred <- as.vector(predict(chla.jack, newdata = as.matrix(val.jack), ncomp = ncomp_chla_CVmodel)[, , 1])
  chla.jack.val.fit <- lm(chla.jack.val.pred~meta(val.jack)$chla)
  chla.jack.stats[[i]] <- c(R2 = summary(chla.jack.val.fit)$r.squared,
                          RMSE = RMSD(meta(val.jack)$chla, chla.jack.val.pred),
                          perRMSE = percentRMSD(meta(val.jack)$chla, chla.jack.val.pred, 0.025, 0.975),
                          bias = mean(chla.jack.val.pred, na.rm = T) - mean(meta(val.jack)$chla, na.rm = T))
  
  chlb.jack.val.pred <- as.vector(predict(chlb.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlb_CVmodel)[,,1])
  chlb.jack.val.fit <- lm(chlb.jack.val.pred~meta(val.jack)$chlb)
  chlb.jack.stats[[i]] <- c(R2 = summary(chlb.jack.val.fit)$r.squared,
                          RMSE = RMSD(meta(val.jack)$chlb,chlb.jack.val.pred),
                          perRMSE = percentRMSD(meta(val.jack)$chlb, chlb.jack.val.pred, 0.025, 0.975),
                          bias = mean(chlb.jack.val.pred, na.rm = T) - mean(meta(val.jack)$chlb, na.rm = T))
  
  carot.jack.val.pred <- as.vector(predict(carot.jack,newdata = as.matrix(val.jack), ncomp = ncomp_carot_CVmodel)[, , 1])
  carot.jack.val.fit <- lm(carot.jack.val.pred~meta(val.jack)$carotenoids)
  carot.jack.stats[[i]] <- c(R2 = summary(carot.jack.val.fit)$r.squared,
                           RMSE = RMSD(meta(val.jack)$carotenoids,carot.jack.val.pred),
                           perRMSE = percentRMSD(meta(val.jack)$carotenoids,carot.jack.val.pred, 0.025, 0.975),
                           bias = mean(carot.jack.val.pred, na.rm = T) - mean(meta(val.jack)$carotenoids, na.rm = T))
  
  carbon.jack.val.pred <- as.vector(predict(carbon.jack, newdata = as.matrix(val.jack), ncomp = ncomp_carbon_CVmodel)[, , 1])
  carbon.jack.val.fit <- lm(carbon.jack.val.pred~meta(val.jack)$carbon)
  carbon.jack.stats[[i]] <- c(R2 = summary(carbon.jack.val.fit)$r.squared,
                            RMSE = RMSD(meta(val.jack)$carbon,carbon.jack.val.pred),
                            perRMSE = percentRMSD(meta(val.jack)$carbon, carbon.jack.val.pred, 0.025, 0.975),
                            bias = mean(carbon.jack.val.pred, na.rm = T) - mean(meta(val.jack)$carbon, na.rm = T))
  
  nitrogen.jack.val.pred <- as.vector(predict(nitrogen.jack, newdata = as.matrix(val.jack), ncomp = ncomp_nitrogen_CVmodel)[, , 1])
  nitrogen.jack.val.fit <- lm(nitrogen.jack.val.pred~meta(val.jack)$nitrogen)
  nitrogen.jack.stats[[i]] <- c(R2 = summary(nitrogen.jack.val.fit)$r.squared,
                              RMSE = RMSD(meta(val.jack)$nitrogen, nitrogen.jack.val.pred),
                              perRMSE = percentRMSD(meta(val.jack)$nitrogen, nitrogen.jack.val.pred, 0.025, 0.975),
                              bias = mean(nitrogen.jack.val.pred, na.rm = T) - mean(meta(val.jack)$nitrogen, na.rm = T))
  
  SLA.jack.coefs[[i]] <- as.vector(coef(SLA.jack, ncomp = ncomp_SLA_CVmodel, intercept = T))
  LMA.jack.coefs[[i]] <- as.vector(coef(LMA.jack, ncomp = ncomp_LMA_CVmodel, intercept = T))
  LDMC.jack.coefs[[i]] <- as.vector(coef(LDMC.jack, ncomp = ncomp_LDMC_CVmodel, intercept = T))
  EWT.jack.coefs[[i]] <- as.vector(coef(EWT.jack, ncomp = ncomp_EWT_CVmodel, intercept = T))
  hemi.jack.coefs[[i]] <- as.vector(coef(hemi.jack, ncomp = ncomp_hemi_CVmodel, intercept = T))
  cell.jack.coefs[[i]] <- as.vector(coef(cell.jack, ncomp = ncomp_cell_CVmodel, intercept = T))
  lignin.jack.coefs[[i]] <- as.vector(coef(lignin.jack, ncomp = ncomp_lignin_CVmodel, intercept = T))
  chla.jack.coefs[[i]] <- as.vector(coef(chla.jack, ncomp = ncomp_chla_CVmodel, intercept = T))
  chlb.jack.coefs[[i]] <- as.vector(coef(chlb.jack,ncomp = ncomp_chlb_CVmodel,intercept = T))
  carot.jack.coefs[[i]] <- as.vector(coef(carot.jack, ncomp = ncomp_carot_CVmodel, intercept = T))
  carbon.jack.coefs[[i]] <- as.vector(coef(carbon.jack, ncomp = ncomp_carbon_CVmodel, intercept = T))
  nitrogen.jack.coefs[[i]] <- as.vector(coef(nitrogen.jack, ncomp=ncomp_nitrogen_CVmodel, intercept = T))
}

SLA.jack.pred <- apply.coefs(SLA.jack.coefs, as.matrix(ref.test))
SLA.jack.stat <- t(apply(SLA.jack.pred ,1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)))))
SLA.jack.df <- data.frame(pred.mean = SLA.jack.stat[, 1],
                        pred.low = SLA.jack.stat[, 2],
                        pred.high = SLA.jack.stat[, 3],
                        Measured = meta(ref.test)$SLA,
                        ncomp = ncomp_SLA_CVmodel,
                        Species = meta(ref.test)$specie,
                        ID = meta(ref.test)$sample_id)

LMA.jack.pred <- apply.coefs(LMA.jack.coefs, as.matrix(ref.test))
LMA.jack.stat <- t(apply(LMA.jack.pred, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)))))
LMA.jack.df <- data.frame(pred.mean = LMA.jack.stat[, 1],
                        pred.low = LMA.jack.stat[, 2],
                        pred.high = LMA.jack.stat[, 3],
                        Measured = meta(ref.test)$LMA,
                        ncomp = ncomp_LMA_CVmodel,
                        Species = meta(ref.test)$specie,
                        ID = meta(ref.test)$sample_id)

LDMC.jack.pred <- apply.coefs(LDMC.jack.coefs, as.matrix(ref.test))
LDMC.jack.stat <- t(apply(LDMC.jack.pred, 1,
                        function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)))))
LDMC.jack.df <- data.frame(pred.mean = LDMC.jack.stat[, 1],
                         pred.low = LDMC.jack.stat[, 2],
                         pred.high = LDMC.jack.stat[, 3],
                         Measured = meta(ref.test)$LDMC,
                         ncomp = ncomp_LDMC_CVmodel,
                         Species = meta(ref.test)$specie,
                         ID = meta(ref.test)$sample_id)

EWT.jack.pred <- apply.coefs(EWT.jack.coefs, as.matrix(ref.test))
EWT.jack.stat <- t(apply(EWT.jack.pred, 1,
                       function(obs) c(mean(obs),quantile(obs, probs = c(0.025, 0.975)))))
EWT.jack.df <- data.frame(pred.mean = EWT.jack.stat[, 1],
                        pred.low = EWT.jack.stat[, 2],
                        pred.high = EWT.jack.stat[,3],
                        Measured = meta(ref.test)$EWT,
                        ncomp = ncomp_EWT_CVmodel,
                        Species = meta(ref.test)$specie,
                        ID = meta(ref.test)$sample_id)

hemi.jack.pred <- apply.coefs(hemi.jack.coefs, as.matrix(ref.test))
hemi.jack.stat <- t(apply(hemi.jack.pred, 1,
                          function(obs) c(mean(obs),quantile(obs, probs = c(0.025, 0.975)))))
hemi.jack.df <- data.frame(pred.mean = hemi.jack.stat[, 1],
                           pred.low = hemi.jack.stat[, 2],
                           pred.high = hemi.jack.stat[, 3],
                           Measured = meta(ref.test)$hemi,
                           ncomp = ncomp_hemi_CVmodel,
                           Species = meta(ref.test)$specie,
                           ID = meta(ref.test)$sample_id)

cell.jack.pred <- apply.coefs(cell.jack.coefs, as.matrix(ref.test))
cell.jack.stat <- t(apply(cell.jack.pred, 1,
                        function(obs) c(mean(obs),quantile(obs, probs = c(0.025, 0.975)))))
cell.jack.df <- data.frame(pred.mean = cell.jack.stat[, 1],
                         pred.low = cell.jack.stat[, 2],
                         pred.high = cell.jack.stat[, 3],
                         Measured = meta(ref.test)$cellulose,
                         ncomp = ncomp_cell_CVmodel,
                         Species = meta(ref.test)$specie,
                         ID = meta(ref.test)$sample_id)

lignin.jack.pred <- apply.coefs(lignin.jack.coefs, as.matrix(ref.test))
lignin.jack.stat <- t(apply(lignin.jack.pred, 1,
                          function(obs) c(mean(obs),quantile(obs, probs = c(0.025, 0.975)))))
lignin.jack.df <- data.frame(pred.mean = lignin.jack.stat[, 1],
                           pred.low = lignin.jack.stat[, 2],
                           pred.high = lignin.jack.stat[, 3],
                           Measured = meta(ref.test)$lignin,
                           ncomp = ncomp_lignin_CVmodel,
                           Species = meta(ref.test)$specie,
                           ID = meta(ref.test)$sample_id)

chla.jack.pred <- apply.coefs(chla.jack.coefs, as.matrix(ref.test))
chla.jack.stat <- t(apply(chla.jack.pred, 1,
                        function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chla.jack.df <- data.frame(pred.mean = chla.jack.stat[, 1],
                         pred.low = chla.jack.stat[, 2],
                         pred.high = chla.jack.stat[, 3],
                         Measured = meta(ref.test)$chla,
                         ncomp = ncomp_chla_CVmodel,
                         Species = meta(ref.test)$specie,
                         ID = meta(ref.test)$sample_id)

chlb.jack.pred <- apply.coefs(chlb.jack.coefs, as.matrix(ref.test))
chlb.jack.stat <- t(apply(chlb.jack.pred, 1,
                        function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)))))
chlb.jack.df <- data.frame(pred.mean = chlb.jack.stat[, 1],
                         pred.low = chlb.jack.stat[, 2],
                         pred.high = chlb.jack.stat[, 3],
                         Measured = meta(ref.test)$chlb,
                         ncomp = ncomp_chlb_CVmodel,
                         Species = meta(ref.test)$specie,
                         ID = meta(ref.test)$sample_id)

carot.jack.pred <- apply.coefs(carot.jack.coefs, as.matrix(ref.test))
carot.jack.stat <- t(apply(carot.jack.pred, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)))))
carot.jack.df <- data.frame(pred.mean = carot.jack.stat[, 1],
                          pred.low = carot.jack.stat[, 2],
                          pred.high = carot.jack.stat[, 3],
                          Measured = meta(ref.test)$carotenoids,
                          ncomp = ncomp_carot_CVmodel,
                          Species = meta(ref.test)$specie,
                          ID = meta(ref.test)$sample_id)

carbon.jack.pred <- apply.coefs(carbon.jack.coefs, as.matrix(ref.test))
carbon.jack.stat <- t(apply(carbon.jack.pred, 1,
                          function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)))))
carbon.jack.df <- data.frame(pred.mean = carbon.jack.stat[, 1],
                           pred.low = carbon.jack.stat[, 2],
                           pred.high = carbon.jack.stat[, 3],
                           Measured = meta(ref.test)$carbon,
                           ncomp = ncomp_carbon_CVmodel,
                           Species = meta(ref.test)$specie,
                           ID = meta(ref.test)$sample_id)

nitrogen.jack.pred <- apply.coefs(nitrogen.jack.coefs, as.matrix(ref.test))
nitrogen.jack.stat <- t(apply(nitrogen.jack.pred, 1,
                            function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)))))
nitrogen.jack.df <- data.frame(pred.mean = nitrogen.jack.stat[, 1],
                             pred.low = nitrogen.jack.stat[, 2],
                             pred.high = nitrogen.jack.stat[, 3],
                             Measured = meta(ref.test)$nitrogen,
                             ncomp = ncomp_nitrogen_CVmodel,
                             Species = meta(ref.test)$specie,
                             ID = meta(ref.test)$sample_id)

all.jack.coef.list <- list(SLA = SLA.jack.coefs,
                         LMA = LMA.jack.coefs,
                         LDMC = LDMC.jack.coefs,
                         EWT = EWT.jack.coefs,
                         Hemi = hemi.jack.coefs,
                         Cell = cell.jack.coefs,
                         Lignin = lignin.jack.coefs,
                         "Chl a" = chla.jack.coefs,
                         "Chl b" = chlb.jack.coefs,
                         Car = carot.jack.coefs,
                         Carbon = carbon.jack.coefs,
                         Nitrogen = nitrogen.jack.coefs)

saveRDS(all.jack.coef.list, "output_plsr/all_jack_coefs_list.rds")

all.jack.df.list <- list(SLA = SLA.jack.df,
                       LMA = LMA.jack.df,
                       LDMC = LDMC.jack.df,
                       EWT = EWT.jack.df,
                       Hemi = hemi.jack.df,
                       Cell = cell.jack.df,
                       Lignin = lignin.jack.df,
                       "Chl a" = chla.jack.df,
                       "Chl b" = chlb.jack.df,
                       Car = carot.jack.df,
                       Carbon = carbon.jack.df,
                       Nitrogen = nitrogen.jack.df)

saveRDS(all.jack.df.list, "output_plsr/all_jack_df_list.rds")

all.jack.stats.list <- list(SLA = SLA.jack.stats,
                          LMA = LMA.jack.stats,
                          LDMC = LDMC.jack.stats,
                          EWT = EWT.jack.stats,
                          Hemi = hemi.jack.stats,
                          Cell = cell.jack.stats,
                          Lignin = lignin.jack.stats,
                          "Chl a" = chla.jack.stats,
                          "Chl b" = chlb.jack.stats,
                          Car = carot.jack.stats,
                          Carbon = carbon.jack.stats,
                          Nitrogen = nitrogen.jack.stats)

saveRDS(all.jack.stats.list, "output_plsr/all_jack_stats_list.rds")

## Barplots----
val_summary <- data.frame(
  variable = names(all.jack.df.list),
  perRMSE = unlist(lapply(all.jack.df.list, function(x) percentRMSD(x$Measured, x$pred.mean, 0.025, 0.975) * 100)),
  RMSE = unlist(lapply(all.jack.df.list, function(x) RMSD(x$Measured, x$pred.mean))),
  R2 = unlist(lapply(all.jack.df.list, function(x) summary(lm(Measured ~ pred.mean, data = x))$r.squared)),
  bias = unlist(lapply(all.jack.df.list, function(x) mean(x$pred.mean - x$Measured, na.rm = T))),
  ncomp = c(ncomp_SLA_CVmodel, ncomp_LMA_CVmodel, ncomp_LDMC_CVmodel, ncomp_EWT_CVmodel,
            ncomp_hemi_CVmodel, ncomp_cell_CVmodel, ncomp_lignin_CVmodel, ncomp_chla_CVmodel,
            ncomp_chlb_CVmodel, ncomp_carot_CVmodel, ncomp_carbon_CVmodel, ncomp_nitrogen_CVmodel))

val_summary <- val_summary[order(-val_summary$R2),] #sort val_summary in descending order based on R2 values

saveRDS(val_summary, "output_plsr/val_summary.rds")
write.csv(val_summary, "output_plsr/val_summary.csv", row.names = F)

# Plot
val_summary$variable <- factor(val_summary$variable, levels = unique(val_summary$variable)) #set the levels of variables

val_R2 <- ggplot(val_summary, aes(y = R2, x = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = expression(italic("R"^2))) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_x_discrete(labels = c("LMA", "EWT", "SLA", "Car", expression(Chl ~ italic(a)), "nitrogen", expression(Chl ~ italic(b)), "Cell", "Lignin", "Carbon", "LDMC", "Hemi"))
val_R2

val_perRMSE <- ggplot(val_summary, aes(y = perRMSE, x = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "%RMSE") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "red", linewidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +
  scale_x_discrete(labels = c("LMA", "EWT", "SLA", "Car", expression(Chl ~ italic(a)), "nitrogen", expression(Chl ~ italic(b)), "Cell", "Lignin", "Carbon", "LDMC", "Hemi"))
val_perRMSE

# Arrange all plots together
barplots <- egg::ggarrange(plots = list(val_R2, val_perRMSE),
               nrow = 2, ncol = 1, labels = c("(A)", "(B)"))

# Save plots
pdf("output_plsr/barplots.pdf", width = 8, height = 8, onefile = F)
barplots
dev.off()

png("output_plsr/barplots.png", width = 8, height = 8, units = "in", res = 300)
barplots
dev.off()

