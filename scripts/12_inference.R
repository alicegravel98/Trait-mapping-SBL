### Inference : apply PLSR models to all imagery

## Load libraries----
library(tidyverse)
library(spectrolab)
library(terra)
library(sf)

# Set working directory
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

source("scripts/00_useful_functions.R") #import useful functions

## Import coefficients----
all_jack_coefs_list <- readRDS("output_plsr/all_jack_coefs_list.rds")

#### Flightline 01####----
# Import spectra
spectra01 <- readRDS("output_allpixels/smooth_spectra01.rds")

spec01 <- as_spectra(spectra01, name_idx = 1, meta_idxs = c(2:3)) #change to spectra object (this takes a long time)

# Normalize spectra
spec01_bn <- t(apply(as.matrix(spec01), 1, function(x) x/sqrt(sum(x^2))))
spec01_bn <- spectra(spec01_bn, bands = bands(spec01), names = names(spec01))
meta(spec01_bn) <- meta(spec01)

rm(spec01)
gc()

saveRDS(spec01_bn, "output_allpixels/spec01_bn.rds") #to load object : spec01_bn <- readRDS("output_allpixels/spec01_bn.rds")

# Apply PLSR models
SLA_pred01 <- apply.coefs(all_jack_coefs_list$SLA, val.spec = spec01_bn, intercept = T)
SLA_stat01 <- t(apply(SLA_pred01, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
SLA_pred_df01 <- data.frame(pred.mean = SLA_stat01[, 1],
                            pred.low = SLA_stat01[, 2],
                            pred.high = SLA_stat01[, 3])
SLA_all01 <- cbind(SLA_pred_df01, sd = SLA_stat01[, 4])

min(SLA_all01$pred.mean)
SLA_neg <- length(which(SLA_all01$pred.mean < 0))
SLA_neg #41 599 pixels are negative values
saveRDS(SLA_all01, "inference/flg01/SLA_all01.rds") #to load object : SLA_all01 <- readRDS("inference/flg01/SLA_all01.rds")

rm(SLA_pred01, SLA_pred_df01, SLA_stat01)
gc()

LMA_pred01 <- apply.coefs(all_jack_coefs_list$LMA, val.spec = spec01_bn, intercept = T)
LMA_stat01 <- t(apply(LMA_pred01, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LMA_pred_df01 <- data.frame(pred.mean = LMA_stat01[, 1],
                            pred.low = LMA_stat01[, 2],
                            pred.high = LMA_stat01[, 3])
LMA_all01 <- cbind(LMA_pred_df01, sd = LMA_stat01[, 4])

min(LMA_all01$pred.mean)
LMA_neg <- length(which(LMA_all01$pred.mean < 0))
LMA_neg #23 597 pixels have negative values
saveRDS(LMA_all01, "inference/flg01/LMA_all01.rds") #to load object : LMA_all01 <- readRDS("inference/flg01/LMA_all01.rds")

rm(LMA_pred01, LMA_pred_df01, LMA_stat01)
gc()

LDMC_pred01 <- apply.coefs(all_jack_coefs_list$LDMC, val.spec = spec01_bn, intercept = T)
LDMC_stat01 <- t(apply(LDMC_pred01, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LDMC_pred_df01 <- data.frame(pred.mean = LDMC_stat01[, 1],
                             pred.low = LDMC_stat01[, 2],
                             pred.high = LDMC_stat01[, 3])
LDMC_all01 <- cbind(LDMC_pred_df01, sd = LDMC_stat01[, 4])

min(LDMC_all01$pred.mean)
LDMC_neg <- length(which(LDMC_all01$pred.mean < 0))
LDMC_neg #13 pixels have negative values
saveRDS(LDMC_all01, "inference/flg01/LDMC_all01.rds") #to load object : LDMC_all01 <- readRDS("inference/flg01/LDMC_all01.rds")

rm(LDMC_pred01, LDMC_pred_df01, LDMC_stat01)
gc()

EWT_pred01 <- apply.coefs(all_jack_coefs_list$EWT, val.spec = spec01_bn, intercept = T)
EWT_stat01 <- t(apply(EWT_pred01, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
EWT_pred_df01 <- data.frame(pred.mean = EWT_stat01[, 1],
                            pred.low = EWT_stat01[, 2],
                            pred.high = EWT_stat01[, 3])
EWT_all01 <- cbind(EWT_pred_df01, sd = EWT_stat01[, 4])

min(EWT_all01$pred.mean)
EWT_neg <- length(which(EWT_all01$pred.mean < 0))
EWT_neg #13 310 pixels have negative values
saveRDS(EWT_all01, "inference/flg01/EWT_all01.rds") #to load object : EWT_all01 <- readRDS("inference/flg01/EWT_all01.rds")

rm(EWT_pred01, EWT_pred_df01, EWT_stat01)
gc()

hemi_pred01 <- apply.coefs(all_jack_coefs_list$hemi, val.spec = spec01_bn, intercept = T)
hemi_stat01 <- t(apply(hemi_pred01, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
hemi_pred_df01 <- data.frame(pred.mean = hemi_stat01[, 1],
                             pred.low = hemi_stat01[, 2],
                             pred.high = hemi_stat01[, 3])
hemi_all01 <- cbind(hemi_pred_df01, sd = hemi_stat01[, 4])

min(hemi_all01$pred.mean)
hemi_neg <- length(which(hemi_all01$pred.mean < 0))
hemi_neg #5540 pixels have negative values
saveRDS(hemi_all01, "inference/flg01/hemi_all01.rds") #to load object : hemi_all01 <- readRDS("inference/flg01/hemi_all01.rds")

rm(hemi_pred01, hemi_pred_df01, hemi_stat01)
gc()

cell_pred01 <- apply.coefs(all_jack_coefs_list$cell, val.spec = spec01_bn, intercept = T)
cell_stat01 <- t(apply(cell_pred01, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
cell_pred_df01 <- data.frame(pred.mean = cell_stat01[, 1],
                             pred.low = cell_stat01[, 2],
                             pred.high = cell_stat01[, 3])
cell_all01 <- cbind(cell_pred_df01, sd = cell_stat01[, 4])

min(cell_all01$pred.mean)
cell_neg <- length(which(cell_all01$pred.mean < 0))
cell_neg #80 pixels have negative values
saveRDS(cell_all01, "inference/flg01/cell_all01.rds") #to load object : cell_all01 <- readRDS("inference/flg01/cell_all01.rds")

rm(cell_pred01, cell_pred_df01, cell_stat01)
gc()

lignin_pred01 <- apply.coefs(all_jack_coefs_list$lignin, val.spec = spec01_bn, intercept = T)
lignin_stat01 <- t(apply(lignin_pred01, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
lignin_pred_df01 <- data.frame(pred.mean = lignin_stat01[, 1],
                               pred.low = lignin_stat01[, 2],
                               pred.high = lignin_stat01[, 3])
lignin_all01 <- cbind(lignin_pred_df01, sd = lignin_stat01[, 4])

min(lignin_all01$pred.mean)
lignin_neg <- length(which(lignin_all01$pred.mean < 0))
lignin_neg #16 pixels have negative values
saveRDS(lignin_all01, "inference/flg01/lignin_all01.rds") #to load object : lignin_all01 <- readRDS("inference/flg01/lignin_all01.rds")

rm(lignin_pred01, lignin_pred_df01, lignin_stat01)
gc()

chla_pred01 <- apply.coefs(all_jack_coefs_list$chla, val.spec = spec01_bn, intercept = T)
chla_stat01 <- t(apply(chla_pred01, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chla_pred_df01 <- data.frame(pred.mean = chla_stat01[, 1],
                             pred.low = chla_stat01[, 2],
                             pred.high = chla_stat01[, 3])
chla_all01 <- cbind(chla_pred_df01, sd = chla_stat01[, 4])

min(chla_all01$pred.mean)
chla_neg <- length(which(chla_all01$pred.mean < 0))
chla_neg #4178 pixels have negative values
saveRDS(chla_all01, "inference/flg01/chla_all01.rds") #to load object : chla_all01 <- readRDS("inference/flg01/chla_all01.rds")

rm(chla_pred01, chla_pred_df01, chla_stat01)
gc()

chlb_pred01 <- apply.coefs(all_jack_coefs_list$chlb, val.spec = spec01_bn, intercept = T)
chlb_stat01 <- t(apply(chlb_pred01, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025 ,0.975)), sd(obs))))
chlb_pred_df01 <- data.frame(pred.mean = chlb_stat01[, 1],
                             pred.low = chlb_stat01[, 2],
                             pred.high = chlb_stat01[, 3])
chlb_all01 <- cbind(chlb_pred_df01, sd = chlb_stat01[, 4])

min(chlb_all01$pred.mean)
chlb_neg <- length(which(chlb_all01$pred.mean < 0))
chlb_neg #3400 pixels have negative values
saveRDS(chlb_all01, "inference/flg01/chlb_all01.rds") #to load object : chlb_all01 <- readRDS("inference/flg01/chlb_all01.rds")

rm(chlb_pred01, chlb_pred_df01, chlb_stat01)
gc()

carot_pred01 <- apply.coefs(all_jack_coefs_list$carot, val.spec = spec01_bn, intercept = T)
carot_stat01 <- t(apply(carot_pred01, 1,
                        function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
carot_pred_df01 <- data.frame(pred.mean = carot_stat01[, 1],
                              pred.low = carot_stat01[, 2],
                              pred.high = carot_stat01[, 3])
carot_all01 <- cbind(carot_pred_df01, sd = carot_stat01[, 4])

min(carot_all01$pred.mean)
carot_neg <- length(which(carot_all01$pred.mean < 0))
carot_neg #4572 pixels have negative values
saveRDS(carot_all01, "inference/flg01/carot_all01.rds") #to load object : carot_all01 <- readRDS("inference/flg01/carot_all01.rds")

rm(carot_pred01, carot_pred_df01, carot_stat01)
gc()

carbon_pred01 <- apply.coefs(all_jack_coefs_list$carbon, val.spec = spec01_bn, intercept = T)
carbon_stat01 <- t(apply(carbon_pred01, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
carbon_pred_df01 <- data.frame(pred.mean = carbon_stat01[, 1],
                               pred.low = carbon_stat01[, 2],
                               pred.high = carbon_stat01[, 3])
carbon_all01 <- cbind(carbon_pred_df01, sd = carbon_stat01[, 4])

min(carbon_all01$pred.mean)
carbon_neg <- length(which(carbon_all01$pred.mean < 0))
carbon_neg #0 pixels have negative values
saveRDS(carbon_all01, "inference/flg01/carbon_all01.rds") #to load object : carbon_all01 <- readRDS("inference/flg01/carbon_all01.rds")

rm(carbon_pred01, carbon_pred_df01, carbon_stat01)
gc()

nitrogen_pred01 <- apply.coefs(all_jack_coefs_list$nitrogen, val.spec = spec01_bn, intercept = T)
nitrogen_stat01 <- t(apply(nitrogen_pred01, 1,
                           function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
nitrogen_pred_df01 <- data.frame(pred.mean = nitrogen_stat01[, 1],
                                 pred.low = nitrogen_stat01[, 2],
                                 pred.high = nitrogen_stat01[, 3])
nitrogen_all01 <- cbind(nitrogen_pred_df01, sd = nitrogen_stat01[, 4])

min(nitrogen_all01$pred.mean)
nitrogen_neg <- length(which(nitrogen_all01$pred.mean < 0))
nitrogen_neg #15 480 pixels have negative values
saveRDS(nitrogen_all01, "inference/flg01/nitrogen_all01.rds") #to load object : nitrogen_all01 <- readRDS("inference/flg01/nitrogen_all01.rds")

rm(nitrogen_pred01, nitrogen_pred_df01, nitrogen_stat01)
gc()

# Extract coordinates
x_coords01 <- spectra01$x
y_coords01 <- spectra01$y

# Add coordinates to results
SLA_all01$x <- x_coords01
SLA_all01$y <- y_coords01
LMA_all01$x <- x_coords01
LMA_all01$y <- y_coords01
LDMC_all01$x <- x_coords01
LDMC_all01$y <- y_coords01
EWT_all01$x <- x_coords01
EWT_all01$y <- y_coords01
hemi_all01$x <- x_coords01
hemi_all01$y <- y_coords01
cell_all01$x <- x_coords01
cell_all01$y <- y_coords01
lignin_all01$x <- x_coords01
lignin_all01$y <- y_coords01
chla_all01$x <- x_coords01
chla_all01$y <- y_coords01
chlb_all01$x <- x_coords01
chlb_all01$y <- y_coords01
carot_all01$x <- x_coords01
carot_all01$y <- y_coords01
carbon_all01$x <- x_coords01
carbon_all01$y <- y_coords01
nitrogen_all01$x <- x_coords01
nitrogen_all01$y <- y_coords01

# Select columns of interest for mapping
SLA_01 <- SLA_all01 %>%
  dplyr::select(pred.SLA = pred.mean, sd.SLA = sd, x, y) %>%
  `rownames<-`(NULL)
LMA_01 <- LMA_all01 %>%
  dplyr::select(pred.LMA = pred.mean, sd.LMA = sd, x, y) %>%
  `rownames<-`(NULL)
LDMC_01 <- LDMC_all01 %>%
  dplyr::select(pred.LDMC = pred.mean, sd.LDMC = sd, x, y) %>%
  `rownames<-`(NULL)
EWT_01 <- EWT_all01 %>%
  dplyr::select(pred.EWT = pred.mean, sd.EWT = sd, x, y) %>%
  `rownames<-`(NULL)
hemi_01 <- hemi_all01 %>%
  dplyr::select(pred.hemi = pred.mean, sd.hemi = sd, x, y) %>%
  `rownames<-`(NULL)
cell_01 <- cell_all01 %>%
  dplyr::select(pred.cell = pred.mean, sd.cell = sd, x, y) %>%
  `rownames<-`(NULL)
lignin_01 <- lignin_all01 %>%
  dplyr::select(pred.lignin = pred.mean, sd.lignin = sd, x, y) %>%
  `rownames<-`(NULL)
chla_01 <- chla_all01 %>%
  dplyr::select(pred.chla = pred.mean, sd.chla = sd, x, y) %>%
  `rownames<-`(NULL)
chlb_01 <- chlb_all01 %>%
  dplyr::select(pred.chlb = pred.mean, sd.chlb = sd, x, y) %>%
  `rownames<-`(NULL)
carot_01 <- carot_all01 %>%
  dplyr::select(pred.carot = pred.mean, sd.carot = sd, x, y) %>%
  `rownames<-`(NULL)
carbon_01 <- carbon_all01 %>%
  dplyr::select(pred.carbon = pred.mean, sd.carbon = sd, x, y) %>%
  `rownames<-`(NULL)
nitrogen_01 <- nitrogen_all01 %>%
  dplyr::select(pred.nitrogen = pred.mean, sd.nitrogen = sd, x, y) %>%
  `rownames<-`(NULL)

# Join all dataframes
df_list_01 = list(SLA_01, LMA_01, LDMC_01, EWT_01, hemi_01, cell_01, lignin_01, chla_01, chlb_01, carot_01, carbon_01, nitrogen_01) 

df_merged01 <- df_list_01 %>%
  reduce(right_join, by = c("x", "y")) %>%
  relocate("x", "y", .before = "pred.SLA")
saveRDS(df_merged01, "inference/flg01/df_merged01.rds") #to load object : df_merged01 <- readRDS("inference/flg01/df_merged01.rds")

# Add relative uncertainties (CV)
df_merged01 <- df_merged01 %>%
  mutate(CV.SLA = sd.SLA/pred.SLA,
         CV.LMA = sd.LMA/pred.LMA,
         CV.LDMC = sd.LDMC/pred.LDMC,
         CV.EWT = sd.EWT/pred.EWT,
         CV.hemi = sd.hemi/pred.hemi,
         CV.cell = sd.cell/pred.cell,
         CV.lignin = sd.lignin/pred.lignin,
         CV.chla = sd.chla/pred.chla,
         CV.chlb = sd.chlb/pred.chlb,
         CV.carot = sd.carot/pred.carot,
         CV.carbon = sd.carbon/pred.carbon,
         CV.nitrogen = sd.nitrogen/pred.nitrogen)

# Create 3 dataframes (traits, sd and CV)
df_traits01 <- df_merged01 %>%
  dplyr::select(starts_with("pred"), x, y) %>%
  rename_with(~sub("pred\\.", "", .), starts_with("pred"))
df_sd01 <- df_merged01 %>%
  dplyr::select(starts_with("sd"), x, y) %>%
  rename_with(~sub("sd\\.", "", .), starts_with("sd"))
df_CV01 <- df_merged01 %>%
  dplyr::select(starts_with("CV"), x, y) %>%
  rename_with(~sub("CV\\.", "", .), starts_with("CV"))
saveRDS(df_traits01, "inference/flg01/df_traits01.rds")
saveRDS(df_sd01, "inference/flg01/df_sd01.rds")
saveRDS(df_CV01, "inference/flg01/df_CV01.rds")

#### Flightline 02####----
# Import spectra
spectra02 <- readRDS("output_allpixels/smooth_spectra02.rds")

spec02 <- as_spectra(spectra02, name_idx = 1, meta_idxs = c(2:3)) #change to spectra object (this takes a long time)

# Normalize spectra
spec02_bn <- t(apply(as.matrix(spec02), 1, function(x) x/sqrt(sum(x^2))))
spec02_bn <- spectra(spec02_bn, bands = bands(spec02), names = names(spec02))
spectrolab::meta(spec02_bn) <- spectrolab::meta(spec02)

rm(spec02)
gc()

saveRDS(spec02_bn, "output_allpixels/spec02_bn.rds") #to load object : spec02_bn <- readRDS("output_allpixels/spec02_bn.rds")

# Apply PLSR models
SLA_pred02 <- apply.coefs(all_jack_coefs_list$SLA, val.spec = spec02_bn, intercept = T)
SLA_stat02 <- t(apply(SLA_pred02, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
SLA_pred_df02 <- data.frame(pred.mean = SLA_stat02[, 1],
                            pred.low = SLA_stat02[, 2],
                            pred.high = SLA_stat02[, 3])
SLA_all02 <- cbind(SLA_pred_df02, sd = SLA_stat02[, 4])

min(SLA_all02$pred.mean)
SLA_neg <- length(which(SLA_all02$pred.mean < 0))
SLA_neg #23 728 pixels are negative values
saveRDS(SLA_all02, "inference/flg02/SLA_all02.rds") #to load object : SLA_all02 <- readRDS("inference/flg02/SLA_all02.rds")

rm(SLA_pred02, SLA_pred_df02, SLA_stat02)
gc()

LMA_pred02 <- apply.coefs(all_jack_coefs_list$LMA, val.spec = spec02_bn, intercept = T)
LMA_stat02 <- t(apply(LMA_pred02, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LMA_pred_df02 <- data.frame(pred.mean = LMA_stat02[, 1],
                            pred.low = LMA_stat02[, 2],
                            pred.high = LMA_stat02[, 3])
LMA_all02 <- cbind(LMA_pred_df02, sd = LMA_stat02[, 4])

min(LMA_all02$pred.mean)
LMA_neg <- length(which(LMA_all02$pred.mean < 0))
LMA_neg #22 566 pixels have negative values
saveRDS(LMA_all02, "inference/flg02/LMA_all02.rds") #to load object : LMA_all02 <- readRDS("inference/flg02/LMA_all02.rds")

rm(LMA_pred02, LMA_pred_df02, LMA_stat02)
gc()

LDMC_pred02 <- apply.coefs(all_jack_coefs_list$LDMC, val.spec = spec02_bn, intercept = T)
LDMC_stat02 <- t(apply(LDMC_pred02, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LDMC_pred_df02 <- data.frame(pred.mean = LDMC_stat02[, 1],
                             pred.low = LDMC_stat02[, 2],
                             pred.high = LDMC_stat02[, 3])
LDMC_all02 <- cbind(LDMC_pred_df02, sd = LDMC_stat02[, 4])

min(LDMC_all02$pred.mean)
LDMC_neg <- length(which(LDMC_all02$pred.mean < 0))
LDMC_neg #22 pixels have negative values
saveRDS(LDMC_all02, "inference/flg02/LDMC_all02.rds") #to load object : LDMC_all02 <- readRDS("inference/flg02/LDMC_all02.rds")

rm(LDMC_pred02, LDMC_pred_df02, LDMC_stat02)
gc()

EWT_pred02 <- apply.coefs(all_jack_coefs_list$EWT, val.spec = spec02_bn, intercept = T)
EWT_stat02 <- t(apply(EWT_pred02, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
EWT_pred_df02 <- data.frame(pred.mean = EWT_stat02[, 1],
                            pred.low = EWT_stat02[, 2],
                            pred.high = EWT_stat02[, 3])
EWT_all02 <- cbind(EWT_pred_df02, sd = EWT_stat02[, 4])

min(EWT_all02$pred.mean)
EWT_neg <- length(which(EWT_all02$pred.mean < 0))
EWT_neg #16 294 pixels have negative values
saveRDS(EWT_all02, "inference/flg02/EWT_all02.rds") #to load object : EWT_all02 <- readRDS("inference/flg02/EWT_all02.rds")

rm(EWT_pred02, EWT_pred_df02, EWT_stat02)
gc()

hemi_pred02 <- apply.coefs(all_jack_coefs_list$hemi, val.spec = spec02_bn, intercept = T)
hemi_stat02 <- t(apply(hemi_pred02, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
hemi_pred_df02 <- data.frame(pred.mean = hemi_stat02[, 1],
                             pred.low = hemi_stat02[, 2],
                             pred.high = hemi_stat02[, 3])
hemi_all02 <- cbind(hemi_pred_df02, sd = hemi_stat02[, 4])

min(hemi_all02$pred.mean)
hemi_neg <- length(which(hemi_all02$pred.mean < 0))
hemi_neg #1299 pixels have negative values
saveRDS(hemi_all02, "inference/flg02/hemi_all02.rds") #to load object : hemi_all02 <- readRDS("inference/flg02/hemi_all02.rds")

rm(hemi_pred02, hemi_pred_df02, hemi_stat02)
gc()

cell_pred02 <- apply.coefs(all_jack_coefs_list$cell, val.spec = spec02_bn, intercept = T)
cell_stat02 <- t(apply(cell_pred02, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
cell_pred_df02 <- data.frame(pred.mean = cell_stat02[, 1],
                             pred.low = cell_stat02[, 2],
                             pred.high = cell_stat02[, 3])
cell_all02 <- cbind(cell_pred_df02, sd = cell_stat02[, 4])

min(cell_all02$pred.mean)
cell_neg <- length(which(cell_all02$pred.mean < 0))
cell_neg #92 pixels have negative values
saveRDS(cell_all02, "inference/flg02/cell_all02.rds") #to load object : cell_all02 <- readRDS("inference/flg02/cell_all02.rds")

rm(cell_pred02, cell_pred_df02, cell_stat02)
gc()

lignin_pred02 <- apply.coefs(all_jack_coefs_list$lignin, val.spec = spec02_bn, intercept = T)
lignin_stat02 <- t(apply(lignin_pred02, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
lignin_pred_df02 <- data.frame(pred.mean = lignin_stat02[, 1],
                               pred.low = lignin_stat02[, 2],
                               pred.high = lignin_stat02[, 3])
lignin_all02 <- cbind(lignin_pred_df02, sd = lignin_stat02[, 4])

min(lignin_all02$pred.mean)
lignin_neg <- length(which(lignin_all02$pred.mean < 0))
lignin_neg #13 pixels have negative values
saveRDS(lignin_all02, "inference/flg02/lignin_all02.rds") #to load object : lignin_all02 <- readRDS("inference/flg02/lignin_all02.rds")

rm(lignin_pred02, lignin_pred_df02, lignin_stat02)
gc()

chla_pred02 <- apply.coefs(all_jack_coefs_list$chla, val.spec = spec02_bn, intercept = T)
chla_stat02 <- t(apply(chla_pred02, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chla_pred_df02 <- data.frame(pred.mean = chla_stat02[, 1],
                             pred.low = chla_stat02[, 2],
                             pred.high = chla_stat02[, 3])
chla_all02 <- cbind(chla_pred_df02, sd = chla_stat02[, 4])

min(chla_all02$pred.mean)
chla_neg <- length(which(chla_all02$pred.mean < 0))
chla_neg #4686 pixels have negative values
saveRDS(chla_all02, "inference/flg02/chla_all02.rds") #to load object : chla_all02 <- readRDS("inference/flg02/chla_all02.rds")

rm(chla_pred02, chla_pred_df02, chla_stat02)
gc()

chlb_pred02 <- apply.coefs(all_jack_coefs_list$chlb, val.spec = spec02_bn, intercept = T)
chlb_stat02 <- t(apply(chlb_pred02, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chlb_pred_df02 <- data.frame(pred.mean = chlb_stat02[, 1],
                             pred.low = chlb_stat02[, 2],
                             pred.high = chlb_stat02[, 3])
chlb_all02 <- cbind(chlb_pred_df02, sd = chlb_stat02[, 4])

min(chlb_all02$pred.mean)
chlb_neg <- length(which(chlb_all02$pred.mean < 0))
chlb_neg #3728 pixels have negative values
saveRDS(chlb_all02, "inference/flg02/chlb_all02.rds") #to load object : chlb_all02 <- readRDS("inference/flg02/chlb_all02.rds")

rm(chlb_pred02, chlb_pred_df02, chlb_stat02)
gc()

carot_pred02 <- apply.coefs(all_jack_coefs_list$carot, val.spec = spec02_bn, intercept = T)
carot_stat02 <- t(apply(carot_pred02, 1,
                        function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
carot_pred_df02 <- data.frame(pred.mean = carot_stat02[, 1],
                              pred.low = carot_stat02[, 2],
                              pred.high = carot_stat02[, 3])
carot_all02 <- cbind(carot_pred_df02, sd = carot_stat02[, 4])

min(carot_all02$pred.mean)
carot_neg <- length(which(carot_all02$pred.mean < 0))
carot_neg #4247 pixels have negative values
saveRDS(carot_all02, "inference/flg02/carot_all02.rds") # to load object: carot_all02 <- readRDS("inference/flg02/carot_all02.rds")

rm(carot_pred02, carot_pred_df02, carot_stat02)
gc()

carbon_pred02 <- apply.coefs(all_jack_coefs_list$carbon, val.spec = spec02_bn, intercept = T)
carbon_stat02 <- t(apply(carbon_pred02, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
carbon_pred_df02 <- data.frame(pred.mean = carbon_stat02[, 1],
                               pred.low = carbon_stat02[, 2],
                               pred.high = carbon_stat02[, 3])
carbon_all02 <- cbind(carbon_pred_df02, sd = carbon_stat02[, 4])

min(carbon_all02$pred.mean)
carbon_neg <- length(which(carbon_all02$pred.mean < 0))
carbon_neg #0 pixels have negative values
saveRDS(carbon_all02, "inference/flg02/carbon_all02.rds") # to load object: carbon_all02 <- readRDS("inference/flg02/carbon_all02.rds")

rm(carbon_pred02, carbon_pred_df02, carbon_stat02)
gc()

nitrogen_pred02 <- apply.coefs(all_jack_coefs_list$nitrogen, val.spec = spec02_bn, intercept = T)
nitrogen_stat02 <- t(apply(nitrogen_pred02, 1,
                           function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
nitrogen_pred_df02 <- data.frame(pred.mean = nitrogen_stat02[, 1],
                                 pred.low = nitrogen_stat02[, 2],
                                 pred.high = nitrogen_stat02[, 3])
nitrogen_all02 <- cbind(nitrogen_pred_df02, sd = nitrogen_stat02[, 4])

min(nitrogen_all02$pred.mean)
nitrogen_neg <- length(which(nitrogen_all02$pred.mean < 0))
nitrogen_neg #9686 pixels have negative values
saveRDS(nitrogen_all02, "inference/flg02/nitrogen_all02.rds") # to load object: nitrogen_all02 <- readRDS("inference/flg02/nitrogen_all02.rds")

rm(nitrogen_pred02, nitrogen_pred_df02, nitrogen_stat02)
gc()

# Extract coordinates
x_coords02 <- spectra02$x
y_coords02 <- spectra02$y

# Add coordinates to results
SLA_all02$x <- x_coords02
SLA_all02$y <- y_coords02
LMA_all02$x <- x_coords02
LMA_all02$y <- y_coords02
LDMC_all02$x <- x_coords02
LDMC_all02$y <- y_coords02
EWT_all02$x <- x_coords02
EWT_all02$y <- y_coords02
hemi_all02$x <- x_coords02
hemi_all02$y <- y_coords02
cell_all02$x <- x_coords02
cell_all02$y <- y_coords02
lignin_all02$x <- x_coords02
lignin_all02$y <- y_coords02
chla_all02$x <- x_coords02
chla_all02$y <- y_coords02
chlb_all02$x <- x_coords02
chlb_all02$y <- y_coords02
carot_all02$x <- x_coords02
carot_all02$y <- y_coords02
carbon_all02$x <- x_coords02
carbon_all02$y <- y_coords02
nitrogen_all02$x <- x_coords02
nitrogen_all02$y <- y_coords02

# Select columns of interest for mapping
SLA_02 <- SLA_all02 %>%
  dplyr::select(pred.SLA = pred.mean, sd.SLA = sd, x, y) %>%
  `rownames<-`(NULL)
LMA_02 <- LMA_all02 %>%
  dplyr::select(pred.LMA = pred.mean, sd.LMA = sd, x, y) %>%
  `rownames<-`(NULL)
LDMC_02 <- LDMC_all02 %>%
  dplyr::select(pred.LDMC = pred.mean, sd.LDMC = sd, x, y) %>%
  `rownames<-`(NULL)
EWT_02 <- EWT_all02 %>%
  dplyr::select(pred.EWT = pred.mean, sd.EWT = sd, x, y) %>%
  `rownames<-`(NULL)
hemi_02 <- hemi_all02 %>%
  dplyr::select(pred.hemi = pred.mean, sd.hemi = sd, x, y) %>%
  `rownames<-`(NULL)
cell_02 <- cell_all02 %>%
  dplyr::select(pred.cell = pred.mean, sd.cell = sd, x, y) %>%
  `rownames<-`(NULL)
lignin_02 <- lignin_all02 %>%
  dplyr::select(pred.lignin = pred.mean, sd.lignin = sd, x, y) %>%
  `rownames<-`(NULL)
chla_02 <- chla_all02 %>%
  dplyr::select(pred.chla = pred.mean, sd.chla = sd, x, y) %>%
  `rownames<-`(NULL)
chlb_02 <- chlb_all02 %>%
  dplyr::select(pred.chlb = pred.mean, sd.chlb = sd, x, y) %>%
  `rownames<-`(NULL)
carot_02 <- carot_all02 %>%
  dplyr::select(pred.carot = pred.mean, sd.carot = sd, x, y) %>%
  `rownames<-`(NULL)
carbon_02 <- carbon_all02 %>%
  dplyr::select(pred.carbon = pred.mean, sd.carbon = sd, x, y) %>%
  `rownames<-`(NULL)
nitrogen_02 <- nitrogen_all02 %>%
  dplyr::select(pred.nitrogen = pred.mean, sd.nitrogen = sd, x, y) %>%
  `rownames<-`(NULL)

# Join all dataframes
df_list_02 = list(SLA_02, LMA_02, LDMC_02, EWT_02, hemi_02, cell_02, lignin_02, chla_02, chlb_02, carot_02, carbon_02, nitrogen_02) 

df_merged02 <- df_list_02 %>%
  reduce(right_join, by = c("x", "y"))
saveRDS(df_merged02, "inference/flg02/df_merged02.rds") #to load object : df_merged02 <- readRDS("inference/flg02/df_merged02.rds")

# Add relative uncertainties (CV)
df_merged02 <- df_merged02 %>%
  mutate(CV.SLA = sd.SLA/pred.SLA,
         CV.LMA = sd.LMA/pred.LMA,
         CV.LDMC = sd.LDMC/pred.LDMC,
         CV.EWT = sd.EWT/pred.EWT,
         CV.hemi = sd.hemi/pred.hemi,
         CV.cell = sd.cell/pred.cell,
         CV.lignin = sd.lignin/pred.lignin,
         CV.chla = sd.chla/pred.chla,
         CV.chlb = sd.chlb/pred.chlb,
         CV.carot = sd.carot/pred.carot,
         CV.carbon = sd.carbon/pred.carbon,
         CV.nitrogen = sd.nitrogen/pred.nitrogen)

# Create 3 dataframes (traits, sd and CV)
df_traits02 <- df_merged02 %>%
  dplyr::select(starts_with("pred"), x, y) %>%
  rename_with(~sub("pred\\.", "", .), starts_with("pred"))
df_sd02 <- df_merged02 %>%
  dplyr::select(starts_with("sd"), x, y) %>%
  rename_with(~sub("sd\\.", "", .), starts_with("sd"))
df_CV02 <- df_merged02 %>%
  dplyr::select(starts_with("CV"), x, y) %>%
  rename_with(~sub("CV\\.", "", .), starts_with("CV"))
saveRDS(df_traits02, "inference/flg02/df_traits02.rds")
saveRDS(df_sd02, "inference/flg02/df_sd02.rds")
saveRDS(df_CV02, "inference/flg02/df_CV02.rds")

#### Flightline 03####----
# Import spectra
spectra03 <- readRDS("output_allpixels/smooth_spectra03.rds")

spec03 <- as_spectra(spectra03, name_idx = 1, meta_idxs = c(2:3)) #change to spectra object (this takes a long time)

# Normalize spectra
spec03_bn <- t(apply(as.matrix(spec03), 1, function(x) x/sqrt(sum(x^2))))
spec03_bn <- spectra(spec03_bn, bands = bands(spec03), names = names(spec03))
spectrolab::meta(spec03_bn) <- spectrolab::meta(spec03)

rm(spec03)
gc()

saveRDS(spec03_bn, "output_allpixels/spec03_bn.rds") #to load object : spec03_bn <- readRDS("output_allpixels/spec03_bn.rds")

# Apply PLSR models
SLA_pred03 <- apply.coefs(all_jack_coefs_list$SLA, val.spec = spec03_bn, intercept = T)
SLA_stat03 <- t(apply(SLA_pred03, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
SLA_pred_df03 <- data.frame(pred.mean = SLA_stat03[, 1],
                            pred.low = SLA_stat03[, 2],
                            pred.high = SLA_stat03[, 3])
SLA_all03 <- cbind(SLA_pred_df03, sd = SLA_stat03[, 4])

min(SLA_all03$pred.mean)
SLA_neg <- length(which(SLA_all03$pred.mean < 0))
SLA_neg #23 405 pixels are negative values
saveRDS(SLA_all03, "inference/flg03/SLA_all03.rds") #to load object : SLA_all03 <- readRDS("inference/flg03/SLA_all03.rds")

rm(SLA_pred03, SLA_pred_df03, SLA_stat03)
gc()

LMA_pred03 <- apply.coefs(all_jack_coefs_list$LMA, val.spec = spec03_bn, intercept = T)
LMA_stat03 <- t(apply(LMA_pred03, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LMA_pred_df03 <- data.frame(pred.mean = LMA_stat03[, 1],
                            pred.low = LMA_stat03[, 2],
                            pred.high = LMA_stat03[, 3])
LMA_all03 <- cbind(LMA_pred_df03, sd = LMA_stat03[, 4])

min(LMA_all03$pred.mean)
LMA_neg <- length(which(LMA_all03$pred.mean < 0))
LMA_neg #19 689 pixels have negative values
saveRDS(LMA_all03, "inference/flg03/LMA_all03.rds") #to load object : LMA_all03 <- readRDS("inference/flg03/LMA_all03.rds")

rm(LMA_pred03, LMA_pred_df03, LMA_stat03)
gc()

LDMC_pred03 <- apply.coefs(all_jack_coefs_list$LDMC, val.spec = spec03_bn, intercept = T)
LDMC_stat03 <- t(apply(LDMC_pred03, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LDMC_pred_df03 <- data.frame(pred.mean = LDMC_stat03[, 1],
                             pred.low = LDMC_stat03[, 2],
                             pred.high = LDMC_stat03[, 3])
LDMC_all03 <- cbind(LDMC_pred_df03, sd = LDMC_stat03[, 4])

min(LDMC_all03$pred.mean)
LDMC_neg <- length(which(LDMC_all03$pred.mean < 0))
LDMC_neg #2 pixels have negative values
saveRDS(LDMC_all03, "inference/flg03/LDMC_all03.rds") #to load object : LDMC_all03 <- readRDS("inference/flg03/LDMC_all03.rds")

rm(LDMC_pred03, LDMC_pred_df03, LDMC_stat03)
gc()

EWT_pred03 <- apply.coefs(all_jack_coefs_list$EWT, val.spec = spec03_bn, intercept = T)
EWT_stat03 <- t(apply(EWT_pred03, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
EWT_pred_df03 <- data.frame(pred.mean = EWT_stat03[, 1],
                            pred.low = EWT_stat03[, 2],
                            pred.high = EWT_stat03[, 3])
EWT_all03 <- cbind(EWT_pred_df03, sd = EWT_stat03[, 4])

min(EWT_all03$pred.mean)
EWT_neg <- length(which(EWT_all03$pred.mean < 0))
EWT_neg #13 451 pixels have negative values
saveRDS(EWT_all03, "inference/flg03/EWT_all03.rds") #to load object : EWT_all03 <- readRDS("inference/flg03/EWT_all03.rds")

rm(EWT_pred03, EWT_pred_df03, EWT_stat03)
gc()

hemi_pred03 <- apply.coefs(all_jack_coefs_list$hemi, val.spec = spec03_bn, intercept = T)
hemi_stat03 <- t(apply(hemi_pred03, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
hemi_pred_df03 <- data.frame(pred.mean = hemi_stat03[, 1],
                             pred.low = hemi_stat03[, 2],
                             pred.high = hemi_stat03[, 3])
hemi_all03 <- cbind(hemi_pred_df03, sd = hemi_stat03[, 4])

min(hemi_all03$pred.mean)
hemi_neg <- length(which(hemi_all03$pred.mean < 0))
hemi_neg #5086 pixels have negative values
saveRDS(hemi_all03, "inference/flg03/hemi_all03.rds") #to load object : hemi_all03 <- readRDS("inference/flg03/hemi_all03.rds")

rm(hemi_pred03, hemi_pred_df03, hemi_stat03)
gc()

cell_pred03 <- apply.coefs(all_jack_coefs_list$cell, val.spec = spec03_bn, intercept = T)
cell_stat03 <- t(apply(cell_pred03, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
cell_pred_df03 <- data.frame(pred.mean = cell_stat03[, 1],
                             pred.low = cell_stat03[, 2],
                             pred.high = cell_stat03[, 3])
cell_all03 <- cbind(cell_pred_df03, sd = cell_stat03[, 4])

min(cell_all03$pred.mean)
cell_neg <- length(which(cell_all03$pred.mean < 0))
cell_neg #35 pixels have negative values
saveRDS(cell_all03, "inference/flg03/cell_all03.rds") #to load object : cell_all03 <- readRDS("inference/flg03/cell_all03.rds")

rm(cell_pred03, cell_pred_df03, cell_stat03)
gc()

lignin_pred03 <- apply.coefs(all_jack_coefs_list$lignin, val.spec = spec03_bn, intercept = T)
lignin_stat03 <- t(apply(lignin_pred03, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
lignin_pred_df03 <- data.frame(pred.mean = lignin_stat03[, 1],
                               pred.low = lignin_stat03[, 2],
                               pred.high = lignin_stat03[, 3])
lignin_all03 <- cbind(lignin_pred_df03, sd = lignin_stat03[, 4])

min(lignin_all03$pred.mean)
lignin_neg <- length(which(lignin_all03$pred.mean < 0))
lignin_neg #12 pixels have negative values
saveRDS(lignin_all03, "inference/flg03/lignin_all03.rds") #to load object : lignin_all03 <- readRDS("inference/flg03/lignin_all03.rds")

rm(lignin_pred03, lignin_pred_df03, lignin_stat03)
gc()

chla_pred03 <- apply.coefs(all_jack_coefs_list$chla, val.spec = spec03_bn, intercept = T)
chla_stat03 <- t(apply(chla_pred03, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chla_pred_df03 <- data.frame(pred.mean = chla_stat03[, 1],
                             pred.low = chla_stat03[, 2],
                             pred.high = chla_stat03[, 3])
chla_all03 <- cbind(chla_pred_df03, sd = chla_stat03[, 4])

min(chla_all03$pred.mean)
chla_neg <- length(which(chla_all03$pred.mean < 0))
chla_neg #4920 pixels have negative values
saveRDS(chla_all03, "inference/flg03/chla_all03.rds") #to load object : chla_all03 <- readRDS("inference/flg03/chla_all03.rds")

rm(chla_pred03, chla_pred_df03, chla_stat03)
gc()

chlb_pred03 <- apply.coefs(all_jack_coefs_list$chlb, val.spec = spec03_bn, intercept = T)
chlb_stat03 <- t(apply(chlb_pred03, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chlb_pred_df03 <- data.frame(pred.mean = chlb_stat03[, 1],
                             pred.low = chlb_stat03[, 2],
                             pred.high = chlb_stat03[, 3])
chlb_all03 <- cbind(chlb_pred_df03, sd = chlb_stat03[, 4])

min(chlb_all03$pred.mean)
chlb_neg <- length(which(chlb_all03$pred.mean < 0))
chlb_neg #4096 pixels have negative values
saveRDS(chlb_all03, "inference/flg03/chlb_all03.rds") #to load object : chlb_all03 <- readRDS("inference/flg03/chlb_all03.rds")

rm(chlb_pred03, chlb_pred_df03, chlb_stat03)
gc()

carot_pred03 <- apply.coefs(all_jack_coefs_list$carot, val.spec = spec03_bn, intercept = T)
carot_stat03 <- t(apply(carot_pred03, 1,
                        function(obs) c(mean(obs), quantile(obs, probs=c(0.025, 0.975)), sd(obs))))
carot_pred_df03 <- data.frame(pred.mean = carot_stat03[, 1],
                              pred.low = carot_stat03[, 2],
                              pred.high = carot_stat03[, 3])
carot_all03 <- cbind(carot_pred_df03, sd = carot_stat03[, 4])

min(carot_all03$pred.mean)
carot_neg <- length(which(carot_all03$pred.mean < 0))
carot_neg #3750 pixels have negative values
saveRDS(carot_all03, "inference/flg03/carot_all03.rds") #to load object : carot_all03 <- readRDS("inference/flg03/carot_all03.rds")

rm(carot_pred03, carot_pred_df03, carot_stat03)
gc()

carbon_pred03 <- apply.coefs(all_jack_coefs_list$carbon, val.spec = spec03_bn, intercept = T)
carbon_stat03 <- t(apply(carbon_pred03, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
carbon_pred_df03<-data.frame(pred.mean = carbon_stat03[, 1],
                             pred.low = carbon_stat03[, 2],
                             pred.high = carbon_stat03[, 3])
carbon_all03 <- cbind(carbon_pred_df03, sd = carbon_stat03[, 4])

min(carbon_all03$pred.mean)
carbon_neg <- length(which(carbon_all03$pred.mean < 0))
carbon_neg #0 pixels have negative values
saveRDS(carbon_all03, "inference/flg03/carbon_all03.rds") #to load object : carbon_all03 <- readRDS("inference/flg03/carbon_all03.rds")

rm(carbon_pred03, carbon_pred_df03, carbon_stat03)
gc()

nitrogen_pred03 <- apply.coefs(all_jack_coefs_list$nitrogen, val.spec = spec03_bn, intercept = T)
nitrogen_stat03 <- t(apply(nitrogen_pred03, 1,
                           function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
nitrogen_pred_df03 <- data.frame(pred.mean = nitrogen_stat03[, 1],
                                 pred.low = nitrogen_stat03[, 2],
                                 pred.high = nitrogen_stat03[, 3])
nitrogen_all03 <- cbind(nitrogen_pred_df03, sd = nitrogen_stat03[, 4])

min(nitrogen_all03$pred.mean)
nitrogen_neg <- length(which(nitrogen_all03$pred.mean < 0))
nitrogen_neg #8218 pixels have negative values
saveRDS(nitrogen_all03, "inference/flg03/nitrogen_all03.rds") #to load object : nitrogen_all03 <- readRDS("inference/flg03/nitrogen_all03.rds")

rm(nitrogen_pred03, nitrogen_pred_df03, nitrogen_stat03)
gc()

# Extract coordinates
x_coords03 <- spectra03$x
y_coords03 <- spectra03$y

# Add coordinates to results
SLA_all03$x <- x_coords03
SLA_all03$y <- y_coords03
LMA_all03$x <- x_coords03
LMA_all03$y <- y_coords03
LDMC_all03$x <- x_coords03
LDMC_all03$y <- y_coords03
EWT_all03$x <- x_coords03
EWT_all03$y <- y_coords03
hemi_all03$x <- x_coords03
hemi_all03$y <- y_coords03
cell_all03$x <- x_coords03
cell_all03$y <- y_coords03
lignin_all03$x <- x_coords03
lignin_all03$y <- y_coords03
chla_all03$x <- x_coords03
chla_all03$y <- y_coords03
chlb_all03$x <- x_coords03
chlb_all03$y <- y_coords03
carot_all03$x <- x_coords03
carot_all03$y <- y_coords03
carbon_all03$x <- x_coords03
carbon_all03$y <- y_coords03
nitrogen_all03$x <- x_coords03
nitrogen_all03$y <- y_coords03

# Select columns of interest for mapping
SLA_03 <- SLA_all03 %>%
  dplyr::select(pred.SLA = pred.mean, sd.SLA = sd, x, y) %>%
  `rownames<-`(NULL)
LMA_03 <- LMA_all03 %>%
  dplyr::select(pred.LMA = pred.mean, sd.LMA = sd, x, y) %>%
  `rownames<-`(NULL)
LDMC_03 <- LDMC_all03 %>%
  dplyr::select(pred.LDMC = pred.mean, sd.LDMC = sd, x, y) %>%
  `rownames<-`(NULL)
EWT_03 <- EWT_all03 %>%
  dplyr::select(pred.EWT = pred.mean, sd.EWT = sd, x, y) %>%
  `rownames<-`(NULL)
hemi_03 <- hemi_all03 %>%
  dplyr::select(pred.hemi = pred.mean, sd.hemi = sd, x, y) %>%
  `rownames<-`(NULL)
cell_03 <- cell_all03 %>%
  dplyr::select(pred.cell = pred.mean, sd.cell = sd, x, y) %>%
  `rownames<-`(NULL)
lignin_03 <- lignin_all03 %>%
  dplyr::select(pred.lignin = pred.mean, sd.lignin = sd, x, y) %>%
  `rownames<-`(NULL)
chla_03 <- chla_all03 %>%
  dplyr::select(pred.chla = pred.mean, sd.chla = sd, x, y) %>%
  `rownames<-`(NULL)
chlb_03 <- chlb_all03 %>%
  dplyr::select(pred.chlb = pred.mean, sd.chlb = sd, x, y) %>%
  `rownames<-`(NULL)
carot_03 <- carot_all03 %>%
  dplyr::select(pred.carot = pred.mean, sd.carot = sd, x, y) %>%
  `rownames<-`(NULL)
carbon_03 <- carbon_all03 %>%
  dplyr::select(pred.carbon = pred.mean, sd.carbon = sd, x, y) %>%
  `rownames<-`(NULL)
nitrogen_03 <- nitrogen_all03 %>%
  dplyr::select(pred.nitrogen = pred.mean, sd.nitrogen = sd, x, y) %>%
  `rownames<-`(NULL)

# Join all dataframes
df_list_03 = list(SLA_03, LMA_03, LDMC_03, EWT_03, hemi_03, cell_03, lignin_03, chla_03, chlb_03, carot_03, carbon_03, nitrogen_03) 

df_merged03 <- df_list_03 %>%
  reduce(right_join, by = c("x", "y"))
saveRDS(df_merged03, "inference/flg03/df_merged03.rds") #to load object : df_merged03 <- readRDS("inference/flg03/df_merged03.rds")

# Add relative uncertainties (CV)
df_merged03 <- df_merged03 %>%
  mutate(CV.SLA = sd.SLA/pred.SLA,
         CV.LMA = sd.LMA/pred.LMA,
         CV.LDMC = sd.LDMC/pred.LDMC,
         CV.EWT = sd.EWT/pred.EWT,
         CV.hemi = sd.hemi/pred.hemi,
         CV.cell = sd.cell/pred.cell,
         CV.lignin = sd.lignin/pred.lignin,
         CV.chla = sd.chla/pred.chla,
         CV.chlb = sd.chlb/pred.chlb,
         CV.carot = sd.carot/pred.carot,
         CV.carbon = sd.carbon/pred.carbon,
         CV.nitrogen = sd.nitrogen/pred.nitrogen)

# Create 3 dataframes (traits, sd and CV)
df_traits03 <- df_merged03 %>%
  dplyr::select(starts_with("pred"), x, y) %>%
  rename_with(~sub("pred\\.", "", .), starts_with("pred"))
df_sd03 <- df_merged03 %>%
  dplyr::select(starts_with("sd"), x, y) %>%
  rename_with(~sub("sd\\.", "", .), starts_with("sd"))
df_CV03 <- df_merged03 %>%
  dplyr::select(starts_with("CV"), x, y) %>%
  rename_with(~sub("CV\\.", "", .), starts_with("CV"))
saveRDS(df_traits03, "inference/flg03/df_traits03.rds")
saveRDS(df_sd03, "inference/flg03/df_sd03.rds")
saveRDS(df_CV03, "inference/flg03/df_CV03.rds")

#### Flightline 04####----
# Import spectra
spectra04 <- readRDS("output_allpixels/smooth_spectra04.rds")

spec04 <- as_spectra(spectra04, name_idx = 1, meta_idxs = c(2:3)) #change to spectra object (this takes a long time)

# Normalize spectra
spec04_bn <- t(apply(as.matrix(spec04), 1, function(x) x/sqrt(sum(x^2))))
spec04_bn <- spectra(spec04_bn, bands = bands(spec04), names = names(spec04))
spectrolab::meta(spec04_bn) <- spectrolab::meta(spec04)

rm(spec04)
gc()

saveRDS(spec04_bn, "output_allpixels/spec04_bn.rds") #to load object : spec04_bn <- readRDS("output_allpixels/spec04_bn.rds")

# Apply PLSR models
SLA_pred04 <- apply.coefs(all_jack_coefs_list$SLA, val.spec = spec04_bn, intercept = T)
SLA_stat04 <- t(apply(SLA_pred04, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
SLA_pred_df04 <- data.frame(pred.mean = SLA_stat04[, 1],
                            pred.low = SLA_stat04[, 2],
                            pred.high = SLA_stat04[, 3])
SLA_all04 <- cbind(SLA_pred_df04, sd = SLA_stat04[, 4])

min(SLA_all04$pred.mean)
SLA_neg <- length(which(SLA_all04$pred.mean < 0))
SLA_neg #14 236 pixels are negative values
saveRDS(SLA_all04, "inference/flg04/SLA_all04.rds") #to load object : SLA_all04 <- readRDS("inference/flg04/SLA_all04.rds")

rm(SLA_pred04, SLA_pred_df04, SLA_stat04)
gc()

LMA_pred04 <- apply.coefs(all_jack_coefs_list$LMA, val.spec = spec04_bn, intercept = T)
LMA_stat04 <- t(apply(LMA_pred04, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LMA_pred_df04 <- data.frame(pred.mean = LMA_stat04[, 1],
                            pred.low = LMA_stat04[, 2],
                            pred.high = LMA_stat04[, 3])
LMA_all04 <- cbind(LMA_pred_df04, sd = LMA_stat04[, 4])

min(LMA_all04$pred.mean)
LMA_neg <- length(which(LMA_all04$pred.mean < 0))
LMA_neg #25 670 pixels have negative values
saveRDS(LMA_all04, "inference/flg04/LMA_all04.rds") #to load object : LMA_all04 <- readRDS("inference/flg04/LMA_all04.rds")

rm(LMA_pred04, LMA_pred_df04, LMA_stat04)
gc()

LDMC_pred04 <- apply.coefs(all_jack_coefs_list$LDMC, val.spec = spec04_bn, intercept = T)
LDMC_stat04 <- t(apply(LDMC_pred04, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LDMC_pred_df04 <- data.frame(pred.mean = LDMC_stat04[, 1],
                             pred.low = LDMC_stat04[, 2],
                             pred.high = LDMC_stat04[, 3])
LDMC_all04 <- cbind(LDMC_pred_df04, sd = LDMC_stat04[, 4])

min(LDMC_all04$pred.mean)
LDMC_neg <- length(which(LDMC_all04$pred.mean < 0))
LDMC_neg #9 pixels have negative values
saveRDS(LDMC_all04, "inference/flg04/LDMC_all04.rds") #to load object : LDMC_all04 <- readRDS("inference/flg04/LDMC_all04.rds")

rm(LDMC_pred04, LDMC_pred_df04, LDMC_stat04)
gc()

EWT_pred04 <- apply.coefs(all_jack_coefs_list$EWT, val.spec = spec04_bn, intercept = T)
EWT_stat04 <- t(apply(EWT_pred04, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
EWT_pred_df04 <- data.frame(pred.mean = EWT_stat04[, 1],
                            pred.low = EWT_stat04[, 2],
                            pred.high = EWT_stat04[, 3])
EWT_all04 <- cbind(EWT_pred_df04, sd = EWT_stat04[, 4])

min(EWT_all04$pred.mean)
EWT_neg <- length(which(EWT_all04$pred.mean < 0))
EWT_neg #20 340 pixels have negative values
saveRDS(EWT_all04, "inference/flg04/EWT_all04.rds") #to load object : EWT_all04 <- readRDS("inference/flg04/EWT_all04.rds")

rm(EWT_pred04, EWT_pred_df04, EWT_stat04)
gc()

hemi_pred04 <- apply.coefs(all_jack_coefs_list$hemi, val.spec = spec04_bn, intercept = T)
hemi_stat04 <- t(apply(hemi_pred04, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
hemi_pred_df04 <- data.frame(pred.mean = hemi_stat04[, 1],
                             pred.low = hemi_stat04[, 2],
                             pred.high = hemi_stat04[, 3])
hemi_all04 <- cbind(hemi_pred_df04, sd = hemi_stat04[, 4])

min(hemi_all04$pred.mean)
hemi_neg <- length(which(hemi_all04$pred.mean < 0))
hemi_neg #1496 pixels have negative values
saveRDS(hemi_all04, "inference/flg04/hemi_all04.rds") #to load object : hemi_all04 <- readRDS("inference/flg04/hemi_all04.rds")

rm(hemi_pred04, hemi_pred_df04, hemi_stat04)
gc()

cell_pred04 <- apply.coefs(all_jack_coefs_list$cell, val.spec = spec04_bn, intercept = T)
cell_stat04 <- t(apply(cell_pred04, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
cell_pred_df04 <- data.frame(pred.mean = cell_stat04[, 1],
                             pred.low = cell_stat04[, 2],
                             pred.high = cell_stat04[, 3])
cell_all04 <- cbind(cell_pred_df04, sd = cell_stat04[, 4])

min(cell_all04$pred.mean)
cell_neg <- length(which(cell_all04$pred.mean < 0))
cell_neg #82 pixels have negative values
saveRDS(cell_all04, "inference/flg04/cell_all04.rds") #to load object : cell_all04 <- readRDS("inference/flg04/cell_all04.rds")

rm(cell_pred04, cell_pred_df04, cell_stat04)
gc()

lignin_pred04 <- apply.coefs(all_jack_coefs_list$lignin, val.spec = spec04_bn, intercept = T)
lignin_stat04 <- t(apply(lignin_pred04, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
lignin_pred_df04 <- data.frame(pred.mean = lignin_stat04[, 1],
                               pred.low = lignin_stat04[, 2],
                               pred.high = lignin_stat04[, 3])
lignin_all04 <- cbind(lignin_pred_df04, sd = lignin_stat04[, 4])

min(lignin_all04$pred.mean)
lignin_neg <- length(which(lignin_all04$pred.mean < 0))
lignin_neg #9 pixels have negative values
saveRDS(lignin_all04, "inference/flg04/lignin_all04.rds") #to load object : lignin_all04 <- readRDS("inference/flg04/lignin_all04.rds")

rm(lignin_pred04, lignin_pred_df04, lignin_stat04)
gc()

chla_pred04 <- apply.coefs(all_jack_coefs_list$chla, val.spec = spec04_bn, intercept = T)
chla_stat04 <- t(apply(chla_pred04, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chla_pred_df04 <- data.frame(pred.mean = chla_stat04[, 1],
                             pred.low = chla_stat04[, 2],
                             pred.high = chla_stat04[, 3])
chla_all04 <- cbind(chla_pred_df04, sd = chla_stat04[, 4])

min(chla_all04$pred.mean)
chla_neg <- length(which(chla_all04$pred.mean < 0))
chla_neg #3309 pixels have negative values
saveRDS(chla_all04, "inference/flg04/chla_all04.rds") #to load object : chla_all04 <- readRDS("inference/flg04/chla_all04.rds")

rm(chla_pred04, chla_pred_df04, chla_stat04)
gc()

chlb_pred04 <- apply.coefs(all_jack_coefs_list$chlb, val.spec = spec04_bn, intercept = T)
chlb_stat04 <- t(apply(chlb_pred04, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chlb_pred_df04 <- data.frame(pred.mean = chlb_stat04[, 1],
                             pred.low = chlb_stat04[, 2],
                             pred.high = chlb_stat04[, 3])
chlb_all04 <- cbind(chlb_pred_df04, sd = chlb_stat04[, 4])

min(chlb_all04$pred.mean)
chlb_neg <- length(which(chlb_all04$pred.mean < 0))
chlb_neg #2327 pixels have negative values
saveRDS(chlb_all04, "inference/flg04/chlb_all04.rds") #to load object : chlb_all04 <- readRDS("inference/flg04/chlb_all04.rds")

rm(chlb_pred04, chlb_pred_df04, chlb_stat04)
gc()

carot_pred04 <- apply.coefs(all_jack_coefs_list$carot, val.spec = spec04_bn, intercept = T)
carot_stat04 <- t(apply(carot_pred04, 1,
                        function(obs) c(mean(obs), quantile(obs, probs=c(0.025, 0.975)), sd(obs))))
carot_pred_df04 <- data.frame(pred.mean = carot_stat04[, 1],
                              pred.low = carot_stat04[, 2],
                              pred.high = carot_stat04[, 3])
carot_all04 <- cbind(carot_pred_df04, sd = carot_stat04[, 4])

min(carot_all04$pred.mean)
carot_neg <- length(which(carot_all04$pred.mean < 0))
carot_neg #3035 pixels have negative values
saveRDS(carot_all04, "inference/flg04/carot_all04.rds") #to load object : carot_all04 <- readRDS("inference/flg04/carot_all04.rds")

rm(carot_pred04, carot_pred_df04, carot_stat04)
gc()

carbon_pred04 <- apply.coefs(all_jack_coefs_list$carbon, val.spec = spec04_bn, intercept = T)
carbon_stat04 <- t(apply(carbon_pred04, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
carbon_pred_df04 <-data.frame(pred.mean = carbon_stat04[, 1],
                              pred.low = carbon_stat04[, 2],
                              pred.high = carbon_stat04[, 3])
carbon_all04 <- cbind(carbon_pred_df04, sd = carbon_stat04[, 4])

min(carbon_all04$pred.mean)
carbon_neg <- length(which(carbon_all04$pred.mean < 0))
carbon_neg #0 pixels have negative values
saveRDS(carbon_all04, "inference/flg04/carbon_all04.rds") #to load object : carbon_all04 <- readRDS("inference/flg04/carbon_all04.rds")

rm(carbon_pred04, carbon_pred_df04, carbon_stat04)
gc()

nitrogen_pred04 <- apply.coefs(all_jack_coefs_list$nitrogen, val.spec = spec04_bn, intercept = T)
nitrogen_stat04 <- t(apply(nitrogen_pred04, 1,
                           function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
nitrogen_pred_df04 <- data.frame(pred.mean = nitrogen_stat04[, 1],
                                 pred.low = nitrogen_stat04[, 2],
                                 pred.high = nitrogen_stat04[, 3])
nitrogen_all04 <- cbind(nitrogen_pred_df04, sd = nitrogen_stat04[, 4])

min(nitrogen_all04$pred.mean)
nitrogen_neg <- length(which(nitrogen_all04$pred.mean < 0))
nitrogen_neg #4873 pixels have negative values
saveRDS(nitrogen_all04, "inference/flg04/nitrogen_all04.rds") #to load object : nitrogen_all04 <- readRDS("inference/flg04/nitrogen_all04.rds")

rm(nitrogen_pred04, nitrogen_pred_df04, nitrogen_stat04)
gc()

# Extract coordinates
x_coords04 <- spectra04$x
y_coords04 <- spectra04$y

# Add coordinates to results
SLA_all04$x <- x_coords04
SLA_all04$y <- y_coords04
LMA_all04$x <- x_coords04
LMA_all04$y <- y_coords04
LDMC_all04$x <- x_coords04
LDMC_all04$y <- y_coords04
EWT_all04$x <- x_coords04
EWT_all04$y <- y_coords04
hemi_all04$x <- x_coords04
hemi_all04$y <- y_coords04
cell_all04$x <- x_coords04
cell_all04$y <- y_coords04
lignin_all04$x <- x_coords04
lignin_all04$y <- y_coords04
chla_all04$x <- x_coords04
chla_all04$y <- y_coords04
chlb_all04$x <- x_coords04
chlb_all04$y <- y_coords04
carot_all04$x <- x_coords04
carot_all04$y <- y_coords04
carbon_all04$x <- x_coords04
carbon_all04$y <- y_coords04
nitrogen_all04$x <- x_coords04
nitrogen_all04$y <- y_coords04

# Select columns of interest for mapping
SLA_04 <- SLA_all04 %>%
  dplyr::select(pred.SLA = pred.mean, sd.SLA = sd, x, y) %>%
  `rownames<-`(NULL)
LMA_04 <- LMA_all04 %>%
  dplyr::select(pred.LMA = pred.mean, sd.LMA = sd, x, y) %>%
  `rownames<-`(NULL)
LDMC_04 <- LDMC_all04 %>%
  dplyr::select(pred.LDMC = pred.mean, sd.LDMC = sd, x, y) %>%
  `rownames<-`(NULL)
EWT_04 <- EWT_all04 %>%
  dplyr::select(pred.EWT = pred.mean, sd.EWT = sd, x, y) %>%
  `rownames<-`(NULL)
hemi_04 <- hemi_all04 %>%
  dplyr::select(pred.hemi = pred.mean, sd.hemi = sd, x, y) %>%
  `rownames<-`(NULL)
cell_04 <- cell_all04 %>%
  dplyr::select(pred.cell = pred.mean, sd.cell = sd, x, y) %>%
  `rownames<-`(NULL)
lignin_04 <- lignin_all04 %>%
  dplyr::select(pred.lignin = pred.mean, sd.lignin = sd, x, y) %>%
  `rownames<-`(NULL)
chla_04 <- chla_all04 %>%
  dplyr::select(pred.chla = pred.mean, sd.chla = sd, x, y) %>%
  `rownames<-`(NULL)
chlb_04 <- chlb_all04 %>%
  dplyr::select(pred.chlb = pred.mean, sd.chlb = sd, x, y) %>%
  `rownames<-`(NULL)
carot_04 <- carot_all04 %>%
  dplyr::select(pred.carot = pred.mean, sd.carot = sd, x, y) %>%
  `rownames<-`(NULL)
carbon_04 <- carbon_all04 %>%
  dplyr::select(pred.carbon = pred.mean, sd.carbon = sd, x, y) %>%
  `rownames<-`(NULL)
nitrogen_04 <- nitrogen_all04 %>%
  dplyr::select(pred.nitrogen = pred.mean, sd.nitrogen = sd, x, y) %>%
  `rownames<-`(NULL)

# Join all dataframes
df_list_04 = list(SLA_04, LMA_04, LDMC_04, EWT_04, hemi_04, cell_04, lignin_04, chla_04, chlb_04, carot_04, carbon_04, nitrogen_04) 

df_merged04 <- df_list_04 %>%
  reduce(right_join, by = c("x", "y"))
saveRDS(df_merged04, "inference/flg04/df_merged04.rds") #to load object : df_merged04 <- readRDS("inference/flg04/df_merged04.rds")

# Add relative uncertainties (CV)
df_merged04 <- df_merged04 %>%
  mutate(CV.SLA = sd.SLA/pred.SLA,
         CV.LMA = sd.LMA/pred.LMA,
         CV.LDMC = sd.LDMC/pred.LDMC,
         CV.EWT = sd.EWT/pred.EWT,
         CV.hemi = sd.hemi/pred.hemi,
         CV.cell = sd.cell/pred.cell,
         CV.lignin = sd.lignin/pred.lignin,
         CV.chla = sd.chla/pred.chla,
         CV.chlb = sd.chlb/pred.chlb,
         CV.carot = sd.carot/pred.carot,
         CV.carbon = sd.carbon/pred.carbon,
         CV.nitrogen = sd.nitrogen/pred.nitrogen)

# Create 3 dataframes (traits,sd and CV)
df_traits04 <- df_merged04 %>%
  dplyr::select(starts_with("pred"), x, y) %>%
  rename_with(~sub("pred\\.", "", .), starts_with("pred"))
df_sd04 <- df_merged04 %>%
  dplyr::select(starts_with("sd"), x, y) %>%
  rename_with(~sub("sd\\.", "", .), starts_with("sd"))
df_CV04 <- df_merged04 %>%
  dplyr::select(starts_with("CV"), x, y) %>%
  rename_with(~sub("CV\\.", "", .), starts_with("CV"))
saveRDS(df_traits04, "inference/flg04/df_traits04.rds")
saveRDS(df_sd04, "inference/flg04/df_sd04.rds")
saveRDS(df_CV04, "inference/flg04/df_CV04.rds")

#### Flightline 05####----
# Import spectra
spectra05 <- readRDS("output_allpixels/smooth_spectra05.rds")

spec05 <- as_spectra(spectra05, name_idx = 1, meta_idxs = c(2:3)) #change to spectra object (this takes a long time)

# Normalize spectra
spec05_bn <- t(apply(as.matrix(spec05), 1, function(x) x/sqrt(sum(x^2))))
spec05_bn <- spectra(spec05_bn, bands = bands(spec05), names = names(spec05))
spectrolab::meta(spec05_bn) <- spectrolab::meta(spec05)

rm(spec05)
gc()

saveRDS(spec05_bn, "output_allpixels/spec05_bn.rds") #to load object : spec05_bn <- readRDS("output_allpixels/spec05_bn.rds")

# Apply PLSR models
SLA_pred05 <- apply.coefs(all_jack_coefs_list$SLA, val.spec = spec05_bn, intercept = T)
SLA_stat05 <- t(apply(SLA_pred05, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
SLA_pred_df05 <- data.frame(pred.mean = SLA_stat05[, 1],
                            pred.low = SLA_stat05[, 2],
                            pred.high = SLA_stat05[, 3])
SLA_all05 <- cbind(SLA_pred_df05, sd = SLA_stat05[, 4])

min(SLA_all05$pred.mean)
SLA_neg <- length(which(SLA_all05$pred.mean < 0))
SLA_neg #21 418 pixels are negative values
saveRDS(SLA_all05, "inference/flg05/SLA_all05.rds") #to load object : SLA_all05 <- readRDS("inference/flg05/SLA_all05.rds")

rm(SLA_pred05, SLA_pred_df05, SLA_stat05)
gc()

LMA_pred05 <- apply.coefs(all_jack_coefs_list$LMA, val.spec = spec05_bn, intercept = T)
LMA_stat05 <- t(apply(LMA_pred05, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LMA_pred_df05 <- data.frame(pred.mean = LMA_stat05[, 1],
                            pred.low = LMA_stat05[, 2],
                            pred.high = LMA_stat05[, 3])
LMA_all05 <- cbind(LMA_pred_df05, sd = LMA_stat05[, 4])

min(LMA_all05$pred.mean)
LMA_neg <- length(which(LMA_all05$pred.mean < 0))
LMA_neg #13 949 pixels have negative values
saveRDS(LMA_all05, "inference/flg05/LMA_all05.rds") #to load object : LMA_all05 <- readRDS("inference/flg05/LMA_all05.rds")

rm(LMA_pred05, LMA_pred_df05, LMA_stat05)
gc()

LDMC_pred05 <- apply.coefs(all_jack_coefs_list$LDMC, val.spec = spec05_bn, intercept = T)
LDMC_stat05 <- t(apply(LDMC_pred05, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LDMC_pred_df05 <- data.frame(pred.mean = LDMC_stat05[, 1],
                             pred.low = LDMC_stat05[, 2],
                             pred.high = LDMC_stat05[, 3])
LDMC_all05 <- cbind(LDMC_pred_df05, sd = LDMC_stat05[, 4])

min(LDMC_all05$pred.mean)
LDMC_neg <- length(which(LDMC_all05$pred.mean < 0))
LDMC_neg #2 pixels have negative values
saveRDS(LDMC_all05, "inference/flg05/LDMC_all05.rds") #to load object : LDMC_all05 <- readRDS("inference/flg05/LDMC_all05.rds")

rm(LDMC_pred05, LDMC_pred_df05, LDMC_stat05)
gc()

EWT_pred05 <- apply.coefs(all_jack_coefs_list$EWT, val.spec = spec05_bn, intercept = T)
EWT_stat05 <- t(apply(EWT_pred05, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
EWT_pred_df05 <- data.frame(pred.mean = EWT_stat05[, 1],
                            pred.low = EWT_stat05[, 2],
                            pred.high = EWT_stat05[, 3])
EWT_all05 <- cbind(EWT_pred_df05, sd = EWT_stat05[, 4])

min(EWT_all05$pred.mean)
EWT_neg <- length(which(EWT_all05$pred.mean < 0))
EWT_neg #10 487 pixels have negative values
saveRDS(EWT_all05, "inference/flg05/EWT_all05.rds") #to load object : EWT_all05 <- readRDS("inference/flg05/EWT_all05.rds")

rm(EWT_pred05, EWT_pred_df05, EWT_stat05)
gc()

hemi_pred05 <- apply.coefs(all_jack_coefs_list$hemi, val.spec = spec05_bn, intercept = T)
hemi_stat05 <- t(apply(hemi_pred05, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
hemi_pred_df05 <- data.frame(pred.mean = hemi_stat05[, 1],
                             pred.low = hemi_stat05[, 2],
                             pred.high = hemi_stat05[, 3])
hemi_all05 <- cbind(hemi_pred_df05, sd = hemi_stat05[, 4])

min(hemi_all05$pred.mean)
hemi_neg <- length(which(hemi_all05$pred.mean < 0))
hemi_neg #7143 pixels have negative values
saveRDS(hemi_all05, "inference/flg05/hemi_all05.rds") #to load object : hemi_all05 <- readRDS("inference/flg05/hemi_all05.rds")

rm(hemi_pred05, hemi_pred_df05, hemi_stat05)
gc()

cell_pred05 <- apply.coefs(all_jack_coefs_list$cell, val.spec = spec05_bn, intercept = T)
cell_stat05 <- t(apply(cell_pred05, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
cell_pred_df05 <- data.frame(pred.mean = cell_stat05[, 1],
                             pred.low = cell_stat05[, 2],
                             pred.high = cell_stat05[, 3])
cell_all05 <- cbind(cell_pred_df05, sd = cell_stat05[, 4])

min(cell_all05$pred.mean)
cell_neg <- length(which(cell_all05$pred.mean < 0))
cell_neg #34 pixels have negative values
saveRDS(cell_all05, "inference/flg05/cell_all05.rds") #to load object : cell_all05 <- readRDS("inference/flg05/cell_all05.rds")

rm(cell_pred05, cell_pred_df05, cell_stat05)
gc()

lignin_pred05 <- apply.coefs(all_jack_coefs_list$lignin, val.spec = spec05_bn, intercept = T)
lignin_stat05 <- t(apply(lignin_pred05, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
lignin_pred_df05 <- data.frame(pred.mean = lignin_stat05[, 1],
                               pred.low = lignin_stat05[, 2],
                               pred.high = lignin_stat05[, 3])
lignin_all05 <- cbind(lignin_pred_df05, sd = lignin_stat05[, 4])

min(lignin_all05$pred.mean)
lignin_neg <- length(which(lignin_all05$pred.mean < 0))
lignin_neg #9 pixels have negative values
saveRDS(lignin_all05, "inference/flg05/lignin_all05.rds") #to load object : lignin_all05 <- readRDS("inference/flg05/lignin_all05.rds")

rm(lignin_pred05, lignin_pred_df05, lignin_stat05)
gc()

chla_pred05 <- apply.coefs(all_jack_coefs_list$chla, val.spec = spec05_bn, intercept = T)
chla_stat05 <- t(apply(chla_pred05, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chla_pred_df05 <- data.frame(pred.mean = chla_stat05[, 1],
                             pred.low = chla_stat05[, 2],
                             pred.high = chla_stat05[, 3])
chla_all05 <- cbind(chla_pred_df05, sd = chla_stat05[, 4])

min(chla_all05$pred.mean)
chla_neg <- length(which(chla_all05$pred.mean < 0))
chla_neg #3355 pixels have negative values
saveRDS(chla_all05, "inference/flg05/chla_all05.rds") #to load object : chla_all05 <- readRDS("inference/flg05/chla_all05.rds")

rm(chla_pred05, chla_pred_df05, chla_stat05)
gc()

chlb_pred05 <- apply.coefs(all_jack_coefs_list$chlb, val.spec = spec05_bn, intercept = T)
chlb_stat05 <- t(apply(chlb_pred05, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chlb_pred_df05 <- data.frame(pred.mean = chlb_stat05[, 1],
                             pred.low = chlb_stat05[, 2],
                             pred.high = chlb_stat05[, 3])
chlb_all05 <- cbind(chlb_pred_df05, sd = chlb_stat05[, 4])

min(chlb_all05$pred.mean)
chlb_neg <- length(which(chlb_all05$pred.mean < 0))
chlb_neg #2443 pixels have negative values
saveRDS(chlb_all05, "inference/flg05/chlb_all05.rds") #to load object : chlb_all05 <- readRDS("inference/flg05/chlb_all05.rds")

rm(chlb_pred05, chlb_pred_df05, chlb_stat05)
gc()

carot_pred05 <- apply.coefs(all_jack_coefs_list$carot, val.spec = spec05_bn, intercept = T)
carot_stat05 <- t(apply(carot_pred05, 1,
                        function(obs) c(mean(obs), quantile(obs, probs=c(0.025, 0.975)), sd(obs))))
carot_pred_df05 <- data.frame(pred.mean = carot_stat05[, 1],
                              pred.low = carot_stat05[, 2],
                              pred.high = carot_stat05[, 3])
carot_all05 <- cbind(carot_pred_df05, sd = carot_stat05[, 4])

min(carot_all05$pred.mean)
carot_neg <- length(which(carot_all05$pred.mean < 0))
carot_neg #3163 pixels have negative values
saveRDS(carot_all05, "inference/flg05/carot_all05.rds") #to load object : carot_all05 <- readRDS("inference/flg05/carot_all05.rds")

rm(carot_pred05, carot_pred_df05, carot_stat05)
gc()

carbon_pred05 <- apply.coefs(all_jack_coefs_list$carbon, val.spec = spec05_bn, intercept = T)
carbon_stat05 <- t(apply(carbon_pred05, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
carbon_pred_df05 <-data.frame(pred.mean = carbon_stat05[, 1],
                              pred.low = carbon_stat05[, 2],
                              pred.high = carbon_stat05[, 3])
carbon_all05 <- cbind(carbon_pred_df05, sd = carbon_stat05[, 4])

min(carbon_all05$pred.mean)

carbon_neg <- length(which(carbon_all05$pred.mean < 0))
carbon_neg #0 pixels have negative values
saveRDS(carbon_all05, "inference/flg05/carbon_all05.rds") #to load object : carbon_all05 <- readRDS("inference/flg05/carbon_all05.rds")

rm(carbon_pred05, carbon_pred_df05, carbon_stat05)
gc()

nitrogen_pred05 <- apply.coefs(all_jack_coefs_list$nitrogen, val.spec = spec05_bn, intercept = T)
nitrogen_stat05 <- t(apply(nitrogen_pred05, 1,
                           function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
nitrogen_pred_df05 <- data.frame(pred.mean = nitrogen_stat05[, 1],
                                 pred.low = nitrogen_stat05[, 2],
                                 pred.high = nitrogen_stat05[, 3])
nitrogen_all05 <- cbind(nitrogen_pred_df05, sd = nitrogen_stat05[, 4])

min(nitrogen_all05$pred.mean)
nitrogen_neg <- length(which(nitrogen_all05$pred.mean < 0))
nitrogen_neg #7291 pixels have negative values
saveRDS(nitrogen_all05, "inference/flg05/nitrogen_all05.rds") #to load object : nitrogen_all05 <- readRDS("inference/flg05/nitrogen_all05.rds")

rm(nitrogen_pred05, nitrogen_pred_df05, nitrogen_stat05)
gc()

# Extract coordinates
x_coords05 <- spectra05$x
y_coords05 <- spectra05$y

# Add coordinates to results
SLA_all05$x <- x_coords05
SLA_all05$y <- y_coords05
LMA_all05$x <- x_coords05
LMA_all05$y <- y_coords05
LDMC_all05$x <- x_coords05
LDMC_all05$y <- y_coords05
EWT_all05$x <- x_coords05
EWT_all05$y <- y_coords05
hemi_all05$x <- x_coords05
hemi_all05$y <- y_coords05
cell_all05$x <- x_coords05
cell_all05$y <- y_coords05
lignin_all05$x <- x_coords05
lignin_all05$y <- y_coords05
chla_all05$x <- x_coords05
chla_all05$y <- y_coords05
chlb_all05$x <- x_coords05
chlb_all05$y <- y_coords05
carot_all05$x <- x_coords05
carot_all05$y <- y_coords05
carbon_all05$x <- x_coords05
carbon_all05$y <- y_coords05
nitrogen_all05$x <- x_coords05
nitrogen_all05$y <- y_coords05

# Select columns of interest for mapping
SLA_05 <- SLA_all05 %>%
  dplyr::select(pred.SLA = pred.mean, sd.SLA = sd, x, y) %>%
  `rownames<-`(NULL)
LMA_05 <- LMA_all05 %>%
  dplyr::select(pred.LMA = pred.mean, sd.LMA = sd, x, y) %>%
  `rownames<-`(NULL)
LDMC_05 <- LDMC_all05 %>%
  dplyr::select(pred.LDMC = pred.mean, sd.LDMC = sd, x, y) %>%
  `rownames<-`(NULL)
EWT_05 <- EWT_all05 %>%
  dplyr::select(pred.EWT = pred.mean, sd.EWT = sd, x, y) %>%
  `rownames<-`(NULL)
hemi_05 <- hemi_all05 %>%
  dplyr::select(pred.hemi = pred.mean, sd.hemi = sd, x, y) %>%
  `rownames<-`(NULL)
cell_05 <- cell_all05 %>%
  dplyr::select(pred.cell = pred.mean, sd.cell = sd, x, y) %>%
  `rownames<-`(NULL)
lignin_05 <- lignin_all05 %>%
  dplyr::select(pred.lignin = pred.mean, sd.lignin = sd, x, y) %>%
  `rownames<-`(NULL)
chla_05 <- chla_all05 %>%
  dplyr::select(pred.chla = pred.mean, sd.chla = sd, x, y) %>%
  `rownames<-`(NULL)
chlb_05 <- chlb_all05 %>%
  dplyr::select(pred.chlb = pred.mean, sd.chlb = sd, x, y) %>%
  `rownames<-`(NULL)
carot_05 <- carot_all05 %>%
  dplyr::select(pred.carot = pred.mean, sd.carot = sd, x, y) %>%
  `rownames<-`(NULL)
carbon_05 <- carbon_all05 %>%
  dplyr::select(pred.carbon = pred.mean, sd.carbon = sd, x, y) %>%
  `rownames<-`(NULL)
nitrogen_05 <- nitrogen_all05 %>%
  dplyr::select(pred.nitrogen = pred.mean, sd.nitrogen = sd, x, y) %>%
  `rownames<-`(NULL)

# Join all dataframes
df_list_05 = list(SLA_05, LMA_05, LDMC_05, EWT_05, hemi_05, cell_05, lignin_05, chla_05, chlb_05, carot_05, carbon_05, nitrogen_05) 

df_merged05 <- df_list_05 %>%
  reduce(right_join, by = c("x", "y"))
saveRDS(df_merged05, "inference/flg05/df_merged05.rds") #to load object : df_merged05 <- readRDS("inference/flg05/df_merged05.rds")

# Add relative uncertainties (CV)
df_merged05 <- df_merged05 %>%
  mutate(CV.SLA = sd.SLA/pred.SLA,
         CV.LMA = sd.LMA/pred.LMA,
         CV.LDMC = sd.LDMC/pred.LDMC,
         CV.EWT = sd.EWT/pred.EWT,
         CV.hemi = sd.hemi/pred.hemi,
         CV.cell = sd.cell/pred.cell,
         CV.lignin = sd.lignin/pred.lignin,
         CV.chla = sd.chla/pred.chla,
         CV.chlb = sd.chlb/pred.chlb,
         CV.carot = sd.carot/pred.carot,
         CV.carbon = sd.carbon/pred.carbon,
         CV.nitrogen = sd.nitrogen/pred.nitrogen)

# Create 3 dataframes (traits,sd and CV)
df_traits05 <- df_merged05 %>%
  dplyr::select(starts_with("pred"), x, y) %>%
  rename_with(~sub("pred\\.", "", .), starts_with("pred"))
df_sd05 <- df_merged05 %>%
  dplyr::select(starts_with("sd"), x, y) %>%
  rename_with(~sub("sd\\.", "", .), starts_with("sd"))
df_CV05 <- df_merged05 %>%
  dplyr::select(starts_with("CV"), x, y) %>%
  rename_with(~sub("CV\\.", "", .), starts_with("CV"))
saveRDS(df_traits05, "inference/flg05/df_traits05.rds")
saveRDS(df_sd05, "inference/flg05/df_sd05.rds")
saveRDS(df_CV05, "inference/flg05/df_CV05.rds")

#### Flightline 06####----
# Import spectra
spectra06 <- readRDS("output_allpixels/smooth_spectra06.rds")

spec06 <- as_spectra(spectra06, name_idx = 1, meta_idxs = c(2:3)) #change to spectra object (this takes a long time)

# Normalize spectra
spec06_bn <- t(apply(as.matrix(spec06), 1, function(x) x/sqrt(sum(x^2))))
spec06_bn <- spectra(spec06_bn, bands = bands(spec06), names = names(spec06))
spectrolab::meta(spec06_bn) <- spectrolab::meta(spec06)

rm(spec06)
gc()

saveRDS(spec06_bn, "output_allpixels/spec06_bn.rds") #to load object : spec06_bn <- readRDS("output_allpixels/spec06_bn.rds")

# Apply PLSR models
SLA_pred06 <- apply.coefs(all_jack_coefs_list$SLA, val.spec = spec06_bn, intercept = T)
SLA_stat06 <- t(apply(SLA_pred06, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
SLA_pred_df06 <- data.frame(pred.mean = SLA_stat06[, 1],
                            pred.low = SLA_stat06[, 2],
                            pred.high = SLA_stat06[, 3])
SLA_all06 <- cbind(SLA_pred_df06, sd = SLA_stat06[, 4])

min(SLA_all06$pred.mean)
SLA_neg <- length(which(SLA_all06$pred.mean < 0))
SLA_neg #18 760 pixels are negative values
saveRDS(SLA_all06, "inference/flg06/SLA_all06.rds") #to load object : SLA_all06 <- readRDS("inference/flg06/SLA_all06.rds")

rm(SLA_pred06, SLA_pred_df06, SLA_stat06)
gc()

LMA_pred06 <- apply.coefs(all_jack_coefs_list$LMA, val.spec = spec06_bn, intercept = T)
LMA_stat06 <- t(apply(LMA_pred06, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LMA_pred_df06 <- data.frame(pred.mean = LMA_stat06[, 1],
                            pred.low = LMA_stat06[, 2],
                            pred.high = LMA_stat06[, 3])
LMA_all06 <- cbind(LMA_pred_df06, sd = LMA_stat06[, 4])

min(LMA_all06$pred.mean)
LMA_neg <- length(which(LMA_all06$pred.mean < 0))
LMA_neg #13 718 pixels have negative values
saveRDS(LMA_all06, "inference/flg06/LMA_all06.rds") #to load object : LMA_all06 <- readRDS("inference/flg06/LMA_all06.rds")

rm(LMA_pred06, LMA_pred_df06, LMA_stat06)
gc()

LDMC_pred06 <- apply.coefs(all_jack_coefs_list$LDMC, val.spec = spec06_bn, intercept = T)
LDMC_stat06 <- t(apply(LDMC_pred06, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LDMC_pred_df06 <- data.frame(pred.mean = LDMC_stat06[, 1],
                             pred.low = LDMC_stat06[, 2],
                             pred.high = LDMC_stat06[, 3])
LDMC_all06 <- cbind(LDMC_pred_df06, sd = LDMC_stat06[, 4])

min(LDMC_all06$pred.mean)
LDMC_neg <- length(which(LDMC_all06$pred.mean < 0))
LDMC_neg #11 pixels have negative values
saveRDS(LDMC_all06, "inference/flg06/LDMC_all06.rds") #to load object : LDMC_all06 <- readRDS("inference/flg06/LDMC_all06.rds")

rm(LDMC_pred06, LDMC_pred_df06, LDMC_stat06)
gc()

EWT_pred06 <- apply.coefs(all_jack_coefs_list$EWT, val.spec = spec06_bn, intercept = T)
EWT_stat06 <- t(apply(EWT_pred06, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
EWT_pred_df06 <- data.frame(pred.mean = EWT_stat06[, 1],
                            pred.low = EWT_stat06[, 2],
                            pred.high = EWT_stat06[, 3])
EWT_all06 <- cbind(EWT_pred_df06, sd = EWT_stat06[, 4])

min(EWT_all06$pred.mean)
EWT_neg <- length(which(EWT_all06$pred.mean < 0))
EWT_neg #19 291 pixels have negative values
saveRDS(EWT_all06, "inference/flg06/EWT_all06.rds") #to load object : EWT_all06 <- readRDS("inference/flg06/EWT_all06.rds")

rm(EWT_pred06, EWT_pred_df06, EWT_stat06)
gc()

hemi_pred06 <- apply.coefs(all_jack_coefs_list$hemi, val.spec = spec06_bn, intercept = T)
hemi_stat06 <- t(apply(hemi_pred06, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
hemi_pred_df06 <- data.frame(pred.mean = hemi_stat06[, 1],
                             pred.low = hemi_stat06[, 2],
                             pred.high = hemi_stat06[, 3])
hemi_all06 <- cbind(hemi_pred_df06, sd = hemi_stat06[, 4])

min(hemi_all06$pred.mean)
hemi_neg <- length(which(hemi_all06$pred.mean < 0))
hemi_neg #5043 pixels have negative values
saveRDS(hemi_all06, "inference/flg06/hemi_all06.rds") #to load object : hemi_all06 <- readRDS("inference/flg06/hemi_all06.rds")

rm(hemi_pred06, hemi_pred_df06, hemi_stat06)
gc()

cell_pred06 <- apply.coefs(all_jack_coefs_list$cell, val.spec = spec06_bn, intercept = T)
cell_stat06 <- t(apply(cell_pred06, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
cell_pred_df06 <- data.frame(pred.mean = cell_stat06[, 1],
                             pred.low = cell_stat06[, 2],
                             pred.high = cell_stat06[, 3])
cell_all06 <- cbind(cell_pred_df06, sd = cell_stat06[, 4])

min(cell_all06$pred.mean)
cell_neg <- length(which(cell_all06$pred.mean < 0))
cell_neg #50 pixels have negative values
saveRDS(cell_all06, "inference/flg06/cell_all06.rds") #to load object : cell_all06 <- readRDS("inference/flg06/cell_all06.rds")

rm(cell_pred06, cell_pred_df06, cell_stat06)
gc()

lignin_pred06 <- apply.coefs(all_jack_coefs_list$lignin, val.spec = spec06_bn, intercept = T)
lignin_stat06 <- t(apply(lignin_pred06, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
lignin_pred_df06 <- data.frame(pred.mean = lignin_stat06[, 1],
                               pred.low = lignin_stat06[, 2],
                               pred.high = lignin_stat06[, 3])
lignin_all06 <- cbind(lignin_pred_df06, sd = lignin_stat06[, 4])

min(lignin_all06$pred.mean)
lignin_neg <- length(which(lignin_all06$pred.mean < 0))
lignin_neg #58 pixels have negative values
saveRDS(lignin_all06, "inference/flg06/lignin_all06.rds") #to load object : lignin_all06 <- readRDS("inference/flg06/lignin_all06.rds")

rm(lignin_pred06, lignin_pred_df06, lignin_stat06)
gc()

chla_pred06 <- apply.coefs(all_jack_coefs_list$chla, val.spec = spec06_bn, intercept = T)
chla_stat06 <- t(apply(chla_pred06, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chla_pred_df06 <- data.frame(pred.mean = chla_stat06[, 1],
                             pred.low = chla_stat06[, 2],
                             pred.high = chla_stat06[, 3])
chla_all06 <- cbind(chla_pred_df06, sd = chla_stat06[, 4])

min(chla_all06$pred.mean)
chla_neg <- length(which(chla_all06$pred.mean < 0))
chla_neg #4960 pixels have negative values
saveRDS(chla_all06, "inference/flg06/chla_all06.rds") #to load object : chla_all06 <- readRDS("inference/flg06/chla_all06.rds")

rm(chla_pred06, chla_pred_df06, chla_stat06)
gc()

chlb_pred06 <- apply.coefs(all_jack_coefs_list$chlb, val.spec = spec06_bn, intercept = T)
chlb_stat06 <- t(apply(chlb_pred06, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chlb_pred_df06 <- data.frame(pred.mean = chlb_stat06[, 1],
                             pred.low = chlb_stat06[, 2],
                             pred.high = chlb_stat06[, 3])
chlb_all06 <- cbind(chlb_pred_df06, sd = chlb_stat06[, 4])

min(chlb_all06$pred.mean)
chlb_neg <- length(which(chlb_all06$pred.mean < 0))
chlb_neg #3753 pixels have negative values
saveRDS(chlb_all06, "inference/flg06/chlb_all06.rds") #to load object : chlb_all06 <- readRDS("inference/flg06/chlb_all06.rds")

rm(chlb_pred06, chlb_pred_df06, chlb_stat06)
gc()

carot_pred06 <- apply.coefs(all_jack_coefs_list$carot, val.spec = spec06_bn, intercept = T)
carot_stat06 <- t(apply(carot_pred06, 1,
                        function(obs) c(mean(obs), quantile(obs, probs=c(0.025, 0.975)), sd(obs))))
carot_pred_df06 <- data.frame(pred.mean = carot_stat06[, 1],
                              pred.low = carot_stat06[, 2],
                              pred.high = carot_stat06[, 3])
carot_all06 <- cbind(carot_pred_df06, sd = carot_stat06[, 4])

min(carot_all06$pred.mean)
carot_neg <- length(which(carot_all06$pred.mean < 0))
carot_neg #3122 pixels have negative values
saveRDS(carot_all06, "inference/flg06/carot_all06.rds") #to load object : carot_all06 <- readRDS("inference/flg06/carot_all06.rds")

rm(carot_pred06, carot_pred_df06, carot_stat06)
gc()

carbon_pred06 <- apply.coefs(all_jack_coefs_list$carbon, val.spec = spec06_bn, intercept = T)
carbon_stat06 <- t(apply(carbon_pred06, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
carbon_pred_df06 <-data.frame(pred.mean = carbon_stat06[, 1],
                              pred.low = carbon_stat06[, 2],
                              pred.high = carbon_stat06[, 3])
carbon_all06 <- cbind(carbon_pred_df06, sd = carbon_stat06[, 4])

min(carbon_all06$pred.mean)
carbon_neg <- length(which(carbon_all06$pred.mean < 0))
carbon_neg #0 pixels have negative values
saveRDS(carbon_all06, "inference/flg06/carbon_all06.rds") #to load object : carbon_all06 <- readRDS("inference/flg06/carbon_all06.rds")

rm(carbon_pred06, carbon_pred_df06, carbon_stat06)
gc()

nitrogen_pred06 <- apply.coefs(all_jack_coefs_list$nitrogen, val.spec = spec06_bn, intercept = T)
nitrogen_stat06 <- t(apply(nitrogen_pred06, 1,
                           function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
nitrogen_pred_df06 <- data.frame(pred.mean = nitrogen_stat06[, 1],
                                 pred.low = nitrogen_stat06[, 2],
                                 pred.high = nitrogen_stat06[, 3])
nitrogen_all06 <- cbind(nitrogen_pred_df06, sd = nitrogen_stat06[, 4])

min(nitrogen_all06$pred.mean)
nitrogen_neg <- length(which(nitrogen_all06$pred.mean < 0))
nitrogen_neg #5460 pixels have negative values
saveRDS(nitrogen_all06, "inference/flg06/nitrogen_all06.rds") #to load object : nitrogen_all06 <- readRDS("inference/flg06/nitrogen_all06.rds")

rm(nitrogen_pred06, nitrogen_pred_df06, nitrogen_stat06)
gc()

# Extract coordinates
x_coords06 <- spectra06$x
y_coords06 <- spectra06$y

# Add coordinates to results
SLA_all06$x <- x_coords06
SLA_all06$y <- y_coords06
LMA_all06$x <- x_coords06
LMA_all06$y <- y_coords06
LDMC_all06$x <- x_coords06
LDMC_all06$y <- y_coords06
EWT_all06$x <- x_coords06
EWT_all06$y <- y_coords06
hemi_all06$x <- x_coords06
hemi_all06$y <- y_coords06
cell_all06$x <- x_coords06
cell_all06$y <- y_coords06
lignin_all06$x <- x_coords06
lignin_all06$y <- y_coords06
chla_all06$x <- x_coords06
chla_all06$y <- y_coords06
chlb_all06$x <- x_coords06
chlb_all06$y <- y_coords06
carot_all06$x <- x_coords06
carot_all06$y <- y_coords06
carbon_all06$x <- x_coords06
carbon_all06$y <- y_coords06
nitrogen_all06$x <- x_coords06
nitrogen_all06$y <- y_coords06

# Select columns of interest for mapping
SLA_06 <- SLA_all06 %>%
  dplyr::select(pred.SLA = pred.mean, sd.SLA = sd, x, y) %>%
  `rownames<-`(NULL)
LMA_06 <- LMA_all06 %>%
  dplyr::select(pred.LMA = pred.mean, sd.LMA = sd, x, y) %>%
  `rownames<-`(NULL)
LDMC_06 <- LDMC_all06 %>%
  dplyr::select(pred.LDMC = pred.mean, sd.LDMC = sd, x, y) %>%
  `rownames<-`(NULL)
EWT_06 <- EWT_all06 %>%
  dplyr::select(pred.EWT = pred.mean, sd.EWT = sd, x, y) %>%
  `rownames<-`(NULL)
hemi_06 <- hemi_all06 %>%
  dplyr::select(pred.hemi = pred.mean, sd.hemi = sd, x, y) %>%
  `rownames<-`(NULL)
cell_06 <- cell_all06 %>%
  dplyr::select(pred.cell = pred.mean, sd.cell = sd, x, y) %>%
  `rownames<-`(NULL)
lignin_06 <- lignin_all06 %>%
  dplyr::select(pred.lignin = pred.mean, sd.lignin = sd, x, y) %>%
  `rownames<-`(NULL)
chla_06 <- chla_all06 %>%
  dplyr::select(pred.chla = pred.mean, sd.chla = sd, x, y) %>%
  `rownames<-`(NULL)
chlb_06 <- chlb_all06 %>%
  dplyr::select(pred.chlb = pred.mean, sd.chlb = sd, x, y) %>%
  `rownames<-`(NULL)
carot_06 <- carot_all06 %>%
  dplyr::select(pred.carot = pred.mean, sd.carot = sd, x, y) %>%
  `rownames<-`(NULL)
carbon_06 <- carbon_all06 %>%
  dplyr::select(pred.carbon = pred.mean, sd.carbon = sd, x, y) %>%
  `rownames<-`(NULL)
nitrogen_06 <- nitrogen_all06 %>%
  dplyr::select(pred.nitrogen = pred.mean, sd.nitrogen = sd, x, y) %>%
  `rownames<-`(NULL)

# Join all dataframes
df_list_06 = list(SLA_06, LMA_06, LDMC_06, EWT_06, hemi_06, cell_06, lignin_06, chla_06, chlb_06, carot_06, carbon_06, nitrogen_06) 

df_merged06 <- df_list_06 %>%
  reduce(right_join, by = c("x", "y"))
saveRDS(df_merged06, "inference/flg06/df_merged06.rds") #to load object : df_merged06 <- readRDS("inference/flg06/df_merged06.rds")

# Add relative uncertainties (CV)
df_merged06 <- df_merged06 %>%
  mutate(CV.SLA = sd.SLA/pred.SLA,
         CV.LMA = sd.LMA/pred.LMA,
         CV.LDMC = sd.LDMC/pred.LDMC,
         CV.EWT = sd.EWT/pred.EWT,
         CV.hemi = sd.hemi/pred.hemi,
         CV.cell = sd.cell/pred.cell,
         CV.lignin = sd.lignin/pred.lignin,
         CV.chla = sd.chla/pred.chla,
         CV.chlb = sd.chlb/pred.chlb,
         CV.carot = sd.carot/pred.carot,
         CV.carbon = sd.carbon/pred.carbon,
         CV.nitrogen = sd.nitrogen/pred.nitrogen)

# Create two dataframes (traits, sd and CV)
df_traits06 <- df_merged06 %>%
  dplyr::select(starts_with("pred"), x, y) %>%
  rename_with(~sub("pred\\.", "", .), starts_with("pred"))
df_sd06 <- df_merged06 %>%
  dplyr::select(starts_with("sd"), x, y) %>%
  rename_with(~sub("sd\\.", "", .), starts_with("sd"))
df_CV06 <- df_merged06 %>%
  dplyr::select(starts_with("CV"), x, y) %>%
  rename_with(~sub("CV\\.", "", .), starts_with("CV"))
saveRDS(df_traits06, "inference/flg06/df_traits06.rds")
saveRDS(df_sd06, "inference/flg06/df_sd06.rds")
saveRDS(df_CV06, "inference/flg06/df_CV06.rds")

#### Flightline 07####----
# Import spectra
spectra07 <- readRDS("output_allpixels/smooth_spectra07.rds")

spec07 <- as_spectra(spectra07, name_idx = 1, meta_idxs = c(2:3)) #change to spectra object (this takes a long time)

# Normalize spectra
spec07_bn <- t(apply(as.matrix(spec07), 1, function(x) x/sqrt(sum(x^2))))
spec07_bn <- spectra(spec07_bn, bands = bands(spec07), names = names(spec07))
spectrolab::meta(spec07_bn) <- spectrolab::meta(spec07)

rm(spec07)
gc()

saveRDS(spec07_bn, "output_allpixels/spec07_bn.rds") #to load object : spec07_bn <- readRDS("output_allpixels/spec07_bn.rds")

# Apply PLSR models
SLA_pred07 <- apply.coefs(all_jack_coefs_list$SLA, val.spec = spec07_bn, intercept = T)
SLA_stat07 <- t(apply(SLA_pred07, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
SLA_pred_df07 <- data.frame(pred.mean = SLA_stat07[, 1],
                            pred.low = SLA_stat07[, 2],
                            pred.high = SLA_stat07[, 3])
SLA_all07 <- cbind(SLA_pred_df07, sd = SLA_stat07[, 4])

min(SLA_all07$pred.mean)
SLA_neg <- length(which(SLA_all07$pred.mean < 0))
SLA_neg #14 238 pixels are negative values
saveRDS(SLA_all07, "inference/flg07/SLA_all07.rds") #to load object : SLA_all07 <- readRDS("inference/flg07/SLA_all07.rds")

rm(SLA_pred07, SLA_pred_df07, SLA_stat07)
gc()

LMA_pred07 <- apply.coefs(all_jack_coefs_list$LMA, val.spec = spec07_bn, intercept = T)
LMA_stat07 <- t(apply(LMA_pred07, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LMA_pred_df07 <- data.frame(pred.mean = LMA_stat07[, 1],
                            pred.low = LMA_stat07[, 2],
                            pred.high = LMA_stat07[, 3])
LMA_all07 <- cbind(LMA_pred_df07, sd = LMA_stat07[, 4])

min(LMA_all07$pred.mean)
LMA_neg <- length(which(LMA_all07$pred.mean < 0))
LMA_neg #11 268 pixels have negative values
saveRDS(LMA_all07, "inference/flg07/LMA_all07.rds") #to load object : LMA_all07 <- readRDS("inference/flg07/LMA_all07.rds")

rm(LMA_pred07, LMA_pred_df07, LMA_stat07)
gc()

LDMC_pred07 <- apply.coefs(all_jack_coefs_list$LDMC, val.spec = spec07_bn, intercept = T)
LDMC_stat07 <- t(apply(LDMC_pred07, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
LDMC_pred_df07 <- data.frame(pred.mean = LDMC_stat07[, 1],
                             pred.low = LDMC_stat07[, 2],
                             pred.high = LDMC_stat07[, 3])
LDMC_all07 <- cbind(LDMC_pred_df07, sd = LDMC_stat07[, 4])

min(LDMC_all07$pred.mean)
LDMC_neg <- length(which(LDMC_all07$pred.mean < 0))
LDMC_neg #11 pixels have negative values
saveRDS(LDMC_all07, "inference/flg07/LDMC_all07.rds") #to load object : LDMC_all07 <- readRDS("inference/flg07/LDMC_all07.rds")

rm(LDMC_pred07, LDMC_pred_df07, LDMC_stat07)
gc()

EWT_pred07 <- apply.coefs(all_jack_coefs_list$EWT, val.spec = spec07_bn, intercept = T)
EWT_stat07 <- t(apply(EWT_pred07, 1,
                      function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
EWT_pred_df07 <- data.frame(pred.mean = EWT_stat07[, 1],
                            pred.low = EWT_stat07[, 2],
                            pred.high = EWT_stat07[, 3])
EWT_all07 <- cbind(EWT_pred_df07, sd = EWT_stat07[, 4])

min(EWT_all07$pred.mean)
EWT_neg <- length(which(EWT_all07$pred.mean < 0))
EWT_neg #8882 pixels have negative values
saveRDS(EWT_all07, "inference/flg07/EWT_all07.rds") #to load object : EWT_all07 <- readRDS("inference/flg07/EWT_all07.rds")

rm(EWT_pred07, EWT_pred_df07, EWT_stat07)
gc()

hemi_pred07 <- apply.coefs(all_jack_coefs_list$hemi, val.spec = spec07_bn, intercept = T)
hemi_stat07 <- t(apply(hemi_pred07, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
hemi_pred_df07 <- data.frame(pred.mean = hemi_stat07[, 1],
                             pred.low = hemi_stat07[, 2],
                             pred.high = hemi_stat07[, 3])
hemi_all07 <- cbind(hemi_pred_df07, sd = hemi_stat07[, 4])

min(hemi_all07$pred.mean)
hemi_neg <- length(which(hemi_all07$pred.mean < 0))
hemi_neg #8766 pixels have negative values
saveRDS(hemi_all07, "inference/flg07/hemi_all07.rds") #to load object : hemi_all07 <- readRDS("inference/flg07/hemi_all07.rds")

rm(hemi_pred07, hemi_pred_df07, hemi_stat07)
gc()

cell_pred07 <- apply.coefs(all_jack_coefs_list$cell, val.spec = spec07_bn, intercept = T)
cell_stat07 <- t(apply(cell_pred07, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
cell_pred_df07 <- data.frame(pred.mean = cell_stat07[, 1],
                             pred.low = cell_stat07[, 2],
                             pred.high = cell_stat07[, 3])
cell_all07 <- cbind(cell_pred_df07, sd = cell_stat07[, 4])

min(cell_all07$pred.mean)
cell_neg <- length(which(cell_all07$pred.mean < 0))
cell_neg #59 pixels have negative values
saveRDS(cell_all07, "inference/flg07/cell_all07.rds") #to load object : cell_all07 <- readRDS("inference/flg07/cell_all07.rds")

rm(cell_pred07, cell_pred_df07, cell_stat07)
gc()

lignin_pred07 <- apply.coefs(all_jack_coefs_list$lignin, val.spec = spec07_bn, intercept = T)
lignin_stat07 <- t(apply(lignin_pred07, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
lignin_pred_df07 <- data.frame(pred.mean = lignin_stat07[, 1],
                               pred.low = lignin_stat07[, 2],
                               pred.high = lignin_stat07[, 3])
lignin_all07 <- cbind(lignin_pred_df07, sd = lignin_stat07[, 4])

min(lignin_all07$pred.mean)
lignin_neg <- length(which(lignin_all07$pred.mean < 0))
lignin_neg #24 pixels have negative values
saveRDS(lignin_all07, "inference/flg07/lignin_all07.rds") #to load object : lignin_all07 <- readRDS("inference/flg07/lignin_all07.rds")

rm(lignin_pred07, lignin_pred_df07, lignin_stat07)
gc()

chla_pred07 <- apply.coefs(all_jack_coefs_list$chla, val.spec = spec07_bn, intercept = T)
chla_stat07 <- t(apply(chla_pred07, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chla_pred_df07 <- data.frame(pred.mean = chla_stat07[, 1],
                             pred.low = chla_stat07[, 2],
                             pred.high = chla_stat07[, 3])
chla_all07 <- cbind(chla_pred_df07, sd = chla_stat07[, 4])

min(chla_all07$pred.mean)
chla_neg <- length(which(chla_all07$pred.mean < 0))
chla_neg #2814 pixels have negative values
saveRDS(chla_all07, "inference/flg07/chla_all07.rds") #to load object : chla_all07 <- readRDS("inference/flg07/chla_all07.rds")

rm(chla_pred07, chla_pred_df07, chla_stat07)
gc()

chlb_pred07 <- apply.coefs(all_jack_coefs_list$chlb, val.spec = spec07_bn, intercept = T)
chlb_stat07 <- t(apply(chlb_pred07, 1,
                       function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
chlb_pred_df07 <- data.frame(pred.mean = chlb_stat07[, 1],
                             pred.low = chlb_stat07[, 2],
                             pred.high = chlb_stat07[, 3])
chlb_all07 <- cbind(chlb_pred_df07, sd = chlb_stat07[, 4])

min(chlb_all07$pred.mean)
chlb_neg <- length(which(chlb_all07$pred.mean < 0))
chlb_neg #2189 pixels have negative values
saveRDS(chlb_all07, "inference/flg07/chlb_all07.rds") #to load object : chlb_all07 <- readRDS("inference/flg07/chlb_all07.rds")

rm(chlb_pred07, chlb_pred_df07, chlb_stat07)
gc()

carot_pred07 <- apply.coefs(all_jack_coefs_list$carot, val.spec = spec07_bn, intercept = T)
carot_stat07 <- t(apply(carot_pred07, 1,
                        function(obs) c(mean(obs), quantile(obs, probs=c(0.025, 0.975)), sd(obs))))
carot_pred_df07 <- data.frame(pred.mean = carot_stat07[, 1],
                              pred.low = carot_stat07[, 2],
                              pred.high = carot_stat07[, 3])
carot_all07 <- cbind(carot_pred_df07, sd = carot_stat07[, 4])

min(carot_all07$pred.mean)
carot_neg <- length(which(carot_all07$pred.mean < 0))
carot_neg #1926 pixels have negative values
saveRDS(carot_all07, "inference/flg07/carot_all07.rds") #to load object : carot_all07 <- readRDS("inference/flg07/carot_all07.rds")

rm(carot_pred07, carot_pred_df07, carot_stat07)
gc()

carbon_pred07 <- apply.coefs(all_jack_coefs_list$carbon, val.spec = spec07_bn, intercept = T)
carbon_stat07 <- t(apply(carbon_pred07, 1,
                         function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
carbon_pred_df07 <-data.frame(pred.mean = carbon_stat07[, 1],
                              pred.low = carbon_stat07[, 2],
                              pred.high = carbon_stat07[, 3])
carbon_all07 <- cbind(carbon_pred_df07, sd = carbon_stat07[, 4])

min(carbon_all07$pred.mean)
carbon_neg <- length(which(carbon_all07$pred.mean < 0))
carbon_neg #0 pixels have negative values
saveRDS(carbon_all07, "inference/flg07/carbon_all07.rds") #to load object : carbon_all07 <- readRDS("inference/flg07/carbon_all07.rds")

rm(carbon_pred07, carbon_pred_df07, carbon_stat07)
gc()

nitrogen_pred07 <- apply.coefs(all_jack_coefs_list$nitrogen, val.spec = spec07_bn, intercept = T)
nitrogen_stat07 <- t(apply(nitrogen_pred07, 1,
                           function(obs) c(mean(obs), quantile(obs, probs = c(0.025, 0.975)), sd(obs))))
nitrogen_pred_df07 <- data.frame(pred.mean = nitrogen_stat07[, 1],
                                 pred.low = nitrogen_stat07[, 2],
                                 pred.high = nitrogen_stat07[, 3])
nitrogen_all07 <- cbind(nitrogen_pred_df07, sd = nitrogen_stat07[, 4])

min(nitrogen_all07$pred.mean)
nitrogen_neg <- length(which(nitrogen_all07$pred.mean < 0))
nitrogen_neg #4350 pixels have negative values
saveRDS(nitrogen_all07, "inference/flg07/nitrogen_all07.rds") #to load object : nitrogen_all07 <- readRDS("inference/flg07/nitrogen_all07.rds")

rm(nitrogen_pred07, nitrogen_pred_df07, nitrogen_stat07)
gc()

# Extract coordinates
x_coords07 <- spectra07$x
y_coords07 <- spectra07$y

# Add coordinates to results
SLA_all07$x <- x_coords07
SLA_all07$y <- y_coords07
LMA_all07$x <- x_coords07
LMA_all07$y <- y_coords07
LDMC_all07$x <- x_coords07
LDMC_all07$y <- y_coords07
EWT_all07$x <- x_coords07
EWT_all07$y <- y_coords07
hemi_all07$x <- x_coords07
hemi_all07$y <- y_coords07
cell_all07$x <- x_coords07
cell_all07$y <- y_coords07
lignin_all07$x <- x_coords07
lignin_all07$y <- y_coords07
chla_all07$x <- x_coords07
chla_all07$y <- y_coords07
chlb_all07$x <- x_coords07
chlb_all07$y <- y_coords07
carot_all07$x <- x_coords07
carot_all07$y <- y_coords07
carbon_all07$x <- x_coords07
carbon_all07$y <- y_coords07
nitrogen_all07$x <- x_coords07
nitrogen_all07$y <- y_coords07

# Select columns of interest for mapping
SLA_07 <- SLA_all07 %>%
  dplyr::select(pred.SLA = pred.mean, sd.SLA = sd, x, y) %>%
  `rownames<-`(NULL)
LMA_07 <- LMA_all07 %>%
  dplyr::select(pred.LMA = pred.mean, sd.LMA = sd, x, y) %>%
  `rownames<-`(NULL)
LDMC_07 <- LDMC_all07 %>%
  dplyr::select(pred.LDMC = pred.mean, sd.LDMC = sd, x, y) %>%
  `rownames<-`(NULL)
EWT_07 <- EWT_all07 %>%
  dplyr::select(pred.EWT = pred.mean, sd.EWT = sd, x, y) %>%
  `rownames<-`(NULL)
hemi_07 <- hemi_all07 %>%
  dplyr::select(pred.hemi = pred.mean, sd.hemi = sd, x, y) %>%
  `rownames<-`(NULL)
cell_07 <- cell_all07 %>%
  dplyr::select(pred.cell = pred.mean, sd.cell = sd, x, y) %>%
  `rownames<-`(NULL)
lignin_07 <- lignin_all07 %>%
  dplyr::select(pred.lignin = pred.mean, sd.lignin = sd, x, y) %>%
  `rownames<-`(NULL)
chla_07 <- chla_all07 %>%
  dplyr::select(pred.chla = pred.mean, sd.chla = sd, x, y) %>%
  `rownames<-`(NULL)
chlb_07 <- chlb_all07 %>%
  dplyr::select(pred.chlb = pred.mean, sd.chlb = sd, x, y) %>%
  `rownames<-`(NULL)
carot_07 <- carot_all07 %>%
  dplyr::select(pred.carot = pred.mean, sd.carot = sd, x, y) %>%
  `rownames<-`(NULL)
carbon_07 <- carbon_all07 %>%
  dplyr::select(pred.carbon = pred.mean, sd.carbon = sd, x, y) %>%
  `rownames<-`(NULL)
nitrogen_07 <- nitrogen_all07 %>%
  dplyr::select(pred.nitrogen = pred.mean, sd.nitrogen = sd, x, y) %>%
  `rownames<-`(NULL)

# Join all dataframes
df_list_07 = list(SLA_07, LMA_07, LDMC_07, EWT_07, hemi_07, cell_07, lignin_07, chla_07, chlb_07, carot_07, carbon_07, nitrogen_07) 

df_merged07 <- df_list_07 %>%
  reduce(right_join, by = c("x", "y"))
saveRDS(df_merged07, "inference/flg07/df_merged07.rds") #to load object : df_merged07 <- readRDS("inference/flg07/df_merged07.rds")

# Add relative uncertainties (CV)
df_merged07 <- df_merged07 %>%
  mutate(CV.SLA = sd.SLA/pred.SLA,
         CV.LMA = sd.LMA/pred.LMA,
         CV.LDMC = sd.LDMC/pred.LDMC,
         CV.EWT = sd.EWT/pred.EWT,
         CV.hemi = sd.hemi/pred.hemi,
         CV.cell = sd.cell/pred.cell,
         CV.lignin = sd.lignin/pred.lignin,
         CV.chla = sd.chla/pred.chla,
         CV.chlb = sd.chlb/pred.chlb,
         CV.carot = sd.carot/pred.carot,
         CV.carbon = sd.carbon/pred.carbon,
         CV.nitrogen = sd.nitrogen/pred.nitrogen)

# Create 3 dataframes (traits, sd and CV)
df_traits07 <- df_merged07 %>%
  dplyr::select(starts_with("pred"), x, y) %>%
  rename_with(~sub("pred\\.", "", .), starts_with("pred"))
df_sd07 <- df_merged07 %>%
  dplyr::select(starts_with("sd"), x, y) %>%
  rename_with(~sub("sd\\.", "", .), starts_with("sd"))
df_CV07 <- df_merged07 %>%
  dplyr::select(starts_with("CV"), x, y) %>%
  rename_with(~sub("CV\\.", "", .), starts_with("CV"))
saveRDS(df_traits07, "inference/flg07/df_traits07.rds")
saveRDS(df_sd07, "inference/flg07/df_sd07.rds")
saveRDS(df_CV07, "inference/flg07/df_CV07.rds")
