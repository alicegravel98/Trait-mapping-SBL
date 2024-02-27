### Split data into calibration and validation

## Install all packages and load libraries----
install.packages("caret")
library(tidyverse)
library(spectrolab)
library(caret)

## Set working directory----
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Import data----
plsr_data <- read_csv("input_data/plsr_data.csv", col_names = T)

## Convert to spectra object----
wv <- as.numeric(colnames(plsr_data[19:238])) #define wavelengths

spec <- as_spectra(plsr_data, name_idx = 1, meta_idxs = 2:18)

## Create division----
train.sample <- createDataPartition(
  y = meta(spec)$specie,
  p = .70,
  list = FALSE
) #Here we get a warning that a class have a single record

test.sample<-setdiff(1:nrow(as.matrix(spec)),train.sample)

ref.train<-spec[train.sample,]
ref.test<-spec[test.sample,]

# Save as R objects
saveRDS(ref.train,"output_plsr/ref_train.rds")
saveRDS(ref.test,"output_plsr/ref_test.rds")
saveRDS(wv, "output_plsr/wavelength.rds")
