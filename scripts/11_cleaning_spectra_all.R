### Preprocess spectra

## Load libraries----
library(tidyverse)
library(hsdar)
library(spectrolab)

# Set working directory
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

#### Flightline 01 ####----

# Import data
CASI_spectra01 <- readRDS("output_allpixels/mask_CASI_flg01.rds") %>%
dplyr::select(-c(ID, fraction,Shadow, NDVI, Vegetation))
SASI_spectra01 <- readRDS("output_allpixels/mask_SASI_flg01.rds") %>%
  dplyr::select(-c(ID, fraction))

# Join tables
full_spectra01 <- right_join(CASI_spectra01, SASI_spectra01, by = c("x","y")) %>%
  relocate(x, y, .before = Band_400nm)
full_spectra01 <- na.omit(full_spectra01) #remove pixels that have no spectra in SASI (NDVI filter/Shadow threshold)

rm(CASI_spectra01)
rm(SASI_spectra01)
gc()

saveRDS(full_spectra01, file = "output_allpixels/full_spectra01.rds") #To load object : full_spectra01 <- readRDS("output_allpixels/full_spectra01.rds")

# Cleaning spectra (same process than for sampled trees)
# Convert to spectra object
colnames(full_spectra01) <- gsub("Band_|nm", "", colnames(full_spectra01)) #rename wavelength columns to remove "Band_" and "nm"

spec01 <- as_spectra(full_spectra01, meta_idxs = c(1:2))
summary(spec01) #check spectra object

gc()

# Match the reflectance data
splice_bands = 998 #define the splice band
sensors_matched01 = match_sensors(x = spec01, splice_at = splice_bands, fixed_sensor = 2,
                                interpolate_wvl = c(5, 1))

rm(full_spectra01)
rm(spec01)
gc()

# Convert spectra object to a dataframe (wide format)
spec_df01 <- as.data.frame(sensors_matched01, fix_names = "none", metadata = T)

# Create object of class speclib
bands01 <-  spec_df01 %>%
  dplyr::select(4:223) #extract only the bands

wavelengths01 <- as.numeric(names(spec_df01[, 4:223])) #define wavelengths
speclib01 <- speclib(spectra = t(bands01), wavelength = wavelengths01)

rm(sensors_matched01)
gc()

# Smooth spectra
smooth01 <- noiseFiltering(speclib01, method = 'sgolay', n = 3, p = 2)
hsdar::spectra(smooth01)[hsdar::spectra(smooth01) < 0] <- 0

# Add all important supplementary information
Meta01 <- spec_df01 %>% 
  dplyr::select(1:3) #select metadata
SI(smooth01) <- Meta01 # add metadata to the supplementary information of the spectral library
colnames(SI(smooth01)) <- colnames(Meta01)
SI(smooth01)

rm(spec_df01)
gc()

saveRDS(smooth01, "output_allpixels/smooth01_speclib.rds") #To load object : smooth01 <- readRDS("output_allpixels/smooth01_speclib.rds")

plot(smooth01, FUN = 1:1000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth01, FUN = 1001:2000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth01, FUN = 2001:3000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth01, FUN = 3001:4000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth01, FUN = 4001:5000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth01, FUN = 5001:6000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")

# Convert speclib object to dataframe
SI_smooth01 <- do.call(data.frame, smooth01@SI@SI_data) #extract SI
spectra01 <-smooth01@spectra@spectra_ma %>% data.frame() #extract reflectance data
colnames(spectra01) <- as.integer(smooth01@wavelength) # define colnames of Bands
smoothed01 <- cbind(SI_smooth01, spectra01) 
saveRDS(smoothed01, "output_allpixels/smooth_spectra01.rds") #save data

gc()

#### Flightline 02 ####----

# Import data
CASI_spectra02 <- readRDS("output_allpixels/mask_CASI_flg02.rds") %>%
  dplyr::select(-c(ID, fraction,Shadow, NDVI, Vegetation))
SASI_spectra02 <- readRDS("output_allpixels/mask_SASI_flg02.rds") %>%
  dplyr::select(-c(ID, fraction))

# Join tables
full_spectra02 <- right_join(CASI_spectra02, SASI_spectra02, by = c("x","y")) %>%
  relocate(x, y, .before = Band_400nm)
full_spectra02 <- na.omit(full_spectra02) #remove pixels that have no spectra in SASI (NDVI filter/Shadow threshold)

rm(CASI_spectra02)
rm(SASI_spectra02)
gc()

saveRDS(full_spectra02, file = "output_allpixels/full_spectra02.rds") #To load object : full_spectra02 <- readRDS("output_allpixels/full_spectra02.rds")

# Cleaning spectra (same process than for sampled trees)
# Convert to spectra object
colnames(full_spectra02) <- gsub("Band_|nm", "", colnames(full_spectra02)) #rename wavelength columns to remove "Band_" and "nm"

spec02 <- as_spectra(full_spectra02, meta_idxs = c(1:2))
summary(spec02) #check spectra object

gc()

# Match the reflectance data
splice_bands = 998 #define the splice band
sensors_matched02 = match_sensors(x = spec02, splice_at = splice_bands, fixed_sensor = 2,
                                  interpolate_wvl = c(5, 1))

rm(full_spectra02)
rm(spec02)
gc()

# Convert spectra object to a dataframe (wide format)
spec_df02 <- as.data.frame(sensors_matched02, fix_names = "none", metadata = T)

# Create object of class speclib (load Bands into spectral library)
bands02 <-  spec_df02 %>%
  dplyr::select(4:223) #extract only the bands

wavelengths02 <- as.numeric(names(spec_df02[, 4:223])) #define wavelengths
speclib02 <- speclib(spectra = t(bands02), wavelength = wavelengths02)

rm(sensors_matched02)
gc()

# Smooth spectra
smooth02 <- noiseFiltering(speclib02, method = 'sgolay', n = 3, p = 2)
hsdar::spectra(smooth02)[hsdar::spectra(smooth02) < 0] <- 0

# Add all important supplementary information
Meta02 <- spec_df02 %>% 
  dplyr::select(1:3) #select metadata
SI(smooth02) <- Meta02 # add metadata to the supplementary information of the spectral library
colnames(SI(smooth02)) <- colnames(Meta02)
SI(smooth02)

rm(spec_df02)
gc()

saveRDS(smooth02, "output_allpixels/smooth02_speclib.rds") #To load object : smooth02 <- readRDS("output_allpixels/smooth02_speclib.rds")

plot(smooth02, FUN = 1:1000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth02, FUN = 1001:2000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth02, FUN = 2001:3000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth02, FUN = 3001:4000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth02, FUN = 4001:5000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth02, FUN = 5001:6000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")

# Convert speclib object to dataframe
SI_smooth02 <- do.call(data.frame, smooth02@SI@SI_data) #extract SI
spectra02 <-smooth02@spectra@spectra_ma %>% data.frame() #extract reflectance data
colnames(spectra02) <- as.integer(smooth02@wavelength) # define colnames of Bands
smoothed02 <- cbind(SI_smooth02, spectra02) 
saveRDS(smoothed02, "output_allpixels/smooth_spectra02.rds") #save data

gc()

#### Flightline 03 ####----

# Import data
CASI_spectra03 <- readRDS("output_allpixels/mask_CASI_flg03.rds") %>%
  dplyr::select(-c(ID, fraction,Shadow, NDVI, Vegetation))
SASI_spectra03 <- readRDS("output_allpixels/mask_SASI_flg03.rds") %>%
  dplyr::select(-c(ID, fraction))

# Join tables
full_spectra03 <- right_join(CASI_spectra03, SASI_spectra03, by = c("x","y")) %>%
  relocate(x, y, .before = Band_400nm)
full_spectra03 <- na.omit(full_spectra03) #remove pixels that have no spectra in SASI (NDVI filter/Shadow threshold)

rm(CASI_spectra03)
rm(SASI_spectra03)
gc()

saveRDS(full_spectra03, file = "output_allpixels/full_spectra03.rds") #To load object : full_spectra03 <- readRDS("output_allpixels/full_spectra03.rds")

# Cleaning spectra (same process than for sampled trees)
# Convert to spectra object
colnames(full_spectra03) <- gsub("Band_|nm", "", colnames(full_spectra03)) #rename wavelength columns to remove "Band_" and "nm"

spec03 <- as_spectra(full_spectra03, meta_idxs = c(1:2))
summary(spec03) #check spectra object

gc()

# Match the reflectance data
splice_bands = 998 #define the splice band
sensors_matched03 = match_sensors(x = spec03, splice_at = splice_bands, fixed_sensor = 2,
                                  interpolate_wvl = c(5, 1))

rm(full_spectra03)
rm(spec03)
gc()

# Convert spectra object to a dataframe (wide format)
spec_df03 <- as.data.frame(sensors_matched03, fix_names = "none", metadata = T)

# Create object of class speclib (load Bands into spectral library)
bands03 <-  spec_df03 %>%
  dplyr::select(4:223) #extract only the bands

wavelengths03 <- as.numeric(names(spec_df03[, 4:223])) #define wavelengths
speclib03 <- speclib(spectra = t(bands03), wavelength = wavelengths03)

rm(sensors_matched03)
gc()

# Smooth spectra
smooth03 <- noiseFiltering(speclib03, method = 'sgolay', n = 3, p = 2)
hsdar::spectra(smooth03)[hsdar::spectra(smooth03) < 0] <- 0

# Add all important supplementary information
Meta03 <- spec_df03 %>% 
  dplyr::select(1:3) #select metadata
SI(smooth03) <- Meta03 # add metadata to the supplementary information of the spectral library
colnames(SI(smooth03)) <- colnames(Meta03)
SI(smooth03)

rm(spec_df03)
gc()

saveRDS(smooth03, "output_allpixels/smooth03_speclib.rds") #To load object : smooth03 <- readRDS("output_allpixels/smooth03_speclib.rds")

plot(smooth03, FUN = 1:1000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth03, FUN = 1001:2000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth03, FUN = 2001:3000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth03, FUN = 3001:4000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth03, FUN = 4001:5000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth03, FUN = 5001:6000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")

# Convert speclib object to dataframe
SI_smooth03 <- do.call(data.frame, smooth03@SI@SI_data) #extract SI
spectra03 <-smooth03@spectra@spectra_ma %>% data.frame() #extract reflectance data
colnames(spectra03) <- as.integer(smooth03@wavelength) # define colnames of Bands
smoothed03 <- cbind(SI_smooth03, spectra03) 
saveRDS(smoothed03, "output_allpixels/smooth_spectra03.rds") #save data

gc()

#### Flightline 04 ####----

# Import data
CASI_spectra04 <- readRDS("output_allpixels/mask_CASI_flg04.rds") %>%
  dplyr::select(-c(ID, fraction,Shadow, NDVI, Vegetation))
SASI_spectra04 <- readRDS("output_allpixels/mask_SASI_flg04.rds") %>%
  dplyr::select(-c(ID, fraction))

# Join tables
full_spectra04 <- right_join(CASI_spectra04, SASI_spectra04, by = c("x","y")) %>%
  relocate(x, y, .before = Band_400nm)
full_spectra04 <- na.omit(full_spectra04) #remove pixels that have no spectra in SASI (NDVI filter/Shadow threshold)

rm(CASI_spectra04)
rm(SASI_spectra04)
gc()

saveRDS(full_spectra04, file = "output_allpixels/full_spectra04.rds") #To load object : full_spectra04 <- readRDS("output_allpixels/full_spectra04.rds")

# Cleaning spectra (same process than for sampled trees)
# Convert to spectra object
colnames(full_spectra04) <- gsub("Band_|nm", "", colnames(full_spectra04)) #rename wavelength columns to remove "Band_" and "nm"

spec04 <- as_spectra(full_spectra04, meta_idxs = c(1:2))
summary(spec04) #check spectra object

gc()

# Match the reflectance data
splice_bands = 998 #define the splice band
sensors_matched04 = match_sensors(x = spec04, splice_at = splice_bands, fixed_sensor = 2,
                                  interpolate_wvl = c(5, 1))

rm(full_spectra04)
rm(spec04)
gc()

# Convert spectra object to a dataframe (wide format)
spec_df04 <- as.data.frame(sensors_matched04, fix_names = "none", metadata = T)

# Create object of class speclib (load Bands into spectral library)
bands04 <-  spec_df04 %>%
  dplyr::select(4:223) #extract only the bands

wavelengths04 <- as.numeric(names(spec_df04[, 4:223])) #define wavelengths
speclib04 <- speclib(spectra = t(bands04), wavelength = wavelengths04)

rm(sensors_matched04)
gc()

# Smooth spectra
smooth04 <- noiseFiltering(speclib04, method = 'sgolay', n = 3, p = 2)
hsdar::spectra(smooth04)[hsdar::spectra(smooth04) < 0] <- 0

# Add all important supplementary information
Meta04 <- spec_df04 %>% 
  dplyr::select(1:3) #select metadata
SI(smooth04) <- Meta04 # add metadata to the supplementary information of the spectral library
colnames(SI(smooth04)) <- colnames(Meta04)
SI(smooth04)

rm(spec_df04)
gc()

saveRDS(smooth04, "output_allpixels/smooth04_speclib.rds") #To load object : smooth04 <- readRDS("output_allpixels/smooth04_speclib.rds")

plot(smooth04, FUN = 1:1000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth04, FUN = 1001:2000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth04, FUN = 2001:3000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth04, FUN = 3001:4000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth04, FUN = 4001:5000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth04, FUN = 5001:6000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")

# Convert speclib object to dataframe
SI_smooth04 <- do.call(data.frame, smooth04@SI@SI_data) #extract SI
spectra04 <-smooth04@spectra@spectra_ma %>% data.frame() #extract reflectance data
colnames(spectra04) <- as.integer(smooth04@wavelength) # define colnames of Bands
smoothed04 <- cbind(SI_smooth04, spectra04) 
saveRDS(smoothed04, "output_allpixels/smooth_spectra04.rds") #save data

gc()

#### Flightline 05 ####----

# Import data
CASI_spectra05 <- readRDS("output_allpixels/mask_CASI_flg05.rds") %>%
  dplyr::select(-c(ID, fraction,Shadow, NDVI, Vegetation))
SASI_spectra05 <- readRDS("output_allpixels/mask_SASI_flg05.rds") %>%
  dplyr::select(-c(ID, fraction))

# Join tables
full_spectra05 <- right_join(CASI_spectra05, SASI_spectra05, by = c("x","y")) %>%
  relocate(x, y, .before = Band_400nm)
full_spectra05 <- na.omit(full_spectra05) #remove pixels that have no spectra in SASI (NDVI filter/Shadow threshold)

rm(CASI_spectra05)
rm(SASI_spectra05)
gc()

saveRDS(full_spectra05, file = "output_allpixels/full_spectra05.rds") #To load object : full_spectra05 <- readRDS("output_allpixels/full_spectra05.rds")

# Cleaning spectra (same process than for sampled trees)
# Convert to spectra object
colnames(full_spectra05) <- gsub("Band_|nm", "", colnames(full_spectra05)) #rename wavelength columns to remove "Band_" and "nm"

spec05 <- as_spectra(full_spectra05, meta_idxs = c(1:2))
summary(spec05) #check spectra object

gc()

# Match the reflectance data
splice_bands = 998 #define the splice band
sensors_matched05 = match_sensors(x = spec05, splice_at = splice_bands, fixed_sensor = 2,
                                  interpolate_wvl = c(5, 1))

rm(full_spectra05)
rm(spec05)
gc()

# Convert spectra object to a dataframe (wide format)
spec_df05 <- as.data.frame(sensors_matched05, fix_names = "none", metadata = T)

# Create object of class speclib (load Bands into spectral library)
bands05 <-  spec_df05 %>%
  dplyr::select(4:223) #extract only the bands

wavelengths05 <- as.numeric(names(spec_df05[, 4:223])) #define wavelengths
speclib05 <- speclib(spectra = t(bands05), wavelength = wavelengths05)

rm(sensors_matched05)
gc()

# Smooth spectra
smooth05 <- noiseFiltering(speclib05, method = 'sgolay', n = 3, p = 2)
hsdar::spectra(smooth05)[hsdar::spectra(smooth05) < 0] <- 0

# Add all important supplementary information
Meta05 <- spec_df05 %>% 
  dplyr::select(1:3) #select metadata
SI(smooth05) <- Meta05 # add metadata to the supplementary information of the spectral library
colnames(SI(smooth05)) <- colnames(Meta05)
SI(smooth05)

rm(spec_df05)
gc()

saveRDS(smooth05, "output_allpixels/smooth05_speclib.rds") #To load object : smooth05 <- readRDS("output_allpixels/smooth05_speclib.rds")

plot(smooth05, FUN = 1:1000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth05, FUN = 1001:2000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth05, FUN = 2001:3000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth05, FUN = 3001:4000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth05, FUN = 4001:5000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth05, FUN = 5001:6000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")

# Convert speclib object to dataframe
SI_smooth05 <- do.call(data.frame, smooth05@SI@SI_data) #extract SI
spectra05 <-smooth05@spectra@spectra_ma %>% data.frame() #extract reflectance data
colnames(spectra05) <- as.integer(smooth05@wavelength) # define colnames of Bands
smoothed05 <- cbind(SI_smooth05, spectra05) 
saveRDS(smoothed05, "output_allpixels/smooth_spectra05.rds") #save data

gc()

#### Flightline 06 ####----

# Import data
CASI_spectra06 <- readRDS("output_allpixels/mask_CASI_flg06.rds") %>%
  dplyr::select(-c(ID, fraction,Shadow, NDVI, Vegetation))
SASI_spectra06 <- readRDS("output_allpixels/mask_SASI_flg06.rds") %>%
  dplyr::select(-c(ID, fraction))

# Join tables
full_spectra06 <- right_join(CASI_spectra06, SASI_spectra06, by = c("x","y")) %>%
  relocate(x, y, .before = Band_400nm)
full_spectra06 <- na.omit(full_spectra06) #remove pixels that have no spectra in SASI (NDVI filter/Shadow threshold)

rm(CASI_spectra06)
rm(SASI_spectra06)
gc()

saveRDS(full_spectra06, file = "output_allpixels/full_spectra06.rds") #To load object : full_spectra06 <- readRDS("output_allpixels/full_spectra06.rds")

# Cleaning spectra (same process than for sampled trees)
# Convert to spectra object
colnames(full_spectra06) <- gsub("Band_|nm", "", colnames(full_spectra06)) #rename wavelength columns to remove "Band_" and "nm"

spec06 <- as_spectra(full_spectra06, meta_idxs = c(1:2))
summary(spec06) #check spectra object

gc()

# Match the reflectance data
splice_bands = 998 #define the splice band
sensors_matched06 = match_sensors(x = spec06, splice_at = splice_bands, fixed_sensor = 2,
                                  interpolate_wvl = c(5, 1))

rm(full_spectra06)
rm(spec06)
gc()

# Convert spectra object to a dataframe (wide format)
spec_df06 <- as.data.frame(sensors_matched06, fix_names = "none", metadata = T)

# Create object of class speclib (load Bands into spectral library)
bands06 <-  spec_df06 %>%
  dplyr::select(4:223) #extract only the bands

wavelengths06 <- as.numeric(names(spec_df06[, 4:223])) #define wavelengths
speclib06 <- speclib(spectra = t(bands06), wavelength = wavelengths06)

rm(sensors_matched06)
gc()

# Smooth spectra
smooth06 <- noiseFiltering(speclib06, method = 'sgolay', n = 3, p = 2)
hsdar::spectra(smooth06)[hsdar::spectra(smooth06) < 0] <- 0

# Add all important supplementary information
Meta06 <- spec_df06 %>% 
  dplyr::select(1:3) #select metadata
SI(smooth06) <- Meta06 # add metadata to the supplementary information of the spectral library
colnames(SI(smooth06)) <- colnames(Meta06)
SI(smooth06)

rm(spec_df06)
gc()

saveRDS(smooth06, "output_allpixels/smooth06_speclib.rds") #To load object : smooth06 <- readRDS("output_allpixels/smooth06_speclib.rds")

plot(smooth06, FUN = 1:1000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth06, FUN = 1001:2000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth06, FUN = 2001:3000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth06, FUN = 3001:4000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth06, FUN = 4001:5000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth06, FUN = 5001:6000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")

# Convert speclib object to dataframe
SI_smooth06 <- do.call(data.frame, smooth06@SI@SI_data) #extract SI
spectra06 <-smooth06@spectra@spectra_ma %>% data.frame() #extract reflectance data
colnames(spectra06) <- as.integer(smooth06@wavelength) # define colnames of Bands
smoothed06 <- cbind(SI_smooth06, spectra06) 
saveRDS(smoothed06, "output_allpixels/smooth_spectra06.rds") #save data

gc()

#### Flightline 07 ####----

# Import data
CASI_spectra07 <- readRDS("output_allpixels/mask_CASI_flg07.rds") %>%
  dplyr::select(-c(ID, fraction,Shadow, NDVI, Vegetation))
SASI_spectra07 <- readRDS("output_allpixels/mask_SASI_flg07.rds") %>%
  dplyr::select(-c(ID, fraction))

# Join tables
full_spectra07 <- right_join(CASI_spectra07, SASI_spectra07, by = c("x","y")) %>%
  relocate(x, y, .before = Band_400nm)
full_spectra07 <- na.omit(full_spectra07) #remove pixels that have no spectra in SASI (NDVI filter/Shadow threshold)

rm(CASI_spectra07)
rm(SASI_spectra07)
gc()

saveRDS(full_spectra07, file = "output_allpixels/full_spectra07.rds") #To load object : full_spectra07 <- readRDS("output_allpixels/full_spectra07.rds")

# Cleaning spectra (same process than for sampled trees)
# Convert to spectra object
colnames(full_spectra07) <- gsub("Band_|nm", "", colnames(full_spectra07)) #rename wavelength columns to remove "Band_" and "nm"

spec07 <- as_spectra(full_spectra07, meta_idxs = c(1:2))
summary(spec07) #check spectra object

gc()

# Match the reflectance data
splice_bands = 998 #define the splice band
sensors_matched07 = match_sensors(x = spec07, splice_at = splice_bands, fixed_sensor = 2,
                                  interpolate_wvl = c(5, 1))

rm(full_spectra07)
rm(spec07)
gc()

# Convert spectra object to a dataframe (wide format)
spec_df07 <- as.data.frame(sensors_matched07, fix_names = "none", metadata = T)

# Create object of class speclib (load Bands into spectral library)
bands07 <-  spec_df07 %>%
  dplyr::select(4:223) #extract only the bands

wavelengths07 <- as.numeric(names(spec_df07[, 4:223])) #define wavelengths
speclib07 <- speclib(spectra = t(bands07), wavelength = wavelengths07)

rm(sensors_matched07)
gc()

# Smooth spectra
smooth07 <- noiseFiltering(speclib07, method = 'sgolay', n = 3, p = 2)
hsdar::spectra(smooth07)[hsdar::spectra(smooth07) < 0] <- 0

# Add all important supplementary information
Meta07 <- spec_df07 %>% 
  dplyr::select(1:3) #select metadata
SI(smooth07) <- Meta07 # add metadata to the supplementary information of the spectral library
colnames(SI(smooth07)) <- colnames(Meta07)
SI(smooth07)

rm(spec_df07)
gc()

saveRDS(smooth07, "output_allpixels/smooth07_speclib.rds") #To load object : smooth07 <- readRDS("output_allpixels/smooth07_speclib.rds")

plot(smooth07, FUN = 1:1000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth07, FUN = 1001:2000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth07, FUN = 2001:3000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth07, FUN = 3001:4000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth07, FUN = 4001:5000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")
plot(smooth07, FUN = 5001:6000, main = "Canopy-level spectra (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)")

# Convert speclib object to dataframe
SI_smooth07 <- do.call(data.frame, smooth07@SI@SI_data) #extract SI
spectra07 <-smooth07@spectra@spectra_ma %>% data.frame() #extract reflectance data
colnames(spectra07) <- as.integer(smooth07@wavelength) # define colnames of Bands
smoothed07 <- cbind(SI_smooth07, spectra07) 
saveRDS(smoothed07, "output_allpixels/smooth_spectra07.rds") #save data

gc()