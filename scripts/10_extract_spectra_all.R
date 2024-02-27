### Extract spectra from all pixels (one flightline at a time)

## Install packages and load libraries----
library(tidyverse)
library(hsdar)
library(sf)
library(terra)

# Set working directory
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Define wavelenghts----
CASIWavelengths = c(
  376.473999,  381.256989,  386.039001,  390.821014,  395.604004,  400.386993,
  405.170013,  409.953003,  414.735992,  419.519989,  424.303009,  429.087006,
  433.869995,  438.653992,  443.437988,  448.221985,  453.006012,  457.790009,
  462.575012,  467.359009,  472.144012,  476.928009,  481.713013,  486.497986,
  491.282013,  496.066986,  500.851990,  505.636993,  510.421997,  515.208008,
  519.992981,  524.778015,  529.562988,  534.348999,  539.133972,  543.919006,
  548.705017,  553.489990,  558.276001,  563.060974,  567.846985,  572.632019,
  577.418030,  582.203979,  586.989014,  591.775024,  596.560974,  601.346008,
  606.132019,  610.916992,  615.703003,  620.489014,  625.273987,  630.059998,
  634.844971,  639.630981,  644.416016,  649.202026,  653.987000,  658.773010,
  663.557983,  668.343018,  673.127991,  677.914001,  682.698975,  687.484009,
  692.268982,  697.054016,  701.838989,  706.622986,  711.408020,  716.192993,
  720.976990,  725.762024,  730.546021,  735.330017,  740.114990,  744.898987,
  749.682983,  754.466980,  759.250000,  764.033997,  768.817993,  773.601013,
  778.383972,  783.168030,  787.950989,  792.734009,  797.515991,  802.299011,
  807.081970,  811.864014,  816.645996,  821.427979,  826.210022,  830.992004,
  835.773987,  840.554993,  845.335999,  850.117004,  854.898010,  859.679016,
  864.460022,  869.239990,  874.020020,  878.799988,  883.580017,  888.359985,
  893.138977,  897.919006,  902.697998,  907.476013,  912.255005,  917.033020,
  921.812012,  926.590027,  931.367004,  936.145020,  940.921997,  945.698975,
  950.476013,  955.252014,  960.028992,  964.804993,  969.580994,  974.356018,
  979.130981,  983.906006,  988.681030,  993.455994,  998.229980, 1003.004028,
  1007.776978, 1012.551025, 1017.323975, 1022.096985, 1026.869019, 1031.641968,
  1036.413940, 1041.185059, 1045.957031, 1050.728027, 1055.498047, 1060.269043)

SASIWavelengths = c(
  957.5,  972.5,  987.5, 1002.5, 1017.5, 1032.5, 1047.5, 1062.5, 1077.5, 1092.5, 
  1107.5, 1122.5, 1137.5, 1152.5, 1167.5, 1182.5, 1197.5, 1212.5, 1227.5, 1242.5, 
  1257.5, 1272.5, 1287.5, 1302.5, 1317.5, 1332.5, 1347.5, 1362.5, 1377.5, 1392.5, 
  1407.5, 1422.5, 1437.5, 1452.5, 1467.5, 1482.5, 1497.5, 1512.5, 1527.5, 1542.5, 
  1557.5, 1572.5, 1587.5, 1602.5, 1617.5, 1632.5, 1647.5, 1662.5, 1677.5, 1692.5, 
  1707.5, 1722.5, 1737.5, 1752.5, 1767.5, 1782.5, 1797.5, 1812.5, 1827.5, 1842.5, 
  1857.5, 1872.5, 1887.5, 1902.5, 1917.5, 1932.5, 1947.5, 1962.5, 1977.5, 1992.5, 
  2007.5, 2022.5, 2037.5, 2052.5, 2067.5, 2082.5, 2097.5, 2112.5, 2127.5, 2142.5, 
  2157.5, 2172.5, 2187.5, 2202.5, 2217.5, 2232.5, 2247.5, 2262.5, 2277.5, 2292.5, 
  2307.5, 2322.5, 2337.5, 2352.5, 2367.5, 2382.5, 2397.5, 2412.5, 2427.5, 2442.5)

## Import study area----
# Set coordinate reference systems (CRS)
CRS_SBL <- 32618 # epsg WGS84 UTM18N : 32618

# Import study area
SBL <- st_read("input_data/study_area/study_area.shp") %>% 
  st_transform(crs=CRS_SBL)

# Plot polygons
plot(st_geometry(SBL), border = "red")


#### Flightline 01---- ####

## CASI data----
# Set path
RasterCASIpath01 <- "D:/Hyperspectral_data/CASI/Tiles_all/Flg01"

# Create a list of raster
CASI_raster01 <- list.files(path = RasterCASIpath01, recursive = T,
                          pattern = "\\.TIF$", 
                          full.names = T)
CASI_raster01 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", CASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
CASI_df <- lapply(CASI_raster01[1:11], function(file) { 
  df_CASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_CASI)[3:146] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_CASI)
})

# Rbind lists
CASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

CASI_spectra01 <- CASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(CASI_df)
gc() #Free unused memory (do this often!!)

saveRDS(CASI_spectra01, file = "output_allpixels/raw_CASI_flg01.rds") #To load object : CASI_spectra01 <- readRDS("output_allpixels/raw_CASI_flg01.rds")

# Filter data based on shadows and non-vegetation
CASI_filter01 <- CASI_spectra01 %>%
  mutate(Shadow = case_when(
    Band_802nm > as.numeric(0.1) ~ "NO", 
    Band_802nm <= as.numeric(0.1) ~ "YES")) %>%
  mutate(NDVI = (Band_798nm-Band_640nm)/(Band_798nm+Band_640nm)) %>% #Red(639.630981 - Band 56) and NIR(797.516nm - Band 89)
  mutate(Vegetation = case_when(
    NDVI > 0.6 ~ "YES", 
    NDVI <= 0.6 ~ "NO")) %>%
  dplyr::filter(Shadow=="NO") %>% #remove shadowed areas
  dplyr::filter(Vegetation=="YES") #filter areas with low NDVI

rm(CASI_spectra01)
gc()

saveRDS(CASI_filter01, file = "output_allpixels/CASI_filter01.rds") #To load object : CASI_filter01 <- readRDS("output_allpixels/CASI_filter01.rds")

# Extract only the bands
CASI_bands01 <-  CASI_filter01 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
CASI_speclib01 <- speclib(spectra = t(CASI_bands01), wavelength = CASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(CASI_speclib01) <- paste0(rownames(CASI_filter01), "_",CASI_filter01$"ID") # add IDs

# Extract all important supplementary information
Meta01 <- CASI_filter01 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(CASI_speclib01) <- Meta01 
colnames(SI(CASI_speclib01)) <- colnames(Meta01)
CASI_speclib01 # check summary

# Check dimensions of speclib object
dim(CASI_speclib01) #gives nspectra and nbands (144)

#Plot raw spectra
pdf("output_allpixels/raw_CASI_spectra_01.pdf") #open a pdf device and choose file name

plot(CASI_speclib01, main = "Flightline 01 raw CASI spectra (nspectra = 3 264 579, nbands = 144)", FUN = "max", col = "red", ylim = c(0, 1))
plot(CASI_speclib01, FUN = "min", col = "blue", new = F)
plot(CASI_speclib01, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = c(350,1000), up = c(400,1100))
hsdar::mask(CASI_speclib01) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI01 <- do.call(data.frame, CASI_speclib01@SI@SI_data) #extract SI
CASI_spec01 <- CASI_speclib01@spectra@spectra_ma %>% data.frame() # spectra object as dataframe
colnames(CASI_spec01) <- paste0("Band_", as.integer(CASI_speclib01@wavelength), "nm") #define colnames of bands
CASI_spec01 <- cbind(CASI_spec01, SI01) #bind spectra and SI

saveRDS(CASI_spec01, file = "output_allpixels/mask_CASI_flg01.rds")

rm(CASI_speclib01)
rm(CASI_spec01)
gc()

## SASI data----
# Set path
RasterSASIpath01 <- "D:/Hyperspectral_data/SASI/Tiles_all/Flg01"

# Create a list of raster
SASI_raster01 <- list.files(path = RasterSASIpath01, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
SASI_raster01 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", SASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
SASI_df <- lapply(SASI_raster01[1:12], function(file) { 
  df_SASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_SASI)[3:102] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_SASI)
})

# Rbind lists
SASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

SASI_spectra01 <- SASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(SASI_df)
gc()

saveRDS(SASI_spectra01, file = "output_allpixels/raw_SASI_flg01.rds") #To load object : SASI_spectra01 <- readRDS("output_allpixels/raw_SASI_flg01.rds")

# Filter data based on shadows and non-vegetation
coord_CASI01 <- CASI_filter01$geometry #extract CASI's coordinates
SASI_filter01 <- SASI_spectra01 %>%
  dplyr::filter(geometry %in% coord_CASI01) #keep only SASI points that have same coordinates as SASI

rm(SASI_spectra01)
rm(CASI_filter01)
gc()

saveRDS(SASI_filter01, file = "output_allpixels/SASI_filter01.rds") #To load object : SASI_filter01 <- readRDS("output_allpixels/SASI_filter01.rds")

# Extract only the bands
SASI_bands01 <-  SASI_filter01 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
SASI_speclib01 <- speclib(spectra = t(SASI_bands01), wavelength = SASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(SASI_speclib01) <- paste0(rownames(SASI_filter01), "_",SASI_filter01$"ID") # add IDs

# Extract all important supplementary information
Meta01 <- SASI_filter01 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(SASI_speclib01) <- Meta01 
colnames(SI(SASI_speclib01)) <- colnames(Meta01)
SASI_speclib01 # check summary

# Check dimensions of speclib object
dim(SASI_speclib01) #gives nspectra and nbands (100)

#Plot raw spectra
pdf("output_allpixels/raw_SASI_spectra_01.pdf") #open a pdf device and choose file name

plot(SASI_speclib01, main = "Flightline 01 raw SASI spectra (nspectra = 2 978 991 , nbands = 100)", FUN = "max", col = "red", ylim = c(0, 1))
plot(SASI_speclib01, FUN = "min", col = "blue", new = F)
plot(SASI_speclib01, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = 2400, up = 2500)
hsdar::mask(SASI_speclib01) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI01 <- do.call(data.frame, SASI_speclib01@SI@SI_data) #extract SI
SASI_spec01 <- SASI_speclib01@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(SASI_spec01) <- paste0("Band_", as.integer(SASI_speclib01@wavelength), "nm") #define colnames of bands
SASI_spec01 <- cbind(SASI_spec01, SI01) #bind spectra and SI

saveRDS(SASI_spec01, file = "output_allpixels/mask_SASI_flg01.rds")

rm(SASI_speclib01)
rm(SASI_spec01)
rm(SASI_filter01)
gc()

#### Flightline 02---- ####

## CASI data----
# Set path
RasterCASIpath02 <- "D:/Hyperspectral_data/CASI/Tiles_all/Flg02"

# Create a list of raster
CASI_raster02 <- list.files(path = RasterCASIpath02, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
CASI_raster02 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", CASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
CASI_df <- lapply(CASI_raster02[1:12], function(file) { 
  df_CASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_CASI)[3:146] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_CASI)
})

# Rbind lists
CASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

CASI_spectra02 <- CASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(CASI_df)
gc() #Free unused memory (do this often!!)

saveRDS(CASI_spectra02, file = "output_allpixels/raw_CASI_flg02.rds") #To load object : CASI_spectra02 <- readRDS("output_allpixels/raw_CASI_flg02.rds")

# Filter data based on shadows and non-vegetation
CASI_filter02 <- CASI_spectra02 %>%
  mutate(Shadow = case_when(
    Band_802nm > as.numeric(0.1) ~ "NO", 
    Band_802nm <= as.numeric(0.1) ~ "YES")) %>%
  mutate(NDVI = (Band_798nm-Band_640nm)/(Band_798nm+Band_640nm)) %>% #Red(639.630981 - Band 56) and NIR(797.516nm - Band 89)
  mutate(Vegetation = case_when(
    NDVI > 0.6 ~ "YES", 
    NDVI <= 0.6 ~ "NO")) %>%
  dplyr::filter(Shadow=="NO") %>% #remove shadowed areas
  dplyr::filter(Vegetation=="YES") #filter areas with low NDVI

rm(CASI_spectra02)
gc()

saveRDS(CASI_filter02, file = "output_allpixels/CASI_filter02.rds") #To load object : CASI_filter02 <- readRDS("output_allpixels/CASI_filter02.rds")

# Extract only the bands
CASI_bands02 <-  CASI_filter02 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
CASI_speclib02 <- speclib(spectra = t(CASI_bands02), wavelength = CASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(CASI_speclib02) <- paste0(rownames(CASI_filter02), "_",CASI_filter02$"ID") # add IDs

# Extract all important supplementary information
Meta02 <- CASI_filter02 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(CASI_speclib02) <- Meta02 
colnames(SI(CASI_speclib02)) <- colnames(Meta02)
CASI_speclib02 # check summary

# Check dimensions of speclib object
dim(CASI_speclib02) #gives nspectra and nbands (144)

#Plot raw spectra
pdf("output_allpixels/raw_CASI_spectra_02.pdf") #open a pdf device and choose file name

plot(CASI_speclib02, main = "Flightline 02 raw CASI spectra (nspectra = 2 888 071, nbands = 144)", FUN = "max", col = "red", ylim = c(0, 1))
plot(CASI_speclib02, FUN = "min", col = "blue", new = F)
plot(CASI_speclib02, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = c(350,1000), up = c(400,1100))
hsdar::mask(CASI_speclib02) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI02 <- do.call(data.frame, CASI_speclib02@SI@SI_data) #extract SI
CASI_spec02 <- CASI_speclib02@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(CASI_spec02) <- paste0("Band_", as.integer(CASI_speclib02@wavelength), "nm") #define colnames of bands
CASI_spec02 <- cbind(CASI_spec02, SI02) #bind spectra and SI

saveRDS(CASI_spec02, file = "output_allpixels/mask_CASI_flg02.rds")

rm(CASI_speclib02)
rm(CASI_spec02)
gc()

## SASI data----
# Set path
RasterSASIpath02 <- "D:/Hyperspectral_data/SASI/Tiles_all/Flg02"

# Create a list of raster
SASI_raster02 <- list.files(path = RasterSASIpath02, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
SASI_raster02 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", SASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
SASI_df <- lapply(SASI_raster02[1:12], function(file) { 
  df_SASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_SASI)[3:102] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_SASI)
})

# Rbind lists
SASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

SASI_spectra02 <- SASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(SASI_df)
gc()

saveRDS(SASI_spectra02, file = "output_allpixels/raw_SASI_flg02.rds") #To load object : SASI_spectra02 <- readRDS("output_allpixels/raw_SASI_flg02.rds")

# Filter data based on shadows and non-vegetation
coord_CASI02 <- CASI_filter02$geometry #extract CASI's coordinates
SASI_filter02 <- SASI_spectra02 %>%
  dplyr::filter(geometry %in% coord_CASI02) #keep only SASI points that have same coordinates as SASI

rm(SASI_spectra02)
rm(CASI_filter02)
gc()

saveRDS(SASI_filter02, file = "output_allpixels/SASI_filter02.rds") #To load object : SASI_filter02 <- readRDS("output_allpixels/SASI_filter02.rds")

# Extract only the bands
SASI_bands02 <-  SASI_filter02 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
SASI_speclib02 <- speclib(spectra = t(SASI_bands02), wavelength = SASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(SASI_speclib02) <- paste0(rownames(SASI_filter02), "_",SASI_filter02$"ID") # add IDs

# Extract all important supplementary information
Meta02 <- SASI_filter02 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(SASI_speclib02) <- Meta02 
colnames(SI(SASI_speclib02)) <- colnames(Meta02)
SASI_speclib02 # check summary

# Check dimensions of speclib object
dim(SASI_speclib02) #gives nspectra and nbands (100)

#Plot raw spectra
pdf("output_allpixels/raw_SASI_spectra_02.pdf") #open a pdf device and choose file name

plot(SASI_speclib02, main = "Flightline 02 raw SASI spectra (nspectra = 2 763 413, nbands = 100)", FUN = "max", col = "red", ylim = c(0, 1))
plot(SASI_speclib02, FUN = "min", col = "blue", new = F)
plot(SASI_speclib02, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = 2400, up = 2500)
hsdar::mask(SASI_speclib02) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI02 <- do.call(data.frame, SASI_speclib02@SI@SI_data) #extract SI
SASI_spec02 <- SASI_speclib02@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(SASI_spec02) <- paste0("Band_", as.integer(SASI_speclib02@wavelength), "nm") #define colnames of bands
SASI_spec02 <- cbind(SASI_spec02, SI02) #bind spectra and SI

saveRDS(SASI_spec02, file = "output_allpixels/mask_SASI_flg02.rds")

rm(SASI_speclib02)
rm(SASI_spec02)
rm(SASI_filter02)
gc()

#### Flightline 03---- ####

## CASI data----
# Set path
RasterCASIpath03 <- "D:/Hyperspectral_data/CASI/Tiles_all/Flg03"

# Create a list of raster
CASI_raster03 <- list.files(path = RasterCASIpath03, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
CASI_raster03 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", CASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
CASI_df <- lapply(CASI_raster03[1:12], function(file) { 
  df_CASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_CASI)[3:146] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_CASI)
})

# Rbind lists
CASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

CASI_spectra03 <- CASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(CASI_df)
gc() #Free unused memory (do this often!!)

saveRDS(CASI_spectra03, file = "output_allpixels/raw_CASI_flg03.rds") #To load object : CASI_spectra03 <- readRDS("output_allpixels/raw_CASI_flg03.rds")

# Filter data based on shadows and non-vegetation
CASI_filter03 <- CASI_spectra03 %>%
  mutate(Shadow = case_when(
    Band_802nm > as.numeric(0.1) ~ "NO", 
    Band_802nm <= as.numeric(0.1) ~ "YES")) %>%
  mutate(NDVI = (Band_798nm-Band_640nm)/(Band_798nm+Band_640nm)) %>% #Red(639.630981 - Band 56) and NIR(797.516nm - Band 89)
  mutate(Vegetation = case_when(
    NDVI > 0.6 ~ "YES", 
    NDVI <= 0.6 ~ "NO")) %>%
  dplyr::filter(Shadow=="NO") %>% #remove shadowed areas
  dplyr::filter(Vegetation=="YES") #filter areas with low NDVI

rm(CASI_spectra03)
gc()

saveRDS(CASI_filter03, file = "output_allpixels/CASI_filter03.rds") #To load object : CASI_filter03 <- readRDS("output_allpixels/CASI_filter03.rds")

# Extract only the bands
CASI_bands03 <-  CASI_filter03 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
CASI_speclib03 <- speclib(spectra = t(CASI_bands03), wavelength = CASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(CASI_speclib03) <- paste0(rownames(CASI_filter03), "_",CASI_filter03$"ID") # add IDs

# Extract all important supplementary information
Meta03 <- CASI_filter03 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(CASI_speclib03) <- Meta03 
colnames(SI(CASI_speclib03)) <- colnames(Meta03)
CASI_speclib03 # check summary

# Check dimensions of speclib object
dim(CASI_speclib03) #gives nspectra and nbands (144)

#Plot raw spectra
pdf("output_allpixels/raw_CASI_spectra_03.pdf") #open a pdf device and choose file name

plot(CASI_speclib03, main = "Flightline 03 raw CASI spectra (nspectra = 2 757 496 , nbands = 144)", FUN = "max", col = "red", ylim = c(0, 1))
plot(CASI_speclib03, FUN = "min", col = "blue", new = F)
plot(CASI_speclib03, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = c(350,1000), up = c(400,1100))
hsdar::mask(CASI_speclib03) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI03 <- do.call(data.frame, CASI_speclib03@SI@SI_data) #extract SI
CASI_spec03 <- CASI_speclib03@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(CASI_spec03) <- paste0("Band_", as.integer(CASI_speclib03@wavelength), "nm") #define colnames of bands
CASI_spec03 <- cbind(CASI_spec03, SI03) #bind spectra and SI

saveRDS(CASI_spec03, file = "output_allpixels/mask_CASI_flg03.rds")

rm(CASI_speclib03)
rm(CASI_spec03)
gc()

## SASI data----
# Set path
RasterSASIpath03 <- "D:/Hyperspectral_data/SASI/Tiles_all/Flg03"

# Create a list of raster
SASI_raster03 <- list.files(path = RasterSASIpath03, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
SASI_raster03 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", SASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
SASI_df <- lapply(SASI_raster03[1:11], function(file) { 
  df_SASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_SASI)[3:102] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_SASI)
})

# Rbind lists
SASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

SASI_spectra03 <- SASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(SASI_df)
gc()

saveRDS(SASI_spectra03, file = "output_allpixels/raw_SASI_flg03.rds") #To load object : SASI_spectra03 <- readRDS("output_allpixels/raw_SASI_flg03.rds")

# Filter data based on shadows and non-vegetation
coord_CASI03 <- CASI_filter03$geometry #extract CASI's coordinates
SASI_filter03 <- SASI_spectra03 %>%
  dplyr::filter(geometry %in% coord_CASI03) #keep only SASI points that have same coordinates as SASI

rm(SASI_spectra03)
rm(CASI_filter03)
gc()

saveRDS(SASI_filter03, file = "output_allpixels/SASI_filter03.rds") #To load object : SASI_filter03 <- readRDS("output_allpixels/SASI_filter03.rds")

# Extract only the bands
SASI_bands03 <-  SASI_filter03 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
SASI_speclib03 <- speclib(spectra = t(SASI_bands03), wavelength = SASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(SASI_speclib03) <- paste0(rownames(SASI_filter03), "_",SASI_filter03$"ID") # add IDs

# Extract all important supplementary information
Meta03 <- SASI_filter03 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(SASI_speclib03) <- Meta03 
colnames(SI(SASI_speclib03)) <- colnames(Meta03)
SASI_speclib03 # check summary

# Check dimensions of speclib object
dim(SASI_speclib03) #gives nspectra and nbands (100)

#Plot raw spectra
pdf("output_allpixels/raw_SASI_spectra_03.pdf") #open a pdf device and choose file name

plot(SASI_speclib03, main = "Flightline 03 raw SASI spectra (nspectra = 2 649 880, nbands = 100)", FUN = "max", col = "red", ylim = c(0, 1))
plot(SASI_speclib03, FUN = "min", col = "blue", new = F)
plot(SASI_speclib03, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = 2400, up = 2500)
hsdar::mask(SASI_speclib03) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI03 <- do.call(data.frame, SASI_speclib03@SI@SI_data) #extract SI
SASI_spec03 <- SASI_speclib03@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(SASI_spec03) <- paste0("Band_", as.integer(SASI_speclib03@wavelength), "nm") #define colnames of bands
SASI_spec03 <- cbind(SASI_spec03, SI03) #bind spectra and SI

saveRDS(SASI_spec03, file = "output_allpixels/mask_SASI_flg03.rds")

rm(SASI_speclib03)
rm(SASI_spec03)
rm(SASI_filter03)
gc()

#### Flightline 04---- ####

## CASI data----
# Set path
RasterCASIpath04 <- "D:/Hyperspectral_data/CASI/Tiles_all/Flg04"

# Create a list of raster
CASI_raster04 <- list.files(path = RasterCASIpath04, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
CASI_raster04 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", CASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
CASI_df <- lapply(CASI_raster04[1:10], function(file) { 
  df_CASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_CASI)[3:146] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_CASI)
})

# Rbind lists
CASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

CASI_spectra04 <- CASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(CASI_df)
gc() #Free unused memory (do this often!!)

saveRDS(CASI_spectra04, file = "output_allpixels/raw_CASI_flg04.rds") #To load object : CASI_spectra04 <- readRDS("output_allpixels/raw_CASI_flg04.rds")

# Filter data based on shadows and non-vegetation
CASI_filter04 <- CASI_spectra04 %>%
  mutate(Shadow = case_when(
    Band_802nm > as.numeric(0.1) ~ "NO", 
    Band_802nm <= as.numeric(0.1) ~ "YES")) %>%
  mutate(NDVI = (Band_798nm-Band_640nm)/(Band_798nm+Band_640nm)) %>% #Red(639.630981 - Band 56) and NIR(797.516nm - Band 89)
  mutate(Vegetation = case_when(
    NDVI > 0.6 ~ "YES", 
    NDVI <= 0.6 ~ "NO")) %>%
  dplyr::filter(Shadow=="NO") %>% #remove shadowed areas
  dplyr::filter(Vegetation=="YES") #filter areas with low NDVI

rm(CASI_spectra04)
gc()

saveRDS(CASI_filter04, file = "output_allpixels/CASI_filter04.rds") #To load object : CASI_filter04 <- readRDS("output_allpixels/CASI_filter04.rds")

# Extract only the bands
CASI_bands04 <-  CASI_filter04 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
CASI_speclib04 <- speclib(spectra = t(CASI_bands04), wavelength = CASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(CASI_speclib04) <- paste0(rownames(CASI_filter04), "_",CASI_filter04$"ID") # add IDs

# Extract all important supplementary information
Meta04 <- CASI_filter04 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(CASI_speclib04) <- Meta04 
colnames(SI(CASI_speclib04)) <- colnames(Meta04)
CASI_speclib04 # check summary

# Check dimensions of speclib object
dim(CASI_speclib04) #gives nspectra and nbands (144)

#Plot raw spectra
pdf("output_allpixels/raw_CASI_spectra_04.pdf") #open a pdf device and choose file name

plot(CASI_speclib04, main = "Flightline 04 raw CASI spectra (nspectra = 3 209 736 , nbands = 144)", FUN = "max", col = "red", ylim = c(0, 1))
plot(CASI_speclib04, FUN = "min", col = "blue", new = F)
plot(CASI_speclib04, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = c(350,1000), up = c(400,1100))
hsdar::mask(CASI_speclib04) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI04 <- do.call(data.frame, CASI_speclib04@SI@SI_data) #extract SI
CASI_spec04 <- CASI_speclib04@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(CASI_spec04) <- paste0("Band_", as.integer(CASI_speclib04@wavelength), "nm") #define colnames of bands
CASI_spec04 <- cbind(CASI_spec04, SI04) #bind spectra and SI

saveRDS(CASI_spec04, file = "output_allpixels/mask_CASI_flg04.rds")

rm(CASI_speclib04)
rm(CASI_spec04)
gc()

## SASI data----
# Set path
RasterSASIpath04 <- "D:/Hyperspectral_data/SASI/Tiles_all/Flg04"

# Create a list of raster
SASI_raster04 <- list.files(path = RasterSASIpath04, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
SASI_raster04 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", SASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
SASI_df <- lapply(SASI_raster04[1:13], function(file) { 
  df_SASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_SASI)[3:102] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_SASI)
})

# Rbind lists
SASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

SASI_spectra04 <- SASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(SASI_df)
gc()

saveRDS(SASI_spectra04, file = "output_allpixels/raw_SASI_flg04.rds") #To load object : SASI_spectra04 <- readRDS("output_allpixels/raw_SASI_flg04.rds")

# Filter data based on shadows and non-vegetation
coord_CASI04 <- CASI_filter04$geometry #extract CASI's coordinates
SASI_filter04 <- SASI_spectra04 %>%
  dplyr::filter(geometry %in% coord_CASI04) #keep only SASI points that have same coordinates as SASI

rm(SASI_spectra04)
rm(CASI_filter04)
gc()

saveRDS(SASI_filter04, file = "output_allpixels/SASI_filter04.rds") #To load object : SASI_filter04 <- readRDS("output_allpixels/SASI_filter04.rds")

# Extract only the bands
SASI_bands04 <-  SASI_filter04 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
SASI_speclib04 <- speclib(spectra = t(SASI_bands04), wavelength = SASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(SASI_speclib04) <- paste0(rownames(SASI_filter04), "_",SASI_filter04$"ID") # add IDs

# Extract all important supplementary information
Meta04 <- SASI_filter04 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(SASI_speclib04) <- Meta04 
colnames(SI(SASI_speclib04)) <- colnames(Meta04)
SASI_speclib04 # check summary

# Check dimensions of speclib object
dim(SASI_speclib04) #gives nspectra and nbands (100)

#Plot raw spectra
pdf("output_allpixels/raw_SASI_spectra_04.pdf") #open a pdf device and choose file name

plot(SASI_speclib04, main = "Flightline 04 raw SASI spectra (nspectra = 3 180 371, nbands = 100)", FUN = "max", col = "red", ylim = c(0, 1))
plot(SASI_speclib04, FUN = "min", col = "blue", new = F)
plot(SASI_speclib04, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = 2400, up = 2500)
hsdar::mask(SASI_speclib04) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI04 <- do.call(data.frame, SASI_speclib04@SI@SI_data) #extract SI
SASI_spec04 <- SASI_speclib04@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(SASI_spec04) <- paste0("Band_", as.integer(SASI_speclib04@wavelength), "nm") #define colnames of bands
SASI_spec04 <- cbind(SASI_spec04, SI04) #bind spectra and SI

saveRDS(SASI_spec04, file = "output_allpixels/mask_SASI_flg04.rds")

rm(SASI_speclib04)
rm(SASI_spec04)
rm(SASI_filter04)
gc()

#### Flightline 05---- ####

## CASI data----
# Set path
RasterCASIpath05 <- "D:/Hyperspectral_data/CASI/Tiles_all/Flg05"

# Create a list of raster
CASI_raster05 <- list.files(path = RasterCASIpath05, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
CASI_raster05 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", CASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
CASI_df <- lapply(CASI_raster05[1:12], function(file) { 
  df_CASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_CASI)[3:146] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_CASI)
})

# Rbind lists
CASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

CASI_spectra05 <- CASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(CASI_df)
gc() #Free unused memory (do this often!!)

saveRDS(CASI_spectra05, file = "output_allpixels/raw_CASI_flg05.rds") #To load object : CASI_spectra05 <- readRDS("output_allpixels/raw_CASI_flg05.rds")

# Filter data based on shadows and non-vegetation
CASI_filter05 <- CASI_spectra05 %>%
  mutate(Shadow = case_when(
    Band_802nm > as.numeric(0.1) ~ "NO", 
    Band_802nm <= as.numeric(0.1) ~ "YES")) %>%
  mutate(NDVI = (Band_798nm-Band_640nm)/(Band_798nm+Band_640nm)) %>% #Red(639.630981 - Band 56) and NIR(797.516nm - Band 89)
  mutate(Vegetation = case_when(
    NDVI > 0.6 ~ "YES", 
    NDVI <= 0.6 ~ "NO")) %>%
  dplyr::filter(Shadow=="NO") %>% #remove shadowed areas
  dplyr::filter(Vegetation=="YES") #filter areas with low NDVI

rm(CASI_spectra05)
gc()

saveRDS(CASI_filter05, file = "output_allpixels/CASI_filter05.rds") #To load object : CASI_filter05 <- readRDS("output_allpixels/CASI_filter05.rds")

# Extract only the bands
CASI_bands05 <-  CASI_filter05 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
CASI_speclib05 <- speclib(spectra = t(CASI_bands05), wavelength = CASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(CASI_speclib05) <- paste0(rownames(CASI_filter05), "_",CASI_filter05$"ID") # add IDs

# Extract all important supplementary information
Meta05 <- CASI_filter05 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(CASI_speclib05) <- Meta05 
colnames(SI(CASI_speclib05)) <- colnames(Meta05)
CASI_speclib05 # check summary

# Check dimensions of speclib object
dim(CASI_speclib05) #gives nspectra and nbands (144)

#Plot raw spectra
pdf("output_allpixels/raw_CASI_spectra_05.pdf") #open a pdf device and choose file name

plot(CASI_speclib05, main = "Flightline 05 raw CASI spectra (nspectra = 2 933 280, nbands = 144)", FUN = "max", col = "red", ylim = c(0, 1))
plot(CASI_speclib05, FUN = "min", col = "blue", new = F)
plot(CASI_speclib05, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = c(350,1000), up = c(400,1100))
hsdar::mask(CASI_speclib05) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI05 <- do.call(data.frame, CASI_speclib05@SI@SI_data) #extract SI
CASI_spec05 <- CASI_speclib05@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(CASI_spec05) <- paste0("Band_", as.integer(CASI_speclib05@wavelength), "nm") #define colnames of bands
CASI_spec05 <- cbind(CASI_spec05, SI05) #bind spectra and SI

saveRDS(CASI_spec05, file = "output_allpixels/mask_CASI_flg05.rds")

rm(CASI_speclib05)
rm(CASI_spec05)
gc()

## SASI data----
# Set path
RasterSASIpath05 <- "D:/Hyperspectral_data/SASI/Tiles_all/Flg05"

# Create a list of raster
SASI_raster05 <- list.files(path = RasterSASIpath05, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
SASI_raster05 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", SASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
SASI_df <- lapply(SASI_raster05[1:12], function(file) { 
  df_SASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_SASI)[3:102] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_SASI)
})

# Rbind lists
SASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

SASI_spectra05 <- SASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(SASI_df)
gc()

saveRDS(SASI_spectra05, file = "output_allpixels/raw_SASI_flg05.rds") #To load object : SASI_spectra05 <- readRDS("output_allpixels/raw_SASI_flg05.rds")

# Filter data based on shadows and non-vegetation
coord_CASI05 <- CASI_filter05$geometry #extract CASI's coordinates
SASI_filter05 <- SASI_spectra05 %>%
  dplyr::filter(geometry %in% coord_CASI05) #keep only SASI points that have same coordinates as SASI

rm(SASI_spectra05)
rm(CASI_filter05)
gc()

saveRDS(SASI_filter05, file = "output_allpixels/SASI_filter05.rds") #To load object : SASI_filter05 <- readRDS("output_allpixels/SASI_filter05.rds")

# Extract only the bands
SASI_bands05 <-  SASI_filter05 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
SASI_speclib05 <- speclib(spectra = t(SASI_bands05), wavelength = SASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(SASI_speclib05) <- paste0(rownames(SASI_filter05), "_",SASI_filter05$"ID") # add IDs

# Extract all important supplementary information
Meta05 <- SASI_filter05 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(SASI_speclib05) <- Meta05 
colnames(SI(SASI_speclib05)) <- colnames(Meta05)
SASI_speclib05 # check summary

# Check dimensions of speclib object
dim(SASI_speclib05) #gives nspectra and nbands (100)

#Plot raw spectra
pdf("output_allpixels/raw_SASI_spectra_05.pdf") #open a pdf device and choose file name

plot(SASI_speclib05, main = "Flightline 05 raw SASI spectra (nspectra = 2 813 225, nbands = 100)", FUN = "max", col = "red", ylim = c(0, 1))
plot(SASI_speclib05, FUN = "min", col = "blue", new = F)
plot(SASI_speclib05, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = 2400, up = 2500)
hsdar::mask(SASI_speclib05) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI05 <- do.call(data.frame, SASI_speclib05@SI@SI_data) #extract SI
SASI_spec05 <- SASI_speclib05@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(SASI_spec05) <- paste0("Band_", as.integer(SASI_speclib05@wavelength), "nm") #define colnames of bands
SASI_spec05 <- cbind(SASI_spec05, SI05) #bind spectra and SI

saveRDS(SASI_spec05, file = "output_allpixels/mask_SASI_flg05.rds")

rm(SASI_speclib05)
rm(SASI_spec05)
rm(SASI_filter05)
gc()

#### Flightline 06---- ####

## CASI data----
# Set path
RasterCASIpath06 <- "D:/Hyperspectral_data/CASI/Tiles_all/Flg06"

# Create a list of raster
CASI_raster06 <- list.files(path = RasterCASIpath06, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
CASI_raster06 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", CASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
CASI_df <- lapply(CASI_raster06[1:8], function(file) { 
  df_CASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_CASI)[3:146] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_CASI)
})

# Rbind lists
CASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

CASI_spectra06 <- CASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(CASI_df)
gc() #Free unused memory (do this often!!)

saveRDS(CASI_spectra06, file = "output_allpixels/raw_CASI_flg06.rds") #To load object : CASI_spectra06 <- readRDS("output_allpixels/raw_CASI_flg06.rds")

# Filter data based on shadows and non-vegetation
CASI_filter06 <- CASI_spectra06 %>%
  mutate(Shadow = case_when(
    Band_802nm > as.numeric(0.1) ~ "NO", 
    Band_802nm <= as.numeric(0.1) ~ "YES")) %>%
  mutate(NDVI = (Band_798nm-Band_640nm)/(Band_798nm+Band_640nm)) %>% #Red(639.630981 - Band 56) and NIR(797.516nm - Band 89)
  mutate(Vegetation = case_when(
    NDVI > 0.6 ~ "YES", 
    NDVI <= 0.6 ~ "NO")) %>%
  dplyr::filter(Shadow=="NO") %>% #remove shadowed areas
  dplyr::filter(Vegetation=="YES") #filter areas with low NDVI

rm(CASI_spectra06)
gc()

saveRDS(CASI_filter06, file = "output_allpixels/CASI_filter06.rds") #To load object : CASI_filter06 <- readRDS("output_allpixels/CASI_filter06.rds")

# Extract only the bands
CASI_bands06 <-  CASI_filter06 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
CASI_speclib06 <- speclib(spectra = t(CASI_bands06), wavelength = CASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(CASI_speclib06) <- paste0(rownames(CASI_filter06), "_",CASI_filter06$"ID") # add IDs

# Extract all important supplementary information
Meta06 <- CASI_filter06 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(CASI_speclib06) <- Meta06 
colnames(SI(CASI_speclib06)) <- colnames(Meta06)
CASI_speclib06 # check summary

# Check dimensions of speclib object
dim(CASI_speclib06) #gives nspectra and nbands (144)

#Plot raw spectra
pdf("output_allpixels/raw_CASI_spectra_06.pdf") #open a pdf device and choose file name

plot(CASI_speclib06, main = "Flightline 06 raw CASI spectra (nspectra = 2 406 033, nbands = 144)", FUN = "max", col = "red", ylim = c(0, 1))
plot(CASI_speclib06, FUN = "min", col = "blue", new = F)
plot(CASI_speclib06, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = c(350,1000), up = c(400,1100))
hsdar::mask(CASI_speclib06) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI06 <- do.call(data.frame, CASI_speclib06@SI@SI_data) #extract SI
CASI_spec06 <- CASI_speclib06@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(CASI_spec06) <- paste0("Band_", as.integer(CASI_speclib06@wavelength), "nm") #define colnames of bands
CASI_spec06 <- cbind(CASI_spec06, SI06) #bind spectra and SI

saveRDS(CASI_spec06, file = "output_allpixels/mask_CASI_flg06.rds")

rm(CASI_speclib06)
rm(CASI_spec06)
gc()

## SASI data----
# Set path
RasterSASIpath06 <- "D:/Hyperspectral_data/SASI/Tiles_all/Flg06"

# Create a list of raster
SASI_raster06 <- list.files(path = RasterSASIpath06, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
SASI_raster06 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", SASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
SASI_df <- lapply(SASI_raster06[1:9], function(file) { 
  df_SASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_SASI)[3:102] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_SASI)
})

# Rbind lists
SASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

SASI_spectra06 <- SASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(SASI_df)
gc()

saveRDS(SASI_spectra06, file = "output_allpixels/raw_SASI_flg06.rds") #To load object : SASI_spectra06 <- readRDS("output_allpixels/raw_SASI_flg06.rds")

# Filter data based on shadows and non-vegetation
coord_CASI06 <- CASI_filter06$geometry #extract CASI's coordinates
SASI_filter06 <- SASI_spectra06 %>%
  dplyr::filter(geometry %in% coord_CASI06) #keep only SASI points that have same coordinates as SASI

rm(SASI_spectra06)
rm(CASI_filter06)
gc()

saveRDS(SASI_filter06, file = "output_allpixels/SASI_filter06.rds") #To load object : SASI_filter06 <- readRDS("output_allpixels/SASI_filter06.rds")

# Extract only the bands
SASI_bands06 <-  SASI_filter06 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
SASI_speclib06 <- speclib(spectra = t(SASI_bands06), wavelength = SASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(SASI_speclib06) <- paste0(rownames(SASI_filter06), "_",SASI_filter06$"ID") # add IDs

# Extract all important supplementary information
Meta06 <- SASI_filter06 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(SASI_speclib06) <- Meta06 
colnames(SI(SASI_speclib06)) <- colnames(Meta06)
SASI_speclib06 # check summary

# Check dimensions of speclib object
dim(SASI_speclib06) #gives nspectra and nbands (100)

#Plot raw spectra
pdf("output_allpixels/raw_SASI_spectra_06.pdf") #open a pdf device and choose file name

plot(SASI_speclib06, main = "Flightline 06 raw SASI spectra (nspectra = 2 204 716, nbands = 100)", FUN = "max", col = "red", ylim = c(0, 1))
plot(SASI_speclib06, FUN = "min", col = "blue", new = F)
plot(SASI_speclib06, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = 2400, up = 2500)
hsdar::mask(SASI_speclib06) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI06 <- do.call(data.frame, SASI_speclib06@SI@SI_data) #extract SI
SASI_spec06 <- SASI_speclib06@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(SASI_spec06) <- paste0("Band_", as.integer(SASI_speclib06@wavelength), "nm") #define colnames of bands
SASI_spec06 <- cbind(SASI_spec06, SI06) #bind spectra and SI

saveRDS(SASI_spec06, file = "output_allpixels/mask_SASI_flg06.rds")

rm(SASI_speclib06)
rm(SASI_spec06)
rm(SASI_filter06)
gc()

#### Flightline 07---- ####

## CASI data----
# Set path
RasterCASIpath07 <- "D:/Hyperspectral_data/CASI/Tiles_all/Flg07"

# Create a list of raster
CASI_raster07 <- list.files(path = RasterCASIpath07, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
CASI_raster07 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", CASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
CASI_df <- lapply(CASI_raster07[1:11], function(file) { 
  df_CASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_CASI)[3:146] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_CASI)
})

# Rbind lists
CASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

CASI_spectra07 <- CASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(CASI_df)
gc() #Free unused memory (do this often!!)

saveRDS(CASI_spectra07, file = "output_allpixels/raw_CASI_flg07.rds") #To load object : CASI_spectra07 <- readRDS("output_allpixels/raw_CASI_flg07.rds")

# Filter data based on shadows and non-vegetation
CASI_filter07 <- CASI_spectra07 %>%
  mutate(Shadow = case_when(
    Band_802nm > as.numeric(0.1) ~ "NO", 
    Band_802nm <= as.numeric(0.1) ~ "YES")) %>%
  mutate(NDVI = (Band_798nm-Band_640nm)/(Band_798nm+Band_640nm)) %>% #Red(639.630981 - Band 56) and NIR(797.516nm - Band 89)
  mutate(Vegetation = case_when(
    NDVI > 0.6 ~ "YES", 
    NDVI <= 0.6 ~ "NO")) %>%
  dplyr::filter(Shadow=="NO") %>% #remove shadowed areas
  dplyr::filter(Vegetation=="YES") #filter areas with low NDVI

rm(CASI_spectra07)
gc()

saveRDS(CASI_filter07, file = "output_allpixels/CASI_filter07.rds") #To load object : CASI_filter07 <- readRDS("output_allpixels/CASI_filter07.rds")

# Extract only the bands
CASI_bands07 <-  CASI_filter07 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
CASI_speclib07 <- speclib(spectra = t(CASI_bands07), wavelength = CASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(CASI_speclib07) <- paste0(rownames(CASI_filter07), "_",CASI_filter07$"ID") # add IDs

# Extract all important supplementary information
Meta07 <- CASI_filter07 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(CASI_speclib07) <- Meta07 
colnames(SI(CASI_speclib07)) <- colnames(Meta07)
CASI_speclib07 # check summary

# Check dimensions of speclib object
dim(CASI_speclib07) #gives nspectra and nbands (144)

#Plot raw spectra
pdf("output_allpixels/raw_CASI_spectra_07.pdf") #open a pdf device and choose file name

plot(CASI_speclib07, main = "Flightline 07 raw CASI spectra (nspectra = 2 861 160, nbands = 144)", FUN = "max", col = "red", ylim = c(0, 1))
plot(CASI_speclib07, FUN = "min", col = "blue", new = F)
plot(CASI_speclib07, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = c(350,1000), up = c(400,1100))
hsdar::mask(CASI_speclib07) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI07 <- do.call(data.frame, CASI_speclib07@SI@SI_data) #extract SI
CASI_spec07 <- CASI_speclib07@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(CASI_spec07) <- paste0("Band_", as.integer(CASI_speclib07@wavelength), "nm") #define colnames of bands
CASI_spec07 <- cbind(CASI_spec07, SI07) #bind spectra and SI

saveRDS(CASI_spec07, file = "output_allpixels/mask_CASI_flg07.rds")

rm(CASI_speclib07)
rm(CASI_spec07)
gc()

## SASI data----
# Set path
RasterSASIpath07 <- "D:/Hyperspectral_data/SASI/Tiles_all/Flg07"

# Create a list of raster
SASI_raster07 <- list.files(path = RasterSASIpath07, recursive = T,
                            pattern = "\\.TIF$", 
                            full.names = T)
SASI_raster07 #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", SASIWavelengths %>% round(), "nm")
rasternames #To see the result. 

# Apply function to extract all pixels
SASI_df <- lapply(SASI_raster07[1:10], function(file) { 
  df_SASI <- rast(file) %>% 
    terra::extract(SBL, xy = T, exact = T) %>% # extract data within areas of interest
    relocate(c(x, y , fraction), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_SASI)[3:102] <- rasternames #name columns with rasternames (columns 3 to 102 are bands)
  return(df_SASI)
})

# Rbind lists
SASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.))) #remove completely empty columns

SASI_spectra07 <- SASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

rm(SASI_df)
gc()

saveRDS(SASI_spectra07, file = "output_allpixels/raw_SASI_flg07.rds") #To load object : SASI_spectra07 <- readRDS("output_allpixels/raw_SASI_flg07.rds")

# Filter data based on shadows and non-vegetation
coord_CASI07 <- CASI_filter07$geometry #extract CASI's coordinates
SASI_filter07 <- SASI_spectra07 %>%
  dplyr::filter(geometry %in% coord_CASI07) #keep only SASI points that have same coordinates as SASI

rm(SASI_spectra07)
rm(CASI_filter07)
gc()

saveRDS(SASI_filter07, file = "output_allpixels/SASI_filter07.rds") #To load object : SASI_filter07 <- readRDS("output_allpixels/SASI_filter07.rds")

# Extract only the bands
SASI_bands07 <-  SASI_filter07 %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib (load Bands into spectral library)
SASI_speclib07 <- speclib(spectra = t(SASI_bands07), wavelength = SASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Returning and setting ID of spectra in speclib object
idSpeclib(SASI_speclib07) <- paste0(rownames(SASI_filter07), "_",SASI_filter07$"ID") # add IDs

# Extract all important supplementary information
Meta07 <- SASI_filter07 %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(SASI_speclib07) <- Meta07 
colnames(SI(SASI_speclib07)) <- colnames(Meta07)
SASI_speclib07 # check summary

# Check dimensions of speclib object
dim(SASI_speclib07) #gives nspectra and nbands (100)

#Plot raw spectra
pdf("output_allpixels/raw_SASI_spectra_07.pdf") #open a pdf device and choose file name

plot(SASI_speclib07, main = "Flightline 07 raw SASI spectra (nspectra = 2 632 309, nbands = 100)", FUN = "max", col = "red", ylim = c(0, 1))
plot(SASI_speclib07, FUN = "min", col = "blue", new = F)
plot(SASI_speclib07, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Mask edges of the sensor data (noisy)
band_mask <- data.frame(lb = 2400, up = 2500)
hsdar::mask(SASI_speclib07) <- band_mask 

# Convert speclib object to a dataframe with spectra and SI
SI07 <- do.call(data.frame, SASI_speclib07@SI@SI_data) #extract SI
SASI_spec07 <- SASI_speclib07@spectra@spectra_ma %>% data.frame() #spectra object as dataframe
colnames(SASI_spec07) <- paste0("Band_", as.integer(SASI_speclib07@wavelength), "nm") #define colnames of bands
SASI_spec07 <- cbind(SASI_spec07, SI07) #bind spectra and SI

saveRDS(SASI_spec07, file = "output_allpixels/mask_SASI_flg07.rds")

rm(SASI_speclib07)
rm(SASI_spec07)
rm(SASI_filter07)
gc()
