### Extract spectra from raster for polygons

## Install packages and load libraries----
install.packages("hsdar", "sf", "terra")
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

## Read raster format and load only polygons to reduce size of data----

# Set coordinate reference systems (CRS)
CRS_SBL <- 32618 # epsg WGS84 UTM18N : 32618

# Import polygons
polygons <- st_read("input_data/2022-POLYGON-DATA/tree_data.shp") %>% 
  st_transform(crs=CRS_SBL)

# Plot polygons
plot(st_geometry(polygons)) #show plot with only geometry (not all attributes)
plot(polygons["label"]) #show plot with attribute "label"

test <- dplyr::filter(polygons, plant_id == 135134665) #select a sample to test
plot(st_geometry(test))

# Create inner buffer of 0.5 m for each polygons
polygon_buffer <- st_buffer(polygons, dist = -0.5, preserveTopology = T)

test2 <- dplyr::filter(polygon_buffer, plant_id == 135134665) #select same sample
plot(st_geometry(test2)) #plot it to see if buffer worked

# Calculate area of buffer polygons
polygon_buffer$buffer_area <- as.numeric(st_area(polygon_buffer))
polygons$buffer_area <- as.numeric(st_area(polygon_buffer)) #add area buffer to polygon layer

st_write(polygon_buffer, "input_data/2022-POLYGON-DATA/polygon_buffer.geojson", delete_dsn = T, delete_layer = T)

####---------------------------------------------------------------------------------------####

#### CASI data---- ####
# Set path
RasterCASIpath <- "D:/Hyperspectral_data/CASI/L2G"

# Create a list of raster
CASI_list_of_rasters <- list.files(path = RasterCASIpath, recursive = T,
                                 pattern = "Resample.dat", 
                                 full.names = T)
head(CASI_list_of_rasters) #to see results of the list

# Define raster's names
rasternames <- paste0("Band_", CASIWavelengths %>% round(), "nm")
rasternames #check result


# Apply function to extract data of interest (polygons)
CASI_df <- lapply(CASI_list_of_rasters[2:4], function(file) { #polygons are in flightline 2 to 4
  df_CASI <- rast(file) %>% 
    terra::extract(polygon_buffer, xy = T, exact = T) %>% # extract data within area of interest
    mutate(Flightline = substr(basename(file), 2, 14)) %>% #add flightline name as column to table
    relocate(c(x, y, fraction, Flightline), .after = "ID") %>% #fraction is the column from argument "exact" in terra::extract
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_CASI)[4:147] <- rasternames #name columns with rasternames (columns 4 to 147 are bands)
  return(df_CASI)
})

### FOR TERRA::EXTRACT, if exact = TRUE and y has polygons, ###
### the exact fraction of each cell that is covered is ###
### returned as well, for example to compute a weighted mean. ###

# Rbind lists
CASI_df %<>% bind_rows() %>% #operator %<>% is assignment pipe, updates the object.
  select_if(~!all(is.na(.))) #remove completely empty columns

CASI_spectra <- CASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectra completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

# Add polygons metadata by joining points to polygons
# choose "polygons" object or else a tree is missing because points aren't within the buffer
join <- st_join(CASI_spectra, polygons, join = st_within, largest = T, left = F) %>%
  relocate(c(plant_id, specie, label, drainage, long, lat), .before = "ID") %>%
  dplyr::select(-"ID")

length(unique(join$plant_id)) #check if there's 170 trees

# Remove pixels that are not fully in the polygons
## Add a cutoff so that pixels with a fraction > 0.4 will be kept ONLY if ##
## polygon's area < 9m2 (or else we won't have any pixels, these trees are too small) ##
CASI_spectra <- join %>% 
  mutate(cutoff = ifelse(buffer_area < 9, 0.4, 1)) %>% #add a cutoff column
  dplyr::filter(fraction >= cutoff) #keep only pixels that are equal or greater than the cutoff

length(unique(CASI_spectra$plant_id)) #check if there's 170 trees

# Filter data based on shadows and non-vegetation
CASI_filter <- CASI_spectra %>%
  mutate(Shadow = case_when(
    Band_802nm > as.numeric(0.1) ~ "NO", 
    Band_802nm <= as.numeric(0.1) ~ "YES")) %>%
  mutate(NDVI = (Band_798nm-Band_640nm)/(Band_798nm+Band_640nm)) %>% #Red(639.630981 - Band 56) and NIR(797.516nm - Band 89)
  mutate(Vegetation = case_when(
    NDVI > 0.6 ~ "YES", 
    NDVI <= 0.6 ~ "NO")) %>%
  dplyr::filter(Shadow=="NO") %>% #remove shadowed areas
  dplyr::filter(Vegetation=="YES") #filter areas with low NDVI

length(unique(CASI_filter$plant_id)) #check number of trees`= 166

# Check for NA values
apply(CASI_filter, 2, function(x) any(is.na(x))) #FALSE means there are no NAs

# Save data
st_write(CASI_filter, "output_spectra/CASI_filter.geojson", delete_dsn = T, delete_layer = T)

## Convert to a spectral library using hsdar package----
# Extract only the bands
CASI_bands <-  CASI_filter %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib
CASI_speclib <- speclib(spectra = t(CASI_bands), wavelength = CASIWavelengths) #t() function is used to calculate transpose of a matrix or df in R

# Set ID of spectra in speclib object
idSpeclib(CASI_speclib) <- paste0(rownames(CASI_filter), "_",CASI_filter$"plant_id") # add IDs

# Extract all important supplementary information
Meta <- CASI_filter %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(CASI_speclib) <- Meta 
colnames(SI(CASI_speclib)) <- colnames(Meta)
CASI_speclib # check summary

# Check dimensions of speclib object
dim(CASI_speclib) #gives nspectra and nbands (144)

#Plot raw spectra
pdf("figures/raw_CASI_spectra.pdf") #open a pdf device and choose file name

plot(CASI_speclib, main = "Raw CASI spectra (nspectra = 1689, nbands = 144)", FUN = "max", col = "red")
plot(CASI_speclib, FUN = "min", col = "blue", new = F)
plot(CASI_speclib, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Save speclib oject
saveRDS(CASI_speclib, file = "input_data/CASI_speclib_raw.rds") #To load object : CASI_speclib <- readRDS("input_data/CASI_speclib_raw.rds")

## Mask edges of the sensor data (noisy)----
plot(CASI_speclib)
band_mask <- data.frame(lb = c(350, 1000), up = c(400, 1100))
hsdar::mask(CASI_speclib) <- band_mask 
plot(CASI_speclib)
plot(CASI_speclib, FUN = 1:50)

## Divide into flightlines and calculate mean for each polygons----
CASI_02 <- subset(CASI_speclib, Flightline == "220709_SBL-02") #create a subset for each flightline
CASI_03 <- subset(CASI_speclib, Flightline == "220709_SBL-03")
CASI_04 <- subset(CASI_speclib, Flightline == "220709_SBL-04")
mean_02 <- apply(CASI_02, FUN = mean, bySI = "plant_id") #polygons identified by plant_id in SI
mean_03 <- apply(CASI_03, FUN = mean, bySI = "plant_id")
mean_04 <- apply(CASI_04, FUN = mean, bySI = "plant_id")

## Plot all speclibs----
plot(mean_02, col = "orange")
plot(mean_03, col = "green", new = F)
plot(mean_04, col = "red", new = F)

## Check trees that are more in the center of flightline 02 and compare with flightline 03----
plot(subset(mean_02, plant_id == 135133186), col = "orange", ylim = c(0, 1))
plot(subset(mean_03, plant_id == 135133186), col = "green", new = F) #choose 02

plot(subset(mean_02, plant_id == 135133135), col = "orange", ylim = c(0, 1))
plot(subset(mean_03, plant_id == 135133135), col = "green", new = F) #choose 03

plot(subset(mean_02, plant_id == 135141907), col = "orange", ylim = c(0, 1)) #doesn't work (removed with NDVI/shadow filter)
plot(subset(mean_03, plant_id == 135141907), col = "green", new = F) #choose 03

plot(subset(mean_02, plant_id == 135135481), col = "orange", ylim = c(0, 1))
plot(subset(mean_03, plant_id == 135135481), col = "green", new = F) #choose 02

plot(subset(mean_02, plant_id == 135142315), col = "orange", ylim = c(0, 1))
plot(subset(mean_03, plant_id == 135142315), col = "green", new = F) #choose 03 (very similar)

plot(subset(mean_02, plant_id == 135135889), col = "orange", ylim = c(0, 1))
plot(subset(mean_03, plant_id == 135135889), col = "green", new = F) #choose 03 (very similar)

plot(subset(mean_02, plant_id == 135130789), col = "orange", ylim = c(0, 1))
plot(subset(mean_03, plant_id == 135130789), col = "green", new = F) #choose 02 (very similar)

## There are 66 trees in overlap of flightline 02 and 03.
## Except from those specified above, I will choose all trees from flightline 03.

## Check trees that are at limit of flightline 03 and compare spectra with flightline 04----
plot(subset(mean_03, plant_id == 135135685), col = "green", ylim = c(0, 1))
plot(subset(mean_04, plant_id == 135135685), col = "red", new = F) #choose 04

plot(subset(mean_03, plant_id == 135131248), col = "green", ylim = c(0, 1))
plot(subset(mean_04, plant_id == 135131248), col = "red", new = F) #choose 04 (very similar)

plot(subset(mean_03, plant_id == 135128494), col = "green", ylim = c(0, 1)) #doesn't work (removed with NDVI/shadow filter)
plot(subset(mean_04, plant_id == 135128494), col = "red", new = F) #choose 04

plot(subset(mean_03, plant_id == 135136195), col = "green", ylim = c(0, 1))
plot(subset(mean_04, plant_id == 135136195), col = "red", new = F) #choose 03

new_mean_02 <- subset(mean_02, plant_id %in% c(135133186, 135135481, 135130789)) #subset chosen spectra from flg02
new_mean_03 <- subset(mean_03, plant_id != 135133186 & plant_id != 135135481
                      & plant_id != 135130789 & plant_id != 135135685
                      & plant_id != 135131248 & plant_id != 135128494) #subset all spectra from flg03 expect the ones I didn't chose
new_mean_04 <- subset(mean_04, plant_id %in% c(135140989,
                                               135136705,
                                               135131554,
                                               135144661,
                                               135135685,
                                               135131248,
                                               135128494)) #subset chosen spectra from flg04
merge_speclib <- merge(new_mean_02, new_mean_03, new_mean_04) #merge all

# Plot speclib object
plot(merge_speclib, col = "navy")
plot(merge_speclib, FUN = 1:166)

# Check dimensions of speclib object
dim(merge_speclib)

# Plot results
pdf("figures/mean_CASI_spectra.pdf")

plot(merge_speclib, main = "Mean CASI spectra (nspectra = 166, nbands = 126)", FUN = "max", col = "red")
plot(merge_speclib, FUN = "min", col = "blue", new = F)
plot(merge_speclib, FUN = "mean", col = "orange", new = F)

dev.off()

# Save speclib oject
saveRDS(merge_speclib, file = "input_data/CASI_speclib_mean.rds") #To load object : merge_speclib <- readRDS("input_data/CASI_speclib_mean.rds")

## Convert speclib object to a dataframe with spectra and SI----
# Extract SI
SI_mean <- do.call(data.frame, merge_speclib@SI@SI_data)

# Spectra object as dataframe
CASI_mean <- merge_speclib@spectra@spectra_ma %>% data.frame()

# Define colnames of bands
colnames(CASI_mean) <- paste0("Band_", as.integer(merge_speclib@wavelength), "nm")

# Bind spectra and SI
CASI_mean <- cbind(CASI_mean, SI_mean)
CASI_mean

# Join spectra and coordinates of polygon_buffer
polygon_xy <- polygon_buffer %>%
  dplyr::select(plant_id, long,lat, specie) #select only these columns

CASI_mean_xy <- right_join(CASI_mean, polygon_xy, by = "plant_id") 
length(unique(CASI_mean_xy$plant_id)) #count the number of unique plant_id (should be 170 here)

# Convert sf object to dataframe
CASI_mean_large <- CASI_mean_xy %>%
  dplyr::select(-geometry) %>%
  relocate(c(plant_id, long , lat, specie), .before = "Band_400nm")
write_csv(CASI_mean_large, "output_spectra/CASI_mean_large.csv", col_names = T) #save data

## Change to long format to plot----
colnames(CASI_mean_large) <- gsub("Band_|nm", "", colnames(CASI_mean_large)) #rename wavelength columns
CASI_mean_long <- pivot_longer(CASI_mean_large, cols = - c(plant_id,long, lat, specie), names_to = "wavelength", values_to = "val" )
CASI_mean_long$wavelength <- as.numeric(CASI_mean_long$wavelength) #change column wavelength to numeric
write_csv(CASI_mean_long, "output_spectra/CASI_mean_long.csv", col_names = T) #save data

# Visualise all spectra
pdf("figures/CASI_mean_ggplot.pdf")

ggplot(data = CASI_mean_long, aes(x = wavelength,
                                    y = val,
                                    group = plant_id)) + 
  geom_line() +
  ggtitle("Canopy-level reflectance spectra (CASI sensor)") +
  xlab("Wavelength (nm)") +
  ylab("Reflectance")

dev.off()

####---------------------------------------------------------------------------------------####

#### SASI data ---- ####

# Set path
RasterSASIpath <- "D:/Hyperspectral_data/SASI/L2G"

# Create a list of raster
SASI_list_of_rasters <- list.files(path = RasterSASIpath, recursive = T,
                                   pattern = "aligned.dat", 
                                   full.names = T)
head(SASI_list_of_rasters)

# Define raster's names
rasternames1 <- paste0("Band_", SASIWavelengths %>% round(), "nm")
rasternames1


# Apply function to extract data of interest (polygons)
SASI_df <- lapply(SASI_list_of_rasters[2:4], function(file) { 
  df_SASI <- rast(file) %>% 
    terra::extract(polygon_buffer, xy = T, exact = T) %>% # extract data within areas of interest
    mutate(Flightline = substr(basename(file), 1, 13)) %>% #add flightline name as column to table
    relocate(c(x, y , fraction, Flightline), .after = "ID") %>% 
    dplyr::filter(rowSums(is.na(.)) == 0) %>% #remove all rows with completely NA values
    st_as_sf(coords = c("x","y"), crs = CRS_SBL)  #convert to sf point object using the initial crs
  colnames(df_SASI)[4:103] <- rasternames1 #name columns with rasternames (columns 4 to 103 are bands)
  return(df_SASI)
})

# Rbind lists
SASI_df %<>% bind_rows() %>% 
  select_if(~!all(is.na(.)))

SASI_spectra <- SASI_df %>%
  group_by(ID) %>% 
  filter_at(vars(starts_with("Band")), any_vars(. != 0)) %>% #remove all spectral completely zero (black background)
  mutate_at(grep("Band", colnames(.)), list(~./10000))  #scale spectra

# Join points to polygons
join1 <- st_join(SASI_spectra, polygons, join = st_within, largest = T, left = F) %>%
  relocate(c(plant_id, specie, label, drainage, long, lat), .before = "ID") %>%
  dplyr::select(-"ID")

length(unique(join1$plant_id)) #check if there's 170 trees

# Remove pixels that are not fully in the polygons
SASI_spectra <- join1 %>% 
  mutate(cutoff = ifelse(buffer_area < 9, 0.4, 1)) %>%
  dplyr::filter(fraction >= cutoff)

length(unique(SASI_spectra$plant_id))

coord_CASI <- CASI_filter$geometry #extract CASI's coordinates
SASI_filter <- SASI_spectra %>%
  dplyr::filter(geometry %in% coord_CASI) #keep only SASI points that have same coordinates as SASI

length(unique(SASI_filter$plant_id))

# Check for NA values
apply(SASI_filter, 2, function(x) any(is.na(x))) #FALSE means there are no NAs

# Save data
st_write(SASI_filter, "output_spectra/SASI_filter.geojson", delete_dsn = T, delete_layer = T)

## Convert to a spectral library using hsdar package----
# Extract only the bands
SASI_bands <-  SASI_filter %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  dplyr::select(grep("Band", colnames(.)))

# Create object of class speclib
SASI_speclib <- speclib(spectra = t(SASI_bands), wavelength = SASIWavelengths)

# Set ID of spectra in speclib object
idSpeclib(SASI_speclib) <- paste0(rownames(SASI_filter), "_",SASI_filter$"plant_id") # add IDs

# Extract all important supplementary information
Meta1 <- SASI_filter %>% 
  ungroup() %>% 
  dplyr::select(!grep("Band",names(.))) %>% # selects all columns except the bands
  dplyr::mutate(x=sf::st_coordinates(.)[,1]) %>% # add coordinates to Metainformation
  dplyr::mutate(y=sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

# Add metainformation to the supplementary information of the spectral library
SI(SASI_speclib) <- Meta1 
colnames(SI(SASI_speclib)) <- colnames(Meta1)
SASI_speclib # check summary

# Check dimensions of speclib object
dim(SASI_speclib) #gives nspectra and nbands 

#Plot spectra
pdf("figures/raw_SASI_spectra.pdf")

plot(SASI_speclib, main = "Raw SASI spectra (nspectra = 1684 , nbands = 100)", FUN = "max", col = "red", ylim = c(0, 0.8))
plot(SASI_speclib, FUN = "min", col = "blue", new = F)
plot(SASI_speclib, FUN = "mean", col = "orange", new = F)

dev.off()

# Save speclib oject
saveRDS(SASI_speclib, file = "input_data/SASI_speclib_raw.rds") # To load object : SASI_speclib <- readRDS("input_data/SASI_speclib_raw.rds")

## Mask edges of the sensor data (noisy)----
plot(SASI_speclib)
band_mask1 <- data.frame(lb = 2400, up = 2500)
hsdar::mask(SASI_speclib) <- band_mask1 
plot(SASI_speclib)
plot(SASI_speclib, FUN = 1:50)

## Divide into both flightlines and calculate mean for each polygons----
SASI_02 <- subset(SASI_speclib, Flightline == "220709_SBL-02")
SASI_03 <- subset(SASI_speclib, Flightline == "220709_SBL-03")
SASI_04 <- subset(SASI_speclib, Flightline == "220709_SBL-04")
SASImean_02 <- apply(SASI_02, FUN = mean, bySI = "plant_id")
SASImean_03 <- apply(SASI_03, FUN = mean, bySI = "plant_id") 
SASImean_04 <- apply(SASI_04, FUN = mean, bySI = "plant_id")

## Plot all speclibs and choose same flightline as CASI----
plot(SASImean_02, col = "orange", ylim = c(0, 1))
plot(SASImean_03, col = "green", new = F)
plot(SASImean_04, col = "red", new = F)

new_SASImean_02 <- subset(SASImean_02, plant_id %in% c(135133186, 135135481, 135130789))
new_SASImean_03 <- subset(SASImean_03, plant_id != 135133186 & plant_id != 135135481
                          & plant_id != 135130789 & plant_id != 135135685
                          & plant_id != 135131248 & plant_id != 135128494)
new_SASImean_04 <- subset(SASImean_04, plant_id %in% c(135140989,
                                                       135136705,
                                                       135131554,
                                                       135144661,
                                                       135135685,
                                                       135131248,
                                                       135128494))
merge_speclib1 <- merge(new_SASImean_02, new_SASImean_03, new_SASImean_04) #merge all

# Plot speclib object
plot(merge_speclib1, col = "green")
plot(merge_speclib1, FUN = 1:166)

# Check dimensions of speclib object
dim(merge_speclib1)

# Plot results
pdf("figures/mean_SASI_spectra.pdf") #open a pdf device and choose file name

plot(merge_speclib1, main = "Mean SASI spectra (nspectra = 166, nbands = 97)", FUN = "max", col = "red", ylim = c(0, 0.8))
plot(merge_speclib1, FUN = "min", col = "blue", new = F)
plot(merge_speclib1, FUN = "mean", col = "orange", new = F)

dev.off() #close the pdf device to save plot

# Save speclib oject
saveRDS(merge_speclib1, file = "input_data/SASI_speclib_mean.rds") # To load object : merge_speclib1 <- readRDS("input_data/SASI_speclib_mean.rds")

## Convert speclib object to a dataframe with spectra and SI----
# Extract SI
SI_mean1 <- do.call(data.frame, merge_speclib1@SI@SI_data)

# Object .spectra as dataframe
SASI_mean <- merge_speclib1@spectra@spectra_ma %>% data.frame()

# Define colnames of bands
colnames(SASI_mean) <- paste0("Band_", as.integer(merge_speclib1@wavelength), "nm")

# Bind spectra and SI
SASI_mean <- cbind(SASI_mean, SI_mean1)
SASI_mean

# Join spectra and coordinates of polygon_buffer
SASI_mean_xy <- right_join(SASI_mean, polygon_xy, by = "plant_id") 
length(unique(SASI_mean_xy$plant_id)) #count the number of unique plant_id = 170! GREAT!

# Convert sf object to dataframe
SASI_mean_large <- SASI_mean_xy %>%
  dplyr::select(-geometry) %>%
  relocate(c(plant_id, long , lat, specie), .before = "Band_957nm")
write_csv(SASI_mean_large, "output_spectra/SASI_mean_large.csv", col_names = T)

## Change to long format to plot----
colnames(SASI_mean_large) <- gsub("Band_|nm", "", colnames(SASI_mean_large)) #rename wavelength columns
SASI_mean_long <- pivot_longer(SASI_mean_large, cols = - c(plant_id, long, lat, specie), names_to = "wavelength", values_to = "val" )
SASI_mean_long$wavelength <- as.numeric(SASI_mean_long$wavelength) #change column wavelength to numeric
write_csv(SASI_mean_long, "output_spectra/SASI_mean_long.csv", col_names = T) #save data

# Visualise all spectra
pdf("figures/SASI_mean_ggplot.pdf")

ggplot(data = SASI_mean_long, aes(x = wavelength,
                                    y = val,
                                    group = plant_id)) + 
  geom_line() +
  ggtitle("Canopy-level reflectance spectra (SASI sensor)") +
  xlab("Wavelength (nm)") +
  ylab("Reflectance")

dev.off()
