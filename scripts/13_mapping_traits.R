### Mapping traits and uncertainties 

## Load libraries----
library(raster)

# Set working directory
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

#### Flightline 01####----
# Import data
df_traits01 <- readRDS("inference/flg01/df_traits01.rds")
df_sd01 <- readRDS("inference/flg01/df_sd01.rds")
df_CV01 <- readRDS("inference/flg01/df_CV01.rds")

# Create a multi-band raster of mean predictions
raster_all01 <- rasterFromXYZ(df_traits01[, c("x", "y", "carbon", "carot", "cell", "chla",
                                              "chlb", "EWT", "hemi", "LDMC", "lignin",
                                              "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names <- c("carbon", "car", "cell", "chla", "chlb", "EWT",
                "hemi", "LDMC", "lignin", "LMA", "nitrogen", "SLA") #define band names
names(raster_all01) <- band_names
raster_all01 #check summary to see number of layers

raster_all01[raster_all01 < 0] <- NA #change all negative values to NA
negatives <- any(values(raster_all01) < 0) #if NA, it worked.

spat_raster01 <- as(raster_all01, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_raster01, "mapping/flg01/pred_traits01.tif", overwrite = T) #save datacube

# Create a multi-band raster of uncertainties
raster_sd01 <- rasterFromXYZ(df_sd01[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_sd <- c("sd.carbon", "sd.car", "sd.cell", "sd.chla", "sd.chlb", "sd.EWT",
                "sd.hemi", "sd.LDMC", "sd.lignin", "sd.LMA", "sd.nitrogen", "sd.SLA") #define band names
names(raster_sd01) <- band_names_sd
raster_sd01 #check summary to see number of layers

mask01 <- is.na(raster_all01) #create mask to NA values in trait raster
raster_sd01[mask01] <- NA #apply mask to sd raster

spat_rastsd01 <- as(raster_sd01, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastsd01, "mapping/flg01/sd_traits01.tif", overwrite = T) #save datacube

# Create a multi-band raster of relative uncertainties
raster_CV01 <- rasterFromXYZ(df_CV01[, c("x", "y", "SLA", "LMA", "LDMC", "EWT",
                                         "hemi", "cell", "lignin", "chla", "chlb",
                                         "carot", "carbon", "nitrogen")], crs = "EPSG:32618")
band_names_CV <- c("cv.carbon", "cv.car", "cv.cell", "cv.chla", "cv.chlb", "cv.EWT",
                "cv.hemi", "cv.LDMC", "cv.lignin", "cv.LMA", "cv.nitrogen", "cv.SLA") #define band names
names(raster_CV01) <- band_names_CV
raster_CV01 #check summary to see number of layers

mask01 <- is.na(raster_all01) #create mask to NA values in trait raster
raster_CV01[mask01] <- NA #apply mask to sd raster

spat_rastCV01 <- as(raster_CV01, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastCV01, "mapping/flg01/cv_traits01.tif", overwrite = T) #save datacube

#### Flightline 02####----
# Import data
df_traits02 <- readRDS("inference/flg02/df_traits02.rds")
df_sd02 <- readRDS("inference/flg02/df_sd02.rds")
df_CV02 <- readRDS("inference/flg02/df_CV02.rds")

# Create a multi-band raster of mean predictions
raster_all02 <- rasterFromXYZ(df_traits02[, c("x", "y", "carbon", "carot", "cell", "chla",
                                              "chlb", "EWT", "hemi", "LDMC", "lignin",
                                              "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names <- c("carbon", "car", "cell", "chla", "chlb", "EWT",
                "hemi", "LDMC", "lignin", "LMA", "nitrogen", "SLA") #define band names
names(raster_all02) <- band_names
raster_all02 #check summary to see number of layers

raster_all02[raster_all02 < 0] <- NA #change all negative values to NA
negatives <- any(values(raster_all02) < 0) #if NA, it worked.

spat_raster02 <- as(raster_all02, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_raster02, "mapping/flg02/pred_traits02.tif", overwrite = T) #save datacube

# Create a multi-band raster of uncertainties
raster_sd02 <- rasterFromXYZ(df_sd02[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_sd <- c("sd.carbon", "sd.car", "sd.cell", "sd.chla", "sd.chlb", "sd.EWT",
                   "sd.hemi", "sd.LDMC", "sd.lignin", "sd.LMA", "sd.nitrogen", "sd.SLA") #define band names
names(raster_sd02) <- band_names_sd
raster_sd02 #check summary to see number of layers

mask02 <- is.na(raster_all02) #create mask to NA values in trait raster
raster_sd02[mask02] <- NA #apply mask to sd raster

spat_rastsd02 <- as(raster_sd02, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastsd02, "mapping/flg02/sd_traits02.tif", overwrite = T) #save datacube

# Create a multi-band raster of relative uncertainties
raster_CV02 <- rasterFromXYZ(df_CV02[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_CV <- c("cv.carbon", "cv.car", "cv.cell", "cv.chla", "cv.chlb", "cv.EWT",
                   "cv.hemi", "cv.LDMC", "cv.lignin", "cv.LMA", "cv.nitrogen", "cv.SLA") #define band names
names(raster_CV02) <- band_names_CV
raster_CV02 #check summary to see number of layers

mask02 <- is.na(raster_all02) #create mask to NA values in trait raster
raster_CV02[mask02] <- NA #apply mask to sd raster

spat_rastCV02 <- as(raster_CV02, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastCV02, "mapping/flg02/cv_traits02.tif", overwrite = T) #save datacube

#### Flightline 03####----
# Import data
df_traits03 <- readRDS("inference/flg03/df_traits03.rds")
df_sd03 <- readRDS("inference/flg03/df_sd03.rds")
df_CV03 <- readRDS("inference/flg03/df_CV03.rds")

# Create a multi-band raster of mean predictions
raster_all03 <- rasterFromXYZ(df_traits03[, c("x", "y", "carbon", "carot", "cell", "chla",
                                              "chlb", "EWT", "hemi", "LDMC", "lignin",
                                              "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names <- c("carbon", "car", "cell", "chla", "chlb", "EWT",
                "hemi", "LDMC", "lignin", "LMA", "nitrogen", "SLA") #define band names
names(raster_all03) <- band_names
raster_all03 #check summary to see number of layers

raster_all03[raster_all03 < 0] <- NA #change all negative values to NA
negatives <- any(values(raster_all03) < 0) #if NA, it worked.

spat_raster03 <- as(raster_all03, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_raster03, "mapping/flg03/pred_traits03.tif", overwrite = T) #save datacube

# Create a multi-band raster of uncertainties
raster_sd03 <- rasterFromXYZ(df_sd03[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_sd <- c("sd.carbon", "sd.car", "sd.cell", "sd.chla", "sd.chlb", "sd.EWT",
                   "sd.hemi", "sd.LDMC", "sd.lignin", "sd.LMA", "sd.nitrogen", "sd.SLA") #define band names
names(raster_sd03) <- band_names_sd
raster_sd03 #check summary to see number of layers

mask03 <- is.na(raster_all03) #create mask to NA values in trait raster
raster_sd03[mask03] <- NA #apply mask to sd raster

spat_rastsd03 <- as(raster_sd03, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastsd03, "mapping/flg03/sd_traits03.tif", overwrite = T) #save datacube

# Create a multi-band raster of relative uncertainties
raster_CV03 <- rasterFromXYZ(df_CV03[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_CV <- c("cv.carbon", "cv.car", "cv.cell", "cv.chla", "cv.chlb", "cv.EWT",
                   "cv.hemi", "cv.LDMC", "cv.lignin", "cv.LMA", "cv.nitrogen", "cv.SLA") #define band names
names(raster_CV03) <- band_names_CV
raster_CV03 #check summary to see number of layers

mask03 <- is.na(raster_all03) #create mask to NA values in trait raster
raster_CV03[mask03] <- NA #apply mask to sd raster

spat_rastCV03 <- as(raster_CV03, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastCV03, "mapping/flg03/cv_traits03.tif", overwrite = T) #save datacube

#### Flightline 04####----
# Import data
df_traits04 <- readRDS("inference/flg04/df_traits04.rds")
df_sd04 <- readRDS("inference/flg04/df_sd04.rds")
df_CV04 <- readRDS("inference/flg04/df_CV04.rds")

# Create a multi-band raster of mean predictions
raster_all04 <- rasterFromXYZ(df_traits04[, c("x", "y", "carbon", "carot", "cell", "chla",
                                              "chlb", "EWT", "hemi", "LDMC", "lignin",
                                              "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names <- c("carbon", "car", "cell", "chla", "chlb", "EWT",
                "hemi", "LDMC", "lignin", "LMA", "nitrogen", "SLA") #define band names
names(raster_all04) <- band_names
raster_all04 #check summary to see number of layers

raster_all04[raster_all04 < 0] <- NA #change all negative values to NA
negatives <- any(values(raster_all04) < 0) #if NA, it worked.

spat_raster04 <- as(raster_all04, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_raster04, "mapping/flg04/pred_traits04.tif", overwrite = T) #save datacube

# Create a multi-band raster of uncertainties
raster_sd04 <- rasterFromXYZ(df_sd04[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_sd <- c("sd.carbon", "sd.car", "sd.cell", "sd.chla", "sd.chlb", "sd.EWT",
                   "sd.hemi", "sd.LDMC", "sd.lignin", "sd.LMA", "sd.nitrogen", "sd.SLA") #define band names
names(raster_sd04) <- band_names_sd
raster_sd04 #check summary to see number of layers

mask04 <- is.na(raster_all04) #create mask to NA values in trait raster
raster_sd04[mask04] <- NA #apply mask to sd raster

spat_rastsd04 <- as(raster_sd04, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastsd04, "mapping/flg04/sd_traits04.tif", overwrite = T) #save datacube

# Create a multi-band raster of relative uncertainties
raster_CV04 <- rasterFromXYZ(df_CV04[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_CV <- c("cv.carbon", "cv.car", "cv.cell", "cv.chla", "cv.chlb", "cv.EWT",
                   "cv.hemi", "cv.LDMC", "cv.lignin", "cv.LMA", "cv.nitrogen", "cv.SLA") #define band names
names(raster_CV04) <- band_names_CV
raster_CV04 #check summary to see number of layers

mask04 <- is.na(raster_all04) #create mask to NA values in trait raster
raster_CV04[mask04] <- NA #apply mask to sd raster

spat_rastCV04 <- as(raster_CV04, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastCV04, "mapping/flg04/cv_traits04.tif", overwrite = T) #save datacube

#### Flightline 05####----
# Import data
df_traits05 <- readRDS("inference/flg05/df_traits05.rds")
df_sd05 <- readRDS("inference/flg05/df_sd05.rds")
df_CV05 <- readRDS("inference/flg05/df_CV05.rds")

# Create a multi-band raster of mean predictions
raster_all05 <- rasterFromXYZ(df_traits05[, c("x", "y", "carbon", "carot", "cell", "chla",
                                              "chlb", "EWT", "hemi", "LDMC", "lignin",
                                              "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names <- c("carbon", "car", "cell", "chla", "chlb", "EWT",
                "hemi", "LDMC", "lignin", "LMA", "nitrogen", "SLA") #define band names
names(raster_all05) <- band_names
raster_all05 #check summary to see number of layers

raster_all05[raster_all05 < 0] <- NA #change all negative values to NA
negatives <- any(values(raster_all05) < 0) #if NA, it worked.

spat_raster05 <- as(raster_all05, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_raster05, "mapping/flg05/pred_traits05.tif", overwrite = T) #save datacube

# Create a multi-band raster of uncertainties
raster_sd05 <- rasterFromXYZ(df_sd05[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_sd <- c("sd.carbon", "sd.car", "sd.cell", "sd.chla", "sd.chlb", "sd.EWT",
                   "sd.hemi", "sd.LDMC", "sd.lignin", "sd.LMA", "sd.nitrogen", "sd.SLA") #define band names
names(raster_sd05) <- band_names_sd
raster_sd05 #check summary to see number of layers

mask05 <- is.na(raster_all05) #create mask to NA values in trait raster
raster_sd05[mask05] <- NA #apply mask to sd raster

spat_rastsd05 <- as(raster_sd05, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastsd05, "mapping/flg05/sd_traits05.tif", overwrite = T) #save datacube

# Create a multi-band raster of relative uncertainties
raster_CV05 <- rasterFromXYZ(df_CV05[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_CV <- c("cv.carbon", "cv.car", "cv.cell", "cv.chla", "cv.chlb", "cv.EWT",
                   "cv.hemi", "cv.LDMC", "cv.lignin", "cv.LMA", "cv.nitrogen", "cv.SLA") #define band names
names(raster_CV05) <- band_names_CV
raster_CV05 #check summary to see number of layers

mask05 <- is.na(raster_all05) #create mask to NA values in trait raster
raster_CV05[mask05] <- NA #apply mask to sd raster

spat_rastCV05 <- as(raster_CV05, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastCV05, "mapping/flg05/cv_traits05.tif", overwrite = T) #save datacube

#### Flightline 06####----
# Import data
df_traits06 <- readRDS("inference/flg06/df_traits06.rds")
df_sd06 <- readRDS("inference/flg06/df_sd06.rds")
df_CV06 <- readRDS("inference/flg06/df_CV06.rds")

# Create a multi-band raster of mean predictions
raster_all06 <- rasterFromXYZ(df_traits06[, c("x", "y", "carbon", "carot", "cell", "chla",
                                              "chlb", "EWT", "hemi", "LDMC", "lignin",
                                              "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names <- c("carbon", "car", "cell", "chla", "chlb", "EWT",
                "hemi", "LDMC", "lignin", "LMA", "nitrogen", "SLA") #define band names
names(raster_all06) <- band_names
raster_all06 #check summary to see number of layers

raster_all06[raster_all06 < 0] <- NA #change all negative values to NA
negatives <- any(values(raster_all06) < 0) #if NA, it worked.

spat_raster06 <- as(raster_all06, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_raster06, "mapping/flg06/pred_traits06.tif", overwrite = T) #save datacube

# Create a multi-band raster of uncertainties
raster_sd06 <- rasterFromXYZ(df_sd06[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_sd <- c("sd.carbon", "sd.car", "sd.cell", "sd.chla", "sd.chlb", "sd.EWT",
                   "sd.hemi", "sd.LDMC", "sd.lignin", "sd.LMA", "sd.nitrogen", "sd.SLA") #define band names
names(raster_sd06) <- band_names_sd
raster_sd06 #check summary to see number of layers

mask06 <- is.na(raster_all06) #create mask to NA values in trait raster
raster_sd06[mask06] <- NA #apply mask to sd raster

spat_rastsd06 <- as(raster_sd06, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastsd06, "mapping/flg06/sd_traits06.tif", overwrite = T) #save datacube

# Create a multi-band raster of relative uncertainties
raster_CV06 <- rasterFromXYZ(df_CV06[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_CV <- c("cv.carbon", "cv.car", "cv.cell", "cv.chla", "cv.chlb", "cv.EWT",
                   "cv.hemi", "cv.LDMC", "cv.lignin", "cv.LMA", "cv.nitrogen", "cv.SLA") #define band names
names(raster_CV06) <- band_names_CV
raster_CV06 #check summary to see number of layers

mask06 <- is.na(raster_all06) #create mask to NA values in trait raster
raster_CV06[mask06] <- NA #apply mask to sd raster

spat_rastCV06 <- as(raster_CV06, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastCV06, "mapping/flg06/cv_traits06.tif", overwrite = T) #save datacube

#### Flightline 07####----
# Import data
df_traits07 <- readRDS("inference/flg07/df_traits07.rds")
df_sd07 <- readRDS("inference/flg07/df_sd07.rds")
df_CV07 <- readRDS("inference/flg07/df_CV07.rds")

# Create a multi-band raster of mean predictions
raster_all07 <- rasterFromXYZ(df_traits07[, c("x", "y", "carbon", "carot", "cell", "chla",
                                              "chlb", "EWT", "hemi", "LDMC", "lignin",
                                              "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names <- c("carbon", "car", "cell", "chla", "chlb", "EWT",
                "hemi", "LDMC", "lignin", "LMA", "nitrogen", "SLA") #define band names
names(raster_all07) <- band_names
raster_all07 #check summary to see number of layers

raster_all07[raster_all07 < 0] <- NA #change all negative values to NA
negatives <- any(values(raster_all07) < 0) #if NA, it worked.

spat_raster07 <- as(raster_all07, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_raster07, "mapping/flg07/pred_traits07.tif", overwrite = T) #save datacube

# Create a multi-band raster of uncertainties
raster_sd07 <- rasterFromXYZ(df_sd07[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_sd <- c("sd.carbon", "sd.car", "sd.cell", "sd.chla", "sd.chlb", "sd.EWT",
                   "sd.hemi", "sd.LDMC", "sd.lignin", "sd.LMA", "sd.nitrogen", "sd.SLA") #define band names
names(raster_sd07) <- band_names_sd
raster_sd07 #check summary to see number of layers

mask07 <- is.na(raster_all07) #create mask to NA values in trait raster
raster_sd07[mask07] <- NA #apply mask to sd raster

spat_rastsd07 <- as(raster_sd07, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastsd07, "mapping/flg07/sd_traits07.tif", overwrite = T) #save datacube

# Create a multi-band raster of relative uncertainties
raster_CV07 <- rasterFromXYZ(df_CV07[, c("x", "y", "carbon", "carot", "cell", "chla",
                                         "chlb", "EWT", "hemi", "LDMC", "lignin",
                                         "LMA", "nitrogen", "SLA")], crs = "EPSG:32618")
band_names_CV <- c("cv.carbon", "cv.car", "cv.cell", "cv.chla", "cv.chlb", "cv.EWT",
                   "cv.hemi", "cv.LDMC", "cv.lignin", "cv.LMA", "cv.nitrogen", "cv.SLA") #define band names
names(raster_CV07) <- band_names_CV
raster_CV07 #check summary to see number of layers

mask07 <- is.na(raster_all07) #create mask to NA values in trait raster
raster_CV07[mask07] <- NA #apply mask to sd raster

spat_rastCV07 <- as(raster_CV07, "SpatRaster") #change to spatraster object (for band names)
writeRaster(spat_rastCV07, "mapping/flg07/cv_traits07.tif", overwrite = T) #save datacube