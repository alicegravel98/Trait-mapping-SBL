### Predicted vs measured trait

## Load libraries----
library(tidyverse)
library(sf)
library(raster)
library(terra)
library(egg)

## Set working directory----
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Import data----
CRS_SBL <- 32618 # epsg WGS84 UTM18N : 32618

# Import measured traits
traits <- c("SLA", "LMA", "LDMC", "EWT", "hemicellulose", "cellulose",
            "lignin", "chla", "chlb", "carotenoids", "carbon", "nitrogen")

measured_traits <- read.csv("input_data/plsr_data.csv", header = T, sep = ",", check.names = F) %>%
  dplyr::select(all_of(traits), "specie", "plant_id")

unique(measured_traits$specie)
length(unique(measured_traits$specie)) #16 classes

# Import trait datacube
mosaic_pred <- terra::rast("mapping/mosaic_pred_traits.tif") #trait predictions
names(mosaic_pred) <- traits #define layers' names
mosaic_pred #check result

# Import polygons
Z1 <- st_read("input_data/polygons_Cloutier/Z1_polygons.shp") %>% 
  st_transform(crs = CRS_SBL)
plot(Z1["Label"])

Z2 <- st_read("input_data/polygons_Cloutier/Z2_polygons.shp") %>% 
  st_transform(crs = CRS_SBL)
plot(Z2["Label"])

Z3 <- st_read("input_data/polygons_Cloutier/Z3_polygons.shp") %>% 
  st_transform(crs = CRS_SBL)
plot(Z3["Label"])

# Merge zones
merged_polygons <- rbind(Z1, Z2, Z3)
plot(merged_polygons["Label"])

# Define unwanted labels
labels <- unique(merged_polygons$Label)
labels #check result
remove_labels <- c("Acer", "Mort", "Feuillus",
                   "Conifere", "FRNI", "Betula", "Populus",
                   "PRPE", "POBA", "BEPO")
labels_to_keep <- labels[!labels %in% remove_labels]
labels_to_keep #check result

# Remove unwanted labels from merged_polygons
merged_polygons <- merged_polygons %>%
  dplyr::filter(Label %in% labels_to_keep) %>%
  dplyr::rename(label = Label)
unique(merged_polygons$label) #check result

merged_polygons$label <- ifelse(merged_polygons$label == "PIMA" |
                                  merged_polygons$label == "PIRU" |
                                  merged_polygons$label == "PIGL",
                         "Picea", merged_polygons$label) #Put all Picea together
unique(merged_polygons$label) #check result

## To add labels related to species to measured traits ##----
# Import species and labels
label_sp <- read_sf("input_data/2022-POLYGON-DATA/tree_data.shp") %>%
  dplyr::select(label, specie, plant_id) %>%
  st_drop_geometry()
head(label_sp) #check result

# Put all Picea together (labels and species)
label_sp$label <- ifelse(label_sp$label == "PIMA" | label_sp$label == "PIRU",
                          "Picea", label_sp$label)
label_sp$specie <- ifelse(label_sp$specie == "Picea mariana" | label_sp$specie == "Picea rubens",
                         "Picea sp.", label_sp$specie)

# Check labels
length(unique(label_sp$label)) #16 classes
unique(label_sp$specie)

# Add labels to measured_traits
label_sp$plant_id <- as.numeric(label_sp$plant_id)
measured_traits <- dplyr::right_join(measured_traits, label_sp, by = c("plant_id", "specie")) %>%
  drop_na()

## Extract values and keep attributes of polygons----
estimated <- raster::extract(mosaic_pred, merged_polygons, exact = T) %>%
  dplyr::filter(fraction == 1)
head(estimated)

saveRDS(estimated, "input_data/estimated_traits_polygons.rds") #to read data estimated <- readRDS("input_data/estimated_traits_polygons.rds")

## Plot----
# SLA
# sd_estimated_SLA <- sd(estimated$SLA, na.rm = TRUE) #compare sd
# sd_measured_SLA <- sd(measured_traits$SLA, na.rm = TRUE)

SLA <- ggplot() +
  geom_density(data = measured_traits, aes(x = SLA, color = "Measured"), fill = NA, linewidth = 0.8) +
  geom_density(data = estimated, aes(x = SLA, color = "Estimated"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = expression("SLA (m"^2*" kg"^-1*")"), y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue")) #+
  #annotate("text", x = 0, y = Inf, hjust = 0, vjust = 1.5, 
           #label = paste("SD (Estimated):", round(sd_estimated_SLA, 2), "\nSD (Measured):", round(sd_measured_SLA, 2)),
           #size = 5, color = "black")
print(SLA)

# LMA
LMA <- ggplot() +
  geom_density(data = estimated, aes(x = LMA, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = LMA, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = expression("LMA (g m"^-2*")"), y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(LMA)

# LDMC
LDMC <- ggplot() +
  geom_density(data = estimated, aes(x = LDMC, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = LDMC, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = expression("LDMC (mg g"^-1*")"), y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(LDMC)

# EWT
EWT <- ggplot() +
  geom_density(data = estimated, aes(x = EWT, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = EWT, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = "EWT (mm)", y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(EWT)

# Hemicellulose
hemi <- ggplot() +
  geom_density(data = estimated, aes(x = hemicellulose, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = hemicellulose, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = "Hemicellulose (%)", y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(hemi)

# Cellulose
cell <- ggplot() +
  geom_density(data = estimated, aes(x = cellulose, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = cellulose, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = "Cellulose (%)", y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(cell)

# Lignin
lignin <- ggplot() +
  geom_density(data = estimated, aes(x = lignin, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = lignin, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = "Lignin (%)", y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(lignin)

# Chlorophyll a
chla <- ggplot() +
  geom_density(data = estimated, aes(x = chla, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = chla, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = expression("Chl" ~ italic(a) ~ "(mg g"^-1*")"), y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(chla)

# Chlorophyll b
chlb <- ggplot() +
  geom_density(data = estimated, aes(x = chlb, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = chlb, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = expression("Chl" ~ italic(b) ~ "(mg g"^-1*")"), y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(chlb)

# Carotenoids
car <- ggplot() +
  geom_density(data = estimated, aes(x = carotenoids, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = carotenoids, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = expression("Car (mg g"^-1*")"), y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(car)

# Carbon
carbon <- ggplot() +
  geom_density(data = estimated, aes(x = carbon, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = carbon, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = "Carbon (%)", y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(carbon)

# Nitrogen
nitrogen <- ggplot() +
  geom_density(data = estimated, aes(x = nitrogen, color = "Estimated"), fill = NA, linewidth = 0.8) +
  geom_density(data = measured_traits, aes(x = nitrogen, color = "Measured"), fill = NA, linewidth = 0.8) +
  labs(title = NULL, x = "Nitrogen (%)", y = NULL) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(size = 22),
    legend.position = c(0.85,0.85),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("Estimated" = "limegreen", "Measured" = "royalblue"))

print(nitrogen)

## Save plots----
all_plots <- egg::ggarrange(plots = list(carbon, car, cell, chla,
                                         chlb, EWT, hemi, LDMC,
                                         lignin, LMA, nitrogen, SLA),
                            ncol = 3, nrow = 4)
pdf("figures/trait_distribution_plots.pdf", width = 18, height = 18, onefile = F)
all_plots
dev.off()

png("figures/trait_distribution_plots.png", width = 18, height = 18, units = "in", res = 300)
all_plots
dev.off()