### Preprocess spectra

## Install all packages and load libraries----
install.packages("egg")

library(tidyverse)
library(hsdar)
library(spectrolab)
library(egg) #ggarrange function

# Set working directory
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Import data----
CASI_spectra <- read_csv("output_spectra/CASI_mean_large.csv", col_names = T)
SASI_spectra <- read_csv("output_spectra/SASI_mean_large.csv", col_names = T)

# Join tables
full_spectra <- right_join(CASI_spectra, SASI_spectra, by = c("plant_id", "long", "lat", "specie"))
full_spectra #check result
full_spectra <- na.omit(full_spectra) #remove 4 trees that have no spectra
write_csv(full_spectra, "output_spectra/full_spectra.csv", col_names = T) #save data

## Cleaning spectra----
# Convert to spectra object
colnames(full_spectra) <- gsub("Band_|nm", "", colnames(full_spectra)) #rename wavelength columns to remove "Band_" and "nm"

spec <- as_spectra(full_spectra, name_idx = 1, meta_idxs = c(2:4))
summary(spec) #check spectra object
plot_interactive(spec)

# Match the reflectance data
guess_splice_at(spec) #check the guess splice band (bound between sensors)

splice_bands = 998 #define the splice band
sensors_matched = match_sensors(x = spec, splice_at = splice_bands, fixed_sensor = 2,
                                interpolate_wvl = c(5, 1))
plot_interactive(sensors_matched)

# Plot to see the change
lwd = 0.5
cex = 0.8

plot(spec, main = "Reflectance", 
     lwd = lwd, cex.main = cex, cex.lab = cex, cex.axis = cex, ylim = (c(0,1)))

plot(sensors_matched, col = "red", add = T, 
     lwd = lwd, cex.main = cex, cex.lab = cex, cex.axis = cex)

pdf("figures/BEFORE_match_sensors.pdf", width = 10, height = 6)
plot(spec, main = "Reflectance BEFORE matching", 
     lwd = lwd, cex.main = cex, cex.lab = cex, cex.axis = cex)
dev.off()

pdf("figures/AFTER_match_sensors.pdf", width = 10, height = 6)
plot(sensors_matched, main = "Reflectance AFTER matching", col = "red", add = F, 
     lwd = lwd, cex.main = cex, cex.lab = cex, cex.axis = cex)
dev.off()

# Convert spectra object to a dataframe (wide format)
spec_df <- as.data.frame(sensors_matched, fix_names = "none", metadata = T) %>%
  dplyr::rename(plant_id = sample_name) %>%
  dplyr::mutate(plant_id = gsub("spec_", "", plant_id))

## Create object of class speclib ----
bands <-  spec_df %>%
  dplyr::select(5:224) #extract only the bands

wavelengths <- as.numeric(names(spec_df[, 5:224])) #define wavelengths
speclib <- speclib(spectra = t(bands), wavelength = wavelengths)
plot(speclib)

## Smooth spectra----
smooth <- noiseFiltering(speclib, method = 'sgolay', n = 3, p = 2)
hsdar::spectra(smooth)[hsdar::spectra(smooth) < 0] <- 0
plot(smooth)

# Add all important supplementary information
Meta <- spec_df %>% 
  dplyr::select(1:4) #select metadata
SI(smooth) <- Meta # add metadata to the supplementary information of the spectral library
colnames(SI(smooth)) <- colnames(Meta)
SI(smooth)

colors <- rainbow(length(unique(smooth@SI@SI_data$plant_id))) #set colors
pdf("figures/smooth_spectra.pdf", width = 10, height = 6)
plot(smooth, FUN = 1:166, main = "Canopy-level spectra of 166 trees (smoothed 
     with hsdar, n = 3, p = 2, and matched sensors at 998nm)", col = colors)
dev.off()

# Convert speclib object to dataframe
SI_smooth <- do.call(data.frame, smooth@SI@SI_data) #extract SI
spectra <-smooth@spectra@spectra_ma %>% data.frame() #extract reflectance data
colnames(spectra) <- as.integer(smooth@wavelength) # define colnames of Bands
smoothed <- cbind(SI_smooth, spectra) 

# Put all Picea together
smoothed$specie <- ifelse(smoothed$specie == "Picea mariana" | smoothed$specie == "Picea rubens",
                              "Picea sp.", smoothed$specie)
unique(smoothed$specie) #check results

# Save data
write_csv(smoothed, "output_spectra/smooth_spectra.csv", col_names = T) # to load data : smoothed <- read_csv("output_spectra/smooth_spectra.csv", col_names = T)

## Plot spectra by vegetation type----
# Change to long format
spectra_long <- pivot_longer(smoothed, cols = - c(1:4), names_to = "wavelength", values_to = "reflectance")
spectra_long$wavelength <- as.numeric(spectra_long$wavelength)

# Define broadleaf and conifer
spectra_long$type <- ifelse(spectra_long$specie %in% c("Acer pensylvanicum", "Acer rubrum", "Acer saccharum",
                                                        "Betula alleghaniensis", "Populus grandidentata",
                                                        "Fagus grandifolia", "Ostrya virginiana", "Quercus rubra",
                                                        "Populus tremuloides", "Betula papyrifera"), "broadleaf", "conifer")
# Group by vegetation type
type_spectra <- group_by(spectra_long, type, wavelength)
mean_spectra_type <- summarize(type_spectra, mean_reflectance = mean(reflectance, na.rm = T),
                               sd_reflectance = sd(reflectance, na.rm = T))

# Plot
plot <- ggplot(data = mean_spectra_type, aes(x = wavelength, y = mean_reflectance, group = type, col = type)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("broadleaf" = "limegreen", "conifer" = "royalblue")) +
  scale_fill_manual(values = c("broadleaf" = "limegreen", "conifer" = "royalblue")) +
  xlab("Wavelength (nm)") +
  ylab("Reflectance") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 15,),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1, 'cm')) +
  theme(legend.position = c(0.8, 0.8),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

plot

## Plot normalized spectra----
# Normalize spectra
smoothed <- as_spectra(smoothed, name_idx = 1, meta_idxs = c(2:4))
norm <- normalize(smoothed)
plot_interactive(norm)

# Change to dataframe and long format
norm_df <- as.data.frame(norm)
norm_dflong <- pivot_longer(norm_df, cols = - c(1:5), names_to = "wavelength", values_to = "reflectance")
norm_dflong$wavelength <- as.numeric(norm_dflong$wavelength)

# Define broadleaf and conifer
norm_dflong$type <- ifelse(norm_dflong$specie %in% c("Acer pensylvanicum", "Acer rubrum", "Acer saccharum",
                                                       "Betula alleghaniensis", "Populus grandidentata",
                                                       "Fagus grandifolia", "Ostrya virginiana", "Quercus rubra",
                                                       "Populus tremuloides", "Betula papyrifera"), "broadleaf", "conifer")
# Group by vegetation type
type_specBN <- group_by(norm_dflong, type, wavelength)
mean_specBN_type <- summarize(type_specBN, mean_reflectance = mean(reflectance, na.rm = T),
                               sd_reflectance = sd(reflectance, na.rm = T))
# Plot
plot_norm <- ggplot(data = mean_specBN_type, aes(x = wavelength, y = mean_reflectance, group = type, col = type)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("broadleaf" = "limegreen", "conifer" = "royalblue")) +
  scale_fill_manual(values = c("broadleaf" = "limegreen", "conifer" = "royalblue")) +
  xlab("Wavelength (nm)") +
  ylab("Brightness normalized reflectance") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))


plot_norm

# Save figure
veg_type <- egg::ggarrange(plots = list(plot, plot_norm), nrow = 2, labels = c("(A)", "(B)"))

pdf("figures/mean_spectra_vtype.pdf", width = 8, height = 9)
veg_type
dev.off()

png("figures/mean_spectra_vtype.png", width = 8, height = 9, units = "in", res = 300)
veg_type
dev.off()

## Plot spectra by specie----
# Add labels
species_labels <- data.frame(
  specie = c("Abies balsamea", "Acer pensylvanicum", "Acer rubrum", "Acer saccharum",
              "Betula alleghaniensis", "Betula papyrifera", "Fagus grandifolia",
              "Larix laricina", "Ostrya virginiana", "Picea sp.", "Pinus strobus",
              "Populus grandidentata", "Populus tremuloides", "Quercus rubra",
              "Thuja occidentalis", "Tsuga canadensis"),
  Label = c("ABBA", "ACPE", "ACRU", "ACSA", "BEAL", "BEPA", "FAGR", "LALA", "OSVI",
            "Picea", "PIST", "POGR", "POTR", "QURU", "THOC", "TSCA")
)

spectra_long <- spectra_long %>% 
  left_join(species_labels, by = "specie")

# Group by label
sp_spectra <- group_by(spectra_long, Label, wavelength)
mean_spectra_sp <- summarize(sp_spectra, mean_reflectance = mean(reflectance, na.rm = T))

# Change order of legend (max spectra to min spectra)
reflectance_754nm <- mean_spectra_sp %>% #calculate ref at 754nm
  dplyr::filter(wavelength == 754) %>%
  arrange(desc(mean_reflectance)) #sort lines of df in descending order based on reflectance at 754nm

label_order <- reflectance_754nm$Label #extract Label column
mean_spectra_sp$Label <- factor(mean_spectra_sp$Label, levels = label_order) #reorder 'Label' factor levels based on label_order

# Plot
colors <- c("#cab2d6", "#33a02c", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f","#b15928",
            "#FF6666", "#999933", "#66CCCC", "#1f78b4", "#a6cee3","#ffff00", "#ff7f00",
            "#6a3d9a","#CC66CC") 

plot_sp <- ggplot(data = mean_spectra_sp, aes(x = wavelength,
                                   y = mean_reflectance,
                                   group = Label, col = Label)) + 
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors) +
  xlab("Wavelength (nm)") +
  ylab("Reflectance") +
  theme_classic() +
  labs(col = "Label") +
theme(panel.border = element_rect(color = "black",
                                  fill = NA,
                                  linewidth = 1),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15))
plot_sp

# Save figure
pdf("figures/mean_spectra_specie_colors.pdf", width = 8, height = 5)
plot_sp
dev.off()

png("figures/mean_spectra_specie_colors.png", width = 8, height = 5, units = "in", res = 300)
plot_sp
dev.off()
