### Model interpretation : important wavelengths

## Load libraries----
install.packages("splus2R")
library(tidyverse)
library(reshape2)
library(splus2R) #function peaks()

# Set working directory
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Import VIP----
VIP1.df <- readRDS("output_plsr/VIP1_df.rds")
VIP2.df <- readRDS("output_plsr/VIP2_df.rds")
VIP3.df <- readRDS("output_plsr/VIP3_df.rds")
VIP4.df <- readRDS("output_plsr/VIP4_df.rds")

VIP1.long <- melt(VIP1.df, id.vars = "wavelength")
VIP2.long <- melt(VIP2.df, id.vars = "wavelength")
VIP3.long <- melt(VIP3.df, id.vars = "wavelength")
VIP4.long <- melt(VIP4.df, id.vars = "wavelength")

## Find local maxima VIP1 : EWT, LDMC, LMA, SLA.----
#EWT
EWT <- VIP1.df %>%
  dplyr::select(EWT, wavelength)

peaks_EWT <- peaks(EWT, span = 11) #return booleans for each columns

peaks_EWT <- peaks_EWT %>%
  dplyr::select(EWT) %>% #select only EWT
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(EWT$wavelength) #create a vector for wavelengths to fill the column in df
peaks_EWT$wavelength <- wavelengths #add wavelength to columns
max_EWT <- peaks_EWT %>%
  dplyr::filter(EWT == TRUE) %>%
  dplyr::pull(wavelength)
print(max_EWT)

ggplot(EWT, aes(x = wavelength, y = EWT)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_EWT, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_EWT <- max_EWT[-2]

#LDMC
LDMC <- VIP1.df %>%
  dplyr::select(LDMC, wavelength)

peaks_LDMC <- peaks(LDMC, span = 9, endbehavior = 1) #endbehavior to consider end as a peak

peaks_LDMC <- peaks_LDMC %>%
  dplyr::select(LDMC) %>% #select only LDMC
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(LDMC$wavelength) #create a vector for wavelengths to fill the column in df
peaks_LDMC$wavelength <- wavelengths #add wavelength to columns
max_LDMC <- peaks_LDMC %>%
  dplyr::filter(LDMC == TRUE) %>%
  dplyr::pull(wavelength)
print(max_LDMC)

ggplot(LDMC, aes(x = wavelength, y = LDMC)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_LDMC, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_LDMC <- max_LDMC[-c(1:7, 11)]

#LMA
LMA <- VIP1.df %>%
  dplyr::select(LMA, wavelength)

peaks_LMA <- peaks(LMA, span = 11) #return booleans for each columns

peaks_LMA <- peaks_LMA %>%
  dplyr::select(LMA) %>% #select only LMA
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(LMA$wavelength) #create a vector for wavelengths to fill the column in df
peaks_LMA$wavelength <- wavelengths #add wavelength to columns
max_LMA <- peaks_LMA %>%
  dplyr::filter(LMA == TRUE) %>%
  dplyr::pull(wavelength)
print(max_LMA)

ggplot(LMA, aes(x = wavelength, y = LMA)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_LMA, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_LMA <- max_LMA[-c(2, 3, 6, 7)]

#SLA
SLA <- VIP1.df %>%
  dplyr::select(SLA, wavelength)

peaks_SLA <- peaks(SLA, span = 11) #return booleans for each columns

peaks_SLA <- peaks_SLA %>%
  dplyr::select(SLA) %>% #select only SLA
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(SLA$wavelength) #create a vector for wavelengths to fill the column in df
peaks_SLA$wavelength <- wavelengths #add wavelength to columns
max_SLA <- peaks_SLA %>%
  dplyr::filter(SLA == TRUE) %>%
  dplyr::pull(wavelength)
print(max_SLA)

ggplot(SLA, aes(x = wavelength, y = SLA)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_SLA, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_SLA <- max_SLA[-c(2, 3)]

## Plot to visualize VIP1----
focal_palette <- c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd")

VIP1_plot <- ggplot(VIP1.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_EWT, linetype = "dashed", color = "#1f77b4", linewidth = 0.8) +
  geom_vline(xintercept = max_LDMC, linetype = "dashed", color = "#2ca02c", linewidth = 0.8) +
  geom_vline(xintercept = max_LMA, linetype = "dashed", color = "#ff7f0e", linewidth = 0.8) +
  geom_vline(xintercept = max_SLA, linetype = "dashed", color = "#d62728", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4), xlim = c(400, 2500)) +
  scale_color_manual(values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

print(VIP1_plot)

VIP1_annotated <- VIP1_plot +
  geom_text(aes(x = 567, y = 2.8, label = "567"), color = "black", size = 5.5, angle = 90, vjust = 0) + 
  geom_text(aes(x = 749, y = 2.8, label = "749"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1077, y = 2.8, label = "1077"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1212, y = 2.8, label = "1212"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1317, y = 2.8, label = "1317"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1422, y = 2.8, label = "1422"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1587, y = 2.8, label = "1587"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1722, y = 2.8, label = "1722"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1827, y = 2.8, label = "1827"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2067, y = 2.4, label = "2067"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2082, y = 3.5, label = "2082"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2382, y = 2.8, label = "2382"), color = "black", size = 5.5, angle = 90, vjust = 0)

print(VIP1_annotated)

## Find local maxima VIP2 : Cellulose, hemicellulose, lignin.----
#Cellulose
cell <- VIP2.df %>%
  dplyr::select(cell, wavelength)

peaks_cell <- peaks(cell, span = 11) #return booleans for each columns

peaks_cell <- peaks_cell %>%
  dplyr::select(cell) %>% #select only cell
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(cell$wavelength) #create a vector for wavelengths to fill the column in df
peaks_cell$wavelength <- wavelengths #add wavelength to columns
max_cell <- peaks_cell %>%
  dplyr::filter(cell == TRUE) %>%
  dplyr::pull(wavelength)
print(max_cell)

ggplot(cell, aes(x = wavelength, y = cell)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_cell, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_cell <- max_cell[-c(1,2, 6:10)]

#Hemicellulose
hemi <- VIP2.df %>%
  dplyr::select(hemi, wavelength)

peaks_hemi <- peaks(hemi, span = 11) #return booleans for each columns

peaks_hemi <- peaks_hemi %>%
  dplyr::select(hemi) %>% #select only hemi
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(hemi$wavelength) #create a vector for wavelengths to fill the column in df
peaks_hemi$wavelength <- wavelengths #add wavelength to columns
max_hemi <- peaks_hemi %>%
  dplyr::filter(hemi == TRUE) %>%
  dplyr::pull(wavelength)
print(max_hemi)

ggplot(hemi, aes(x = wavelength, y = hemi)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_hemi, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_hemi <- max_hemi[-c(3, 4, 9)]

#Lignin
lignin <- VIP2.df %>%
  dplyr::select(lignin, wavelength)

peaks_lignin <- peaks(lignin, span = 11) #return booleans for each columns

peaks_lignin <- peaks_lignin %>%
  dplyr::select(lignin) %>% #select only lignin
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(lignin$wavelength) #create a vector for wavelengths to fill the column in df
peaks_lignin$wavelength <- wavelengths #add wavelength to columns
max_lignin <- peaks_lignin %>%
  dplyr::filter(lignin == TRUE) %>%
  dplyr::pull(wavelength)
print(max_lignin)

ggplot(lignin, aes(x = wavelength, y = lignin)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_lignin, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_lignin <- max_lignin[-c(1:3, 9, 11)]

## Plot to visualize VIP2----
focal_palette <- c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd")

VIP2_plot <- ggplot(VIP2.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_cell, linetype = "dashed", color = "#1f77b4", linewidth = 0.8) +
  geom_vline(xintercept = max_hemi, linetype = "dashed", color = "#2ca02c", linewidth = 0.8) +
  geom_vline(xintercept = max_lignin, linetype = "dashed", color = "#ff7f0e", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4), xlim = c(400, 2500)) +
  scale_color_manual(labels = c("Cell", "Hemi", "Lignin"),
                     values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

print(VIP2_plot)

VIP2_annotated <- VIP2_plot +
  geom_text(aes(x = 563, y = 2.8, label = "563"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 639, y = 2.8, label = "639"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 749, y = 2.8, label = "749"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 893, y = 2.8, label = "893"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1032, y = 2.8, label = "1032"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1077, y = 2.8, label = "1077"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1287, y = 3.5, label = "1287"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1317, y = 2.4, label = "1317"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1587, y = 2.8, label = "1587"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1677, y = 2.8, label = "1677"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1737, y = 2.8, label = "1737"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2217, y = 2.8, label = "2217"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2322, y = 2.8, label = "2322"), color = "black", size = 5.5, angle = 90, vjust = 0)

print(VIP2_annotated)

## Find local maxima VIP3 : Carotenoids, chla, chlb.----
#Carotenoids
car <- VIP3.df %>%
  dplyr::select(car, wavelength)

peaks_car <- peaks(car, span = 11) #return booleans for each columns

peaks_car <- peaks_car %>%
  dplyr::select(car) %>% #select only carot
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(car$wavelength) #create a vector for wavelengths to fill the column in df
peaks_car$wavelength <- wavelengths #add wavelength to columns
max_car <- peaks_car %>%
  dplyr::filter(car == TRUE) %>%
  dplyr::pull(wavelength)
print(max_car)

ggplot(car, aes(x = wavelength, y = car)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_car, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_car <- max_car[-c(2, 6, 7)]

#Chla
chla <- VIP3.df %>%
  dplyr::select(chla, wavelength)

peaks_chla <- peaks(chla, span = 11) #return booleans for each columns

peaks_chla <- peaks_chla %>%
  dplyr::select(chla) %>% #select only chla
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(chla$wavelength) #create a vector for wavelengths to fill the column in df
peaks_chla$wavelength <- wavelengths #add wavelength to columns
max_chla <- peaks_chla %>%
  dplyr::filter(chla == TRUE) %>%
  dplyr::pull(wavelength)
print(max_chla)

ggplot(chla, aes(x = wavelength, y = chla)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_chla, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

#Chlb
chlb <- VIP3.df %>%
  dplyr::select(chlb, wavelength)

peaks_chlb <- peaks(chlb, span = 11) #return booleans for each columns

peaks_chlb <- peaks_chlb %>%
  dplyr::select(chlb) %>% #select only chlb
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(chlb$wavelength) #create a vector for wavelengths to fill the column in df
peaks_chlb$wavelength <- wavelengths #add wavelength to columns
max_chlb <- peaks_chlb %>%
  dplyr::filter(chlb == TRUE) %>%
  dplyr::pull(wavelength)
print(max_chlb)

ggplot(chlb, aes(x = wavelength, y = chlb)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_chlb, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

## Plot to visualize VIP3----
focal_palette <- c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd")

VIP3_plot <- ggplot(VIP3.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_car, linetype = "dashed", color = "#1f77b4", linewidth = 0.8) +
  geom_vline(xintercept = max_chla, linetype = "dashed", color = "#2ca02c", linewidth = 0.8) +
  geom_vline(xintercept = max_chlb, linetype = "dashed", color = "#ff7f0e", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4), xlim = c(400, 2500)) +
  scale_color_manual(labels = expression("Car", Chl ~ italic(a), Chl ~ italic(b)),
                     values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

print(VIP3_plot)

VIP3_annotated <- VIP3_plot +
  geom_text(aes(x = 558, y = 2.8, label = "558"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 740, y = 2.8, label = "740"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1077, y = 2.8, label = "1077"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1317, y = 2.8, label = "1317"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1677, y = 2.8, label = "1677"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2067, y = 2.8, label = "2067"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2202, y = 2.8, label = "2202"), color = "black", size = 5.5, angle = 90, vjust = 0)

print(VIP3_annotated)

## Find local maxima for VIP4 : Carbon, nitrogen.----
#Carbon
carbon <- VIP4.df %>%
  dplyr::select(carbon, wavelength)

peaks_carbon <- peaks(carbon, span = 11, endbehavior = 1) #return booleans for each columns

peaks_carbon <- peaks_carbon %>%
  dplyr::select(carbon) %>% #select only carbon
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(carbon$wavelength) #create a vector for wavelengths to fill the column in df
peaks_carbon$wavelength <- wavelengths #add wavelength to columns
max_carbon <- peaks_carbon %>%
  dplyr::filter(carbon == TRUE) %>%
  dplyr::pull(wavelength)
print(max_carbon)

ggplot(carbon, aes(x = wavelength, y = carbon)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_carbon, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_carbon <- max_carbon[-c(1:2, 4:6)]

#Nitrogen
nitrogen <- VIP4.df %>%
  dplyr::select(nitrogen, wavelength)

peaks_nitrogen <- peaks(nitrogen, span = 11, endbehavior = 1) #return booleans for each columns

peaks_nitrogen <- peaks_nitrogen %>%
  dplyr::select(nitrogen) %>% #select only nitrogen
  mutate(wavelength = as.numeric(0)) #create a column wavelength with zeros

wavelengths <- as.numeric(nitrogen$wavelength) #create a vector for wavelengths to fill the column in df
peaks_nitrogen$wavelength <- wavelengths #add wavelength to columns
max_nitrogen <- peaks_nitrogen %>%
  dplyr::filter(nitrogen == TRUE) %>%
  dplyr::pull(wavelength)
print(max_nitrogen)

ggplot(nitrogen, aes(x = wavelength, y = nitrogen)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_nitrogen, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength") +
  coord_cartesian(ylim = c(0, 4)) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

max_nitrogen <- max_nitrogen[-c(1, 3, 9, 10, 11)]

## Plot to visualize VIP4----
focal_palette <- c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd")

VIP4_plot <- ggplot(VIP4.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = max_carbon, linetype = "dashed", color = "#1f77b4", linewidth = 0.8) +
  geom_vline(xintercept = max_nitrogen, linetype = "dashed", color = "#2ca02c", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength (nm)", color = "Trait") +
  coord_cartesian(ylim = c(0, 4), xlim = c(400, 2500)) +
  scale_color_manual(labels = c("Carbon", "Nitrogen"),
                     values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

print(VIP4_plot)

VIP4_annotated <- VIP4_plot +
  geom_text(aes(x = 558, y = 2.8, label = "558"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 740, y = 2.8, label = "740"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1077, y = 2.8, label = "1077"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1212, y = 2.8, label = "1212"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1287, y = 2.4, label = "1287"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1317, y = 3.5, label = "1317"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1422, y = 2.8, label = "1422"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1617, y = 2.8, label = "1617"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 1722, y = 2.8, label = "1722"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2067, y = 2.8, label = "2067"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2217, y = 2.8, label = "2217"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2367, y = 2.4, label = "2367"), color = "black", size = 5.5, angle = 90, vjust = 0) +
  geom_text(aes(x = 2382, y = 3.5, label = "2382"), color = "black", size = 5.5, angle = 90, vjust = 0)

print(VIP4_annotated)

## Arrange all plots together----
VIP_all <- egg::ggarrange(plots = list(VIP1_annotated, VIP2_annotated, VIP3_annotated, VIP4_annotated),
                          ncol = 1, nrow = 4, labels = c("(A)", "(B)", "(C)", "(D)"))

# Save plots
pdf("output_plsr/VIPs_all_max.pdf", height = 15, width = 10, onefile = F)
VIP_all
dev.off()

png("output_plsr/VIPs_all_max.png", width = 12, height = 10, units = "in", res = 300)
VIP_all
dev.off()

## Water absorption bands----
water1 <- ggplot(VIP1.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = 1332, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1467, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1782, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1977, linetype = "dashed", color = "blue", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4)) +
  scale_color_manual(values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")
  
water2 <- ggplot(VIP2.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = 1332, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1467, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1782, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1977, linetype = "dashed", color = "blue", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4)) +
  scale_color_manual(labels = c("Cell", "Hemi", "Lignin"),
                     values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")
  
water3 <- ggplot(VIP3.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = 1332, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1467, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1782, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1977, linetype = "dashed", color = "blue", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4)) +
  scale_color_manual(labels = expression("Car", Chl ~ italic(a), Chl ~ italic(b)),
                     values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

water4 <- ggplot(VIP4.long, aes(x = wavelength, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_bw() +
  geom_vline(xintercept = 1332, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1467, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1782, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = 1977, linetype = "dashed", color = "blue", linewidth = 0.8) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  labs(y = "VIP", x = "Wavelength", color = "Trait") +
  coord_cartesian(ylim = c(0, 4)) +
  scale_color_manual(labels = c("Carbon", "Nitrogen"),
                     values = focal_palette) +
  geom_abline(slope = 0, intercept = 0.8, linetype = "dashed", linewidth = 1, color = "grey55")

# Arrange all plots together
VIP_water <- egg::ggarrange(plots = list(water1, water2, water3, water4),
                          ncol = 1, nrow = 4)

# Save plots
pdf("output_plsr/water_abs_bands.pdf", height = 15, width = 10, onefile = F)
VIP_water
dev.off()

png("output_plsr/water_abs_bands.png", width = 12, height = 10, units = "in", res = 300)
VIP_water
dev.off()