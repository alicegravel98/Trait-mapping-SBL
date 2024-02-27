### Verification of trait data

## Load libraries----
library(tidyverse)
library(sf)

## Set working directory----
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Import data----
trait_data <- read_csv("input_data/all_traits.csv", col_names = T)
tree_data <- read_sf("input_data/2022-POLYGON-DATA/tree_data.shp") %>%
  dplyr::select(plant_id, label) %>%
  st_drop_geometry()

# Put all Picea together
tree_data$label <- ifelse(tree_data$label == "PIMA" | tree_data$label == "PIRU",
                            "Picea", tree_data$label)
unique(tree_data$label) #check results

tree_data$plant_id <- as.numeric(tree_data$plant_id)
all_traits <- left_join(trait_data, tree_data, by = "plant_id")

## Outlier function----
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

## Verification structural and water-related traits----
# Specific leaf area
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(SLA), round(SLA, 2), "")) #create a column Outlier to label on plot

plot_SLA <- ggplot(data = all_traits,
                   aes(x = SLA, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier),
            position = position_jitter(width = -0.5),
            vjust = -0.8
  )
print(plot_SLA)

# Leaf mass per area
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(LMA), round(LMA, 2), ""))

plot_LMA <- ggplot(data = all_traits,
                   aes(x = LMA, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier),
            position = position_jitter(width = 0.5),
            vjust = -0.8
  )
print(plot_LMA)

# Leaf dry matter content
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(LDMC), round(LDMC, 2), ""))

plot_LDMC <- ggplot(data = all_traits,
                   aes(x = LDMC, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier),
            position = position_jitter(width = 0.5),
            vjust = -0.8
  )
print(plot_LDMC) #No ouliers.

# Equivalent water thickness
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(EWT), round(EWT, 4), ""))

plot_EWT <- ggplot(data = all_traits,
                    aes(x = EWT, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier), vjust = -0.8
  )
print(plot_EWT)

## Verification Cfractions----
# Hemicellulose
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(hemicellulose), round(hemicellulose, 2), ""))

plot_hemi <- ggplot(data = all_traits,
                   aes(x = hemicellulose, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier), vjust = -0.8
  )
print(plot_hemi) #They are in the range of CABO data.

# Cellulose
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(cellulose), round(cellulose, 2), ""))

plot_cell <- ggplot(data = all_traits,
                    aes(x = cellulose, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier), vjust = -0.8
  )
print(plot_cell) #They are in the range of CABO data.

# Lignin
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(lignin), round(lignin, 2), ""))

plot_lignin <- ggplot(data = all_traits,
                    aes(x = lignin, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier),
            position = position_jitter(width = -0.5),
            vjust = 0.8
  )
print(plot_lignin) #They are in the range of CABO data.

## Verification pigments----
# Chlorophyll a
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(chla), round(chla, 2), ""))

plot_chla <- ggplot(data = all_traits,
                      aes(x = chla, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier),
            position = position_jitter(width = -0.5),
            vjust = 0.8
  )
print(plot_chla)

# Chlorophyll b
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(chlb), round(chlb, 2), ""))

plot_chlb <- ggplot(data = all_traits,
                    aes(x = chlb, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier), vjust = 1
  )
print(plot_chlb)

# Carotenoids
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(carotenoid), round(carotenoid, 2), ""))

plot_carot <- ggplot(data = all_traits,
                    aes(x = carotenoid, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier),
            position = position_jitter(width = -0.5), vjust = 1
  )
print(plot_carot)

## Verification Carbone----
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(carbon), round(carbon, 2), ""))

plot_carbon <- ggplot(data = all_traits,
                     aes(x = carbon, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier), vjust = 1
  )
print(plot_carbon)

## Verification Nitrogen----
all_traits <- all_traits %>%
  group_by(label) %>%
  mutate(Outlier = ifelse(is_outlier(nitrogen), round(nitrogen, 2), ""))

plot_nitrogen <- ggplot(data = all_traits,
                      aes(x = nitrogen, y = label, fill = label)) +
  geom_boxplot(show.legend = F, outlier.colour = "red") +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  geom_text(aes(label = all_traits$Outlier),
            position = position_jitter(width = -0.5),
            vjust = -0.5
  )
print(plot_nitrogen)

