### Prepare PLSR data

## Load libraries----
library(tidyverse)
library(sf)
library(gridExtra)

## Set working directory----
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Import csv files (spectra + traits)---
spectra_data <- read.csv("output_spectra/smooth_spectra.csv", header = T, sep = ",", check.names = F)
all_traits <- read.csv("input_data/all_traits.csv", header = T, sep = ",")

## Merge spectra df and traits df----
# Remove missing trees
missing_plant_ids <- setdiff(unique(all_traits$plant_id), unique(spectra_data$plant_id)) #find missing trees
print(missing_plant_ids)
all_traits <- subset(all_traits, !(plant_id %in% missing_plant_ids)) #remove those plant_id for merging
all_traits$plant_id<-as.numeric(all_traits$plant_id) #change plant_id to numeric for merge

# Merge spectra df and trait df
plsr_data <- list(all_traits, spectra_data) %>% #put all data frames into list
  reduce(full_join, by = c("plant_id", "specie")) %>% #merge all data frames together
  relocate(long, lat, .after = specie) #relocate coordinates column

## Save the data (csv file)----
write.csv(x = plsr_data, file = "input_data/plsr_data.csv", row.names = F)

## Create a table for trait stats----
traits_only <- all_traits[, 5:16] #select only traits

min <- round(sapply(traits_only, min), 2) #calculate stats and round to 2 decimals
max <- round(sapply(traits_only, max), 2)
mean <- round(sapply(traits_only, mean), 2)
SD <- round(sapply(traits_only, sd), 2)
CV <- round(SD/mean*100, 2)

unit <- c("m2/kg", "g/m2", "mg/g", "mm", "%", "%", "%", "mg/g", "mg/g", "mg/g", "%", "%") #add units

trait_stats <- rownames_to_column(data.frame(min, max, mean, SD, CV), var = "Traits") %>%
  mutate(unit) %>% #add unit column
  relocate(unit, .after = "Traits") #relocate unit after trait column

write.csv(x = trait_stats, file = "input_data/trait_stats.csv", row.names = F)
