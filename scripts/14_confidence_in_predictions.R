### Confidence in trait predictions

## Load libraries----
library(tidyverse)
library(viridis) #for colors in histograms
library(stringr) #str_to_title()
library(ggpubr) #annotate_figure()
library(grid) #textGrob()

## Set working directory----
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Import data----
# Trait predictions
df_traits01 <- readRDS("inference/flg01/df_traits01.rds")
df_traits02 <- readRDS("inference/flg02/df_traits02.rds")
df_traits03 <- readRDS("inference/flg03/df_traits03.rds")
df_traits04 <- readRDS("inference/flg04/df_traits04.rds")
df_traits05 <- readRDS("inference/flg05/df_traits05.rds")
df_traits06 <- readRDS("inference/flg06/df_traits06.rds")
df_traits07 <- readRDS("inference/flg07/df_traits07.rds")

# Relative uncertainties (CV)
df_CV01 <- readRDS("inference/flg01/df_CV01.rds")
colnames(df_CV01)[1:12] <- paste("CV", colnames(df_CV01)[1:12], sep = ".") #add CV to column names
df_CV02 <- readRDS("inference/flg02/df_CV02.rds")
colnames(df_CV02)[1:12] <- paste("CV", colnames(df_CV02)[1:12], sep = ".")
df_CV03 <- readRDS("inference/flg03/df_CV03.rds")
colnames(df_CV03)[1:12] <- paste("CV", colnames(df_CV03)[1:12], sep = ".") 
df_CV04 <- readRDS("inference/flg04/df_CV04.rds")
colnames(df_CV04)[1:12] <- paste("CV", colnames(df_CV04)[1:12], sep = ".") 
df_CV05 <- readRDS("inference/flg05/df_CV05.rds")
colnames(df_CV05)[1:12] <- paste("CV", colnames(df_CV05)[1:12], sep = ".") 
df_CV06 <- readRDS("inference/flg06/df_CV06.rds")
colnames(df_CV06)[1:12] <- paste("CV", colnames(df_CV06)[1:12], sep = ".") 
df_CV07 <- readRDS("inference/flg07/df_CV07.rds")
colnames(df_CV07)[1:12] <- paste("CV", colnames(df_CV07)[1:12], sep = ".") 

## Bind rows (all flightlines)----
df_traits <- bind_rows(df_traits01, df_traits02, df_traits03, df_traits04, df_traits05, df_traits06, df_traits07)
df_CV <- bind_rows(df_CV01, df_CV02, df_CV03, df_CV04, df_CV05, df_CV06, df_CV07)

rm(df_traits01, df_traits02, df_traits03, df_traits04, df_traits05, df_traits06, df_traits07)
rm(df_CV01, df_CV02, df_CV03, df_CV04, df_CV05, df_CV06, df_CV07)

## Join all----
join <- right_join(df_traits, df_CV, by = c("x", "y"))

## Create dataframe for each trait (with prediction and CV)----
traits <- c("SLA", "LMA", "LDMC", "EWT", "hemi", "cell", "lignin", "chla", "chlb", "carot", "carbon", "nitrogen")

list_dataframes <- list()
for(trait in traits) {
  df_trait <- join[, c(trait, paste0("CV.", trait))]
  list_dataframes[[trait]] <- df_trait
}

df_SLA <- list_dataframes[["SLA"]] %>%
  filter_all(all_vars(. >= 0))
df_LMA <- list_dataframes[["LMA"]] %>%
  filter_all(all_vars(. >= 0))
df_LDMC <- list_dataframes[["LDMC"]] %>%
  filter_all(all_vars(. >= 0))
df_EWT <- list_dataframes[["EWT"]] %>%
  filter_all(all_vars(. >= 0))
df_hemi <- list_dataframes[["hemi"]] %>%
  filter_all(all_vars(. >= 0))
df_cell <- list_dataframes[["cell"]] %>%
  filter_all(all_vars(. >= 0))
df_lignin <- list_dataframes[["lignin"]] %>%
  filter_all(all_vars(. >= 0))
df_chla <- list_dataframes[["chla"]] %>%
  filter_all(all_vars(. >= 0))
df_chlb <- list_dataframes[["chlb"]] %>%
  filter_all(all_vars(. >= 0))
df_carot <- list_dataframes[["carot"]] %>%
  filter_all(all_vars(. >= 0))
df_carbon <- list_dataframes[["carbon"]] %>%
  filter_all(all_vars(. >= 0))
df_nitrogen <- list_dataframes[["nitrogen"]] %>%
  filter_all(all_vars(. >= 0))

## Define range trait----
# [traitmin - | traitmin | ✕ 25%, traitmax - | traitmax | ✕ 25%]
# Create a list of all traits
traits <- list(
  list(trait = "SLA", min = 2.31, max = 20.05),
  list(trait = "LMA", min = 49.87, max = 433.44),
  list(trait = "LDMC", min = 327.20, max = 546.68),
  list(trait = "EWT", min = 0.06, max = 0.35),
  list(trait = "hemi", min = 6.15, max = 28.33),
  list(trait = "cell", min = 6.19, max = 18.73),
  list(trait = "lignin", min = 4.51, max = 19.37),
  list(trait = "chla", min = 0.38, max = 3.58),
  list(trait = "chlb", min = 0.12, max = 1),
  list(trait = "carot", min = 0.1, max = 0.83),
  list(trait = "carbon", min = 45.58, max = 51.19),
  list(trait = "nitrogen", min = 0.69, max = 3.06)
)

## Calculate high confidence pixels----
## What % of pixels are within the range of measured traits?
# Create an empty dataframe to stock values
result_range <- data.frame(trait = character(), perc_in_range = numeric())

# Loop for all traits
for (trait_info in traits) {
  trait_df <- get(paste0("df_", trait_info$trait))  
  in_range <- which(trait_df[[trait_info$trait]] >= trait_info$min & 
                      trait_df[[trait_info$trait]] <= trait_info$max)
  percentage_in_range <- 100 * length(in_range) / nrow(trait_df)
  result_range <- rbind(result_range, data.frame(trait = trait_info$trait, perc_in_range = percentage_in_range))
}

print(result_range)

## What % of pixels (within the range of measured traits) have a CV lower than 25%? = high confidence pixels
# Create an empty dataframe to stock values
result_df <- data.frame(trait = character(), perc_confidence = numeric())

# Loop for all traits
for (trait_info in traits) {
  trait_df <- get(paste0("df_", trait_info$trait))
  in_range <- which(trait_df[[trait_info$trait]] >= trait_info$min & 
                      trait_df[[trait_info$trait]] <= trait_info$max)
  low_cv <- which(trait_df[[paste0("CV.", trait_info$trait)]] < 0.25)
  confident_indices <- intersect(in_range, low_cv)
  percentage_confidence <- 100 * length(confident_indices) / length(in_range)
  result_df <- rbind(result_df, data.frame(trait = trait_info$trait, perc_confidence = percentage_confidence))
}

print(result_df)

# Change trait names for uniformity
result_range$trait <- ifelse(result_range$trait %in% c("SLA", "LMA", "EWT", "LDMC"),
                          result_range$trait, 
                          str_to_title(result_range$trait)) #Capitalize first letter for some traits
result_range$trait <- sub("Carot", "Carotenoids", result_range$trait) #change carotenoids'name
result_range$trait <- sub("Chla", "Chl a", result_range$trait)
result_range$trait <- sub("Chlb", "Chl b", result_range$trait)

print(result_range)

result_df$trait <- ifelse(result_df$trait %in% c("SLA", "LMA", "EWT", "LDMC"),
                            result_df$trait, 
                            str_to_title(result_df$trait)) #Capitalize first letter for some traits
result_df$trait <- sub("Carot", "Carotenoids", result_df$trait) #change carotenoids'name
result_df$trait <- sub("Chla", "Chl a", result_df$trait)
result_df$trait <- sub("Chlb", "Chl b", result_df$trait)

print(result_df)

## Histogram----
plot_range <- ggplot(result_range, aes(x = trait, y = perc_in_range, fill = trait)) +
  geom_bar(stat = "identity", fill = "grey55") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"))+
  guides(fill = "none") + #remove legend
  labs(y = "Pixels within the range of measured traits (%)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110)) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "red", linewidth = 1)
print(plot_range)

plot_conf <- ggplot(result_df, aes(x = trait, y = perc_confidence, fill = trait)) +
  geom_bar(stat = "identity", fill = "grey55") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"))+
  guides(fill = "none") + #remove legend
  labs(y = "High confidence pixels (%)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110)) +
  scale_x_discrete(labels = c("Carbon", "Carotenoids", "Cell", expression(Chl ~ italic(a)), expression(Chl ~ italic(b)), "EWT", "Hemi", "LDMC", "Lignin", "LMA", "Nitrogen", "SLA")) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "red", linewidth = 1)
print(plot_conf)

# Save plot
pdf("figures/barplots_confidence_all.pdf", width = 7, height = 4, onefile = F)
plot_conf
dev.off()

png("figures/barplots_confidence_all.png", width = 7, height = 4, units = "in", res = 300)
plot_conf
dev.off()
