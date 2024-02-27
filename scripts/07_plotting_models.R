### Plotting PLSR models

## Load libraries----
library(tidyverse)
library(egg)

## Set working directory----
setwd("C:/Users/alice/Desktop/DossieR/Trait-mapping-SBL")

## Import data----
all.jack.df.list <-readRDS("output_plsr/all_jack_df_list.rds")
val_summary<-readRDS("output_plsr/val_summary.rds")

## Select R2 and %RMSE----
SLA_R2 <- val_summary %>%
  dplyr::filter(variable == "SLA") %>%
  dplyr::select(R2)
LMA_R2 <- val_summary %>%
  dplyr::filter(variable == "LMA") %>%
  dplyr::select(R2)
LDMC_R2 <- val_summary %>%
  dplyr::filter(variable == "LDMC") %>%
  dplyr::select(R2)
EWT_R2 <- val_summary %>%
  dplyr::filter(variable == "EWT") %>%
  dplyr::select(R2)
hemi_R2 <- val_summary %>%
  dplyr::filter(variable == "Hemi") %>%
  dplyr::select(R2)
cell_R2 <- val_summary %>%
  dplyr::filter(variable == "Cell") %>%
  dplyr::select(R2)
lignin_R2 <- val_summary %>%
  dplyr::filter(variable == "Lignin") %>%
  dplyr::select(R2)
chla_R2 <- val_summary %>%
  dplyr::filter(variable == "Chl a") %>%
  dplyr::select(R2)
chlb_R2 <- val_summary %>%
  dplyr::filter(variable == "Chl b") %>%
  dplyr::select(R2)
carot_R2 <- val_summary %>%
  dplyr::filter(variable == "Car") %>%
  dplyr::select(R2)
carbon_R2 <- val_summary %>%
  dplyr::filter(variable == "Carbon") %>%
  dplyr::select(R2)
nitrogen_R2 <- val_summary %>%
  dplyr::filter(variable == "Nitrogen") %>%
  dplyr::select(R2)
SLA_perRMSE <- val_summary %>%
  dplyr::filter(variable == "SLA") %>%
  dplyr::select(perRMSE)
LMA_perRMSE <- val_summary %>%
  dplyr::filter(variable == "LMA") %>%
  dplyr::select(perRMSE)
LDMC_perRMSE <- val_summary %>%
  dplyr::filter(variable == "LDMC") %>%
  dplyr::select(perRMSE)
EWT_perRMSE <- val_summary %>%
  dplyr::filter(variable == "EWT") %>%
  dplyr::select(perRMSE)
hemi_perRMSE <- val_summary %>%
  dplyr::filter(variable == "Hemi") %>%
  dplyr::select(perRMSE)
cell_perRMSE <- val_summary %>%
  dplyr::filter(variable == "Cell") %>%
  dplyr::select(perRMSE)
lignin_perRMSE <- val_summary %>%
  dplyr::filter(variable == "Lignin") %>%
  dplyr::select(perRMSE)
chla_perRMSE <- val_summary %>%
  dplyr::filter(variable == "Chl a") %>%
  dplyr::select(perRMSE)
chlb_perRMSE <- val_summary %>%
  dplyr::filter(variable == "Chl b") %>%
  dplyr::select(perRMSE)
carot_perRMSE <- val_summary %>%
  dplyr::filter(variable == "Car") %>%
  dplyr::select(perRMSE)
carbon_perRMSE <- val_summary %>%
  dplyr::filter(variable == "Carbon") %>%
  dplyr::select(perRMSE)
nitrogen_perRMSE <- val_summary %>%
  dplyr::filter(variable == "Nitrogen") %>%
  dplyr::select(perRMSE)

## Define axis limits----
SLA_all<-with(all.jack.df.list$SLA,c(pred.low,pred.high,Measured))
SLA_upper<-max(SLA_all,na.rm=T)+1
SLA_lower<-min(SLA_all,na.rm=T)-1    

LMA_all<-with(all.jack.df.list$LMA,c(pred.low,pred.high,Measured))
LMA_upper<-max(LMA_all,na.rm=T)+0.2
LMA_lower<-min(LMA_all,na.rm=T)-0.2  

LDMC_all<-with(all.jack.df.list$LDMC,c(pred.low,pred.high,Measured))
LDMC_upper<-max(LDMC_all,na.rm=T)+10
LDMC_lower<-min(LDMC_all,na.rm=T)-10

EWT_all<-with(all.jack.df.list$EWT,c(pred.low,pred.high,Measured))
EWT_upper<-max(EWT_all,na.rm=T)
EWT_lower<-min(EWT_all,na.rm=T)

hemi_all<-with(all.jack.df.list$Hemi,c(pred.low,pred.high,Measured))
hemi_upper<-max(hemi_all,na.rm=T)+2
hemi_lower<-min(hemi_all,na.rm=T)-2

cell_all<-with(all.jack.df.list$Cell,c(pred.low,pred.high,Measured))
cell_upper<-max(cell_all,na.rm=T)+1
cell_lower<-min(cell_all,na.rm=T)-1  

lignin_all<-with(all.jack.df.list$Lignin,c(pred.low,pred.high,Measured))
lignin_upper<-max(lignin_all,na.rm=T)+0.5
lignin_lower<-min(lignin_all,na.rm=T)-0.5

chla_all<-with(all.jack.df.list$"Chl a",c(pred.low,pred.high,Measured))
chla_upper<-max(chla_all,na.rm=T)
chla_lower<-min(chla_all,na.rm=T)  

chlb_all<-with(all.jack.df.list$"Chl b",c(pred.low,pred.high,Measured))
chlb_upper<-max(chlb_all,na.rm=T)
chlb_lower<-min(chlb_all,na.rm=T)

carot_all<-with(all.jack.df.list$Car,c(pred.low,pred.high,Measured))
carot_upper<-max(carot_all,na.rm=T)
carot_lower<-min(carot_all,na.rm=T)

carbon_all<-with(all.jack.df.list$Carbon,c(pred.low,pred.high,Measured))
carbon_upper<-max(as.numeric(carbon_all),na.rm=T)+1
carbon_lower<-min(as.numeric(carbon_all),na.rm=T)-1 

nitrogen_all<-with(all.jack.df.list$Nitrogen,c(pred.low,pred.high,Measured))
nitrogen_upper<-max(as.numeric(nitrogen_all),na.rm=T)
nitrogen_lower<-min(as.numeric(nitrogen_all),na.rm=T) 


## Plot models----
SLA.val.plot <- ggplot(all.jack.df.list$SLA, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(SLA_lower, SLA_upper), ylim = c(SLA_lower, SLA_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = expression("Measured SLA (m"^2*" kg"^-1*")"), x = expression("Predicted SLA (m"^2*" kg"^-1*")")) +
  annotate("text", x = SLA_lower, y = SLA_upper, label = bquote(R^2 == .(sprintf("%.2f", SLA_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = SLA_lower, y = SLA_upper - (SLA_upper - SLA_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", SLA_perRMSE)), hjust = 0, vjust = 1, size = 6)

SLA.val.plot

LMA.val.plot <- ggplot(all.jack.df.list$LMA, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(LMA_lower, LMA_upper), ylim = c(LMA_lower, LMA_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = expression("Measured LMA (g m"^-2*")"), x = expression("Predicted LMA (g m"^-2*")")) +
  annotate("text", x = LMA_lower, y = LMA_upper, label = bquote(R^2 == .(sprintf("%.2f", LMA_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = LMA_lower, y = LMA_upper - (LMA_upper - LMA_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", LMA_perRMSE)), hjust = 0, vjust = 1, size = 6)

LMA.val.plot

LDMC.val.plot <- ggplot(all.jack.df.list$LDMC, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(LDMC_lower, LDMC_upper), ylim = c(LDMC_lower, LDMC_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = expression("Measured LDMC (mg g"^-1*")"), x = expression("Predicted LDMC (mg g"^-1*")")) +
  annotate("text", x = LDMC_lower, y = LDMC_upper, label = bquote(R^2 == .(sprintf("%.2f", LDMC_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = LDMC_lower, y = LDMC_upper - (LDMC_upper - LDMC_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", LDMC_perRMSE)), hjust = 0, vjust = 1, size = 6)

LDMC.val.plot

EWT.val.plot <- ggplot(all.jack.df.list$EWT, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(EWT_lower, EWT_upper), ylim = c(EWT_lower, EWT_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = "Measured EWT (mm)", x = "Predicted EWT (mm)") +
  annotate("text", x = EWT_lower, y = EWT_upper, label = bquote(R^2 == .(sprintf("%.2f", EWT_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = EWT_lower, y = EWT_upper - (EWT_upper - EWT_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", EWT_perRMSE)), hjust = 0, vjust = 1, size = 6)

EWT.val.plot

hemi.val.plot <- ggplot(all.jack.df.list$Hemi, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(hemi_lower, hemi_upper), ylim = c(hemi_lower, hemi_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = "Measured Hemicellulose (%)", x = "Predicted Hemicellulose (%)") +
  annotate("text", x = hemi_lower, y = hemi_upper, label = bquote(R^2 == .(sprintf("%.2f", hemi_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = hemi_lower, y = hemi_upper - (hemi_upper - hemi_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", hemi_perRMSE)), hjust = 0, vjust = 1, size = 6)

hemi.val.plot

cell.val.plot <- ggplot(all.jack.df.list$Cell, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(cell_lower, cell_upper), ylim = c(cell_lower, cell_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = "Measured Cellulose (%)", x = "Predicted Cellulose (%)") +
  annotate("text", x = cell_lower, y = cell_upper, label = bquote(R^2 == .(sprintf("%.2f", cell_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = cell_lower, y = cell_upper - (cell_upper - cell_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", cell_perRMSE)), hjust = 0, vjust = 1, size = 6)

cell.val.plot

lignin.val.plot <- ggplot(all.jack.df.list$Lignin, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(lignin_lower, lignin_upper), ylim = c(lignin_lower, lignin_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = "Measured Lignin (%)", x = "Predicted Lignin (%)") +
  annotate("text", x = lignin_lower, y = lignin_upper, label = bquote(R^2 == .(sprintf("%.2f", lignin_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = lignin_lower, y = lignin_upper - (lignin_upper - lignin_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", lignin_perRMSE)), hjust = 0, vjust = 1, size = 6)

lignin.val.plot

chla.val.plot <- ggplot(all.jack.df.list$"Chl a", aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(chla_lower, chla_upper), ylim = c(chla_lower, chla_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = ("Measured Chl" ~ italic(a) ~ "(mg g"^-1*")"), x = expression("Predicted Chl" ~ italic(a) ~ "(mg g"^-1*")")) +
  annotate("text", x = chla_lower, y = chla_upper, label = bquote(R^2 == .(sprintf("%.2f", chla_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = chla_lower, y = chla_upper - (chla_upper - chla_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", chla_perRMSE)), hjust = 0, vjust = 1, size = 6)

chla.val.plot

chlb.val.plot <- ggplot(all.jack.df.list$"Chl b", aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(chlb_lower, chlb_upper), ylim = c(chlb_lower, chlb_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = ("Measured Chl" ~ italic(b) ~ "(mg g"^-1*")"), x = expression("Predicted Chl" ~ italic(b) ~ "(mg g"^-1*")")) +
  annotate("text", x = chlb_lower, y = chlb_upper, label = bquote(R^2 == .(sprintf("%.2f", chlb_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = chlb_lower, y = chlb_upper - (chlb_upper - chlb_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", chlb_perRMSE)), hjust = 0, vjust = 1, size = 6)

chlb.val.plot

carot.val.plot <- ggplot(all.jack.df.list$Car, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(carot_lower, carot_upper), ylim = c(carot_lower, carot_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = expression("Measured Carotenoids (mg g"^-1*")"), x = expression("Predicted Carotenoids (mg g"^-1*")")) +
  annotate("text", x = carot_lower, y = carot_upper, label = bquote(R^2 == .(sprintf("%.2f", carot_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = carot_lower, y = carot_upper - (carot_upper - carot_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", carot_perRMSE)), hjust = 0, vjust = 1, size = 6)

carot.val.plot

carbon.val.plot <- ggplot(all.jack.df.list$Carbon, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(carbon_lower, carbon_upper), ylim = c(carbon_lower, carbon_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = "Measured Carbon (%)", x = "Predicted Carbon (%)") +
  annotate("text", x = carbon_lower, y = carbon_upper, label = bquote(R^2 == .(sprintf("%.2f", carbon_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = carbon_lower, y = carbon_upper - (carbon_upper - carbon_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", carbon_perRMSE)), hjust = 0, vjust = 1, size = 6)

carbon.val.plot

nitrogen.val.plot <- ggplot(all.jack.df.list$Nitrogen, aes(y = Measured, x = pred.mean)) +
  geom_errorbarh(aes(y = Measured, xmin = pred.low, xmax = pred.high), color = "gray") +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "red1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 2) +
  coord_cartesian(xlim = c(nitrogen_lower, nitrogen_upper), ylim = c(nitrogen_lower, nitrogen_upper)) +
  theme(text = element_text(size = 22), plot.title = element_text(size = 22)) +
  labs(y = "Measured Nitrogen (%)", x = "Predicted Nitrogen (%)") +
  annotate("text", x = nitrogen_lower, y = nitrogen_upper, label = bquote(R^2 == .(sprintf("%.2f", nitrogen_R2))), hjust = 0, vjust = 1, size = 6) +
  annotate("text", x = nitrogen_lower, y = nitrogen_upper - (nitrogen_upper - nitrogen_lower) * 0.09, label = paste0("%RMSE = ", sprintf("%.2f", nitrogen_perRMSE)), hjust = 0, vjust = 1, size = 6)

nitrogen.val.plot

## Save plots----
all_plots <- egg::ggarrange(plots = list(LMA.val.plot, EWT.val.plot, SLA.val.plot,
                            carot.val.plot, chla.val.plot, nitrogen.val.plot,
                            chlb.val.plot, cell.val.plot, lignin.val.plot,
                            carbon.val.plot, LDMC.val.plot, hemi.val.plot), ncol = 3, nrow = 4)

pdf("output_plsr/validation_plots.pdf", width = 18, height = 18, onefile = F)
all_plots
dev.off()

png("output_plsr/validation_plots.png", width = 18, height = 18, units = "in", res = 300)
all_plots
dev.off()
