################################################################################
# This Script generates the Variogram plots of Field-Y for the paper. The      #
# figure will have three panels corresponding to three unit area sizes:        #
# 1.2, 8.6 and 60 m2                                                           #
################################################################################

################################################################################
# Load required packages                                                       #  
################################################################################
library(here)
library(readxl)
library(tidyverse)
library(sp)
library(spatstat)
library(ggpubr)
library(rgdal)
library(maptools)
library(sf)
library(gstat)
library(geoR)
library(GWmodel)
library(raster)
library(rgeos)
library(automap)
library(GGally)
library(viridis)
library(RColorBrewer)
library(ggtext)
library(rlang)
library(cowplot)
library(patchwork)
################################################################################
# Source utility functions                                                     #
################################################################################
source(here::here("RCode", "Utils_ForPaper.R"))
################################################################################
################################################################################
#                                                                              #        
#                              Field-Y                                         #
#                         -----------------                                    #
#            Semi-variances and Variograms of fungicide resistance             #
# frequency as functions of distance from differing unit-area-sizes of the     #
# nested transect sampling design.                                             #
# Plots included: RAKED (r) strategy for (a) 1.2 m2                            #
#                                        (b) 8.6 m2                            #
#                                        (c) 60 m2                             #
#                                                                              #
################################################################################
################################################################################
# Read in the Spatial Dependency data --- sampled using a nested scheme        #
################################################################################
spat_dep_data <- readRDS(here::here("DerivedData", "spat_dep_data.rds"))
spat_dep_sf_projected <- readRDS(here::here("DerivedData", 
                                            "spat_dep_sf_projected.rds"))
##################################################################################
# Add the Projected Coordinates to the Spatial Dependency Dataset                #
##################################################################################
spat_dep_data2 <-  spat_dep_data %>% 
  bind_cols(st_coordinates(spat_dep_sf_projected))
# Define the plotting attributes
point_size <- 4 
axis_text_size <- 16
axis_text_face <- "bold"
axis_title_size <- 16
axis_title_face <- "bold"
################################################################################
#                                                                              #  
# a. Paddock = Y; variable = dPCR; Sampling design = "b" (rake and 1.2 m2)     #
#                                                                              #  
################################################################################
# Filter out the necessary data
padd_Y_1.2sqm <- spat_dep_data2 %>% 
  dplyr::filter(PaddockFac == "Y", 
                NestedSubsample == "b")
# Compute semi-variances
padd_Y_1.2sqm_vgram_cloud <- compute_svgm(padd_Y_1.2sqm, Raw_dPCR)
################################################################################
# Create the semi-variance point cloud                                         #
################################################################################
padd_Y_1.2sqm_vgram_cloud_plt <- padd_Y_1.2sqm_vgram_cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha = 0.1, size= point_size, 
             fill="black", shape = 21, color="black") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = axis_text_size, face = axis_text_face),
        axis.text.y = element_text(size = axis_text_size, face = axis_text_face),
        axis.title = element_text(size = axis_title_size, face = axis_title_face)) +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2)) +
  scale_x_continuous(breaks = c(0, 100, 250, 400, 500), limits = c(0, 500))

padd_Y_1.2sqm_vgram_cloud_plt2 <- padd_Y_1.2sqm_vgram_cloud_plt +
  ggtitle(bquote(bold(atop("(a) Field-SD-b, 1.2"~m^ 2~"(r)")))) +
  theme(plot.title = element_text(size = 15, face = "bold"))
################################################################################
# Construct the variogram                                                      #
################################################################################
################################################################################
padd_Y_1.2sqm_vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data2, 
                                                           var=Raw_dPCR, 
                                                           paddock="Y", 
                                                           samp_design="b",
                                        dist_brks=c(0, 75, 150, 450, 550),
                                                           boot_sim=1000, 
                                                           boot_seed=100)
################################################################################
# Save the data                                                                #
################################################################################
saveRDS(padd_Y_1.2sqm_vgram_data,
        here::here("DerivedData", "padd_Y_1.2sqm_vgram_data.rds"))
################################################################################
padd_Y_1.2sqm_lowlim_data <- padd_Y_1.2sqm_vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  tidyr::pivot_longer(cols = -c(Dist, SemiVariance), 
                      names_to="LeftQuantile", 
                      values_to="VarLowLim")

padd_Y_1.2sqm_upplim_data <- padd_Y_1.2sqm_vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  tidyr::pivot_longer(cols = -Dist, 
                      names_to="RightQuantile", 
                      values_to="VarUppLim") %>%
  dplyr::select(-Dist)

padd_Y_1.2sqm_low_and_upplim_data <- cbind(padd_Y_1.2sqm_lowlim_data, 
                                           padd_Y_1.2sqm_upplim_data) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile, 
                                      labels=c("80%", "90%", "95%")))

padd_Y_1.2sqm_dpcr_vgram_plt <- padd_Y_1.2sqm_low_and_upplim_data %>% 
  ggplot(aes(Dist, SemiVariance)) +
  geom_point(shape = 19, size=3) +
  geom_line(linetype = "dashed") +
  geom_ribbon(aes(ymin = VarLowLim, ymax=VarUppLim, 
                  group= LeftQuantile, 
                  fill = LeftQuantile),
              alpha=0.20) +
  theme_classic() +
  labs(fill = "") +
  scale_fill_manual(values = c("gray30", "gray40", "gray50")) +
  theme(legend.position = c(0.68, 0.92),
        legend.direction="horizontal",
        legend.margin = margin(0,0,0,0 ,"mm"),
        legend.key.size = unit(7, "mm"),
        legend.box.background = element_rect(fill = "white", color="white"),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size = axis_text_size, face = axis_text_face),
        axis.text.y = element_text(size = axis_text_size, face = axis_text_face),
        axis.title = element_text(size = axis_title_size, face = axis_title_face)) +
  guides(fill = guide_legend(override.aes = 
                               list(alpha=0.7,fill = c("gray40", "gray60", "gray80"))),
         size = guide_legend(order=10)) +
  labs(x = "Distance, h (meter)",
       y = expression(gamma(h))) +
  scale_x_continuous(breaks = c(0, 100, 250, 400, 500), limits = c(0, 500)) +
  ylim(0, 0.006)
################################################################################
padd_Y_1.2sqm_dPCR_vgram_jointplt <- (padd_Y_1.2sqm_vgram_cloud_plt2 / 
                                        padd_Y_1.2sqm_dpcr_vgram_plt)
################################################################################


################################################################################
#                                                                              #  
# b. Paddock = Y; variable = dPCR; Sampling design = "c" (rake and 8.6 m2)     #
#                                                                              #  
################################################################################
# Filter out the necessary data
padd_Y_8.6sqm <- spat_dep_data2 %>% 
  dplyr::filter(PaddockFac == "Y", 
                NestedSubsample == "c")
# Compute semi-variances
padd_Y_8.6sqm_vgram_cloud <- compute_svgm(padd_Y_8.6sqm, Raw_dPCR)
################################################################################
# Create the semi-variance point cloud                                         #
################################################################################
padd_Y_8.6sqm_vgram_cloud_plt <- padd_Y_8.6sqm_vgram_cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha = 0.1, size= point_size, 
             fill="black", shape = 21, color="black") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = axis_text_size, face = axis_text_face),
        axis.text.y = element_text(size = axis_text_size, face = axis_text_face),
  axis.title = element_text(size = axis_title_size, face = axis_title_face)) +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2)) +
  scale_x_continuous(breaks = c(0, 100, 250, 400, 500), limits = c(0, 550))

padd_Y_8.6sqm_vgram_cloud_plt2 <- padd_Y_8.6sqm_vgram_cloud_plt +
  ggtitle(bquote(bold(atop("(b) Field-SD-b, 8.6"~m^ 2~"(r)")))) +
  theme(plot.title = element_text(size = 15, face = "bold"))
################################################################################
# Construct the variogram                                                      #
################################################################################
################################################################################
padd_Y_8.6sqm_vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data2, 
                                                            var=Raw_dPCR, 
                                                            paddock="Y", 
                                                            samp_design="c",
                                            dist_brks=c(0, 75, 250, 450, 550),
                                                            boot_sim=1000, 
                                                            boot_seed=100)
################################################################################
# Save the data                                                                #
################################################################################
saveRDS(padd_Y_8.6sqm_vgram_data,
        here::here("DerivedData", "padd_Y_8.6sqm_vgram_data.rds"))
################################################################################
################################################################################
padd_Y_8.6sqm_vgram_data <- readRDS(here::here("DerivedData", 
                                               "padd_Y_8.6sqm_vgram_data.rds"))
################################################################################
padd_Y_8.6sqm_lowlim_data <- padd_Y_8.6sqm_vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  tidyr::pivot_longer(cols = -c(Dist, SemiVariance), 
                      names_to="LeftQuantile", 
                      values_to="VarLowLim")

padd_Y_8.6sqm_upplim_data <- padd_Y_8.6sqm_vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  tidyr::pivot_longer(cols = -Dist, 
                      names_to="RightQuantile", 
                      values_to="VarUppLim") %>%
  dplyr::select(-Dist)

padd_Y_8.6sqm_low_and_upplim_data <- cbind(padd_Y_8.6sqm_lowlim_data, 
                                           padd_Y_8.6sqm_upplim_data) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile, 
                                      labels=c("80%", "90%", "95%")))

padd_Y_8.6sqm_dpcr_vgram_plt <- padd_Y_8.6sqm_low_and_upplim_data %>% 
  ggplot(aes(Dist, SemiVariance)) +
  geom_point(shape = 19, size=3) +
  geom_line(linetype = "dashed") +
  geom_ribbon(aes(ymin = VarLowLim, ymax=VarUppLim, 
                  group= LeftQuantile, 
                  fill = LeftQuantile),
              alpha=0.20) +
  theme_classic() +
  labs(fill = "") +
  scale_fill_manual(values = c("gray30", "gray40", "gray50")) +
  theme(legend.position = c(0.68, 0.92),
        legend.direction="horizontal",
        legend.margin = margin(0,0,0,0 ,"mm"),
        legend.key.size = unit(7, "mm"),
        legend.box.background = element_rect(fill = "white", color="white"),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size = axis_text_size, face = axis_text_face),
        axis.text.y = element_text(size = axis_text_size, face = axis_text_face),
        axis.title = element_text(size = axis_title_size, face = axis_title_face)) +
  guides(fill = guide_legend(override.aes = 
                               list(alpha=0.7,fill = c("gray40", "gray60", "gray80"))),
         size = guide_legend(order=10)) +
  labs(x = "Distance, h (meter)",
       y = expression(gamma(h))) +
  ylim(0, 0.001) +
  scale_x_continuous(breaks = c(0, 100, 250, 400, 500), limits = c(0, 550))
################################################################################
padd_Y_8.6sqm_dPCR_vgram_jointplt <- (padd_Y_8.6sqm_vgram_cloud_plt2 / 
                                        padd_Y_8.6sqm_dpcr_vgram_plt)
################################################################################


################################################################################
#                                                                              #  
# c. Paddock = Y; variable = dPCR; Sampling design = "d" (rake and 60 m2)      #
#                                                                              #  
################################################################################
# Filter out the necessary data
padd_Y_60sqm <- spat_dep_data2 %>% 
  dplyr::filter(PaddockFac == "Y", 
                NestedSubsample == "d")
# Compute semi-variances
padd_Y_60sqm_vgram_cloud <- compute_svgm(padd_Y_60sqm, Raw_dPCR)
################################################################################
# Create the semi-variance point cloud                                         #
################################################################################
padd_Y_60sqm_vgram_cloud_plt <- padd_Y_60sqm_vgram_cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha = 0.1, size= point_size, 
             fill="black", shape = 21, color="black") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = axis_text_size, face = axis_text_face),
        axis.text.y = element_text(size = axis_text_size, face = axis_text_face),
        axis.title = element_text(size = axis_title_size, face = axis_title_face)) +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2))  +
  xlim(0, 650)

padd_Y_60sqm_vgram_cloud_plt2 <- padd_Y_60sqm_vgram_cloud_plt +
  ggtitle(bquote(bold(atop("(c) Field-SD-b, 60"~m^ 2~"(r)")))) +
  theme(plot.title = element_text(size = 15, face = "bold"))
################################################################################
# Construct the variogram                                                      #
################################################################################
################################################################################
padd_Y_60sqm_vgram_data <- readRDS(here::here("DerivedData", 
                                              "padd_Y_60sqm_vgram_data.rds"))
################################################################################
padd_Y_60sqm_lowlim_data <- padd_Y_60sqm_vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  tidyr::pivot_longer(cols = -c(Dist, SemiVariance), 
                      names_to="LeftQuantile", 
                      values_to="VarLowLim")

padd_Y_60sqm_upplim_data <- padd_Y_60sqm_vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  tidyr::pivot_longer(cols = -Dist, 
                      names_to="RightQuantile", 
                      values_to="VarUppLim") %>%
  dplyr::select(-Dist)

padd_Y_60sqm_low_and_upplim_data <- cbind(padd_Y_60sqm_lowlim_data, 
                                          padd_Y_60sqm_upplim_data) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile, 
                                      labels=c("80%", "90%", "95%")))

options(scipen=10000)
padd_Y_60sqm_dpcr_vgram_plt <- padd_Y_60sqm_low_and_upplim_data %>% 
  ggplot(aes(Dist, SemiVariance)) +
  geom_point(shape = 19, size=3) +
  geom_line(linetype = "dashed") +
  geom_ribbon(aes(ymin = VarLowLim, ymax=VarUppLim, 
                  group= LeftQuantile, 
                  fill = LeftQuantile),
              alpha=0.20) +
  theme_classic() +
  labs(fill = "") +
  scale_fill_manual(values = c("gray30", "gray40", "gray50")) +
  theme(legend.position = c(0.68, 0.92),
        legend.direction="horizontal",
        legend.margin = margin(0,0,0,0 ,"mm"),
        legend.key.size = unit(7, "mm"),
        legend.box.background = element_rect(fill = "white", color="white"),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size = axis_text_size, face = axis_text_face),
        axis.text.y = element_text(size = axis_text_size, face = axis_text_face),
        axis.title = element_text(size = axis_title_size, face = axis_title_face)) +
  guides(fill = guide_legend(override.aes = 
                               list(alpha=0.7,fill = c("gray40", "gray60", "gray80"))),
         size = guide_legend(order=10)) +
  labs(x = "Distance, h (meter)",
       y = expression(gamma(h))) +
  ylim(0, 0.0004) +
  xlim(0, 650)
################################################################################
padd_Y_60sqm_dPCR_vgram_jointplt <- (padd_Y_60sqm_vgram_cloud_plt2 / 
                                       padd_Y_60sqm_dpcr_vgram_plt)
################################################################################

################################################################################
#                                                                              #
#       Field-Y Variogram Plot for Paper                                       #
#                                                                              #
################################################################################
padd_Y_dPCR_vgram_jointplt <- (padd_Y_1.2sqm_dPCR_vgram_jointplt |
                                 padd_Y_8.6sqm_dPCR_vgram_jointplt |
                                 padd_Y_60sqm_dPCR_vgram_jointplt)

ggsave(here::here("Figures", 
                  "Field_Y_VariogramPlots.jpeg"),
       padd_Y_dPCR_vgram_jointplt, height = 9, width = 15)

ggsave(here::here("Figures", 
                  "Field_Y_VariogramPlots.png"),
       padd_Y_dPCR_vgram_jointplt, height = 9, width = 15)

ggsave(here::here("Figures", 
                  "Field_Y_VariogramPlots.tiff"),
       padd_Y_dPCR_vgram_jointplt, height = 9, width = 15)
























