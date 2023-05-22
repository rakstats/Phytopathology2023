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
# Analysis of Unit Size data                                                   #
################################################################################
unit_data <- readRDS(here::here("DerivedData", "unit_data.rds"))
glimpse(unit_data)
################################################################################
#                                                                              #
# Figure 7 in the paper: Spatial Distribution of fungicide resistance frequency# 
#       in stubble-borne Ptt (Ptt-FR) for different sampling strategies and    #
#       unit area sizes for the fields Field-M, Field-P, and Field-V           #
#                                                                              #
################################################################################
################################################################################
# We should remove all data points with Cyp51A_copies <= 29                    #    
################################################################################
################################################################################
# Create sampling strategy labels for faceting                                 #
################################################################################
lab_vals2  <- c(bquote(bold(atop("1.2"~m^ 2~"(hp)"))),
                bquote(bold(atop("1.2"~m^ 2~"(r)"))),
                bquote(bold(atop("8.6"~m^ 2~"(r)"))),
                bquote(bold(atop("60"~m^ 2~"(r)"))))
################################################################################
# Create new labels for Paddocks/Fields for the paper                          #
################################################################################
unit_data2 <- unit_data %>% 
              dplyr::filter(Cyp51A_copies > 29) %>%
              dplyr::mutate(PaddockFac = factor(PaddockFac,
                    levels = c("M", "P", "V"),
                    labels = c("SS-i", "SS-ii", "SS-iii"))) %>%
              dplyr::mutate(NestedSubsample = factor(NestedSubsample,
                                                     labels = lab_vals2))
################################################################################
# Create sf object using lon, lat and crs = 4326                               #
################################################################################
unit_sf <- sf::st_as_sf(unit_data2, coords = c("lon", "lat"), crs = 4326)
################################################################################
# Project Lon and Lat to planar coordinates (in meters) using crs=3112         #
################################################################################
unit_sf_projected <- sf::st_transform(unit_sf, 3112)
################################################################################
# To maintain anonymity, apply an affine transformation to the coordinates     #
# so that the origin of the transformed coordinates lies at (0,0)              #
################################################################################
unit_coords0 <- st_geometry(unit_sf_projected)
unit_bbox0 <- st_bbox(unit_coords0)
unit_coords1 <- unit_coords0 - c(unit_bbox0['xmin'], unit_bbox0['ymin'])
# Check the transformed coordinates
unit_coords1 %>% ggplot() + geom_sf() # check origin is at (0,0)
################################################################################
# Set this new coordinate system with origin at (0,0) to the projected data    #
################################################################################
unit_sf_projected2 <- st_set_geometry(unit_sf_projected, unit_coords1)
################################################################################
# Get a sense of the data to choose colour scheme and colour labels            #
################################################################################
summary(unit_sf_projected2$Raw_dPCR)
freq_col_lims <- c(0, 0.50)
freq_col_brks <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
freq_col_palette <- brewer.pal(6, name = "YlGnBu")
freq_col_labels <- c("0%", "10%", "20%", "30%", "40%", "50%")
################################################################################
# The following construction is OPTIONAL for using geom_sf_text()              #
# To use geom_sf_text, we can create an sf object that follows the same        #
# coordinate projection system (i.e., in meters with the origin at (0,0)) with #
# the desired labels to be inserted at the specified coordinates in the plot.  #
# However, we have not used this; instead we have used annotation(), which     #
# worked well for our figure. We have included a commented line of code        #
# if anyone wants to use the geom_sf_text() function for annotating their plots#
################################################################################
field_text_df <- data.frame(lon = c(1000, 5200, 6000),
                            lat = c(3000, 5000, 500),
                            text = c("SS-i", "SS-ii", "SS-iii"))
field_text_sf <- st_as_sf(field_text_df, coords = c("lon", "lat"))
################################################################################
#              Create the plot                                                 #
################################################################################
dPCR_spplot <- unit_sf_projected2 %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=2) +
  geom_sf(shape = 1, size=2.1) +
  coord_sf() +
  facet_wrap(vars(NestedSubsample), 
             labeller = label_parsed, nrow=2) +
  scale_color_gradientn(colours = freq_col_palette,
                        breaks = freq_col_brks,
                        limits = freq_col_lims,
                        labels = freq_col_labels,
                        guide = guide_colorbar(barwidth = unit(8, "cm"))) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size= rel(1), margin = margin(6,0,0,0)),
        axis.text = element_text(face="bold"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(face="bold"),
        axis.title = element_text(face="bold")) +
  labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")),
       x = "Relative distance (m)",
       y = "Relative distance (m)") +
  annotate("text", x = 1000 , y = 3000, size = 3.5, label = "SS-i", fontface = "bold") +
  annotate("text", x = 5200 , y = 5000, size = 3.5, label = "SS-ii", fontface = "bold") +
  annotate("text", x = 6000 , y = 500, size = 3.5, label = "SS-iii", fontface = "bold")
  # geom_sf_text(data = field_text_sf, aes(label = text))
################################################################################
# Save the figure                                                              #
################################################################################
ggsave(here::here("Figures", "Fig7_Spmap_Ptt_FRfreq.jpeg"),
       dPCR_spplot, width = 8.25, height = 6)
ggsave(here::here("Figures", "Fig7_Spmap_Ptt_FRfreq.tiff"),
       dPCR_spplot, width = 8.25, height = 6)
###############################  END  ##########################################

################################################################################
#                                                                              #
# Figure 2 in the paper: Semi-variances and variograms of fungicide resistance #
# frequency as functions of distance from differing unit-area-sizes of the     #
# nested transect sampling design.                                             #
# Plots included: RAKED (r) strategy for (a) 60 m2 for Field-W and (b) Field-Y #
#               and (c) 1.2 m2 for Field-W and (d) 8.6 m2 for Field-W          #
#                                                                              #
################################################################################
################################################################################
# Read in the Spatial Dependency data --- sampled using a nested scheme        #
################################################################################
spat_dep_data <- readRDS(here::here("DerivedData", "spat_dep_data.rds"))
spat_dep_sf_projected <- readRDS(here::here("DerivedData", "spat_dep_sf_projected.rds"))
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
# a. Paddock = W; variable = dPCR; Sampling design = "d" (rake and 60 m2)      #
#                                                                              #  
################################################################################
# Filter out the necessary data
padd_W_60sqm <- spat_dep_data2 %>% 
                dplyr::filter(PaddockFac == "W", 
                              NestedSubsample == "d")
# Compute semi-variances
padd_W_60sqm_vgram_cloud <- compute_svgm(padd_W_60sqm, Raw_dPCR)
################################################################################
# Create the semi-variance point cloud                                         #
################################################################################
padd_W_60sqm_vgram_cloud_plt <- padd_W_60sqm_vgram_cloud %>% 
                                ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha = 0.1, size= point_size, 
             fill="black", shape = 21, color="black") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = axis_text_size, face = axis_text_face),
        axis.text.y = element_text(size = axis_text_size, face = axis_text_face),
        axis.title = element_text(size = axis_title_size, face = axis_title_face)) +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2)) +
  xlim(0, 450)

padd_W_60sqm_vgram_cloud_plt2 <- padd_W_60sqm_vgram_cloud_plt +
      ggtitle(bquote(bold(atop("(a) Field-SD-a, 60"~m^ 2~"(r)")))) +
      theme(plot.title = element_text(size = 15, face = "bold"))
################################################################################
# Construct the variogram                                                      #
################################################################################
################################################################################
#May not need to run the following code as the data were saved and loaded below#
################################################################################
padd_W_60sqm_vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data2, 
                                                         var=Raw_dPCR, 
                                                         paddock="W", 
                                                         samp_design="d",
                                         dist_brks=c(0, 75, 150, 250, 350, 500),
                                                         boot_sim=1000, 
                                                         boot_seed=100)
########################## Already Saved #######################################
# saveRDS(padd_W_60sqm_vgram_data,                                             #
#         here::here("DerivedData", "padd_W_60sqm_vgram_data.rds"))            #
################################################################################
padd_W_60sqm_vgram_data <- readRDS(here::here("DerivedData", 
                                              "padd_W_60sqm_vgram_data.rds"))
################################################################################

padd_W_60sqm_lowlim_data <- padd_W_60sqm_vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  tidyr::pivot_longer(cols = -c(Dist, SemiVariance), 
                      names_to="LeftQuantile", 
                      values_to="VarLowLim")

padd_W_60sqm_upplim_data <- padd_W_60sqm_vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  tidyr::pivot_longer(cols = -Dist, 
                      names_to="RightQuantile", 
                      values_to="VarUppLim") %>%
  dplyr::select(-Dist)

padd_W_60sqm_low_and_upplim_data <- cbind(padd_W_60sqm_lowlim_data, 
                                          padd_W_60sqm_upplim_data) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile, 
                                      labels=c("80%", "90%", "95%")))

padd_W_60sqm_dpcr_vgram_plt <- padd_W_60sqm_low_and_upplim_data %>% 
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
  theme(legend.position = c(0.72, 0.92),
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
  ylim(c(0, 0.07)) + xlim(0, 450)
################################################################################
padd_W_60sqm_dPCR_vgram_jointplt <- (padd_W_60sqm_vgram_cloud_plt2 / padd_W_60sqm_dpcr_vgram_plt)
# plot_annotation(tag_levels = list(c("Field: W, 60 m2 (r)")))
################################################################################
################################################################################
#                                                                              #
# b. Paddock = Y; variable = dPCR; Sampling design = "d" (rake and 60 m2)      #
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
       y = expression(0.5*(z[s+h] - z[s])^2)) +
  xlim(0, 650)

padd_Y_60sqm_vgram_cloud_plt2 <- padd_Y_60sqm_vgram_cloud_plt +
  ggtitle(bquote(bold(atop("(b) Field-SD-b, 60"~m^ 2~"(r)")))) +
  theme(plot.title = element_text(size = 15, face = "bold"))
################################################################################
# Construct the variogram                                                      #
################################################################################
################################################################################
#May not need to run the following code as the data were saved and loaded below#
################################################################################
options(scipen=10000)
padd_Y_60sqm_vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data2, 
                                                           var=Raw_dPCR, 
                                                           paddock="Y", 
                                                           samp_design="d",
                                  dist_brks=c(0, 75, 150, 250, 350, 450, 650), 
                                                           boot_sim=1000, 
                                                           boot_seed=100)
#######################Already Saved ###########################################
# saveRDS(padd_Y_60sqm_vgram_data,                                             #
#         here::here("DerivedData", "padd_Y_60sqm_vgram_data.rds"))            #
################################################################################
padd_Y_60sqm_vgram_data <- readRDS(here::here("DerivedData", 
                                              "padd_Y_60sqm_vgram_data"))
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
  theme(legend.position = c(0.72, 0.92),
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
padd_Y_60sqm_dPCR_vgram_jointplt <- (padd_Y_60sqm_vgram_cloud_plt2 / padd_Y_60sqm_dpcr_vgram_plt)
################################################################################
figure2_top_panel <- (padd_W_60sqm_dPCR_vgram_jointplt | padd_Y_60sqm_dPCR_vgram_jointplt)

ggsave(here::here("Figures", "Figure2_a_b_Variogram_Field_W_and_Y_60m2_rake.jpeg"),
       figure2_top_panel, height = 11, width = 11)
################################################################################


################################################################################
#                                                                              #  
# c. Paddock = W; variable = dPCR; Sampling design = "b" (rake and 1.2 m2)     #
#                                                                              #  
################################################################################
# Filter out the necessary data
padd_W_1.2sqm <- spat_dep_data2 %>% 
  dplyr::filter(PaddockFac == "W", 
                NestedSubsample == "b")
# Compute semi-variances
padd_W_1.2sqm_vgram_cloud <- compute_svgm(padd_W_1.2sqm, Raw_dPCR)
################################################################################
# Create the semi-variance point cloud                                         #
################################################################################
padd_W_1.2sqm_vgram_cloud_plt <- padd_W_1.2sqm_vgram_cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha = 0.1, size= point_size, 
             fill="black", shape = 21, color="black") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = axis_text_size, face = axis_text_face),
        axis.text.y = element_text(size = axis_text_size, face = axis_text_face),
        axis.title = element_text(size = axis_title_size, face = axis_title_face)) +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2)) +
  xlim(0, 350)

padd_W_1.2sqm_vgram_cloud_plt2 <- padd_W_1.2sqm_vgram_cloud_plt +
  ggtitle(bquote(bold(atop("(c) Field-SD-a, 1.2"~m^ 2~"(r)")))) +
  theme(plot.title = element_text(size = 15, face = "bold"))
################################################################################
# Construct the variogram                                                      #
################################################################################
################################################################################
#May not need to run the following code as the data were saved and loaded below#
################################################################################
padd_W_1.2sqm_vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data2, 
                                                            var=Raw_dPCR, 
                                                            paddock="W", 
                                                            samp_design="b",
                                              dist_brks=c(0, 60, 150, 250, 350),
                                                            boot_sim=1000, 
                                                            boot_seed=100)
#######################Already Saved ###########################################
# saveRDS(padd_W_1.2sqm_vgram_data,                                            #                                            
#         here::here("DerivedData", "padd_W_1.2sqm_vgram_data.rds"))           #       
################################################################################

padd_W_1.2sqm_vgram_data <- readRDS(here::here("DerivedData", 
                                              "padd_W_1.2sqm_vgram_data.rds"))
################################################################################

padd_W_1.2sqm_lowlim_data <- padd_W_1.2sqm_vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  tidyr::pivot_longer(cols = -c(Dist, SemiVariance), 
                      names_to="LeftQuantile", 
                      values_to="VarLowLim")

padd_W_1.2sqm_upplim_data <- padd_W_1.2sqm_vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  tidyr::pivot_longer(cols = -Dist, 
                      names_to="RightQuantile", 
                      values_to="VarUppLim") %>%
  dplyr::select(-Dist)

padd_W_1.2sqm_low_and_upplim_data <- cbind(padd_W_1.2sqm_lowlim_data, 
                                           padd_W_1.2sqm_upplim_data) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile, 
                                      labels=c("80%", "90%", "95%")))

padd_W_1.2sqm_dpcr_vgram_plt <- padd_W_1.2sqm_low_and_upplim_data %>% 
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
  theme(legend.position = c(0.72, 0.92),
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
  ylim(0, 0.004) + xlim(0, 350)
################################################################################
padd_W_1.2sqm_dPCR_vgram_jointplt <- (padd_W_1.2sqm_vgram_cloud_plt2 / padd_W_1.2sqm_dpcr_vgram_plt)
################################################################################

################################################################################
#                                                                              #  
# d. Paddock = W; variable = dPCR; Sampling design = "c" (rake and 8.6 m2)     #
#                                                                              #  
################################################################################
# Filter out the necessary data
padd_W_8.6sqm <- spat_dep_data2 %>% 
  dplyr::filter(PaddockFac == "W", 
                NestedSubsample == "c")
# Compute semi-variances
padd_W_8.6sqm_vgram_cloud <- compute_svgm(padd_W_8.6sqm, Raw_dPCR)
################################################################################
# Create the semi-variance point cloud                                         #
################################################################################
padd_W_8.6sqm_vgram_cloud_plt <- padd_W_8.6sqm_vgram_cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha = 0.1, size= point_size, 
             fill="black", shape = 21, color="black") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = axis_text_size, face = axis_text_face),
        axis.text.y = element_text(size = axis_text_size, face = axis_text_face),
        axis.title = element_text(size = axis_title_size, face = axis_title_face)) +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2)) +
  xlim(0, 350)

padd_W_8.6sqm_vgram_cloud_plt2 <- padd_W_8.6sqm_vgram_cloud_plt +
  ggtitle(bquote(bold(atop("(d) Field-SD-a, 8.6"~m^ 2~"(r)")))) +
  theme(plot.title = element_text(size = 15, face = "bold"))
################################################################################
# Construct the variogram                                                      #
################################################################################
################################################################################
#May not need to run the following code as the data were saved and loaded below#
################################################################################
padd_W_8.6sqm_vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data2, 
                                                            var=Raw_dPCR, 
                                                            paddock="W", 
                                                            samp_design="c",
                                            dist_brks=c(0, 70, 150, 300, 350),
                                                            boot_sim=1000, 
                                                            boot_seed=100)
#######################Already Saved ###########################################
# saveRDS(padd_W_8.6sqm_vgram_data,                                            #                                          
#        here::here("DerivedData", "padd_W_8.6sqm_vgram_data.rds"))            #         
################################################################################
padd_W_8.6sqm_vgram_data <- readRDS(here::here("DerivedData", 
                                              "padd_W_8.6sqm_vgram_data.rds"))
################################################################################
padd_W_8.6sqm_lowlim_data <- padd_W_8.6sqm_vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  tidyr::pivot_longer(cols = -c(Dist, SemiVariance), 
                      names_to="LeftQuantile", 
                      values_to="VarLowLim")

padd_W_8.6sqm_upplim_data <- padd_W_8.6sqm_vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  tidyr::pivot_longer(cols = -Dist, 
                      names_to="RightQuantile", 
                      values_to="VarUppLim") %>%
  dplyr::select(-Dist)

padd_W_8.6sqm_low_and_upplim_data <- cbind(padd_W_8.6sqm_lowlim_data, 
                                           padd_W_8.6sqm_upplim_data) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile, 
                                      labels=c("80%", "90%", "95%")))

padd_W_8.6sqm_dpcr_vgram_plt <- padd_W_8.6sqm_low_and_upplim_data %>% 
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
  theme(legend.position = c(0.72, 0.92),
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
  ylim(0, 0.003) +
  xlim(0, 350)
################################################################################
padd_W_8.6sqm_dPCR_vgram_jointplt <- (padd_W_8.6sqm_vgram_cloud_plt2 / padd_W_8.6sqm_dpcr_vgram_plt)
################################################################################
################################################################################
figure2_bottom_panel <- (padd_W_1.2sqm_dPCR_vgram_jointplt | padd_W_8.6sqm_dPCR_vgram_jointplt)

ggsave(here::here("Figures", 
                  "Figure2_c_d_Variogram_Field_W_1pt2m2_and_8pt6m2_rake.jpeg"),
       figure2_bottom_panel, height = 11, width = 11)
################################################################################










