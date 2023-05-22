################################################################################
# Load required packages                                                       #  
################################################################################
library(here)
library(readxl)
library(tidyverse)
library(sp)
library(spatstat)
library(ggpubr)
library(sf)
library(geoR)
library(viridis)
library(RColorBrewer)
library(ggtext)
library(patchwork)
library(gridExtra)
################################################################################
# Source utility functions                                                     #
################################################################################
source(here::here("RCode", "Utils_ForPaper.R"))
################################################################################
# Define labels for sampling strategies                                        #
################################################################################
lab_vals2  <- c(bquote(bold(atop("1.2"~m^ 2~"(hp)"))),
                bquote(bold(atop("1.2"~m^ 2~"(r)"))),
                bquote(bold(atop("8.6"~m^ 2~"(r)"))),
                bquote(bold(atop("60"~m^ 2~"(r)"))))
################################################################################
# Load data recorded using spatially-nested sampling strategies                #
################################################################################
spat_dep_data <- readRDS(here::here("DerivedData",
                                    "spat_dep_data.rds"))
spat_dep_sf_projected <- readRDS(here::here("DerivedData",
                                    "spat_dep_sf_projected.rds"))
################################################################################
# Apply an affine-transformation to shift the origin at (0,0)                  #
################################################################################
# Obtain the projected x-y coordinates
sp_dep_coords0 <- st_geometry(spat_dep_sf_projected)
# Obtain the extent of bounding box for the study region
sp_dep_bbox0 <- st_bbox(sp_dep_coords0)
# Shift origin to (0,0)
sp_dep_coords1 <- sp_dep_coords0 - c(sp_dep_bbox0['xmin'], sp_dep_bbox0['ymin'])
plot(sp_dep_coords1)
sp_dep_coords1 %>% ggplot() + geom_sf()
################################################################################
# Apply the new coordinate system to the sf object                             #
################################################################################
spat_dep_sf_projected2 <- st_set_geometry(spat_dep_sf_projected,
                                          sp_dep_coords1)
################################################################################
# Apply appropriate labels for Nested subsamples                               #
################################################################################
spat_dep_sf_projected3 <- spat_dep_sf_projected2 %>%
  dplyr::mutate(NestedSubsample = factor(NestedSubsample,
                                         labels = lab_vals2[c(2,3,4)]))
################################################################################
# Define color attributes for plotting Ptt FR frequency                        #
################################################################################
freq_col_lims <- c(0, 0.45)
freq_col_brks <- c(0.0, 0.15, 0.30, 0.45)
freq_col_palette <- brewer.pal(6, name = "YlGnBu")
freq_col_labels <- c("0%", "15%", "30%", "45%")
axis_text_size <- 15
################################################################################
#                                                                              #
# Field-W 60 m2 sampling positions plot                                        #
#                                                                              #
################################################################################
nested_levels <- levels(spat_dep_sf_projected3$NestedSubsample)

dPCR_paddW_60m2 <- spat_dep_sf_projected3 %>% 
  dplyr::filter(PaddockFac == "W", NestedSubsample == nested_levels[3]) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=10) +
  geom_sf(shape = 1, size=10.1) +
  coord_sf() +
  theme_bw() +
  facet_grid(. ~ NestedSubsample, labeller = label_parsed) + # switch = "y"
  theme(legend.position = "none",
        strip.text = element_text(size=16, margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=(axis_text_size+1)),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 20, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_gradientn(colours = freq_col_palette,
                        breaks = freq_col_brks,
                        limits = freq_col_lims,
                        labels = freq_col_labels,
                        guide = guide_colorbar(barwidth = unit(6, "cm"))) +
  labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")),
       x = "Relative position distance (m)",
       y = "Relative position distance (m)")

dPCR_paddW_60m2_annotated <- dPCR_paddW_60m2 + 
                   annotate("rect", xmin = 215, xmax = 245,
                            ymin = 340, ymax = 360,
                            alpha = 0, color = "gray20", linewidth=1.2) +
                  annotate("rect", xmin = 150, xmax = 180,
                            ymin = 100, ymax = 120,
                            alpha = 0, color = "gray20", linewidth=1.2) +
                  annotate("rect", xmin = 180, xmax = 210,
                                  ymin = 15, ymax = 35,
                                  alpha = 0, color = "gray20", linewidth=1.2)

dPCR_paddW_60m2_annotated
################################################################################
#                                                                              #
# Field-W 8.6 m2 sampling positions for Rep-1                                  #
#                                                                              #
################################################################################
dPCR_paddW_8.6m2_rep1_v0 <- spat_dep_sf_projected3 %>% 
  bind_cols(st_coordinates(spat_dep_sf_projected3)) %>%
  filter(PaddockFac == "W", 
         NestedSubsample == nested_levels[2],
         ReplicationNo == 1) %>%
  # ggplot() +
  ggplot(aes(X, Y, fill = Raw_dPCR), size=9) +
  geom_point(shape = 21, size=9.1, stroke=1.5) +
  facet_grid(. ~ NestedSubsample, labeller = label_parsed) +
  theme_bw() +
  theme(strip.text = element_text(size=16, 
                                  margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=axis_text_size),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  xlim(205, 265) +
  ylim(346, 356) +
  scale_fill_gradientn(colours = freq_col_palette,
                        breaks = freq_col_brks,
                        limits = freq_col_lims,
                        labels = freq_col_labels) +
  coord_fixed()

dPCR_paddW_8.6m2_rep1 <- spat_dep_sf_projected3 %>% 
  bind_cols(st_coordinates(spat_dep_sf_projected3)) %>%
  filter(PaddockFac == "W", 
         NestedSubsample == nested_levels[2],
         ReplicationNo == 1) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=8) +
  geom_sf(shape = 1, size=8.1, stroke=1.5) +
  coord_sf(xlim = c(205, 265), ylim = c(346, 356),
           label_axes = list(bottom = "E", right = "N")) +
  facet_grid(. ~ NestedSubsample, labeller = label_parsed) +
  theme_bw() +
  theme(strip.text = element_text(size=16, 
                                  margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=axis_text_size),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  scale_color_gradientn(colours = freq_col_palette,
                       breaks = freq_col_brks,
                       limits = freq_col_lims,
                       labels = freq_col_labels) 

dPCR_paddW_8.6m2_rep1_annotated <- dPCR_paddW_8.6m2_rep1 +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.placement = "outside") +
    annotate("rect", xmin = 232, xmax = 238,
           ymin = 350, ymax = 352,
           alpha = 0, color = "gray20", linewidth=1.2)
################################################################################
#                                                                              #
# Field-W 8.6 m2 sampling positions for Rep-2                                  #
#                                                                              #
################################################################################
dPCR_paddW_8.6m2_rep2 <- spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W", 
         NestedSubsample == nested_levels[2],
         ReplicationNo == 2) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=8) +
  geom_sf(shape = 1, size=8.1, stroke=1.5) +
  coord_sf(xlim = c(140, 200), ylim = c(108.5, 118.5),
           label_axes = list(bottom = "E", right = "N")) +
  theme_bw() +
  theme(strip.text = element_text(size=16, 
                                  margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=axis_text_size),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  scale_color_gradientn(colours = freq_col_palette,
                        breaks = freq_col_brks,
                        limits = freq_col_lims,
                        labels = freq_col_labels)

dPCR_paddW_8.6m2_rep2_annotated <- dPCR_paddW_8.6m2_rep2 +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  annotate("rect", xmin = 168, xmax = 174,
           ymin = 112, ymax = 114,
           alpha = 0, color = "gray20", linewidth=1.2)

dPCR_paddW_8.6m2_rep2_annotated
################################################################################
#                                                                              #
# Field-W 8.6 m2 sampling positions for Rep-3                                  #
#                                                                              #
################################################################################
dPCR_paddW_8.6m2_rep3 <- spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W", 
         NestedSubsample == nested_levels[2],
         ReplicationNo == 3) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=8) +
  geom_sf(shape = 1, size=8.1, stroke=1.5) +
  coord_sf(xlim = c(165, 225), ylim = c(20, 30),
           label_axes = list(bottom = "E", right = "N")) +
  theme_bw() +
  theme(strip.text = element_text(size=16, 
                                  margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=axis_text_size),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  scale_color_gradientn(colours = freq_col_palette,
                        breaks = freq_col_brks,
                        limits = freq_col_lims,
                        labels = freq_col_labels,
                        guide = guide_colorbar(barwidth = unit(8, "cm")))



dPCR_paddW_8.6m2_rep3_with_legend <- dPCR_paddW_8.6m2_rep3 +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 14, face="bold"),
        legend.title = element_text(size = 18, face="bold")) +
  labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")))
  

dPCR_paddW_8.6m2_rep3_annotated <-  dPCR_paddW_8.6m2_rep3 +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  annotate("rect", xmin = 195, xmax = 201,
           ymin = 24.2, ymax = 26.2,
           alpha = 0, color = "gray20", linewidth=1.2)


dPCR_paddW_8.6m2_rep3_annotated

################################################################################
# Extract Legend                                                               #
################################################################################
dPCR_legend <- extract_legend(dPCR_paddW_8.6m2_rep3_with_legend)
plot(dPCR_legend)
################################################################################
################################################################################
#                                                                              #
# Field-W 8.6 m2 sampling positions for all replications                       #
#                                                                              #
################################################################################
dPCR_paddW_8.6m2_all_rep <- (dPCR_paddW_8.6m2_rep1_annotated/
   dPCR_paddW_8.6m2_rep2_annotated/
    dPCR_paddW_8.6m2_rep3_annotated)
################################################################################
#                                                                              #
# Field-W 1.2 m2 sampling positions for Rep-1                                  #
#                                                                              #
################################################################################
spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W",
         NestedSubsample == nested_levels[1],
         ReplicationNo == 1) %>% st_coordinates() %>% summary()


dPCR_paddW_1pt2m2_rep1 <- spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W",
         NestedSubsample == nested_levels[1],
         ReplicationNo == 1) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=6) +
  geom_sf(shape = 1, size=6.1, stroke=1.5) +
  coord_sf(xlim = c(230, 240), ylim = c(350, 352),
           label_axes = list(bottom = "E", right = "N")) +
  facet_grid(. ~ NestedSubsample, labeller = label_parsed) +
  theme_bw() +
  theme(strip.text = element_text(size=16, 
                                  margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=axis_text_size),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  scale_color_gradientn(colours = freq_col_palette,
                        breaks = freq_col_brks,
                        limits = freq_col_lims,
                        labels = freq_col_labels)


################################################################################
#                                                                              #
# Field-W 1.2 m2 sampling positions for Rep-2                                  #
#                                                                              #
################################################################################
spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W",
         NestedSubsample == nested_levels[1],
         ReplicationNo == 2) %>% st_coordinates() %>% summary()


dPCR_paddW_1pt2m2_rep2 <- spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W",
         NestedSubsample == nested_levels[1],
         ReplicationNo == 2) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=6) +
  geom_sf(shape = 1, size=6.1, stroke=1.5) +
  coord_sf(xlim = c(165, 175), ylim = c(113, 115),
           label_axes = list(bottom = "E", right = "N")) +
  theme_bw() +
  theme(strip.text = element_text(size=16, 
                                  margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=axis_text_size),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  scale_color_gradientn(colours = freq_col_palette,
                        breaks = freq_col_brks,
                        limits = freq_col_lims,
                        labels = freq_col_labels)


################################################################################
#                                                                              #
# Field-W 1.2 m2 sampling positions for Rep-3                                  #
#                                                                              #
################################################################################
spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W",
         NestedSubsample == nested_levels[1],
         ReplicationNo == 3) %>% st_coordinates() %>% summary()


dPCR_paddW_1pt2m2_rep3 <- spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W",
         NestedSubsample == nested_levels[1],
         ReplicationNo == 3) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=6) +
  geom_sf(shape = 1, size=6.1, stroke=1.5) +
  coord_sf(xlim = c(192, 202), ylim = c(24, 26),
           label_axes = list(bottom = "E", right = "N")) +
  theme_bw() +
  theme(strip.text = element_text(size=16, 
                                  margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=axis_text_size),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  scale_color_gradientn(colours = freq_col_palette,
                        breaks = freq_col_brks,
                        limits = freq_col_lims,
                        labels = freq_col_labels)
################################################################################
dPCR_paddW_1pt2m2_all_rep <- (dPCR_paddW_1pt2m2_rep1/
                                dPCR_paddW_1pt2m2_rep2/
                               dPCR_paddW_1pt2m2_rep3)

################################################################################
# First combine the 1.2 and 8.6 m2                                             #
################################################################################

dPCR_paddW_plt0 <- ggarrange(dPCR_paddW_8.6m2_all_rep,
          dPCR_paddW_1pt2m2_all_rep,
          nrow=1,
          widths = c(2,1.2))

# Add legend
dPCR_paddW_plt1 <- ggarrange(dPCR_paddW_plt0,
                             dPCR_legend, 
                             nrow = 2,
                             heights = c(5,1))
################################################################################
# First combine all sampling strategies                                        #
################################################################################
dPCR_paddW_final_plt <- ggarrange(dPCR_paddW_60m2_annotated,
                             dPCR_paddW_plt1, nrow=1,
                             widths = c(2.5, 3.2))






ggsave(here::here("Figures", "dPCR_paddW_plt.jpeg"),
       dPCR_paddW_final_plt, height = 8, width= 16)

ggsave(here::here("Figures", "dPCR_paddW_plt.tiff"),
       dPCR_paddW_final_plt, height = 8, width= 16)

################################################################################
################################################################################
dPCR_paddW_8pt6m2 <- spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W", NestedSubsample == nested_levels[2]) %>%
  ggplot() +
  #geom_sf(data = rep_poly_data, color = "red", fill = "white") +
  geom_sf(aes(color = Raw_dPCR), size=5) +
  geom_sf(shape = 1, size=5.1) +
  coord_sf() +
  theme_bw() +
  facet_grid(.~ NestedSubsample, labeller = label_parsed) + # switch = "y"
  theme(legend.position = "none",
        strip.text = element_text(size=16, margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 14, face="bold")
        # legend.text = element_text(face="bold", size=12)
  ) +
  scale_color_gradientn(colours = freq_col_palette,
                        breaks = c(0.0, 0.15, 0.30, 0.45),
                        limits = c(0, 0.45),
                        labels = c("0%", "15%", "30%", "45%"),
                        guide = guide_colorbar(barwidth = unit(6, "cm"))) +
  labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")),
       x = "Relative position distance (m)",
       y = "Relative position distance (m)")

