library(here)
library(tidyverse)
library(sf)
library(gstat)
library(sp)
library(rlang)
library(ggpubr)
#################################################################################
# Load the file that contains utility functions for bootstrapped envelopes      #
#################################################################################
source(here::here("FilesForLeon_9June2022", "RCode", "Utils_ForPaper.R"))
################################################################################
# Prepare the Spatial Dependency data                                          #
################################################################################
spat_dep_data <- readRDS(here::here("DerivedData", "spat_dep_data.rds"))
spat_dep_sf_projected <- readRDS(here::here("DerivedData", "spat_dep_sf_projected.rds"))
##################################################################################
# Add the Projected Coordinates to the Spatial Dependency Dataset                #
##################################################################################
spat_dep_data <-  spat_dep_data %>% 
  bind_cols(st_coordinates(spat_dep_sf_projected))
################################################################################
#1.1a Paddock = W; variable = dPCR; Sampling design = b                        #
################################################################################
# Create the Variogram Cloud
padd_W_1.2sqm <- spat_dep_data %>% 
                 dplyr::filter(PaddockFac == "W", 
                               NestedSubsample == "b")

paddW_1.2sqm_Vgram_Cloud <- compute_svgm(padd_W_1.2sqm, Raw_dPCR)

paddW_1.2sqm_Vgram_Cloud_plt <- paddW_1.2sqm_Vgram_Cloud %>% ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha=0.1, size=3, fill="blue", shape = 21, color="black") +
  theme_minimal() +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2))

paddW_1.2sqm_Vgram_Cloud_plt

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddW_1pt2sqm_dPCR_Vgram_Cloud.png"),
       plot = paddW_1.2sqm_Vgram_Cloud_plt,
       height = 3, width=4.5)
############################################################
paddW_1.2sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                    var=Raw_dPCR, 
                                                    paddock="W", 
                                                    samp_design="b",
                                              dist_brks=c(0, 60, 200, 300, 350), 
                                                         boot_sim=1000, 
                                                         boot_seed=100)

lowLim_data <- paddW_1.2sqm_Vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_data <- paddW_1.2sqm_Vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  dplyr::select(-Dist)

low_and_uppLim_data <- cbind(lowLim_data, uppLim_data) %>%
                       mutate(LeftQuantile = factor(LeftQuantile,
                                  labels=c("80%", "90%", "95%")))

paddW_1.2sqm_dPCR_Vgram_Plt <- low_and_uppLim_data %>% ggplot(aes(Dist, SemiVariance)) +
  geom_point(shape = 19, size=3) +
  geom_line(linetype = "dashed") +
  geom_ribbon(aes(ymin = VarLowLim, ymax=VarUppLim, 
                  group= LeftQuantile, 
                  fill = LeftQuantile),
              alpha=0.20) +
  theme_classic() +
  labs(fill = "") +
  scale_fill_manual(values = c("gray30", "gray40", "gray50")) +
  theme(legend.position = c(0.75, 0.92),
        legend.direction="horizontal",
        legend.margin = margin(0,0,0,0 ,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.box.background = element_rect(fill = "white", color="white")) +
  guides(fill = guide_legend(override.aes = 
         list(alpha=0.7,fill = c("gray40", "gray60", "gray80"))),
         size = guide_legend(order=10)) +
  labs(x = "Distance, h (meter)",
       y = expression(gamma(h))) +
  ylim(c(0, 0.0035))

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddW_1pt2sqm_dPCR_Variogram.png"),
       plot = paddW_1.2sqm_dPCR_Vgram_Plt,
       height = 3, width=4.5)

paddW_1.2sqm_dPCR_Vgram_JointPlt <- ggarrange(paddW_1.2sqm_Vgram_Cloud_plt,
                                              paddW_1.2sqm_dPCR_Vgram_Plt)

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddW_1pt2sqm_dPCR_Cloud_and_Variogram.png"),
       plot = paddW_1.2sqm_dPCR_Vgram_JointPlt,
       height = 3, width=9)
##########################################################
################################################################################
#1.2a Paddock = W; variable = dPCR; Sampling design = c                        #
################################################################################
# Create the Variogram Cloud
padd_W_7sqm <- spat_dep_data %>% 
  dplyr::filter(PaddockFac == "W", 
                NestedSubsample == "c")

paddW_7sqm_Vgram_Cloud <- compute_svgm(padd_W_7sqm, Raw_dPCR)

paddW_7sqm_Vgram_Cloud_plt <- paddW_7sqm_Vgram_Cloud %>% 
                              ggplot(aes(x=Dist, y=Svar)) +
                              geom_point(alpha=0.1, size=3, 
                                        fill="blue", shape = 21, 
                                        color="black") +
                              theme_minimal() +
                              labs(x = "Distance, h (meter)", 
                                   y = expression(0.5*(z[s+h] - z[s])^2))

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddW_7sqm_dPCR_Vgram_Cloud.png"),
       plot = paddW_7sqm_Vgram_Cloud_plt,
       height = 3, width=4.5)
############################################################
paddW_7sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                         var=Raw_dPCR, 
                                                         paddock="W", 
                                                         samp_design="c",
                                                         dist_brks=c(0, 60, 150, 300, 350), 
                                                         boot_sim=1000, 
                                                         boot_seed=100)

lowLim_data <- paddW_7sqm_Vgram_data %>%
  select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_data <- paddW_7sqm_Vgram_data %>%
  select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  select(-Dist)

low_and_uppLim_data <- cbind(lowLim_data, uppLim_data) %>%
  mutate(LeftQuantile = factor(LeftQuantile,
                               labels=c("80%", "90%", "95%")))

paddW_7sqm_dPCR_Vgram_Plt <- low_and_uppLim_data %>% ggplot(aes(Dist, SemiVariance)) +
  geom_point(shape = 19, size=3) +
  geom_line(linetype = "dashed") +
  geom_ribbon(aes(ymin = VarLowLim, ymax=VarUppLim, 
                  group= LeftQuantile, 
                  fill = LeftQuantile),
              alpha=0.20) +
  theme_classic() +
  labs(fill = "") +
  scale_fill_manual(values = c("gray30", "gray40", "gray50")) +
  theme(legend.position = c(0.65, 0.92),
        legend.direction="horizontal",
        legend.margin = margin(0,0,0,0 ,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.box.background = element_rect(fill = "white", color="white")) +
  guides(fill = guide_legend(override.aes = 
                               list(alpha=0.7,fill = c("gray40", "gray60", "gray80"))),
         size = guide_legend(order=10)) +
  labs(x = "Distance, h (meter)",
       y = expression(gamma(h))) +
  ylim(c(0, 0.0025))

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddW_7sqm_dPCR_Variogram.png"),
       plot = paddW_7sqm_dPCR_Vgram_Plt,
       height = 3, width=4.5)

paddW_7sqm_dPCR_Vgram_JointPlt <- ggarrange(paddW_7sqm_Vgram_Cloud_plt,
                                            paddW_7sqm_dPCR_Vgram_Plt)



ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddW_7sqm_dPCR_Cloud_and_Variogram.png"),
       plot = paddW_7sqm_dPCR_Vgram_JointPlt,
       height = 3, width=9)

paddW_7sqm_dPCR_Vgram_JointPlt2 <- ggarrange(paddW_7sqm_dPCR_Vgram_JointPlt,
                                             paddW_7sqm_dPCR_Vgram_JointPlt)
ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddW_7sqm_dPCR_Cloud_and_Variogram.png"),
       plot = paddW_7sqm_dPCR_Vgram_JointPlt2,
       height = 3, width=18)
