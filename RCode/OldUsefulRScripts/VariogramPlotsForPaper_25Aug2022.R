################################################################################
# Load required packages                                                       #
################################################################################
library(here)
library(tidyverse)
library(sf)
library(gstat)
library(sp)
library(rlang)
library(ggpubr)
################################################################################
# Load the utility script                                                      #
################################################################################
source(here::here("FilesForLeon_9June2022", "RCode", "Utils_ForPaper.R"))
################################################################################
# 1. Plot 50 sqm data for paddocks W and Y                                     #
################################################################################
################################################################################
# Prepare the Spatial Dependency data                                          #
################################################################################
spat_dep_data <- readRDS(here::here("DerivedData", "spat_dep_data.rds"))
spat_dep_sf_projected <- readRDS(here::here("DerivedData", "spat_dep_sf_projected.rds"))
##################################################################################
# Add the Projected Coordinates to the Spatial Dependency Dataset                #
##################################################################################
spat_dep_data <-  spat_dep_data %>% bind_cols(st_coordinates(spat_dep_sf_projected))
################################################################################
#1.1a Paddock = W; variable = dPCR; Sampling design = d                        #
################################################################################
# Create the Variogram Cloud
padd_W_50sqm <- spat_dep_data %>% 
  dplyr::filter(PaddockFac == "W", NestedSubsample == "d")

paddW_50sqm_Vgram_Cloud <- compute_svgm(padd_W_50sqm, Raw_dPCR)

paddW_50sqm_Vgram_Cloud_plt <- paddW_50sqm_Vgram_Cloud %>% 
                               ggplot(aes(x=Dist, y=Svar)) +
                               geom_point(alpha=0.1, size=3, 
                                          fill="blue", shape = 21, 
                                          color="black") +
                               theme_minimal() +
                               labs(x = "Distance, h (meter)", 
                                    y = expression(0.5*(z[s+h] - z[s])^2),
                                    subtitle = "Field: W") +
                               xlim(0, 450)

paddW_50sqm_Vgram_Cloud_plt
################################################################################
# Compute the semivariogram and associated simulation envelopes                #
################################################################################
paddW_50sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                          var=Raw_dPCR, 
                                                          paddock="W", 
                                                          samp_design="d",
                                                          dist_brks=c(0, 75, 150, 250, 350, 500), 
                                                          boot_sim=1000, 
                                                          boot_seed=100)

lowLim_data_50sqm <- paddW_50sqm_Vgram_data %>%
                     dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_data_50sqm <- paddW_50sqm_Vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  select(-Dist)

low_and_uppLim_50sqmdata <- cbind(lowLim_data_50sqm, uppLim_data_50sqm) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile,
                               labels=c("80%", "90%", "95%")))

paddW_50sqm_dPCR_Vgram_Plt <- low_and_uppLim_50sqmdata %>% ggplot(aes(Dist, SemiVariance)) +
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
  ylim(c(0, 0.07)) + xlim(0, 450)

paddW_50sqm_Vgram_Cloud_plt

paddW_50sqm_dPCR_Vgram_JointPlt <- ggarrange(paddW_50sqm_Vgram_Cloud_plt,
                                             paddW_50sqm_dPCR_Vgram_Plt, 
                                             nrow = 2)
#########################################################################################
################################################################################
#1.2 Paddock = Y; variable = dPCR; Sampling design = d                         #
################################################################################
# Create the Variogram Cloud
padd_Y_50sqm <- spat_dep_data %>% 
  dplyr::filter(PaddockFac == "Y", NestedSubsample == "d")

paddY_50sqm_Vgram_Cloud <- compute_svgm(padd_Y_50sqm, Raw_dPCR)

paddY_50sqm_Vgram_Cloud_plt <- paddY_50sqm_Vgram_Cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha=0.1, size=3, 
             fill="blue", shape = 21, 
             color="black") +
  theme_minimal() +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2),
       subtitle = "Field: Y") +
  xlim(0, 650)

paddY_50sqm_Vgram_Cloud_plt
################################################################################
# Compute the semivariogram and associated simulation envelopes                #
################################################################################
options(scipen=10000)
paddY_50sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                          var=Raw_dPCR, 
                                                          paddock="Y", 
                                                          samp_design="d",
                                              dist_brks=c(0, 75, 150, 250, 350, 450, 650), 
                                                          boot_sim=1000, 
                                                          boot_seed=100)

lowLim_data_50sqm <- paddY_50sqm_Vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_data_50sqm <- paddY_50sqm_Vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  select(-Dist)

low_and_uppLim_50sqmdata <- cbind(lowLim_data_50sqm, uppLim_data_50sqm) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile,
                                      labels=c("80%", "90%", "95%")))

paddY_50sqm_dPCR_Vgram_Plt <- low_and_uppLim_50sqmdata %>% ggplot(aes(Dist, SemiVariance)) +
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
  ylim(0, 0.0004) +
  xlim(0, 650)

paddY_50sqm_Vgram_Cloud_plt
paddY_50sqm_dPCR_Vgram_Plt

paddY_50sqm_dPCR_Vgram_JointPlt <- ggarrange(paddY_50sqm_Vgram_Cloud_plt,
                                             paddY_50sqm_dPCR_Vgram_Plt, 
                                             nrow = 2)

########################################################################################
paddW_50sqm_dPCR_Vgram_FinalPlot <- ggarrange(paddW_50sqm_dPCR_Vgram_JointPlt,
                                              paddY_50sqm_dPCR_Vgram_JointPlt,
                                              ncol = 2)

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddock_W_and_Y_50sqm_dPCR_Cloud_and_Variogram.png"),
       plot = paddW_50sqm_dPCR_Vgram_FinalPlot,
       height = 8, width=9)

##############################################################################################
################################################################################
# 2. Plot 1.2 sqm and 7 sqm data for paddock: W                                #
################################################################################
# Paddock: W --- 1.2 sqm
# Create the Variogram Cloud
padd_W_1pt2sqm <- spat_dep_data %>% 
  dplyr::filter(PaddockFac == "W", NestedSubsample == "b")

paddW_1pt2sqm_Vgram_Cloud <- compute_svgm(padd_W_1pt2sqm, Raw_dPCR)

paddW_1pt2sqm_Vgram_Cloud_plt <- paddW_1pt2sqm_Vgram_Cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha=0.1, size=3, 
             fill="blue", shape = 21, 
             color="black") +
  theme_minimal() +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2),
       subtitle = bquote(atop("Field: W, 1.2"~m^ 2~"(r)"))) +
  xlim(0, 350)

paddW_1pt2sqm_Vgram_Cloud_plt


paddW_1pt2sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                            var=Raw_dPCR, 
                                                            paddock="W", 
                                                            samp_design="b",
                                                dist_brks=c(0, 60, 150, 250, 350), 
                                                            boot_sim=1000, 
                                                            boot_seed=100)

lowLim_data_1pt2sqm <- paddW_1pt2sqm_Vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_data_1pt2sqm <- paddW_1pt2sqm_Vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  select(-Dist)

low_and_uppLim_1pt2sqmdata <- cbind(lowLim_data_1pt2sqm, uppLim_data_1pt2sqm) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile,
                                      labels=c("80%", "90%", "95%")))

paddW_1pt2sqm_dPCR_Vgram_Plt <- low_and_uppLim_1pt2sqmdata %>% 
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
  ylim(0, 0.004) +
  xlim(0, 350)

paddW_1pt2sqm_Vgram_Cloud_plt
paddW_1pt2sqm_dPCR_Vgram_Plt

paddW_1pt2sqm_dPCR_Vgram_JointPlt <- ggarrange(paddW_1pt2sqm_Vgram_Cloud_plt,
                                               paddW_1pt2sqm_dPCR_Vgram_Plt, 
                                               nrow = 2)
paddW_1pt2sqm_dPCR_Vgram_JointPlt
##################################################################################
# Paddock: W --- 8.6 sqm                                                         #
# Create the Variogram Cloud                                                     #
##################################################################################
padd_W_8pt6sqm <- spat_dep_data %>% 
  dplyr::filter(PaddockFac == "W", NestedSubsample == "c")

paddW_8pt6sqm_Vgram_Cloud <- compute_svgm(padd_W_8pt6sqm, Raw_dPCR)

paddW_8pt6sqm_Vgram_Cloud_plt <- paddW_8pt6sqm_Vgram_Cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha=0.1, size=3, 
             fill="blue", shape = 21, 
             color="black") +
  theme_minimal() +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2),
       subtitle = bquote(atop("Field: W, 8.6"~m^ 2~"(r)"))) +
  xlim(0, 350)

paddW_8pt6sqm_Vgram_Cloud_plt


paddW_8pt6sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                            var=Raw_dPCR, 
                                                            paddock="W", 
                                                            samp_design="c",
                                                            dist_brks=c(0, 70, 150, 300, 350), 
                                                            boot_sim=1000, 
                                                            boot_seed=100)

lowLim_data_8pt6sqm <- paddW_8pt6sqm_Vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_data_8pt6sqm <- paddW_8pt6sqm_Vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  select(-Dist)

low_and_uppLim_8pt6sqmdata <- cbind(lowLim_data_8pt6sqm, uppLim_data_8pt6sqm) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile,
                                      labels=c("80%", "90%", "95%")))

paddW_8pt6sqm_dPCR_Vgram_Plt <- low_and_uppLim_8pt6sqmdata %>% 
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
  ylim(0, 0.003) +
  xlim(0, 350)

paddW_8pt6sqm_Vgram_Cloud_plt
paddW_8pt6sqm_dPCR_Vgram_Plt

paddW_8pt6sqm_dPCR_Vgram_JointPlt <- ggarrange(paddW_8pt6sqm_Vgram_Cloud_plt,
                                               paddW_8pt6sqm_dPCR_Vgram_Plt, 
                                               nrow = 2)
paddW_8pt6sqm_dPCR_Vgram_JointPlt
##################################################################################
paddW_1pt2sqm_and_8pt6sqm_dPCR_Vgram_FinalPlot <- ggarrange(paddW_1pt2sqm_dPCR_Vgram_JointPlt,
                                                paddW_8pt6sqm_dPCR_Vgram_JointPlt,
                                                ncol = 2)

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddock_W_1pt2_and_8pt6_sqm_dPCR_Cloud_and_Variogram.png"),
       plot = paddW_1pt2sqm_and_8pt6sqm_dPCR_Vgram_FinalPlot,
       height = 8, width=9)
####################################################################################
################################################################################
# 3. Plot 1.2 sqm and 7 sqm data for paddock: Y                                #
################################################################################
# Paddock: Y --- 1.2 sqm
# Create the Variogram Cloud
padd_Y_1pt2sqm <- spat_dep_data %>% 
  dplyr::filter(PaddockFac == "Y", NestedSubsample == "b")

paddY_1pt2sqm_Vgram_Cloud <- compute_svgm(padd_Y_1pt2sqm, Raw_dPCR)

paddY_1pt2sqm_Vgram_Cloud_plt <- paddY_1pt2sqm_Vgram_Cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha=0.1, size=3, 
             fill="blue", shape = 21, 
             color="black") +
  theme_minimal() +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2),
       subtitle = bquote(atop("Field: Y, 1.2"~m^ 2~"(r)"))) +
  xlim(0, 550)

paddY_1pt2sqm_Vgram_Cloud_plt


paddY_1pt2sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                            var=Raw_dPCR, 
                                                            paddock="Y", 
                                                            samp_design="b",
                                                            dist_brks=c(0, 60, 200, 425, 550), 
                                                            boot_sim=1000, 
                                                            boot_seed=100)

lowLim_data_1pt2sqm <- paddY_1pt2sqm_Vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_data_1pt2sqm <- paddY_1pt2sqm_Vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  select(-Dist)

low_and_uppLim_1pt2sqmdata <- cbind(lowLim_data_1pt2sqm, uppLim_data_1pt2sqm) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile,
                                      labels=c("80%", "90%", "95%")))

paddY_1pt2sqm_dPCR_Vgram_Plt <- low_and_uppLim_1pt2sqmdata %>% 
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
  ylim(0, 0.005) +
  xlim(0, 550)

paddY_1pt2sqm_Vgram_Cloud_plt
paddY_1pt2sqm_dPCR_Vgram_Plt

paddY_1pt2sqm_dPCR_Vgram_JointPlt <- ggarrange(paddY_1pt2sqm_Vgram_Cloud_plt,
                                               paddY_1pt2sqm_dPCR_Vgram_Plt, 
                                               nrow = 2)
paddY_1pt2sqm_dPCR_Vgram_JointPlt
##################################################################################
# Paddock: W --- 8.6 sqm                                                         #
# Create the Variogram Cloud                                                     #
##################################################################################
padd_Y_8pt6sqm <- spat_dep_data %>% 
  dplyr::filter(PaddockFac == "Y", NestedSubsample == "c")

paddY_8pt6sqm_Vgram_Cloud <- compute_svgm(padd_Y_8pt6sqm, Raw_dPCR)

paddY_8pt6sqm_Vgram_Cloud_plt <- paddY_8pt6sqm_Vgram_Cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha=0.1, size=3, 
             fill="blue", shape = 21, 
             color="black") +
  theme_minimal() +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2),
       subtitle = bquote(atop("Field: Y, 8.6"~m^ 2~"(r)"))) +
  xlim(0, 550)

paddY_8pt6sqm_Vgram_Cloud_plt


paddY_8pt6sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                            var=Raw_dPCR, 
                                                            paddock="Y", 
                                                            samp_design="c",
                                                            dist_brks=c(0, 90, 200, 430, 550), 
                                                            boot_sim=1000, 
                                                            boot_seed=100)

lowLim_data_8pt6sqm <- paddY_8pt6sqm_Vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_data_8pt6sqm <- paddY_8pt6sqm_Vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  select(-Dist)

low_and_uppLim_8pt6sqmdata <- cbind(lowLim_data_8pt6sqm, uppLim_data_8pt6sqm) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile,
                                      labels=c("80%", "90%", "95%")))

paddY_8pt6sqm_dPCR_Vgram_Plt <- low_and_uppLim_8pt6sqmdata %>% 
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
  ylim(0, 0.0008) +
  xlim(0, 550)

paddY_8pt6sqm_Vgram_Cloud_plt
paddY_8pt6sqm_dPCR_Vgram_Plt

paddY_8pt6sqm_dPCR_Vgram_JointPlt <- ggarrange(paddY_8pt6sqm_Vgram_Cloud_plt,
                                               paddY_8pt6sqm_dPCR_Vgram_Plt, 
                                               nrow = 2)
paddY_8pt6sqm_dPCR_Vgram_JointPlt
##################################################################################
paddY_1pt2sqm_and_8pt6sqm_dPCR_Vgram_FinalPlot <- ggarrange(paddY_1pt2sqm_dPCR_Vgram_JointPlt,
                                                            paddY_8pt6sqm_dPCR_Vgram_JointPlt,
                                                            ncol = 2)

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddock_Y_1pt2_and_8pt6_sqm_dPCR_Cloud_and_Variogram.png"),
       plot = paddY_1pt2sqm_and_8pt6sqm_dPCR_Vgram_FinalPlot,
       height = 8, width=9)

