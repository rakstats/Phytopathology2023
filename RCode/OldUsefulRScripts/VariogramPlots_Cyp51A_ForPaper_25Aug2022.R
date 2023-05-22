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
library(scales)
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
#1.1a Paddock = W; variable = Cyp51A; Sampling design = d                        #
################################################################################
# Create the Variogram Cloud
padd_W_60sqm <- spat_dep_data %>% 
  dplyr::filter(PaddockFac == "W", NestedSubsample == "d")

paddW_Cyp51A_60sqm_Vgram_Cloud <- compute_svgm(padd_W_60sqm, Cyp51ACopies)

paddW_Cyp51A_60sqm_Vgram_Cloud_plt <- paddW_Cyp51A_60sqm_Vgram_Cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha=0.1, size=3, 
             fill="blue", shape = 21, 
             color="black") +
  theme_minimal() +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2),
       subtitle = bquote(atop("Field: W, 60"~m^ 2~"(r)"))) +
  xlim(0, 450)

paddW_Cyp51A_60sqm_Vgram_Cloud_plt
################################################################################
# Compute the semivariogram and associated simulation envelopes                #
################################################################################
paddW_Cyp51A_60sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                          var=Cyp51ACopies, 
                                                          paddock="W", 
                                                          samp_design="d",
                                                          dist_brks=c(0, 75, 150, 250, 350, 500), 
                                                          boot_sim=1000, 
                                                          boot_seed=100)

lowLim_Cyp51A_paddW_60sqm <- paddW_Cyp51A_60sqm_Vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_Cyp51A_paddW_60sqm <- paddW_Cyp51A_60sqm_Vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  select(-Dist)

low_and_uppLim_Cyp51A_paddW_60sqmdata <- cbind(lowLim_Cyp51A_paddW_60sqm, uppLim_Cyp51A_paddW_60sqm) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile,
                                      labels=c("80%", "90%", "95%")))

paddW_60sqm_Cyp51A_Vgram_Plt <- low_and_uppLim_Cyp51A_paddW_60sqmdata %>%
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
       y = expression(gamma(h)),
       subtitle = bquote(atop("Field: W, 60"~m^ 2~"(r)"))) +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400), limits = c(0, 400)) +
  scale_y_continuous(labels = scales::comma, limits = c(150000, 1250000))


paddW_60sqm_Cyp51A_Vgram_Plt
########################################################################################
################################################################################
#1.2a Paddock = Y; variable = Cyp51A; Sampling design = d                        #
################################################################################
# Create the Variogram Cloud
padd_Y_60sqm <- spat_dep_data %>% 
  dplyr::filter(PaddockFac == "Y", NestedSubsample == "d")

paddY_Cyp51A_60sqm_Vgram_Cloud <- compute_svgm(padd_Y_60sqm, Cyp51ACopies)

paddY_Cyp51A_60sqm_Vgram_Cloud_plt <- paddY_Cyp51A_60sqm_Vgram_Cloud %>% 
  ggplot(aes(x=Dist, y=Svar)) +
  geom_point(alpha=0.1, size=3, 
             fill="blue", shape = 21, 
             color="black") +
  theme_minimal() +
  labs(x = "Distance, h (meter)", 
       y = expression(0.5*(z[s+h] - z[s])^2),
       subtitle = bquote(atop("Field: Y, 60"~m^ 2~"(r)"))) +
  xlim(0, 650)

paddY_Cyp51A_60sqm_Vgram_Cloud_plt
################################################################################
# Compute the semivariogram and associated simulation envelopes                #
################################################################################
paddY_Cyp51A_60sqm_Vgram_data <-  compute_semivariogram_envelope(data=spat_dep_data, 
                                                                 var=Cyp51ACopies, 
                                                                 paddock="Y", 
                                                                 samp_design="d",
                                                          dist_brks=c(0, 75, 150, 250, 350, 450, 650), 
                                                                 boot_sim=1000, 
                                                                 boot_seed=100)

lowLim_Cyp51A_paddY_60sqm <- paddY_Cyp51A_60sqm_Vgram_data %>%
  dplyr::select(Dist, SemiVariance, contains("Low")) %>%
  pivot_longer(cols = -c(Dist, SemiVariance), names_to="LeftQuantile", values_to="VarLowLim")

uppLim_Cyp51A_paddY_60sqm <- paddY_Cyp51A_60sqm_Vgram_data %>%
  dplyr::select(Dist, contains("Upp")) %>%
  pivot_longer(cols = -Dist, names_to="RightQuantile", values_to="VarUppLim") %>%
  select(-Dist)

low_and_uppLim_Cyp51A_paddY_60sqmdata <- cbind(lowLim_Cyp51A_paddY_60sqm, uppLim_Cyp51A_paddY_60sqm) %>%
  dplyr::mutate(LeftQuantile = factor(LeftQuantile,
                                      labels=c("80%", "90%", "95%")))

paddY_60sqm_Cyp51A_Vgram_Plt <- low_and_uppLim_Cyp51A_paddY_60sqmdata %>%
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
       y = expression(gamma(h)),
       subtitle = bquote(atop("Field: Y, 60"~m^ 2~"(r)"))) +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500), limits = c(0, 550)) +
  scale_y_continuous(labels = scales::comma, limits = c(25000, 450000))


paddY_60sqm_Cyp51A_Vgram_Plt
########################################################################################
paddW_and_Y_60sqm_Cyp51A_Vgram_Plt <- ggarrange(paddW_60sqm_Cyp51A_Vgram_Plt,
                                                paddY_60sqm_Cyp51A_Vgram_Plt,
                                                ncol = 2)
saveRDS(paddW_and_Y_60sqm_Cyp51A_Vgram_Plt,
        here::here("FilesForLeon_9June2022",
                   "SpatDep", "Figures",
                   "paddW_and_Y_60sqm_Cyp51A_Vgram_Plt.rds"))
############################################################################################
# Final Plot
#############################################################################################
paddW_and_Y_Cyp51A_Vgram_FinalPlt <- ggarrange(paddW_and_Y_1pt2sqm_Cyp51A_Vgram_Plt,
                                               paddW_and_Y_8pt6sqm_Cyp51A_Vgram_Plt,
                                               paddW_and_Y_60sqm_Cyp51A_Vgram_Plt,
                                               nrow = 3)
saveRDS(paddW_and_Y_Cyp51A_Vgram_FinalPlt,
        here::here("FilesForLeon_9June2022",
                   "SpatDep", "Figures",
                   "paddW_and_Y_Cyp51A_Vgram_FinalPlt.rds"))

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "paddock_W_and_Y_Cyp51A_Variogram_FinalPlt.png"),
       plot = paddW_and_Y_Cyp51A_Vgram_FinalPlt,
       height = 12, width=9)
