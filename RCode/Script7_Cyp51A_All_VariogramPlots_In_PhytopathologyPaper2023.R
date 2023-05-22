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
library(patchwork)
################################################################################

# (a)
paddW_1pt2sqm_Cyp51A_Vgram_Plt <- readRDS(here::here("DerivedData", 
                                        "paddW_1pt2sqm_Cyp51A_Vgram_Plt.rds"))

# (b)
paddY_1pt2sqm_Cyp51A_Vgram_Plt <- readRDS(here::here("DerivedData", 
                                    "paddY_1pt2sqm_Cyp51A_Vgram_Plt.rds"))

# (c)
paddW_8pt6sqm_Cyp51A_Vgram_Plt <- readRDS(here::here("DerivedData", 
                                       "paddW_8pt6sqm_Cyp51A_Vgram_Plt.rds"))

# (d)
paddY_8pt6sqm_Cyp51A_Vgram_Plt <- readRDS(here::here("DerivedData", 
                                "paddY_8pt6sqm_Cyp51A_Vgram_Plt.rds"))

# (e)
paddW_60sqm_Cyp51A_Vgram_Plt <- readRDS(here::here("DerivedData", 
                                        "paddW_60sqm_Cyp51A_Vgram_Plt.rds"))

# (f)
paddY_60sqm_Cyp51A_Vgram_Plt <- readRDS(here::here("DerivedData", 
                                    "paddY_60sqm_Cyp51A_Vgram_Plt.rds"))
################################################################################


Cyp51A_Vario_Plots <- ((paddW_1pt2sqm_Cyp51A_Vgram_Plt|paddY_1pt2sqm_Cyp51A_Vgram_Plt)/
          (paddW_8pt6sqm_Cyp51A_Vgram_Plt|paddY_8pt6sqm_Cyp51A_Vgram_Plt)/
           (paddW_60sqm_Cyp51A_Vgram_Plt|paddY_60sqm_Cyp51A_Vgram_Plt))

ggsave(here::here("Figures", "Cyp51A_Variogram_Plots.png"),
       Cyp51A_Vario_Plots, height = 12, width= 12)

ggsave(here::here("Figures", "Cyp51A_Variogram_Plots.jpeg"),
       Cyp51A_Vario_Plots, height = 12, width= 12)

ggsave(here::here("Figures", "Cyp51A_Variogram_Plots.tiff"),
       Cyp51A_Vario_Plots, height = 12, width= 12)
