################################################################################
# Load required packages                                                       #  
################################################################################
library(here)
library(readxl)
library(tidyverse)
library(sp)
library(spatstat)
library(plotly)
library(ggpubr)
library(rgdal)
library(maptools)
library(sf)
library(asreml)
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
################################################################################
# Source utility functions                                                     =
################################################################################
source(here::here("FilesForLeon_9June2022", "RCode", "Utils_ForPaper.R"))      #
################################################################################
# Analysis of Unit Size data                                                   #
################################################################################
unit_data <- readRDS(here::here("DerivedData", "unit_data.rds"))
glimpse(unit_data)
################################################################################
# We should remove all data points with Cyp51A_copies <= 29                    #    
################################################################################
unit_data <- unit_data %>% dplyr::filter(Cyp51A_copies > 29)
unit_data2 <- unit_data %>%
  mutate(PaddockFac = factor(PaddockFac,
                             levels = c("M", "P", "V"),
                             labels = paste0("Field-",
                                             c("M", "P", "V"))))
################################################################################
# 1. Boxplot: sqrt(dPCR) against NestedSample                                  #
################################################################################
# Hexadecimal color specification 
display.brewer.pal(n=4, name="Dark2")
color_spec <- brewer.pal(n = 4, name = "Dark2")

color_vals <- c("a" = color_spec[1], "b"= color_spec[2],
                "c" = color_spec[3], "d" = color_spec[4])
lab_vals  <- c(bquote(bold(atop("1.2"~m^ 2, "(hp)"))),
               bquote(bold(atop("1.2"~m^ 2, "(r)"))),
               bquote(bold(atop("8.6"~m^ 2, "(r)"))),
               bquote(bold(atop("60"~m^ 2, "(r)"))))

brk_vals <- c("a", "b", "c", "d")
################################################################################
sqrt_dPCR_boxplt <- plot_dPcr_boxplt(data=unit_data2, 
                                     varname=sqrt_dPCR, 
                                     label="sqrt(dPCR)")
sqrt_dPCR_boxplt2 <- sqrt_dPCR_boxplt + 
                     scale_color_manual(values = color_vals) +
                     scale_x_discrete(breaks = brk_vals, labels = lab_vals) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6), labels = c("0%", "20%", "40%", "60%")) +
                     theme(axis.text.x = element_text(size=12, face="bold"),
                           axis.text.y = element_text(size=12, face="bold"),
                           axis.title = element_text(size=13),
                           axis.title.x = element_text(margin=margin(t = 15, r = 0, b = 0, l = 0)),
                           axis.title.y = element_text(margin=margin(t = 0, r = 13, b = 0, l = 0))) +
                    labs(y = expression(paste("sqrt(", italic(Ptt), "-FR frequency (%))")))

ggsave(here::here("FilesForLeon_9June2022",
                 "UnitSize", "Figures",
                 "boxplt_sqrt_dPCR.png"), plot=sqrt_dPCR_boxplt2,
      width=8, height=5)  
###############################################################################
################################################################################
# 2. Boxplot: log(Cyp51A copies/mg) against NestedSample                       #
################################################################################
log_Cyp51_boxplt <- plot_dPcr_boxplt(data=unit_data2, 
                                     varname=log_Cyp51A_copies, 
                                     label="log(Cyp51A copies/mg)")

log_Cyp51_boxplt2 <- log_Cyp51_boxplt + 
  scale_color_manual(values = color_vals) +
  scale_x_discrete(breaks = brk_vals,
                   labels = lab_vals) +
  theme(axis.text.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12, face="bold"),
        axis.title = element_text(size=13),
        axis.title.x = element_markdown(margin=margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_markdown()) +
        labs(y = "<i>Ptt</i> density log(<i>Cyp51A</i> copies mg<sup>-1</sup>)")

ggsave(here::here("FilesForLeon_9June2022",
                  "UnitSize", "Figures",
                  "boxplt_log_Cyp51.png"), plot=log_Cyp51_boxplt2,
       width=8, height=5)  
################################################################################
################################################################################
# 3. Boxplot: Raw dPCR against NestedSample                                    #
################################################################################
dPCR_boxplt <- plot_dPcr_boxplt(data=unit_data2, 
                                varname=Raw_dPCR, 
                                label="Raw dPCR")

dPCR_boxplt2 <- dPCR_boxplt + 
  scale_color_manual(values = color_vals) +
  scale_x_discrete(breaks = brk_vals, 
                   labels = lab_vals) +
  scale_y_continuous(limits = c(0, 0.50),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), 
                     labels = c("0%", "10%", "20%", "30%", "40%", "50%")) +
  theme(axis.text.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12, face="bold"),
        axis.title = element_text(size=13),
        axis.title.x = element_text(margin=margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin=margin(t = 0, r = 13, b = 0, l = 0))) +
  labs(y = expression(paste(italic(Ptt), "-FR frequency (%)")))

ggsave(here::here("FilesForLeon_9June2022",
                  "UnitSize", "Figures",
                  "boxplt_dPCR.png"), plot=dPCR_boxplt2,
       width=8, height=5)  
################################################################################
# 4. Boxplot: Cyp51A copies/mg against NestedSample                            #
################################################################################
cyp51A_boxplt <- plot_dPcr_boxplt(data=unit_data2, 
                                  varname=Cyp51A_copies, 
                                  label="Cyp51A copies/mg")

cyp51A_boxplt2 <- cyp51A_boxplt + 
  scale_color_manual(values = color_vals) +
  scale_x_discrete(breaks = brk_vals, labels = lab_vals) +
  theme(axis.text.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12, face="bold"),
        axis.title = element_text(size=13),
        axis.title.y = element_markdown()) +
  labs(y = "<i>Ptt</i> density (<i>Cyp51A</i> copies mg<sup>-1</sup>)")

ggsave(here::here("FilesForLeon_9June2022",
                  "UnitSize", "Figures",
                  "boxplt_Cyp51Copies.png"), plot=cyp51A_boxplt2,
       width=8, height=5)  
################################################################################
# 5. Var(Ptt-FR frequency) or Var(dPCR) against Nested Sample                  #                                      
################################################################################
lab_vals2  <- c(bquote(bold(atop("1.2"~m^ 2~"(hp)"))),
                bquote(bold(atop("1.2"~m^ 2~"(r)"))),
                bquote(bold(atop("8.6"~m^ 2~"(r)"))),
                bquote(bold(atop("60"~m^ 2~"(r)"))))

unit_data_summary <- unit_data2 %>%
                     group_by(PaddockFac, NestedSubsample) %>%
                     summarise(count = n(), 
                               SD_RawdPCR = sd(Raw_dPCR),
                               SD_sqrtdPCR = sd(sqrt_dPCR),
                               SD_DNA_copies = sd(Cyp51A_copies),
                               SD_logDNA_copies = sd(log_Cyp51A_copies), .groups = "drop")

sd_plot_rawdPcr <- unit_data_summary %>%
                   ggplot(aes(x = factor(PaddockFac, level = c("Field-V", "Field-M", "Field-P")),
                              y = SD_RawdPCR)) +
                   geom_point(aes(color = NestedSubsample), size=7) +
                   geom_point(shape=1, color = "white", size=7.1) +
                   theme_classic() +
                   theme(legend.position = "top",
                         axis.text.x = element_text(size=12, face="bold"),
                         axis.text.y = element_text(size=12, face="bold"),
                         axis.title.y = element_text(size=13, margin=margin(t = 0, r = 13, b = 0, l = 0)),
                         axis.title.x = element_blank(),
                         legend.text = element_text(size=14)) +
                   scale_color_manual(values = color_vals,
                                      breaks = brk_vals, 
                                      labels = lab_vals2) +
                   guides(color = guide_legend(override.aes = list(size = 4))) +
                   scale_y_continuous(breaks = c(0.04, 0.08, 0.12, 0.16),
                                      labels = c("4%", "8%", "12%", "16%")) +
                   labs(y = expression(paste("Standard deviation of ", italic(Ptt), "-FR frequency (%)")),
                        color = "Sampling strategy")

ggsave(here::here("FilesForLeon_9June2022",
                  "UnitSize", "Figures",
                  "dotplot_sd_dPCR.png"), plot=sd_plot_rawdPcr,
       width=7, height=5)                     
##################################################################################
# 6. Spatial distribution of dPCR or Ptt-FR frequency for the Unit Size data     #
##################################################################################
unit_data3 <- unit_data2 %>%
  mutate(NestedSubsample = factor(NestedSubsample,
                                  labels = lab_vals2))
###################################################################
unit_sf <- st_as_sf(unit_data3, coords = c("lon", "lat"), crs = 4326)
unit_sf_projected <- st_transform(unit_sf, 3112)
# If jittering needed, then run the next two lines:
set.seed(123)
unit_sf_projected <- st_jitter(unit_sf_projected, factor = 0.08)

# Otherwise
unit_coords0 <- st_geometry(unit_sf_projected)
unit_bbox0 <- st_bbox(unit_coords0)



unit_coords1 <- unit_coords0 - c(unit_bbox0['xmin'], unit_bbox0['ymin'])
# Check the transformed coordinates
unit_coords1 %>% ggplot() + geom_sf()

unit_sf_projected2 <- st_set_geometry(unit_sf_projected, unit_coords1)


summary(unit_sf_projected2$Raw_dPCR)

dPCR_spplot <- unit_sf_projected2 %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=3) +
  geom_sf(shape = 1, size=3.1) +
  coord_sf() +
  facet_grid(~NestedSubsample, 
             labeller = label_parsed) +
  # scale_color_viridis_c(option = "viridis", direction = -1) + 
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
                        limits = c(0, 0.50),
                        labels = c("0%", "10%", "20%", "30%", "40%", "50%"),
                        guide = guide_colorbar(barwidth = unit(8, "cm"))) +
  #scale_color_continuous(colours = rev(terrain.colors(126))) +
  theme_classic() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size=12, margin = margin(8,0,0,0)),
        axis.text = element_text(face="bold"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(face="bold")) +
  labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")))

ggsave(here::here("FilesForLeon_9June2022",
                  "UnitSize", "Figures",
                  "dPCR_wojitter_spatial_plot.png"), plot=dPCR_spplot,
       width=10, height=5)
################################################################################
# This is an advanced version where the sizes of the points are varied based   # 
# on the values of the Ptt-FR frequency (dPCR) values                          #
################################################################################
dPCR_spplot2 <- unit_sf_projected2 %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR, size = Raw_dPCR)) +
  geom_sf(aes(size = Raw_dPCR), shape = 1) +
  coord_sf() +
  facet_grid(~NestedSubsample, 
             labeller = label_parsed) +
  # scale_color_viridis_c(option = "viridis", direction = -1) + 
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
                        limits = c(0, 0.50),
                        labels = c("0%", "10%", "20%", "30%", "40%", "50%"),
                        guide = guide_colorbar(barwidth = unit(8, "cm"))) +
  #scale_color_continuous(colours = rev(terrain.colors(126))) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size=12, margin = margin(8,0,0,0)),
        axis.text = element_text(face="bold"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(face="bold")) +
  labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")),
       size = '') +
  scale_size_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4),
                        labels = c("0%", "10%", "20%", "30%", "40%"),
                        guide = guide_legend(override.aes = list(shape = 1)))

ggsave(here::here("FilesForLeon_9June2022",
                  "UnitSize", "Figures",
                  "dPCR_wojitter_spatial_plotV2.png"), plot=dPCR_spplot2,
       width=10, height=5)
################################################################################
# library(grid)
# gt = ggplot_gtable(ggplot_build(d))
# gt$widths[4] = 0.5*gt$widths[4]
# grid.draw(gt)
################################################################################
#  g <- ggplotGrob(d)
#  g$heights[[16]] = unit(0,"cm")
# # 
# # library(grid)
#  grid.newpage()
#  grid.draw(g)
################################################################################
################################################################################
# 7. Spatial distribution of dPCR  for Spatial Dependency Data                 #
################################################################################
spat_dep_data0 <- read_excel(here::here("Excel_Data_sets_Leon",
                            "Master_2021_US_SD_Data_Stubble_Copied.xlsx"),
                            sheet = 1)
glimpse(spat_dep_data0)

spat_dep_data <- spat_dep_data0 %>%
  dplyr::filter(SampleType == "Spatial dependency") %>%
  transmute(lon=longitude,
            lat = latitude,
            PaddockFac = factor(Paddock),
            SampleType = factor(SampleType),
            ReplicationNo = factor(ReplicationNo),
            SampleNumber = factor(SampleNumber),
            NestedSubsample = factor(NestedSubsample, 
                                     levels = c("(1m)", "(7m)", "(50m)"),
                                     labels = c("b", "c", "d")),
            BarleyCultivar = factor(BarleyCultivar),
            SampleUnitSize_m2  = factor(formatC(SampleUnitSize_m2, digits=1, format="f")),
            Raw_dPCR = Raw_dPCR,
            logRaw_dPCR = log(Raw_dPCR + 0.0000001),
            sqrtRaw_dPCR = sqrt(Raw_dPCR),
            Cyp51ACopies = Cyp51ACopies,
            log_Cyp51ACopies =  log(Cyp51ACopies),
            VicCopies = VicCopies,
            FamCopies = FamCopies,
            VicCounts = VicCounts,
            FamAndVicCounts = FamAndVicCounts) %>%
            dplyr::filter(Cyp51ACopies > 29)

glimpse(spat_dep_data)
saveRDS(spat_dep_data, here::here("DerivedData", "spat_dep_data.rds"))
################################################################################
spat_dep_sf <- st_as_sf(spat_dep_data,
                        coords = c("lon", "lat"),
                        crs = 4326)
spat_dep_sf_projected <- st_transform(spat_dep_sf, 3112) 

sp_dep_coords0 <- st_geometry(spat_dep_sf_projected)
sp_dep_bbox0 <- st_bbox(sp_dep_coords0)

sp_dep_coords1 <- sp_dep_coords0 - c(sp_dep_bbox0['xmin'], sp_dep_bbox0['ymin'])
plot(sp_dep_coords1)
sp_dep_coords1 %>% ggplot() + geom_sf()

spat_dep_sf_projected2 <- st_set_geometry(spat_dep_sf_projected,
                                          sp_dep_coords1)
################################################################################
# Create a polygon around the points                                           #
################################################################################
# Rep-1: 1.2 m2 polygons #
##########################
rep1_pol = st_polygon(
  list(
    cbind(
      c(220, 250, 250, 220, 220), 
      c(335, 335, 365, 365, 335))
  )
)
rep1_pol1 <- st_sfc(rep1_pol, crs=st_crs(spat_dep_sf_projected2))
rep1_poly_data <- tibble(NestedSubsample = rep("a", 5),
                         geometry = rep1_pol1)
rep1_poly_data <- st_sf(rep1_poly_data)
#rep1_poly_data %>% ggplot() + geom_sf(color="red")
rep1_poly_data <- rep1_poly_data %>%
                  mutate(NestedSubsample = factor(NestedSubsample,
                                                  labels = lab_vals2[c(2)]))
##########################
# Rep-2: 1.2 m2 polygons #
##########################
rep2_pol = st_polygon(
  list(
    cbind(
      c(155, 185, 185, 155, 155), 
      c(100, 100, 130, 130, 100))
  )
)
rep2_pol2 <- st_sfc(rep2_pol, crs=st_crs(spat_dep_sf_projected2))
rep2_poly_data <- tibble(NestedSubsample = rep("a", 5),
                         geometry = rep2_pol2)
rep2_poly_data <- st_sf(rep2_poly_data)
rep2_poly_data <- rep2_poly_data %>%
  mutate(NestedSubsample = factor(NestedSubsample,
                                  labels = lab_vals2[c(2)]))

##########################
# Rep-3: 1.2 m2 polygons #
##########################
rep3_pol = st_polygon(
  list(
    cbind(
      c(182, 212, 212, 182, 182), 
      c(10,   10, 40,  40, 10))
  )
)
rep3_pol3 <- st_sfc(rep3_pol, crs=st_crs(spat_dep_sf_projected2))
rep3_poly_data <- tibble(NestedSubsample = rep("a", 5),
                         geometry = rep3_pol3)
rep3_poly_data <- st_sf(rep3_poly_data)
rep3_poly_data <- rep3_poly_data %>%
  mutate(NestedSubsample = factor(NestedSubsample,
                                  labels = lab_vals2[c(2)]))

################################################################################
# Create a polygon around the points 8.6 m2                                    #
################################################################################
# Rep-1: 8.6 m2 polygons #
##########################
xmin1 <- 195;xmax1 <- 275;ymin1 <- 325;ymax1 <- 375

rep1_8.6m2_pol = st_polygon(
  list(
    cbind(
      c(xmin1, xmax1, xmax1, xmin1, xmin1), 
      c(ymin1, ymin1, ymax1, ymax1, ymin1))
  )
)
rep1_8.6m2_pol1 <- st_sfc(rep1_8.6m2_pol, crs=st_crs(spat_dep_sf_projected2))
rep1_8.6m2_poly_data <- tibble(NestedSubsample = rep("c", 5),
                               geometry = rep1_8.6m2_pol1)
rep1_8.6m2_poly_data <- st_sf(rep1_8.6m2_poly_data)

rep1_8.6m2_poly_data <- rep1_8.6m2_poly_data %>%
  mutate(NestedSubsample = factor(NestedSubsample,
                                  labels = lab_vals2[c(3)]))

##########################
# Rep-2: 8.6 m2 polygons #
##########################
xmin2 <- 130;xmax2 <- 210;ymin2 <- 88;ymax2 <- 138

rep2_8.6m2_pol = st_polygon(
  list(
    cbind(
      c(xmin2, xmax2, xmax2, xmin2, xmin2), 
      c(ymin2, ymin2, ymax2, ymax2, ymin2))
  )
)
rep2_8.6m2_pol1 <- st_sfc(rep2_8.6m2_pol, crs=st_crs(spat_dep_sf_projected2))
rep2_8.6m2_poly_data <- tibble(NestedSubsample = rep("c", 5),
                               geometry = rep2_8.6m2_pol1)
rep2_8.6m2_poly_data <- st_sf(rep2_8.6m2_poly_data)

rep2_8.6m2_poly_data <- rep2_8.6m2_poly_data %>%
  mutate(NestedSubsample = factor(NestedSubsample,
                                  labels = lab_vals2[c(3)]))

##########################
# Rep-3: 8.6 m2 polygons #
##########################
xmin3 <- 155;xmax3 <- 235;ymin3 <- 0;ymax3<- 50

rep3_8.6m2_pol = st_polygon(
  list(
    cbind(
      c(xmin3, xmax3, xmax3, xmin3, xmin3), 
      c(ymin3, ymin3, ymax3, ymax3, ymin3))
  )
)
rep3_8.6m2_pol1 <- st_sfc(rep3_8.6m2_pol, crs=st_crs(spat_dep_sf_projected2))
rep3_8.6m2_poly_data <- tibble(NestedSubsample = rep("c", 5),
                               geometry = rep3_8.6m2_pol1)
rep3_8.6m2_poly_data <- st_sf(rep3_8.6m2_poly_data)

rep3_8.6m2_poly_data <- rep3_8.6m2_poly_data %>%
  mutate(NestedSubsample = factor(NestedSubsample,
                                  labels = lab_vals2[c(3)]))
###############################################################################
rep_poly_data <- rep1_poly_data %>% 
                    bind_rows(rep2_poly_data) %>%
                      bind_rows(rep3_poly_data) %>%
                       bind_rows(rep1_8.6m2_poly_data) %>%
                         bind_rows(rep2_8.6m2_poly_data) %>%
                            bind_rows(rep3_8.6m2_poly_data) 
################################################################################
spat_dep_sf_projected3 <- spat_dep_sf_projected2 %>%
                          mutate(NestedSubsample = factor(NestedSubsample,
                                                          labels = lab_vals2[c(2,3,4)]))
summary(spat_dep_sf_projected3$Raw_dPCR)

saveRDS(spat_dep_sf_projected2)
##################################################################################
#                                                                                #
# PADDOCK W: SPATIAL DISTN Of dPCR (Ptt-FR freq %) for Spatial-Dependency Data   #
#                                                                                #
##################################################################################
dPCR_paddW_spplot <- spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W") %>%
  ggplot() +
  geom_sf(data = rep_poly_data, color = "red", fill = "white") +
  geom_sf(aes(color = Raw_dPCR), size=4) +
  geom_sf(shape = 1, size=4.1) +
  coord_sf() +
  theme_bw() +
  facet_grid(NestedSubsample ~ ., labeller = label_parsed) + # switch = "y"
  theme(legend.position = "top",
        strip.text.y = element_text(size=13, margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(face="bold", size=12)) +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.15, 0.30, 0.45),
                        limits = c(0, 0.45),
                        labels = c("0%", "15%", "30%", "45%"),
                        guide = guide_colorbar(barwidth = unit(6, "cm"))) +
  labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")))

################################################################################
# Saving dPCR_paddW_spplot                                                     #
################################################################################
ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "dPCR_paddW_spplot.png"), 
       plot=dPCR_paddW_spplot,
       width=5, height=12)

# Above plot without the boxes around the points
dPCR_paddW_spplot_wo_boxes <- spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "W") %>%
  ggplot() +
  #geom_sf(data = rep_poly_data, color = "red", fill = "white") +
  geom_sf(aes(color = Raw_dPCR), size=4) +
  geom_sf(shape = 1, size=4.1) +
  coord_sf() +
  theme_bw() +
  facet_grid(NestedSubsample ~ ., labeller = label_parsed) + # switch = "y"
  theme(legend.position = "top",
        strip.text.y = element_text(size=13, margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(face="bold", size=12)) +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.15, 0.30, 0.45),
                        limits = c(0, 0.45),
                        labels = c("0%", "15%", "30%", "45%"),
                        guide = guide_colorbar(barwidth = unit(6, "cm"))) +
  labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")))

################################################################################
# Saving dPCR_paddW_spplot                                                     #
################################################################################
ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "dPCR_PADDOCK_W_FINAL_SPPLOT_0.png"), 
       plot=dPCR_paddW_spplot_wo_boxes,
       width=5, height=12)
################################################################################
# PADDOCK-W -- 1.2m2 REPLICATIONS                                              #
################################################################################
# Paddock: W, NestedSubsample: 1.2 m2 (rake); REP-1                            #
################################################################################
rep1_dPCR_paddW_1m_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "W" & NestedSubsample == "b" & ReplicationNo == 1) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=9) +
  geom_sf(shape = 1, size=9.1, stroke=1.5) +
  coord_sf(xlim = c(228, 240), ylim = c(350, 352),
           label_axes = list(bottom = "E", right = "N")) +
  #coord_sf(xlim = c(220, 250), ylim = c(335, 365)) +
  theme_bw() +
  theme(axis.text = element_text(face="bold", size=20),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.15, 0.30, 0.45),
                        limits = c(0, 0.45),
                        labels = c("0%", "15%", "30%", "45%")) 
  #labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")))
################################################################################
# Paddock: W, NestedSubsample: 1.2 m2 (rake); REP-2                            #
################################################################################
rep2_dPCR_paddW_1m_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "W" & NestedSubsample == "b" & ReplicationNo == 2) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=9) +
  geom_sf(shape = 1, size=9.1, stroke=1.5) +
  coord_sf(xlim = c(166, 178), ylim = c(112.8, 114.8),
           label_axes = list(bottom = "E", right = "N")) +
  #coord_sf(xlim = c(155, 185), ylim = c(100, 130)) +
  theme_bw() +
  theme(axis.text = element_text(face="bold", size = 20),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.15, 0.30, 0.45),
                        limits = c(0, 0.45),
                        labels = c("0%", "15%", "30%", "45%"))
################################################################################
# Paddock: W, NestedSubsample: 1.2 m2 (rake); REP-3                            #
################################################################################
rep3_dPCR_paddW_1m_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "W" & NestedSubsample == "b" & ReplicationNo == 3) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=9) +
  geom_sf(shape = 1, size=9.1, stroke=1.5) +
  coord_sf(xlim = c(192, 204), ylim = c(24.2, 26.2),
           label_axes = list(bottom = "E", right = "N")) +
  theme_bw() +
  theme(axis.text = element_text(face="bold", size=20),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.15, 0.30, 0.45),
                        limits = c(0, 0.45),
                        labels = c("0%", "15%", "30%", "45%"))
################################################################################
# Paddock: W, NestedSubsample: 1.2 m2 (rake); Combined 3 Reps                  #
################################################################################
allreps_dPCR_paddW_1m_spplot <- ggarrange(rep1_dPCR_paddW_1m_spplot,
                                          rep2_dPCR_paddW_1m_spplot,
                                          rep3_dPCR_paddW_1m_spplot, nrow = 3)
################################################################################
# Saving plot of Combined 3 Reps for Paddock: W, NestedSubsample: 1.2 m2 (rake)#
################################################################################
ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "allreps_dPCR_paddW_1m_spplotV2.png"), plot=allreps_dPCR_paddW_1m_spplot,
       width=10, height=6)
################################################################################
################################################################################
################################################################################
# PADDOCK-W -- 8.6m2 REPLICATIONS                                              #
################################################################################
# Paddock: W, NestedSubsample: 8.6 m2 (rake); REP-1                            #
################################################################################
rep1_dPCR_paddW_8.6m2_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "W" & NestedSubsample == "c" & ReplicationNo == 1) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=9) +
  geom_sf(shape = 1, size=9.1, stroke=1.5) +
  coord_sf(xlim = c(205, 265), ylim = c(346, 356),
           label_axes = list(bottom = "E", right = "N")) +
  theme_bw() +
  theme(axis.text = element_text(face="bold", size=20),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.15, 0.30, 0.45),
                        limits = c(0, 0.45),
                        labels = c("0%", "15%", "30%", "45%")) 
################################################################################
# Paddock: W, NestedSubsample: 8.6 m2 (rake); REP-2                            #
################################################################################
rep2_dPCR_paddW_8.6m2_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "W" & NestedSubsample == "c" & ReplicationNo == 2) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=9) +
  geom_sf(shape = 1, size=9.1, stroke=1.5) +
  coord_sf(xlim = c(140, 200), ylim = c(108.5, 118.5),
           label_axes = list(bottom = "E", right = "N")) +
  theme_bw() +
  theme(axis.text = element_text(face="bold", size = 20),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.15, 0.30, 0.45),
                        limits = c(0, 0.45),
                        labels = c("0%", "15%", "30%", "45%"))
################################################################################
# Paddock: W, NestedSubsample: 8.6 m2 (rake); REP-3                            #
################################################################################
rep3_dPCR_paddW_8.6m2_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "W" & NestedSubsample == "c" & ReplicationNo == 3) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=9) +
  geom_sf(shape = 1, size=9.1, stroke=1.5) +
  coord_sf(xlim = c(165, 225), ylim = c(20, 30),
           label_axes = list(bottom = "E", right = "N")) +
  theme_bw() +
  theme(axis.text = element_text(face="bold", size=20),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.15, 0.30, 0.45),
                        limits = c(0, 0.45),
                        labels = c("0%", "15%", "30%", "45%"))
################################################################################
# Paddock: W, NestedSubsample: 8.6 m2 (rake); Combined 3 Reps                  #
################################################################################
allreps_dPCR_paddW_8.6m2_spplot <- ggarrange(rep1_dPCR_paddW_8.6m2_spplot,
                                             rep2_dPCR_paddW_8.6m2_spplot,
                                             rep3_dPCR_paddW_8.6m2_spplot, nrow = 3)
################################################################################
# Saving plot of Combined 3 Reps for Paddock: W, NestedSubsample: 8.6 m2 (rake)#
################################################################################
ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "allreps_dPCR_paddW_8pt6m2_spplotV2.png"), plot=allreps_dPCR_paddW_8.6m2_spplot,
       width=10, height=6)
################################################################################
################################################################################
##################################################################################
#                                                                                #
# PADDOCK Y: SPATIAL DISTN Of dPCR (Ptt-FR freq %) for Spatial-Dependency Data   #
#                                                                                #
##################################################################################
dPCR_paddY_spplot <- spat_dep_sf_projected3 %>% 
  filter(PaddockFac == "Y") %>%
  ggplot() +
  # geom_sf(data = rep_poly_data, color = "red", fill = "white") +
  geom_sf(aes(color = Raw_dPCR), size=4) +
  geom_sf(shape = 1, size=4.1) +
  coord_sf() +
  theme_bw() +
  facet_grid(NestedSubsample ~ ., labeller = label_parsed) + # switch = "y"
  theme(legend.position = "top",
        strip.text.y = element_text(size=13, margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(face="bold", size=13),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(face="bold", size=12)) +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.05, 0.10, 0.15, 0.20),
                        limits = c(0, 0.20),
                        labels = c("0%", "5%", "10%", "15%", "20%"),
                        guide = guide_colorbar(barwidth = unit(6, "cm"))) +
  labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")))

################################################################################
# Saving dPCR_paddY_spplot                                                     #
################################################################################
ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "dPCR_PADDOCK_Y_FINAL_SPPLOT_0.png"), 
       plot=dPCR_paddY_spplot,
       width=5, height=12)
################################################################################
# Paddock-Y Reps 1.2 m2 REPLICATIONS                                           #
################################################################################
spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "Y", NestedSubsample == "b") %>%
  ggplot() +
  geom_sf(aes(color = ReplicationNo)) +
  coord_sf() +
  scale_x_continuous(breaks = seq(1900, 2400, by=20)) +
  scale_y_continuous(breaks = seq(3600, 3950, by=10))
################################################################################
# PADDOCK-Y -- 1.2m2 REPLICATIONS                                              #
################################################################################
# Paddock: Y, NestedSubsample: 1.2 m2 (rake); REP-1                            #
################################################################################
ax.txt.size <- 28
rep1_dPCR_paddY_1m_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "Y" & NestedSubsample == "b" & ReplicationNo == 1) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=9) +
  geom_sf(shape = 1, size=9.1, stroke=1.5) +
  coord_sf(xlim = c(2215, 2219), ylim = c(3640, 3648), 
           label_axes = list(bottom = "E", right = "N")) +
  scale_x_continuous(breaks = c(2215, 2219)) +
  theme_bw() +
  theme(axis.text.y = element_text(face="bold", size=ax.txt.size),
        axis.text.x = element_text(face="bold", size=ax.txt.size),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.05, 0.10, 0.15, 0.20),
                        limits = c(0, 0.20),
                        labels = c("0%", "5%", "10%", "15%", "20%")) 
#labs(color = expression(paste(italic(Ptt), "-FR frequency (%)")))
################################################################################
# Paddock: Y, NestedSubsample: 1.2 m2 (rake); REP-2                            #
################################################################################
rep2_dPCR_paddY_1m_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "Y" & NestedSubsample == "b" & ReplicationNo == 2) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=9) +
  geom_sf(shape = 1, size=9.1, stroke=1.5) +
  coord_sf(xlim = c(1947, 1951), ylim = c(3910, 3918),
           label_axes = list(bottom = "E", right = "N")) +
  scale_x_continuous(breaks = c(1947, 1951)) +
  #coord_sf(xlim = c(220, 250), ylim = c(335, 365)) +
  theme_bw() +
  theme(axis.text.y = element_text(face="bold", size=ax.txt.size),
        axis.text.x = element_text(face="bold", size=ax.txt.size),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.05, 0.10, 0.15, 0.20),
                        limits = c(0, 0.20),
                        labels = c("0%", "5%", "10%", "15%", "20%")) 
################################################################################
# Paddock: Y, NestedSubsample: 1.2 m2 (rake); REP-3                            #
################################################################################
rep3_dPCR_paddY_1m_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "Y" & NestedSubsample == "b" & ReplicationNo == 3) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=9) +
  geom_sf(shape = 1, size=9.1, stroke=1.5) +
  coord_sf(xlim = c(2330, 2334), ylim = c(3616, 3624), label_axes = list(bottom = "E", right = "N")) +
  scale_x_continuous(breaks = c(2330, 2334)) +
  theme_bw() +
  theme(axis.text.y = element_text(face="bold", size=ax.txt.size),
        axis.text.x = element_text(face="bold", size=ax.txt.size),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.05, 0.10, 0.15, 0.20),
                        limits = c(0, 0.20),
                        labels = c("0%", "5%", "10%", "15%", "20%"))
################################################################################
# Paddock: Y, NestedSubsample: 1.2 m2 (rake); Combined 3 Reps                  #
################################################################################
allreps_dPCR_paddY_1m_spplot <- ggarrange(rep2_dPCR_paddY_1m_spplot,
                                          rep1_dPCR_paddY_1m_spplot,
                                          rep3_dPCR_paddY_1m_spplot, nrow = 3)

allreps_dPCR_paddY_1m_horiz_spplot <- ggarrange(rep2_dPCR_paddY_1m_spplot,
                                                rep1_dPCR_paddY_1m_spplot,
                                                rep3_dPCR_paddY_1m_spplot, ncol = 3)
################################################################################
# Saving plot of Combined 3 Reps for Paddock: Y, NestedSubsample: 1.2 m2 (rake)#
################################################################################
ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "allreps_dPCR_paddY_1m_spplot.png"), plot=allreps_dPCR_paddY_1m_spplot,
       width=6, height=20)

ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "allreps_dPCR_paddY_1m_horiz_spplot.png"), plot=allreps_dPCR_paddY_1m_horiz_spplot,
       width=16, height=8)
################################################################################
################################################################################
# PADDOCK-Y -- 8.6m2 REPLICATIONS                                              #
################################################################################
# PAddock-Y Reps
################################################################################
spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "Y", NestedSubsample == "c") %>%
  ggplot() +
  geom_sf(aes(color = ReplicationNo)) +
  coord_sf() +
  scale_x_continuous(breaks = seq(1900, 2400, by=20)) +
  scale_y_continuous(breaks = seq(3600, 3950, by=10))
################################################################################
# Paddock: Y, NestedSubsample: 8.6 m2 (rake); REP-1                            #
################################################################################
pt_size <- 10.0
rep1_dPCR_paddY_8.6m2_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "Y" & NestedSubsample == "c" & ReplicationNo == 1) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=pt_size) +
  geom_sf(shape = 1, size=(pt_size+0.1), stroke=1.5) +
  coord_sf(xlim = c(2210, 2230), ylim = c(3620, 3670), 
           label_axes = list(bottom = "E", right = "N")) +
  scale_x_continuous(breaks = c(2210, 2230)) +
  theme_bw() +
  theme(axis.text = element_text(face="bold", size=ax.txt.size),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.05, 0.10, 0.15, 0.20),
                        limits = c(0, 0.20),
                        labels = c("0%", "5%", "10%", "15%", "20%")) 
################################################################################
# Paddock: Y, NestedSubsample: 8.6 m2 (rake); REP-2                            #
################################################################################
rep2_dPCR_paddY_8.6m2_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "Y" & NestedSubsample == "c" & ReplicationNo == 2) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=pt_size) +
  geom_sf(shape = 1, size=(pt_size+0.1), stroke=1.5) +
  coord_sf(xlim = c(1940, 1960), ylim = c(3895, 3945),
           label_axes = list(bottom = "E", right = "N")) +
  scale_x_continuous(breaks = c(1940, 1960)) +
  theme_bw() +
  theme(axis.text = element_text(face="bold", size = ax.txt.size),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.05, 0.10, 0.15, 0.20),
                        limits = c(0, 0.20),
                        labels = c("0%", "5%", "10%", "15%", "20%"))
################################################################################
# Paddock: Y, NestedSubsample: 8.6 m2 (rake); REP-3                            #
################################################################################
rep3_dPCR_paddY_8.6m2_spplot <- spat_dep_sf_projected2 %>% 
  filter(PaddockFac == "Y" & NestedSubsample == "c" & ReplicationNo == 3) %>%
  ggplot() +
  geom_sf(aes(color = Raw_dPCR), size=pt_size) +
  geom_sf(shape = 1, size=(pt_size+0.1), stroke=1.5) +
  coord_sf(xlim = c(2320, 2340), ylim = c(3600, 3650),
           label_axes = list(bottom = "E", right = "N")) +
  scale_x_continuous(breaks = c(2320, 2340)) +
  theme_bw() +
  theme(axis.text = element_text(face="bold", size=ax.txt.size),
        legend.position = "none") +
  scale_color_gradientn(colours = rev(terrain.colors(126)),
                        breaks = c(0.0, 0.05, 0.10, 0.15, 0.20),
                        limits = c(0, 0.20),
                        labels = c("0%", "5%", "10%", "15%", "20%"))
################################################################################
# Paddock: Y, NestedSubsample: 8.6 m2 (rake); Combined 3 Reps                  #
################################################################################
allreps_dPCR_paddY_8.6m2_spplot <- ggarrange(rep2_dPCR_paddY_8.6m2_spplot,
                                             rep1_dPCR_paddY_8.6m2_spplot,
                                             rep3_dPCR_paddY_8.6m2_spplot, ncol = 3)
################################################################################
# Saving plot of Combined 3 Reps for Paddock: Y, NestedSubsample: 8.6 m2 (rake)#
################################################################################
ggsave(here::here("FilesForLeon_9June2022",
                  "SpatDep", "Figures",
                  "allreps_dPCR_paddY_8.6m2_horiz_spplot.png"), plot=allreps_dPCR_paddY_8.6m2_spplot,
       width=15, height=8)
################################################################################
################################################################################





























