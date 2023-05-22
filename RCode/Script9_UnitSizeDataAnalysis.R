##############################################################
# Load Required Packages                                     #
##############################################################
library(here)
library(readxl)
library(tidyverse)
library(asreml)
library(sp)
library(sf)
library(spatstat)
library(rgdal)
library(maptools)
library(ggpubr)
library(flextable)
library(viridis)
library(colorspace)
library(ggmap)
library(ggthemes)
library(fitdistrplus)
library(xtable)
library(car)
library(xtable)
library(onewaytests)
##############################################################
# Source utility functions                                   #
##############################################################
source(here::here("RCodeForPaper", "Utils_ForPaper.R"))
##############################################################
# Load the complete Unit Size data                           #
##############################################################
unit_size_data <- readRDS(here::here("DerivedData", "unit_data.rds"))
##############################################################
# To analyse the response variable Ptt fungicide resistance  #
# frequency (measured as percentage of Cyp51A F489L mutation)#
# we exclude observations with the Ptt DNA density <= 29.0   #
#                                                            #
# Note that the response variable resistance frequency is    #
# saved under the variable name *Raw_dPCR* and the DNA       #   
# density is under the name *Cyp51A_copies*. We use the      #
# square-root and logarithm transformations for the response #
# variables Raw_dPCR and Cyp51A_copies, respectively, to meet#
# the normality assumption of the linear model.              #
##############################################################

##############################################################
# Prepare data for the analysis of the response variable     #
# sqrt(dPCR) by removing observations with                   #
#                Cyp51A_copies <= 29.0                       #
##############################################################
frequency_data <- unit_size_data %>%
                  dplyr::filter(Cyp51A_copies > 29.0)
##############################################################
# Check that Cyp51A copies <= 29.0 are removed               #
##############################################################
summary(frequency_data$Cyp51A_copies)
##############################################################
# Prepare data for the analysis of the response variable     #
#              ln(Cyp51A copies).                            #   
# We are going to use all the data without any filtering     #
# for this analysis                                          #
##############################################################
density_data <- unit_size_data
summary(density_data$Cyp51A_copies)
################################################################################
# Variance analysis: the hypothesis from the domain experts is that the        #
#                    variance of the responses should decrease as we increase  #
#                    the size of the sampling unit area size.                  #
#                                                                              #
# We have done some initial exploratory analysis and it has been revealed that #
# there may not be any difference in variances for paddocks M and P, i.e., the #
# homoscedasticity assumption is valid for these two paddocks. In contrast,    #
# for the field V, there may be significant difference between variances of    #
# the response variable sqrt_dPCR for different unit area sizes.               #
#                                                                              #  
# We use the Levene's test to investigate possible variance heteroscedasticity #
# due to the unit area size.                                                   #
################################################################################
# Response: sqrt_dPCR                                                          #
# Perform Levene's test to test for heteroscedasticity                         #
################################################################################
# Comparison of sampling strategies: raking vs hand-picking in 1.2 m2          #
################################################################################
frequency_data_for_strategy_comparison <- frequency_data %>%
                            dplyr::filter(NestedSubsample %in% c("a", "b")) %>%
                        dplyr::mutate(NestedSubsample = factor(NestedSubsample))
frequency_data_for_strategy_comparison$NestedSubsample
# Group the data within each paddock
strategy_data_nested <- frequency_data_for_strategy_comparison %>%
                        dplyr::group_by(PaddockFac) %>%
                        tidyr::nest()
# Perform Levene's test and Brown-Forsythe test within each paddock
strategy_data_levene_test <- strategy_data_nested %>%
                  mutate(LeveneTest = map(data, 
                                          ~ car::leveneTest(sqrt_dPCR ~ NestedSubsample,
                                                            data = .x, 
                                                            center=mean)),
                         BrownForsytheTest = map(data, 
                                                ~ car::leveneTest(sqrt_dPCR ~ NestedSubsample,
                                                                  data = .x, 
                                                                  center=median)))
# Extract F-values and p-values
strategy_data_pvals <- strategy_data_levene_test %>%
  dplyr::mutate(FDist = map_chr(LeveneTest, get_Fdf),
                LeveneFstat = map_dbl(LeveneTest, get_Fvalue),
                LevenePvalue = map_dbl(LeveneTest, get_pvalue),
                BF.Fstat = map_dbl(BrownForsytheTest, get_Fvalue),
                BF.Pvalue = map_dbl(BrownForsytheTest, get_pvalue))
strategy_data_pvals # Levene test indicates there may be heteroscedasticity in Paddock-M
# Let's investigate variances for Paddock-Strategy combinations
frequency_data_for_strategy_comparison %>%
  dplyr::group_by(PaddockFac, NestedSubsample) %>%
  dplyr::summarise(Variance = var(sqrt_dPCR))
########################################################################################
# sqrt(dPCR) ~ sampling strategy                                                       #
########################################################################################
# Separate variances for: Paddock-P and for the two strategies
# variance-1 Paddock: P and strategy: a 
# variance-2 Paddock: P and strategy: b
# variance-3 Paddock: M and V and strategy: a
# variance-4 Paddock: M and V and strategy: b
########################################################################################
frequency_data_for_strategy_comparison_with_section <- frequency_data_for_strategy_comparison %>%
  mutate(Section = case_when(PaddockFac=="P" & NestedSubsample=="a" ~ "Sec-1_P_a",
                             PaddockFac=="P" & NestedSubsample=="b" ~ "Sec-1_P_b",
                             PaddockFac!="P" & NestedSubsample=="a" ~ "Sec-3_MV_a",
                             PaddockFac!="P" & NestedSubsample=="b" ~ "Sec-4_MV_b"),
         Section = factor(Section)) %>%
  arrange(Section)

# sqrt_dPCR ~ NestedSubsample (LMM)
frequency_lmm_for_strategy_comparison <- asreml(data = frequency_data_for_strategy_comparison_with_section,
                                                fixed = sqrt_dPCR ~ NestedSubsample + BarleyCultivar,
                                                random = ~ PaddockFac,
                                                residual = ~ dsum(~id(units) | Section),
                                                na.action = na.method(y = "include"))
# Wald test
wald.asreml(frequency_lmm_for_strategy_comparison)
# Variance components
options(scipen = 999)
summary(frequency_lmm_for_strategy_comparison)$varcomp
# Effects
summary(frequency_lmm_for_strategy_comparison, coef=TRUE)$coef.fixed
########################################################################################
# The sample variance estimates show that the sample strategy-a variance is
# consistently larger than the strategy-b variance for all three paddocks.
# Also, the variance for Paddock-P is quite different than the variances in the other
# two paddocks
########################################################################################
# Levene test to compare variances between two strategies for all paddocks combined    #
########################################################################################
car::leveneTest(sqrt_dPCR ~ NestedSubsample,
                data = frequency_data_for_strategy_comparison,
                center = mean) # p-value is 0.004467 (**); F-value 8.4552 df(1, 102)
# Comparing field-P versus the other two fields
FieldP_Vs_OtherFields <- frequency_data_for_strategy_comparison %>%
                         mutate(PaddockFac2 = case_when(PaddockFac=="P" ~ "Field-P",
                                                        PaddockFac!= "P" ~ "Other"),
                                PaddockFac2 = factor(PaddockFac2)) %>%
                         group_by(NestedSubsample) %>% nest()
FieldP_Vs_OtherFields_LeveneTest <- FieldP_Vs_OtherFields %>%
  mutate(LeveneTest = map(data, 
                          ~ car::leveneTest(sqrt_dPCR ~ PaddockFac,
                                            data = .x, 
                                            center=mean)),
         FDist = map_chr(LeveneTest, get_Fdf),
         LeveneFstat = map_dbl(LeveneTest, get_Fvalue),
         LevenePvalue = map_dbl(LeveneTest, get_pvalue))
FieldP_Vs_OtherFields_LeveneTest
################################################################################
# log(Cyp51A_copy) ~ NestedSubsample                                           #
################################################################################
# Levene's test for comparing the variances of the two sampling strategies
# Perform Levene's test and Brown-Forsythe test within each paddock
density_data_for_strategy_comparison <- density_data %>%
                                dplyr::filter(NestedSubsample%in% c("a","b")) %>%
  dplyr::mutate(NestedSubsample = factor(NestedSubsample))
density_data_for_strategy_comparison$NestedSubsample

# Compute summary for Paddock-specific variances for the two strategies
density_data_for_strategy_comparison %>%
  dplyr::group_by(PaddockFac, NestedSubsample) %>%
  dplyr::summarise(Variance = var(log_Cyp51A_copies))
#######################################################################################
# It seems that each paddock has different variances and there is difference between  #
# variances of the two strategies for some paddocks                                   #
#######################################################################################
# Three paddock homoscedasticity test
density_data_levenetest_stategy <- density_data_for_strategy_comparison %>%
                                  dplyr::group_by(PaddockFac) %>%
                                  tidyr::nest() %>%
                    dplyr::mutate(LeveneTest = map(data, 
                          ~ car::leveneTest(log_Cyp51A_copies ~ NestedSubsample,
                                            data = .x, 
                                            center=mean)),
         FDist = map_chr(LeveneTest, get_Fdf),
         LeveneFstat = map_dbl(LeveneTest, get_Fvalue),
         LevenePvalue = map_dbl(LeveneTest, get_pvalue))
density_data_levenetest_stategy
################################################################################
# Fit the model for the response variable log(Density)                         #
################################################################################
summary(density_data_for_strategy_comparison$Cyp51A_copies)

density_data_for_strategy_comparison2 <- density_data_for_strategy_comparison %>%
   mutate(Section = case_when(PaddockFac=="M" ~ "Sec-1_M",
                              PaddockFac=="V" ~ "Sec-2_V",
                              PaddockFac=="P" & NestedSubsample=="a" ~ "Sec-3_P_a",
                              PaddockFac=="P" & NestedSubsample=="b" ~ "Sec-4_P_b"),
          Section = factor(Section)) %>%
   arrange(Section)

# log_Cyp51A_copies ~ NestedSubsample (LMM)
density_lmm_for_strategy_comparison0 <- asreml(data = density_data_for_strategy_comparison2,
                                              fixed = log_Cyp51A_copies ~ NestedSubsample + BarleyCultivar,
                                              random = ~ PaddockFac,
                                              residual = ~ dsum(~id(units) | Section),
                                              na.action = na.method(y = "include"))
# Should we drop Cultivar term?
wald.asreml(density_lmm_for_strategy_comparison0)
# Final model after dropping the Cultivar
density_lmm_for_strategy_comparison <- asreml(data = density_data_for_strategy_comparison2,
                                               fixed = log_Cyp51A_copies ~ NestedSubsample,
                                               random = ~ PaddockFac,
                                               residual = ~ dsum(~id(units) | Section),
                                               na.action = na.method(y = "include"))
# Wald test
wald.asreml(density_lmm_for_strategy_comparison)
# Variance components
summary(density_lmm_for_strategy_comparison)$varcomp
# Effects
summary(density_lmm_for_strategy_comparison, coef=TRUE)$coef.fixed

################################################################################
# Let's compare the variances between different paddocks                       #
################################################################################

################################################################################
#                                                                              #  
# We interested in comparing unit area sizes 1.2m2 (b), 8.6m2 (c), and 60m2 (d)#
#                                                                              #
################################################################################
frequency_data_for_unit_area_size_comparison <- frequency_data %>%
                                          dplyr::filter(NestedSubsample != "a") %>%
                                          mutate(NestedSubsample = factor(NestedSubsample))

summary(frequency_data_for_unit_area_size_comparison$NestedSubsample)
# Paddock-specific variances for unit area sizes
frequency_data_for_unit_area_size_comparison %>%
  dplyr::group_by(PaddockFac, NestedSubsample) %>%
  dplyr::summarise(Variance = var(sqrt_dPCR))

# Let's look at the Levene's test results for each paddock
frequency_data_for_unit_area_size_comparison %>%
  dplyr::group_by(PaddockFac) %>%
  tidyr::nest() %>%
  mutate(LeveneTest = map(data, ~ car::leveneTest(sqrt_dPCR ~ NestedSubsample,
                                                                 data = .x, 
                                                                 center=mean)),
         FDist = map_chr(LeveneTest, get_Fdf),
         LeveneFstat = map_dbl(LeveneTest, get_Fvalue),
         LevenePvalue = map_dbl(LeveneTest, get_pvalue))
#################################################################################
# Create sections for incorporating difference residual variances
################################################################################
frequency_data_for_unit_area_size_comparison2 <- frequency_data_for_unit_area_size_comparison %>%
   mutate(Section = case_when(PaddockFac=="M" ~ "Sec-1_M",
          PaddockFac=="P" ~ "Sec-2_P",
          PaddockFac=="V" & NestedSubsample!="c" ~ "Sec-3_V_b_and_d",
          PaddockFac=="V" & NestedSubsample=="c" ~ "Sec-4_V_c"),
          Section = factor(Section)) %>%
  arrange(Section)

# sqrt_dPCR ~ NestedSubsample (LMM)
frequency_lmm_for_unitsize_comparison <- asreml(data = frequency_data_for_unit_area_size_comparison2,
                                                fixed = sqrt_dPCR ~ NestedSubsample + BarleyCultivar,
                                                random = ~ PaddockFac,
                                                residual = ~ dsum(~id(units) | Section),
                                                na.action = na.method(y = "include"))
# Wald test
wald.asreml(frequency_lmm_for_unitsize_comparison)
# Variance components
options(scipen = 999)
summary(frequency_lmm_for_unitsize_comparison)$varcomp
# Effects
summary(frequency_lmm_for_unitsize_comparison, coef=TRUE)$coef.fixed
################################################################################
# log(Cyp51A_copy) ~ NestedSubsample  (for comparing unit area sizes)          #
################################################################################
density_data_for_unit_area_comparison <- density_data %>%
  dplyr::filter(NestedSubsample!="a") %>%
  dplyr::mutate(NestedSubsample = factor(NestedSubsample))

density_data_for_unit_area_comparison$NestedSubsample

# Compute summary for Paddock-specific variances for the three unit area sizes
density_data_for_unit_area_comparison %>%
  dplyr::group_by(PaddockFac, NestedSubsample) %>%
  dplyr::summarise(Mean = mean(log_Cyp51A_copies), Variance = var(log_Cyp51A_copies))
################################################################################
# Levene's test for comparing the variances of the three unit area sizes       #  
# Perform Levene's test for each paddock                                       #
################################################################################
# Three paddock homoscedasticity test for unit area size comparisons
density_data_for_unit_area_comparison %>%
  dplyr::group_by(PaddockFac) %>%
  tidyr::nest() %>%
  dplyr::mutate(LeveneTest = map(data, 
                                 ~ car::leveneTest(log_Cyp51A_copies ~ NestedSubsample,
                                                   data = .x, 
                                                   center=mean)),
                FDist = map_chr(LeveneTest, get_Fdf),
                LeveneFstat = map_dbl(LeveneTest, get_Fvalue),
                LevenePvalue = map_dbl(LeveneTest, get_pvalue))

# Because homoscedasticity is rejected for Field-V, we
# perform pairwise comparisons for the three unit area sizes
# b vs c
density_data_for_unit_area_comparison %>%
  filter(PaddockFac == "V", NestedSubsample != "d") %>%
  mutate(NestedSubsample=factor(NestedSubsample)) %>%
  dplyr::group_by(PaddockFac) %>%
  tidyr::nest() %>%
  dplyr::mutate(LeveneTest = map(data, 
                                 ~ car::leveneTest(log_Cyp51A_copies ~ NestedSubsample,
                                                   data = .x, 
                                                   center=mean)),
                FDist = map_chr(LeveneTest, get_Fdf),
                LeveneFstat = map_dbl(LeveneTest, get_Fvalue),
                LevenePvalue = map_dbl(LeveneTest, get_pvalue))

# b vs d
density_data_for_unit_area_comparison %>%
  filter(PaddockFac == "V", NestedSubsample != "c") %>%
  mutate(NestedSubsample=factor(NestedSubsample)) %>%
  dplyr::group_by(PaddockFac) %>%
  tidyr::nest() %>%
  dplyr::mutate(LeveneTest = map(data, 
                                 ~ car::leveneTest(log_Cyp51A_copies ~ NestedSubsample,
                                                   data = .x, 
                                                   center=mean)),
                FDist = map_chr(LeveneTest, get_Fdf),
                LeveneFstat = map_dbl(LeveneTest, get_Fvalue),
                LevenePvalue = map_dbl(LeveneTest, get_pvalue))
# c vs d
density_data_for_unit_area_comparison %>%
  filter(PaddockFac == "V", NestedSubsample != "b") %>%
  mutate(NestedSubsample=factor(NestedSubsample)) %>%
  dplyr::group_by(PaddockFac) %>%
  tidyr::nest() %>%
  dplyr::mutate(LeveneTest = map(data, 
                                 ~ car::leveneTest(log_Cyp51A_copies ~ NestedSubsample,
                                                   data = .x, 
                                                   center=mean)),
                FDist = map_chr(LeveneTest, get_Fdf),
                LeveneFstat = map_dbl(LeveneTest, get_Fvalue),
                LevenePvalue = map_dbl(LeveneTest, get_pvalue))
#########################################################################################
density_data_for_unit_area_comparison2 <- density_data_for_unit_area_comparison %>%
  mutate(Section = case_when(PaddockFac=="M" ~ "Sec-1_M",
                             PaddockFac=="P" ~ "Sec-2_P",
                             PaddockFac=="V" & NestedSubsample=="b" ~ "Sec-3_V_b",
                             PaddockFac=="V" & NestedSubsample!="b" ~ "Sec-4_V_c_and_d"),
         Section = factor(Section)) %>%
  arrange(Section)

# log(Cyp51A copies) ~ NestedSubsample (LMM)
density_lmm_for_unitsize_comparison0 <- asreml(data = density_data_for_unit_area_comparison2,
                                                fixed = log_Cyp51A_copies ~ NestedSubsample + BarleyCultivar,
                                                random = ~ PaddockFac,
                                                residual = ~ dsum(~id(units) | Section),
                                                na.action = na.method(y = "include"))
# Wald test
wald.asreml(density_lmm_for_unitsize_comparison0)
# Drop Cultivar
density_lmm_for_unitsize_comparison <- asreml(data = density_data_for_unit_area_comparison2,
                                              fixed = log_Cyp51A_copies ~ NestedSubsample,
                                              random = ~ PaddockFac,
                                              residual = ~ dsum(~id(units) | Section),
                                              na.action = na.method(y = "include"))
plot(density_lmm_for_unitsize_comparison)
# Wald test
wald.asreml(density_lmm_for_unitsize_comparison)

density_lmm_for_unitsize_comparison2 <- asreml(data = density_data_for_unit_area_comparison2,
                                              fixed = log_Cyp51A_copies ~ NestedSubsample,
                                              random = ~ PaddockFac,
                                              # residual = ~ dsum(~id(units) | PaddockFac),
                                              na.action = na.method(y = "include"))
plot(density_lmm_for_unitsize_comparison2)
# Wald test
wald.asreml(density_lmm_for_unitsize_comparison2)
# Wald test
wald.asreml(density_lmm_for_unitsize_comparison)
####################
summary(density_lmm_for_unitsize_comparison)$aic
summary(density_lmm_for_unitsize_comparison)$bic
# The Area effect is significant
summary(lm(log_Cyp51A_copies ~ NestedSubsample, data = density_data_for_unit_area_comparison2))
plot(lm(log_Cyp51A_copies ~ NestedSubsample, data = density_data_for_unit_area_comparison2))
summary(aov(log_Cyp51A_copies ~ NestedSubsample, data = density_data_for_unit_area_comparison2))
# Variance components
options(scipen = 999)
summary(density_lmm_for_unitsize_comparison)$varcomp
# Effects
summary(density_lmm_for_unitsize_comparison, coef=TRUE)$coef.fixed

density_data_for_unit_area_comparison %>% group_by(NestedSubsample) %>%
  summarise(Average = mean(log_Cyp51A_copies), sd=sd(log_Cyp51A_copies))

density_data_for_unit_area_comparison %>% group_by(PaddockFac, NestedSubsample) %>%
  summarise(Average = mean(log_Cyp51A_copies), sd = sd(log_Cyp51A_copies))
##############################################################################################










