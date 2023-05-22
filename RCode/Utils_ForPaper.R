# This function computes the semi-variogram alongside the bootstrap confidence
# bands for chosen variables 
compute_semivariogram_envelope <- function(data, var, paddock, samp_design,
                                           dist_brks, boot_sim, boot_seed){
  # Subset data for the specified paddock and sampling design 
  data <- data %>%
          dplyr::filter(PaddockFac == paddock, 
                        NestedSubsample == samp_design)
  # 0.5*(zi - zj)^2 for distances ||xi - xj||
  squared_var_diff_and_dist_data <-  compute_svgm(data, {{var}})
  squared_var_diff_and_dist_data <- squared_var_diff_and_dist_data %>%
                        dplyr::mutate(Dist_fac = cut(Dist, breaks = dist_brks))
  # semi-variance and distances for the data
  svar_and_distances <- squared_var_diff_and_dist_data %>%
                        dplyr::group_by(Dist_fac) %>%
                        dplyr::summarise(Dist = round(mean(Dist)),
                                         SemiVariance = mean(Svar))
  # Simulation Envelope Calculation
  squared_var_diff_from_bootstrap <- compute_bootstrap_svgm(data, {{var}}, 
                                                nsim=boot_sim, 
                                                seedval=boot_seed)
  # Create distance bins
  squared_var_diff_from_bootstrap <- squared_var_diff_from_bootstrap %>%
                            dplyr::mutate(Dist_fac = cut(Dist, breaks=dist_brks))
  
  # Semi-variances for all bootstrap samples
  svar_from_bootstrap <- squared_var_diff_from_bootstrap %>%
                         group_by(Dist_fac) %>%
                         summarise_all(.funs = mean) %>%
                         dplyr::select(-Dist, -Svar)
  
  # Get the mean distance in each bin
  avg_dist <- svar_and_distances %>% dplyr::select(Dist_fac, Dist)
  
  # Include mean distances in the semi-variance datasets 
  svar_from_bootstrap <- svar_from_bootstrap %>%
                         left_join(avg_dist,  by = "Dist_fac")
  
  # Compute pointwise envelopes (95%, 90%, and 80%)
    
    svar_from_bootstrap <- svar_from_bootstrap %>% 
    dplyr::select(Dist, contains("Boot")) %>%
    pivot_longer(cols = -Dist, names_to = "BootId", values_to = rlang::as_name(enexpr(var))) %>%
    group_by(Dist) %>%
    summarise(LowLim95  = quantile({{var}}, probs = c(0.025)),
              UppLim95 = quantile({{var}}, probs = c(0.975)),
              LowLim90 = quantile({{var}}, probs = c(0.05)),
              UppLim90 = quantile({{var}}, probs = c(0.95)),
              LowLim80 = quantile({{var}}, probs = c(0.10)),
              UppLim80 = quantile({{var}}, probs = c(0.90))) %>%
    left_join(svar_and_distances, by="Dist")
    
    return(svar_from_bootstrap)
}

#############################################################################
# This function computes for each pair (zi, zj) at a location (si, sj),     #
# 0.5*(zi - zj)^2; the function returns these values alongside the          # 
# Euclidean distance between the points si and sj.                          #
#############################################################################
compute_svgm <- function(data, var){
  val <-  data %>% 
          dplyr::select({{var}}) %>% pull()
  
  xcoord <- expand.grid(x1 = data$X, x2 = data$X)
  ycoord <- expand.grid(y1 = data$Y, y2 = data$Y)
  
  vals <- expand.grid(val1 = val, val2 = val)
  
  dist <- sqrt((xcoord$x1 - xcoord$x2)^2 + (ycoord$y1 - ycoord$y2)^2)
  svar <- 0.5 * (vals$val1 - vals$val2)^2
  
  ret_data <- data.frame(Dist = dist, Svar = svar) %>%
              bind_cols(xcoord, ycoord) %>%
              dplyr::select(x1, y1, x2, y2, Dist, Svar) %>%
              dplyr::filter(Dist > 0)
    
  return(ret_data)
}

###########################################################################
# This function computes bootstrapped versions of the 0.5*(zi - zj)^2     #
# based on the bootstrap samples under the null distribution of no        #
# auto-correlation.                                                        #
###########################################################################
compute_bootstrap_svgm <- function(data, var, 
                                   nsim=1000, 
                                   seedval=100){
  var_vals <- data %>% dplyr::select({{var}}) %>% pull()
  
  xcoord <- expand.grid(x1 = data$X, x2 = data$X)
  ycoord <- expand.grid(y1 = data$Y, y2 = data$Y)
  vals <- expand.grid(val1 = var_vals, val2 = var_vals)
  
  dist <- sqrt((xcoord$x1 - xcoord$x2)^2 + (ycoord$y1 - ycoord$y2)^2)
  svar <- 0.5 * (vals$val1 - vals$val2)^2
  
  nRow <- length(dist)
  
  set.seed(seedval)
  seed_vals <- sample(1:1e7, size=nsim, replace=FALSE)
  
  ret_mat <- matrix(NA, nrow = nRow, ncol = nsim)
  
  for(j in 1:nsim){
    set.seed(seed_vals[j])
    var_vals_boot <- sample(var_vals, replace=TRUE)
    var_expanded_grid <- expand.grid(val1=var_vals_boot, val2=var_vals_boot)
    svar_boot <- 0.5 * (var_expanded_grid$val1 - var_expanded_grid$val2)^2
    ret_mat[,j] <- svar_boot
  }
  
  ret_df0 <- as.data.frame(ret_mat)
  names(ret_df0) <- paste0("BootSvar", 1:nsim)
  ret_df <- data.frame(Dist = dist,
                       Svar = svar) %>%
            bind_cols(ret_df0) %>%
            dplyr::filter(Dist > 0)
  
  return(ret_df)
}



plot_dPcr_density <- function(data, varname, label, adj=2){
    data %>%
    ggplot(aes(x={{varname}})) +
    geom_density(aes(fill = NestedSubsample),
                 adjust = adj) +
    facet_grid(PaddockFac~NestedSubsample) +
    labs(x = label) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill = "white"))
}

plot_dPcr_boxplt <- function(data, 
                             varname, 
                             label="Raw dPCR",
                             strip.text.size = 12,
                             axis.text.size = 12,
                             title.text.size = 15,
                             pt.size = 1.5,
                             out.size=1.5){
     unit_sub_data_means <- data %>%
     group_by(PaddockFac, NestedSubsample) %>%
     summarise({{varname}} := mean({{varname}}, na.rm = TRUE),
               .groups = "drop")
  
    data %>%
    ggplot(aes(x=NestedSubsample, y={{varname}})) +
    geom_boxplot(aes(color = NestedSubsample), fill = "white",
                 outlier.size = out.size) +
    geom_point(aes(color = NestedSubsample), size=pt.size) +
    geom_point(data = unit_sub_data_means, shape = 4, size=round(pt.size*1.5), stroke=1.2) +
    facet_grid(~ PaddockFac) +
    labs(x = "Sampling strategy",
         y = label) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold", size = strip.text.size),
          axis.text = element_text(size = axis.text.size),
          axis.title = element_text(size = title.text.size)) 
}

compute_summary <- function(data, grp_by_var, var, label){
    data %>% group_by({{grp_by_var}}) %>%
    summarise(MIN = min({{var}}, na.rm = TRUE),
              MAX = max({{var}}, na.rm = TRUE),
              MEDIAN = median({{var}}, na.rm = TRUE),
              MEAN = mean({{var}}, na.rm = TRUE),
              SD = sd({{var}}, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Var = label) %>%
    relocate(Var, {{grp_by_var}})
}

get_wald_test_result <- function(obj){
  stopifnot(inherits(obj, "wald"))
  
  ret_df <- as.data.frame(obj)
  ret_df$Component <- row.names(ret_df)
  names(ret_df) <- c("df", "Sum of square", "Wald statistic",
                     "p-value", "Component")
  ret_df <- ret_df[c("Component", "df", "Sum of square", "Wald statistic",
                     "p-value")]
  return(ret_df)
}
################################################################################
################################################################################
get_predictions <- function(Obj, 
                            classify = "NestedSubsample",
                            transformation = "identity"){
  stopifnot(inherits(Obj, "asreml"))
  invisible(capture.output(pred <- predict(Obj, classify = classify, sed = TRUE)$pvals))
  
  if(transformation == "identity"){
    pred$Response.BLUP <- pred$predicted.value
    pred$Approx.StdError <- pred$std.error
  }
  if(transformation == "sqrt"){
    pred$Response.BLUP <- (pred$predicted.value)^2
    pred$Approx.StdError <-  2*(pred$predicted.value)*(pred$std.error)
  }
  if(transformation == "log"){
    pred$Response.BLUP <- exp(pred$predicted.value)
    pred$Approx.StdError <- with(pred, exp(predicted.value)*std.error)
  }
  pred$predicted.value <- NULL
  pred$std.error <- NULL
  pred$status <- NULL
  return(pred)
}
################################################################################
# Extract Legend from the plot                                                 #
################################################################################
extract_legend <- function(plt){
  gtab.obj <- ggplot_gtable(ggplot_build(plt))
  legend.id <- which(sapply(gtab.obj$grob, function(z){z$name}) == "guide-box")
  my.legend <- gtab.obj$grobs[[legend.id]]
  
  return(my.legend)
}