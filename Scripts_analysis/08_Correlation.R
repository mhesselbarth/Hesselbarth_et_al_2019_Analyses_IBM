###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Import libraries and data ####

# load packages #
source("Helper_functions/helper_functions_setup.R")

plot_area <- readr::read_rds("Data/Raw/plot_area_sp.rds")

# density used for model fitting
abiotic_cond_real_fit <- readr::read_rds("Data/Input/abiotic_cond_real_fit.rds") %>% 
  magrittr::extract2("scaled")

# density used for model run
abiotic_cond_real_model <- readr::read_rds("Data/Input/abiotic_cond_real_model.rds") %>% 
  magrittr::extract2("scaled")

# get environmental data
environmetal_data_names <- list.files("Data/Raw/environmental", full.names = TRUE) %>% 
  stringr::str_subset(pattern = "DEM.rds", negate = TRUE)

# import environmental data as raster
environmental_data <- purrr::map(environmetal_data_names, function(x) {
  readr::read_rds(x) %>% raster::rasterFromXYZ()})

# rename list elements
names(environmental_data) <- stringr::str_sub(environmetal_data_names, start = 24, 
                                              end = -5)

# read DEM
dem <- readr::read_rds("Data/Raw/environmental/DEM.rds") 

#### Preprocess data ###

# mask environmental data to plot area
environmental_data <- raster::stack(environmental_data) %>% 
  raster::crop(y = plot_area) %>% 
  raster::mask(mask = plot_area)

# mask DEM to plot area
dem <- raster::crop(dem, y = plot_area) %>% 
  raster::mask(mask = plot_area)

# resample density values to environmental values
abiotic_cond_real_fit_rsmenv <- raster::resample(x = abiotic_cond_real_fit, 
                                                 y = environmental_data) %>% 
  raster::as.data.frame(xy = TRUE)

abiotic_cond_real_model_rsmenv <- raster::resample(x = abiotic_cond_real_model, 
                                                   y = environmental_data) %>% 
  raster::as.data.frame(xy = TRUE)

# resample density values to DEM
abiotic_cond_real_fit_rsmdem <- raster::resample(x = abiotic_cond_real_fit, 
                                                 y = dem) %>% 
  raster::as.data.frame(xy = TRUE)

abiotic_cond_real_model_rsmdem <- raster::resample(x = abiotic_cond_real_model, 
                                                   y = dem) %>% 
  raster::as.data.frame(xy = TRUE)

method <- "spearman"

#### Calculate correlations using coords ####

# get coordinates of cell centers for density data
coords_fit <- raster::coordinates(abiotic_cond_real_fit)

coords_model <- raster::coordinates(abiotic_cond_real_model)

# should be identical
all(coords_fit == coords_fit)

# correlations between density and environmental values
correlations_environmental_coords <- purrr::map_dfr(1:raster::nlayers(environmental_data), function(i) {
  
  # get name of raster layer
  data_name <- names(environmental_data[[i]]) %>% 
    stringr::str_to_lower() %>% 
    stringr::str_replace_all(pattern = "_", replacement = " ")
  
  # get cell id of locations at density cells
  cells_fit <- raster::cellFromXY(object = environmental_data[[i]], 
                                  xy = coords_fit)
  
  cells_model <- raster::cellFromXY(object = environmental_data[[i]],
                                    xy = coords_model)

  # get environmental values at cells 
  values_fit <- raster::values(environmental_data[[i]])[cells_fit]
  
  values_model <- raster::values(environmental_data[[i]])[cells_model]
  
  # calculate correlation
  corr_fit <- cor(raster::values(abiotic_cond_real_fit), values_fit, 
                  method = method, use = "complete.obs")
  
  corr_model <- cor(raster::values(abiotic_cond_real_model), values_model, 
                    method = method, use = "complete.obs")
  
  tibble::tibble(data = data_name, model_fitting = corr_fit, model_run = corr_model, 
                 method = "Cell coords")
}) 

# get correlations between density and DEM
correlations_dem_coords <- purrr::map_dfr(1:raster::nlayers(dem), function(i) {
  
  # get name of raster layers
  data_name <- names(dem[[i]]) %>% 
    stringr::str_to_lower() %>% 
    stringr::str_replace_all(pattern = "_", replacement = " ")
  
  # get cell id of locations at density cells
  cells_fit <- raster::cellFromXY(object = dem[[i]], 
                                  xy = coords_fit)
  
  cells_model <- raster::cellFromXY(object = dem[[i]],
                                    xy = coords_model)
  
  # get DEM values at cell coords
  values_fit <- raster::values(dem[[i]])[cells_fit]
  
  values_model <- raster::values(dem[[i]])[cells_model]
  
  # calculate correlation
  corr_fit <- cor(raster::values(abiotic_cond_real_fit), values_fit, 
                  method = "pearson", use = "complete.obs")
  
  corr_model <- cor(raster::values(abiotic_cond_real_model), values_model, 
                    method = "pearson", use = "complete.obs")
  
  tibble::tibble(data = data_name, model_fitting = corr_fit, model_run = corr_model, 
                 method = "Cell coords")
}) 

#### Calculate correlations using resample ####

# correlations between density and environmental values
correlations_environmental_resample <- purrr::map_dfr(1:raster::nlayers(environmental_data), function(i) {
  
  # get name of raster layer
  data_name <- names(environmental_data[[i]]) %>% 
    stringr::str_to_lower() %>% 
    stringr::str_replace_all(pattern = "_", replacement = " ")
  
  data_temp <- raster::as.data.frame(environmental_data[[i]], xy = TRUE)
  
  values_fit <- dplyr::full_join(x = data_temp, y = abiotic_cond_real_fit_rsmenv, 
                                 by = c("x", "y")) %>% 
    dplyr::filter(complete.cases(.))
  
  values_model <- dplyr::full_join(x = data_temp, y = abiotic_cond_real_model_rsmenv, 
                                   by = c("x", "y")) %>% 
    dplyr::filter(complete.cases(.))
  
  # calculate correlation
  corr_fit <- cor(values_fit[, 3], values_fit[, 4], method = method)
  
  corr_model <- cor(values_model[, 3], values_model[, 4], method = method)
  
  tibble::tibble(data = data_name, model_fitting = corr_fit, model_run = corr_model, 
                 method = "Resample values")
}) 

# correlations between density and environmental values
correlations_dem_resample <- purrr::map_dfr(1:raster::nlayers(dem), function(i) {
  
  # get name of raster layer
  data_name <- names(dem[[i]]) %>% 
    stringr::str_to_lower() %>% 
    stringr::str_replace_all(pattern = "_", replacement = " ")
  
  data_temp <- raster::as.data.frame(dem[[i]], xy = TRUE)
  
  values_fit <- dplyr::full_join(x = data_temp, y = abiotic_cond_real_fit_rsmdem, 
                                 by = c("x", "y")) %>% 
    dplyr::filter(complete.cases(.))
  
  values_model <- dplyr::full_join(x = data_temp, y = abiotic_cond_real_model_rsmdem, 
                                   by = c("x", "y")) %>% 
    dplyr::filter(complete.cases(.))
  
  # calculate correlation
  corr_fit <- cor(values_fit[, 3], values_fit[, 4], method = method)
  
  corr_model <- cor(values_model[, 3], values_model[, 4], method = method)
  
  tibble::tibble(data = data_name, model_fitting = corr_fit, model_run = corr_model, 
                 method = "Resample values")
}) 

#### Results ###

# combine to one dataframe
correlation_full <- dplyr::bind_rows(correlations_environmental_coords, 
                                     correlations_environmental_resample,
                                     correlations_dem_resample, 
                                     correlations_dem_coords) %>% 
  tidyr::pivot_longer(!c(data, method)) %>% 
  dplyr::mutate(name = stringr::str_replace(name, pattern = "_", replacement = " "))
  
# order results by model run and coords method
levels <- correlation_full %>% 
  dplyr::filter(method == "Cell coords", name == "model run") %>% 
  dplyr::arrange(-value) %>% 
  dplyr::pull(data)
  
# make sure correct order for plotitng
correlation_full <- dplyr::mutate(correlation_full, 
                                  data = factor(correlation_full$data, levels = levels),
                                  method = factor(method, levels = c("Resample values", "Cell coords")), 
                                  name = factor(name, levels = c("model fitting", "model run")))
  
(ggplot_correlation <- ggplot(data = correlation_full) +
    geom_col(aes(x = data, y = value, fill = method), position = "dodge", 
             width = 0.75) + 
    geom_text(aes(x = data, y = value, label = round(value, 2), 
                  group = method, col = method), 
              position = position_dodge(0.75)) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(name = "Value matching", values = c("grey75", "black")) + 
    scale_color_manual(name = "Value matching", values = c("black", "grey75")) + 
    facet_wrap(~ name) +
    labs(x = "Environmental data", y = "Correlation value") +
    theme_classic() + 
    coord_flip())

suppoRt::save_ggplot(plot = ggplot_correlation,
                     filename = "ggplot_correlation.png",
                     path = "Figures/Appendix",
                     dpi = dpi, 
                     width = height_full, height = width_full, units = units,
                     overwrite = overwrite)  
