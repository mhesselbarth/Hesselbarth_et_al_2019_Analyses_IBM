###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Create abiotic conditions ####

#### Import libraries and data ####

# load packages #
source("Helper_functions/helper_functions_setup.R")

# import data  #
beech_1999_ppp <- readr::read_rds("Data/Input/beech_1999_ppp.rds")

beech_1999_rec_ppp <- readr::read_rds("Data/Input/beech_1999_rec_ppp.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

plot_area_df <- as.data.frame(plot_area)

############################
####                    ####
#### Field data fitting ####
####                    ####
############################

# # dbh threshold for habitat characterisation #
# dbh_threshold <- quantile(beech_1999_ppp$marks$dbh_99, probs = 0.65)

# filter data using threshold sapling/adult #
beech_1999_ppp <- spatstat::subset.ppp(beech_1999_ppp, dbh_99 > 10)

# get intensity
habitat_im <- spatstat::density.ppp(beech_1999_ppp, at = "pixel", 
                                    weights = beech_1999_ppp$marks$dbh_99, 
                                    dimyx = c(645, 609),
                                    kernel = "epanechnikov", sigma = 5)

# number of rows added at edges #
n_rows <- 3

# convert to raster and add padding #
habitat_ras <- tibble::as_tibble(habitat_im) %>% 
  raster::rasterFromXYZ() %>% 
  landscapemetrics::pad_raster(pad_raster_value = NA, pad_raster_cells = n_rows) %>%
  magrittr::extract2(1)

# add subsequently number of rows at edge #
for (i in 1:n_rows) {
  
  message("\r> Progress: ", i, "/", n_rows, appendLF = "FALSE")
  
  # get all NA cells #
  cells_na <- raster::Which(is.na(habitat_ras),
                            cells = TRUE)
  
  # get neigbors of NA cells #
  neighbours <- raster::adjacent(x = habitat_ras,
                                 cells = cells_na,
                                 directions = 8)
  
  # get mean of all neighboring cells
  neighbours_value <- tibble::tibble(from = neighbours[, 1],
                                     to = neighbours[, 2],
                                     x_from = habitat_ras[neighbours[, 1]], 
                                     x_to = habitat_ras[neighbours[, 2]]) %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(x = mean(x_to, na.rm = TRUE)) %>% 
    dplyr::filter(!is.na(x))
  
  # add values
  habitat_ras[neighbours_value$from] <- neighbours_value$x
}

# scale value to -1 to 1 #
habitat_ras$scaled <- scales::rescale(raster::values(habitat_ras), 
                                      to = c(-1, 1), na.rm = TRUE)

# set names #
names(habitat_ras) <- c("absolute", "scaled")

# ggplot(data = raster::as.data.frame(habitat_ras)) +
#   geom_density(aes(scaled), fill = "#440154FF", alpha = 0.3) + 
#   geom_vline(aes(xintercept = median(scaled, na.rm = TRUE)),
#              linetype = "dashed") +
#   geom_vline(aes(xintercept = quantile(scaled, probs = 0.1, na.rm = TRUE)),
#              linetype = "dashed") +
#   geom_vline(aes(xintercept = quantile(scaled, probs = 0.9, na.rm = TRUE)),
#              linetype = "dashed") +
#   scale_x_continuous(limits = c(-1, 1)) +
#   labs(x = "Scaled intensity value", y = "Density") +
#   theme_classic()

ggplot_abiotic_cond <- ggplot(data = raster::as.data.frame(habitat_ras, xy = TRUE)) + 
  geom_raster(aes(x = x, y = y, fill = scaled)) + 
  geom_polygon(data = plot_area_df, aes(x = x, y = y), fill = NA, col = "black") + 
  scale_fill_viridis_c(name = "Intensity", na.value = "white") + 
  coord_equal() + 
  guides(size = FALSE) + 
  theme_void(base_size = 10) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"))

#### Save data ####
suppoRt::save_rds(object = habitat_ras, 
                  filename = "abiotic_cond_real_fit.rds", 
                  path = "Data/Input/", overwrite = overwrite)

suppoRt::save_ggplot(plot = ggplot_abiotic_cond, 
                     filename = "ggplot_abiotic_cond_real_fit.png", 
                     path = "Figures/", overwrite = overwrite,
                     dpi = dpi, 
                     width = width_small, height = height_small, units = units)

############################
####                    ####
#### Reconstructed data ####
####                    ####
############################

# # dbh threshold for habitat characterisation #
# dbh_threshold <- quantile(beech_1999_rec_ppp$marks$dbh, probs = 0.65)

# filter data using threshold sapling-adult #
beech_1999_rec_ppp <- spatstat::subset.ppp(beech_1999_rec_ppp, dbh > 10)

# get intensity
habitat_im <- spatstat::density.ppp(beech_1999_rec_ppp, at = "pixel", 
                                    weights = beech_1999_rec_ppp$marks$dbh, 
                                    dimyx = c(645, 609),
                                    kernel = "epanechnikov", sigma = 5)

# number of rows added at edges #
n_rows <- 3

# convert to raster and add padding #
habitat_ras <- tibble::as_tibble(habitat_im) %>% 
  raster::rasterFromXYZ() %>% 
  landscapemetrics::pad_raster(pad_raster_value = NA, pad_raster_cells = n_rows) %>%
  magrittr::extract2(1)

# add subsequently number of rows at edge #
for (i in 1:n_rows) {
  
  message("\r> Progress: ", i, "/", n_rows, appendLF = "FALSE")
  
  # get all NA cells #
  cells_na <- raster::Which(is.na(habitat_ras),
                            cells = TRUE)
  
  # get neigbors of NA cells #
  neighbours <- raster::adjacent(x = habitat_ras,
                                 cells = cells_na,
                                 directions = 8)
  
  # get mean of all neighboring cells
  neighbours_value <- tibble::tibble(from = neighbours[, 1],
                                     to = neighbours[, 2],
                                     x_from = habitat_ras[neighbours[, 1]], 
                                     x_to = habitat_ras[neighbours[, 2]]) %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(x = mean(x_to, na.rm = TRUE)) %>% 
    dplyr::filter(!is.na(x))
  
  # add values
  habitat_ras[neighbours_value$from] <- neighbours_value$x
}

# scale value to -1 to 1 #
habitat_ras$scaled <- scales::rescale(raster::values(habitat_ras), 
                                      to = c(-1, 1), na.rm = TRUE)

# set names #
names(habitat_ras) <- c("absolute", "scaled")

# plot results #
# ggplot(data = raster::as.data.frame(habitat_ras)) +
#   geom_density(aes(scaled), fill = "#440154FF", alpha = 0.3) + 
#   geom_vline(aes(xintercept = median(scaled, na.rm = TRUE)),
#              linetype = "dashed") +
#   geom_vline(aes(xintercept = quantile(scaled, probs = 0.1, na.rm = TRUE)),
#              linetype = "dashed") +
#   geom_vline(aes(xintercept = quantile(scaled, probs = 0.9, na.rm = TRUE)),
#              linetype = "dashed") +
#   scale_x_continuous(limits = c(-1, 1)) + 
#   labs(x = "Scaled intensity value", y = "Density") +
#   theme_classic()

ggplot_abiotic_cond <- ggplot(data = raster::as.data.frame(habitat_ras, xy = TRUE)) + 
  geom_raster(aes(x = x, y = y, fill = scaled)) + 
  geom_polygon(data = plot_area_df, aes(x = x, y = y), fill = NA, col = "black") + 
  scale_fill_viridis_c(name = "Intensity", na.value = "white") + 
  coord_equal() + 
  guides(size = FALSE) + 
  theme_void(base_size = 15) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"))

#### Save data ####

suppoRt::save_rds(object = habitat_ras, 
                  filename = "abiotic_cond_reco.rds", 
                  path = "Data/Input/", overwrite = overwrite)

# suppoRt::save_ggplot(plot = ggplot_abiotic_cond,
#                      filename = "ggplot_abiotic_cond_reco.png",
#                      path = "Figures/", overwrite = overwrite,
#                      dpi = dpi, 
#                      width = width_small, height = height_small, units = units)

##############################
####                      ####
#### Field data model run ####
####                      ####
##############################

# # dbh threshold for habitat characterisation #
# dbh_threshold <- quantile(beech_1999_ppp$marks$dbh_99, probs = 0.65)

# filter data using threshold sapling/adult #
beech_1999_ppp <- spatstat::subset.ppp(beech_1999_ppp, dbh_99 > 10)

# get intensity
habitat_im <- spatstat::density.ppp(beech_1999_ppp, at = "pixel", 
                                    weights = beech_1999_ppp$marks$dbh_99, 
                                    dimyx = c(645, 609),
                                    kernel = "epanechnikov", sigma = 15)

# number of rows added at edges #
n_rows <- 3

# convert to raster and add padding #
habitat_ras <- tibble::as_tibble(habitat_im) %>% 
  raster::rasterFromXYZ() %>% 
  landscapemetrics::pad_raster(pad_raster_value = NA, pad_raster_cells = n_rows) %>%
  magrittr::extract2(1)

# add subsequently number of rows at edge #
for (i in 1:n_rows) {
  
  message("\r> Progress: ", i, "/", n_rows, appendLF = "FALSE")
  
  # get all NA cells #
  cells_na <- raster::Which(is.na(habitat_ras),
                            cells = TRUE)
  
  # get neigbors of NA cells #
  neighbours <- raster::adjacent(x = habitat_ras,
                                 cells = cells_na,
                                 directions = 8)
  
  # get mean of all neighboring cells
  neighbours_value <- tibble::tibble(from = neighbours[, 1],
                                     to = neighbours[, 2],
                                     x_from = habitat_ras[neighbours[, 1]], 
                                     x_to = habitat_ras[neighbours[, 2]]) %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(x = mean(x_to, na.rm = TRUE)) %>% 
    dplyr::filter(!is.na(x))
  
  # add values
  habitat_ras[neighbours_value$from] <- neighbours_value$x
}

# scale value to -1 to 1 #
habitat_ras$scaled <- scales::rescale(raster::values(habitat_ras), 
                                      to = c(-1, 1), na.rm = TRUE)

# set names #
names(habitat_ras) <- c("absolute", "scaled")

# ggplot(data = raster::as.data.frame(habitat_ras)) +
#   geom_density(aes(scaled), fill = "#440154FF", alpha = 0.3) + 
#   geom_vline(aes(xintercept = median(scaled, na.rm = TRUE)),
#              linetype = "dashed") +
#   geom_vline(aes(xintercept = quantile(scaled, probs = 0.1, na.rm = TRUE)),
#              linetype = "dashed") +
#   geom_vline(aes(xintercept = quantile(scaled, probs = 0.9, na.rm = TRUE)),
#              linetype = "dashed") +
#   scale_x_continuous(limits = c(-1, 1)) +
#   labs(x = "Scaled intensity value", y = "Density") +
#   theme_classic()

ggplot_abiotic_cond <- ggplot(data = raster::as.data.frame(habitat_ras, xy = TRUE)) + 
  geom_raster(aes(x = x, y = y, fill = scaled)) + 
  geom_polygon(data = plot_area_df, aes(x = x, y = y), fill = NA, col = "black") + 
  scale_fill_viridis_c(name = "Intensity", na.value = "white") + 
  coord_equal() + 
  guides(size = FALSE) + 
  theme_void(base_size = 10) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"))

#### Save data ####
suppoRt::save_rds(object = habitat_ras, 
                  filename = "abiotic_cond_real_model.rds", 
                  path = "Data/Input/", overwrite = overwrite)

suppoRt::save_ggplot(plot = ggplot_abiotic_cond, 
                     filename = "ggplot_abiotic_cond_real_model.png", 
                     path = "Figures/", overwrite = overwrite,
                     dpi = dpi, 
                     width = width_small, height = height_small, units = units)
