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
habitat_im_fit <- spatstat::density.ppp(beech_1999_ppp, at = "pixel", 
                                        weights = beech_1999_ppp$marks$dbh_99, 
                                        dimyx = c(645, 609),
                                        kernel = "epanechnikov", sigma = 5)

# number of rows added at edges #
n_rows <- 3

# convert to raster and add padding #
habitat_ras_fit <- tibble::as_tibble(habitat_im_fit) %>% 
  raster::rasterFromXYZ() %>% 
  landscapemetrics::pad_raster(pad_raster_value = NA, pad_raster_cells = n_rows) %>%
  magrittr::extract2(1)

# add subsequently number of rows at edge #
for (i in 1:n_rows) {
  
  message("\r> Progress: ", i, "/", n_rows, appendLF = "FALSE")
  
  # get all NA cells #
  cells_na <- raster::Which(is.na(habitat_ras_fit),
                            cells = TRUE)
  
  # get neigbors of NA cells #
  neighbours <- raster::adjacent(x = habitat_ras_fit,
                                 cells = cells_na,
                                 directions = 8)
  
  # get mean of all neighboring cells
  neighbours_value <- tibble::tibble(from = neighbours[, 1],
                                     to = neighbours[, 2],
                                     x_from = habitat_ras_fit[neighbours[, 1]], 
                                     x_to = habitat_ras_fit[neighbours[, 2]]) %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(x = mean(x_to, na.rm = TRUE)) %>% 
    dplyr::filter(!is.na(x))
  
  # add values
  habitat_ras_fit[neighbours_value$from] <- neighbours_value$x
}

# scale value to -1 to 1 #
habitat_ras_fit$scaled <- scales::rescale(raster::values(habitat_ras_fit), 
                                          to = c(-1, 1), na.rm = TRUE)

# set names #
names(habitat_ras_fit) <- c("absolute", "scaled")

#### Save data ####
suppoRt::save_rds(object = habitat_ras_fit, 
                  filename = "abiotic_cond_real_fit.rds", 
                  path = "Data/Input/", overwrite = overwrite)

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
habitat_im_model <- spatstat::density.ppp(beech_1999_ppp, at = "pixel", 
                                          weights = beech_1999_ppp$marks$dbh_99, 
                                          dimyx = c(645, 609),
                                          kernel = "epanechnikov", sigma = 15)

# number of rows added at edges #
n_rows <- 3

# convert to raster and add padding #
habitat_ras_model <- tibble::as_tibble(habitat_im_model) %>% 
  raster::rasterFromXYZ() %>% 
  landscapemetrics::pad_raster(pad_raster_value = NA, pad_raster_cells = n_rows) %>%
  magrittr::extract2(1)

# add subsequently number of rows at edge #
for (i in 1:n_rows) {
  
  message("\r> Progress: ", i, "/", n_rows, appendLF = "FALSE")
  
  # get all NA cells #
  cells_na <- raster::Which(is.na(habitat_ras_model),
                            cells = TRUE)
  
  # get neigbors of NA cells #
  neighbours <- raster::adjacent(x = habitat_ras_model,
                                 cells = cells_na,
                                 directions = 8)
  
  # get mean of all neighboring cells
  neighbours_value <- tibble::tibble(from = neighbours[, 1],
                                     to = neighbours[, 2],
                                     x_from = habitat_ras_model[neighbours[, 1]], 
                                     x_to = habitat_ras_model[neighbours[, 2]]) %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(x = mean(x_to, na.rm = TRUE)) %>% 
    dplyr::filter(!is.na(x))
  
  # add values
  habitat_ras_model[neighbours_value$from] <- neighbours_value$x
}

# scale value to -1 to 1 #
habitat_ras_model$scaled <- scales::rescale(raster::values(habitat_ras_model), 
                                            to = c(-1, 1), na.rm = TRUE)

# set names #
names(habitat_ras_model) <- c("absolute", "scaled")

#### Save data ####
suppoRt::save_rds(object = habitat_ras, 
                  filename = "abiotic_cond_real_model.rds", 
                  path = "Data/Input/", overwrite = overwrite)

##################################
####                          ####
#### Field data plot habitats ####
####                          ####
##################################

habitat_ras_full <- raster::as.data.frame(habitat_ras_fit, xy = TRUE) %>% 
  dplyr::mutate(data = "fitting")

habitat_ras_full <- raster::as.data.frame(habitat_ras_model, xy = TRUE) %>% 
  dplyr::mutate(data = "model") %>% 
  dplyr::bind_rows(habitat_ras_full, .) %>% 
  dplyr::mutate(data = factor(data, 
                              levels = c("fitting", "model"), 
                              labels = c("Competition and growth parameters", 
                                         "Seed establishment and mortality parameters")))

ggplot_abiotic_cond <- ggplot(data = habitat_ras_full) + 
  geom_raster(aes(x = x, y = y, fill = scaled)) + 
  geom_polygon(data = plot_area_df, aes(x = x, y = y), fill = NA, col = "black") + 
  facet_wrap(~ data) + 
  scale_fill_viridis_c(name = "Intensity", na.value = "white") + 
  coord_equal() + 
  guides(size = FALSE) + 
  theme_void(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"))

suppoRt::save_ggplot(plot = ggplot_abiotic_cond, 
                     file = "ggplot_abiotic_cond.png", 
                     path = "Figures/Appendix/", 
                     dpi = 300, units = units, 
                     width = width_full, height = height_small, 
                     overwrite = overwrite)

######################################
####                              ####
#### Reconstructed data fit data  ####
####                              ####
######################################

# # dbh threshold for habitat characterisation #
# dbh_threshold <- quantile(beech_1999_rec_ppp$marks$dbh, probs = 0.65)

# filter data using threshold sapling-adult #
beech_1999_rec_ppp <- spatstat::subset.ppp(beech_1999_rec_ppp, dbh > 10)

# get intensity
habitat_im_fit <- spatstat::density.ppp(beech_1999_rec_ppp, at = "pixel",
                                        weights = beech_1999_rec_ppp$marks$dbh,
                                        dimyx = c(645, 609),
                                        kernel = "epanechnikov", sigma = 5)

# number of rows added at edges #
n_rows <- 3

# convert to raster and add padding #
habitat_ras_fit <- tibble::as_tibble(habitat_im_fit) %>%
  raster::rasterFromXYZ() %>%
  landscapemetrics::pad_raster(pad_raster_value = NA, pad_raster_cells = n_rows) %>%
  magrittr::extract2(1)

# add subsequently number of rows at edge #
for (i in 1:n_rows) {

  message("\r> Progress: ", i, "/", n_rows, appendLF = "FALSE")

  # get all NA cells #
  cells_na <- raster::Which(is.na(habitat_ras_fit),
                            cells = TRUE)

  # get neigbors of NA cells #
  neighbours <- raster::adjacent(x = habitat_ras_fit,
                                 cells = cells_na,
                                 directions = 8)

  # get mean of all neighboring cells
  neighbours_value <- tibble::tibble(from = neighbours[, 1],
                                     to = neighbours[, 2],
                                     x_from = habitat_ras_fit[neighbours[, 1]],
                                     x_to = habitat_ras_fit[neighbours[, 2]]) %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(x = mean(x_to, na.rm = TRUE)) %>%
    dplyr::filter(!is.na(x))

  # add values
  habitat_ras_fit[neighbours_value$from] <- neighbours_value$x
}

# scale value to -1 to 1 #
habitat_ras_fit$scaled <- scales::rescale(raster::values(habitat_ras_fit),
                                      to = c(-1, 1), na.rm = TRUE)

# set names #
names(habitat_ras_fit) <- c("absolute", "scaled")

#### Save data ####

suppoRt::save_rds(object = habitat_ras_fit,
                  filename = "abiotic_cond_reco_fit.rds",
                  path = "Data/Input/", overwrite = overwrite)

######################################
####                              ####
#### Reconstructed data model run ####
####                              ####
######################################

# # dbh threshold for habitat characterisation #
# dbh_threshold <- quantile(beech_1999_rec_ppp$marks$dbh, probs = 0.65)

# filter data using threshold sapling-adult #
beech_1999_rec_ppp <- spatstat::subset.ppp(beech_1999_rec_ppp, dbh > 10)

# get intensity
habitat_im_model <- spatstat::density.ppp(beech_1999_rec_ppp, at = "pixel",
                                          weights = beech_1999_rec_ppp$marks$dbh,
                                          dimyx = c(645, 609),
                                          kernel = "epanechnikov", sigma = 15)

# number of rows added at edges #
n_rows <- 3

# convert to raster and add padding #
habitat_ras_model <- tibble::as_tibble(habitat_im_model) %>%
  raster::rasterFromXYZ() %>%
  landscapemetrics::pad_raster(pad_raster_value = NA, pad_raster_cells = n_rows) %>%
  magrittr::extract2(1)

# add subsequently number of rows at edge #
for (i in 1:n_rows) {

  message("\r> Progress: ", i, "/", n_rows, appendLF = "FALSE")

  # get all NA cells #
  cells_na <- raster::Which(is.na(habitat_ras_model),
                            cells = TRUE)

  # get neigbors of NA cells #
  neighbours <- raster::adjacent(x = habitat_ras_model,
                                 cells = cells_na,
                                 directions = 8)

  # get mean of all neighboring cells
  neighbours_value <- tibble::tibble(from = neighbours[, 1],
                                     to = neighbours[, 2],
                                     x_from = habitat_ras_model[neighbours[, 1]],
                                     x_to = habitat_ras_model[neighbours[, 2]]) %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(x = mean(x_to, na.rm = TRUE)) %>%
    dplyr::filter(!is.na(x))

  # add values
  habitat_ras_model[neighbours_value$from] <- neighbours_value$x
}

# scale value to -1 to 1 #
habitat_ras_model$scaled <- scales::rescale(raster::values(habitat_ras_model),
                                      to = c(-1, 1), na.rm = TRUE)

# set names #
names(habitat_ras_model) <- c("absolute", "scaled")

suppoRt::save_rds(object = habitat_ras_model,
                  filename = "abiotic_cond_reco_model.rds",
                  path = "Data/Input/", overwrite = overwrite)

##########################################
####                                  ####
#### Reconstructed data plot habitats ####
####                                  ####
##########################################

habitat_ras_full <- raster::as.data.frame(habitat_ras_fit, xy = TRUE) %>% 
  dplyr::mutate(data = "fitting")

habitat_ras_full <- raster::as.data.frame(habitat_ras_model, xy = TRUE) %>% 
  dplyr::mutate(data = "model") %>% 
  dplyr::bind_rows(habitat_ras_full, .) %>% 
  dplyr::mutate(data = factor(data, 
                              levels = c("fitting", "model"), 
                              labels = c("Competition and growth parameters", 
                                         "Seed establishment and mortality parameters")))

ggplot_abiotic_cond <- ggplot(data = habitat_ras_full) + 
  geom_raster(aes(x = x, y = y, fill = scaled)) + 
  geom_polygon(data = plot_area_df, aes(x = x, y = y), fill = NA, col = "black") + 
  facet_wrap(~ data) + 
  scale_fill_viridis_c(name = "Intensity", na.value = "white") + 
  coord_equal() + 
  guides(size = FALSE) + 
  theme_void(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"))

suppoRt::save_ggplot(plot = ggplot_abiotic_cond, 
                     file = "ggplot_abiotic_cond_reco.png", 
                     path = "Figures/Appendix/", 
                     dpi = 300, units = units, 
                     width = width_full, height = height_small, 
                     overwrite = overwrite)
