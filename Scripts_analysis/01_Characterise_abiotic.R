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
library(magrittr)
library(landscapemetrics)
library(rabmp)
library(raster)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(spatstat)
library(tidyverse)

source("Helper_functions/helper_functions_abiotic_conditions.R")

# import data  #
pattern_1999_ppp <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

beech_1999_rec_ppp <- readr::read_rds("Data/Input/beech_1999_rec_ppp.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

plot_area_df <- as.data.frame(plot_area)

####################
####            ####
#### Field data ####
####            ####
####################

# filter data #
beech_1999_ppp <- spatstat::subset.ppp(pattern_1999_ppp, 
                                       species == "beech" &  type != "dead")

# dbh threshold for habitat characterisation #
dbh_threshold <- quantile(beech_1999_ppp$marks$dbh_99, probs = 0.95)

# filter data #
beech_1999_ppp <- spatstat::subset.ppp(beech_1999_ppp, dbh_99 > dbh_threshold)

# get intensity
habitat_im <- spatstat::density.ppp(beech_1999_ppp, at = "pixel", 
                                    weights = beech_1999_ppp$marks$dbh_99, 
                                    dimyx = c(645, 609),
                                    kernel = "epanechnikov", sigma = 75)

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

# scale value to 0 - 1 #
habitat_ras$scaled <- habitat_ras[] / max(habitat_ras[], na.rm = TRUE)

# set names #
names(habitat_ras) <- c("absolute", "scaled")

ggplot(data = raster::as.data.frame(habitat_ras)) +
  geom_density(aes(scaled), fill = "#440154FF", alpha = 0.3) + 
  geom_vline(aes(xintercept = mean(scaled, na.rm = TRUE)),
             linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1)) + 
  labs(x = "Scaled intensity value", y = "Density") +
  theme_classic()

ggplot_abiotic_cond <- ggplot(data = raster::as.data.frame(habitat_ras, xy = TRUE)) + 
  geom_raster(aes(x = x, y = y, fill = scaled)) + 
  geom_polygon(data = plot_area_df, aes(x = x, y = y), fill = NA, col = "black") + 
  geom_point(data = tibble::as_tibble(beech_1999_ppp),
             aes(x = x, y = y, size = dbh_99), pch = 1) +
  scale_fill_viridis_c(name = "Intensity", na.value = "white") + 
  coord_equal() + 
  guides(size = FALSE) + 
  theme_void(base_size = 15) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"))

#### Save data ####
overwrite <- FALSE

suppoRt::save_rds(object = habitat_ras, 
                  filename = "abiotic_cond_real.rds", 
                  path = "Data/Input/", overwrite = overwrite)

suppoRt::save_ggplot(plot = ggplot_abiotic_cond, 
                     filename = "ggplot_abiotic_cond_real.png", 
                     path = "Figures/", overwrite = overwrite,
                     dpi = 300, width = 15, height = 15, units = "cm")


############################
####                    ####
#### Reconstructed data ####
####                    ####
############################

# filter data #
beech_1999_rec_ppp <- spatstat::subset.ppp(beech_1999_rec_ppp, 
                                           species == "beech" &  type != "dead")

# dbh threshold for habitat characterisation #
dbh_threshold <- quantile(beech_1999_rec_ppp$marks$dbh, probs = 0.95)

# filter data #
beech_1999_rec_ppp <- spatstat::subset.ppp(beech_1999_rec_ppp, dbh > dbh_threshold)

# get intensity
habitat_im <- spatstat::density.ppp(beech_1999_rec_ppp, at = "pixel", 
                                    weights = beech_1999_rec_ppp$marks$dbh, 
                                    dimyx = c(645, 609),
                                    kernel = "epanechnikov", sigma = 75)

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

# scale value to 0 - 1 #
habitat_ras$scaled <- habitat_ras[] / max(habitat_ras[], na.rm = TRUE)

# set names #
names(habitat_ras) <- c("absolute", "scaled")

ggplot(data = raster::as.data.frame(habitat_ras)) +
  geom_density(aes(scaled), fill = "#440154FF", alpha = 0.3) + 
  geom_vline(aes(xintercept = mean(scaled, na.rm = TRUE)),
             linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1)) + 
  labs(x = "Scaled intensity value", y = "Density") +
  theme_classic()

ggplot_abiotic_cond <- ggplot(data = raster::as.data.frame(habitat_ras, xy = TRUE)) + 
  geom_raster(aes(x = x, y = y, fill = scaled)) + 
  geom_polygon(data = plot_area_df, aes(x = x, y = y), fill = NA, col = "black") + 
  geom_point(data = tibble::as_tibble(beech_1999_rec_ppp),
             aes(x = x, y = y, size = dbh), pch = 1) +
  scale_fill_viridis_c(name = "Intensity", na.value = "white") + 
  coord_equal() + 
  guides(size = FALSE) + 
  theme_void(base_size = 15) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"))

#### Save data ####
overwrite <- FALSE

suppoRt::save_rds(object = habitat_ras, 
                  filename = "abiotic_cond_reco.rds", 
                  path = "Data/Input/", overwrite = overwrite)

suppoRt::save_ggplot(plot = ggplot_abiotic_cond, 
                     filename = "ggplot_abiotic_cond_reco.png", 
                     path = "Figures/", overwrite = overwrite,
                     dpi = 300, width = 15, height = 15, units = "cm")