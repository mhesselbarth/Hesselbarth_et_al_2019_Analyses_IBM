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
library(rabmp)
library(raster)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(spatstat)
library(tidyverse)

source("Helper_functions/helper_functions_abiotic_conditions.R")

# import data  #
pattern_2013_df <- readr::read_rds("Data/Raw/pattern_2013_df.rds")

pattern_1999_ppp <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

plot_area_df <- as.data.frame(plot_area)

# filter data and calculate mean dbh growth #
beech_2013_df <- dplyr::filter(pattern_2013_df, 
                               species == "beech", 
                               !is.na(dbh_99), 
                               !is.na(dbh_13),
                               type == "living", 
                               inside_fence == 0) %>% 
  dplyr::mutate(growth_mean = (growth_07 + growth_13) / 2,
                growth_full = (dbh_13 - dbh_99) / 14) %>% 
  dplyr::filter(growth_mean >= 0,
                growth_full >= 0)

beech_2013_dt <- dplyr::select(beech_2013_df, x, y, dbh_99, type) %>% 
  dplyr::mutate(type = "adult") %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh_99")


# filter data #
beech_1999_ppp <- spatstat::subset.ppp(pattern_1999_ppp, 
                                       species == "beech")

# dbh threshold for habitat characterisation#
dbh_threshold <- quantile(beech_1999_ppp$marks$dbh_99, probs = 0.95)

# filter data #
beech_1999_ppp <- spatstat::subset.ppp(beech_1999_ppp, type != "dead" & 
                                         dbh_99 > dbh_threshold)

# get abiotic conditions 
habitats_im <- spatstat::density.ppp(beech_1999_ppp, at = "pixel", 
                                     weights = beech_1999_ppp$marks$dbh_99, 
                                     dimyx = c(645, 609),
                                     kernel = "epanechnikov", sigma = 75)

habitats_ras <- tibble::as_tibble(habitats_im) %>% 
  raster::rasterFromXYZ()

habitats_ras <- add_padding(habitats_ras)

ggplot_abiotic_cond <- ggplot(data = raster::as.data.frame(habitats_ras, xy = TRUE)) + 
  geom_raster(aes(x = x, y = y, fill = value)) + 
  geom_polygon(data = plot_area_df, aes(x = x, y = y), fill = NA, col = "black") + 
  scale_fill_viridis_c(name = "Intensity", na.value = "white") + 
  coord_equal() + 
  theme_void(base_size = 15) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"))

#### Save data ####
suppoRt::save_rds(object = habitats_ras, 
                  filename = "abiotic_cond.rds", 
                  path = "Data/Input/")

suppoRt::save_ggplot(plot = ggplot_abiotic_cond, 
                     filename = "ggplot_abiotic_cond.png", 
                     path = "Figures/", 
                     dpi = 300, width = 15, height = 15, units = "cm")
