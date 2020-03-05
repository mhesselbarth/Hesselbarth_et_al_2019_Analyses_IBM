###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Model structure Mortality ####

#### Import libraries and data ####

# load packages
library(extrafont)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(spatstat)
library(tidyverse)

# import parameters
parameters_beech_fitted <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt", 
                                                  sep = ";")

pattern_1999_recon <- readr::read_rds("Data/Input/beech_1999_ppp_rec.rds")

plot_area <- tibble::as_tibble(pattern_1999_recon$window)

#### Pre-processing of input data ####
set.seed(42)
sample_id <- sample(1:pattern_1999_recon$n, size = pattern_1999_recon$n)

input_data <- tibble::as_tibble(pattern_1999_recon) %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  dplyr::filter(species == "beech", id %in% sample_id) %>%
  dplyr::select(-id, -species) %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh")

rm(pattern_1999_recon)

#### Mortality probs ####
mort_prob <- purrr::map_dbl(0:120, function(x) 
  rabmp:::rcpp_calculate_mortality_probs(dbh = x, 
                                         int_early = parameters_beech_fitted$mort_int_early,
                                         dbh_early = parameters_beech_fitted$mort_dbh_early,
                                         int_late = parameters_beech_fitted$mort_int_late,
                                         dbh_late = parameters_beech_fitted$mort_dbh_late,
                                         dinc = parameters_beech_fitted$mort_dinc))

plot_mort_prob <- ggplot() + 
  geom_line(aes(x = 0:120, y = mort_prob)) + 
  labs(x = "DBH [cm]", y = "Mortality probability") + 
  theme_classic(base_size = 15) +
  theme(text = element_text(family = "Calibri Light"))

#### Save plots ####
overwrite <- FALSE

suppoRt::save_ggplot(plot = plot_mort_prob, 
                    filename = "ggplot_structure_mortality.png", 
                    path = "C:/Users/Maximilian/ownCloud/13_Disputation/Figures/",
                    dpi = 300, units = "mm",
                    width = 200, height = 125, 
                    overwrite = FALSE)
