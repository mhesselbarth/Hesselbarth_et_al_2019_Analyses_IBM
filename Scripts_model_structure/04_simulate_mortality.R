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
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(spatstat)
library(tidyverse)

# import parameters
parameters_beech_fitted <- rabmp::read_parameters("Data/Input/parameters_beech_fitted.txt")

pattern_1999_recon <- readr::read_rds("Data/Input/beech_1999_rec_ppp.rds")

plot_area <- tibble::as_tibble(pattern_1999_recon$window)

#### Pre-processing of input data ####
set.seed(42)
sample_id <- sample(1:pattern_1999_recon$n, size = pattern_1999_recon$n)

input_data <- tibble::as_tibble(pattern_1999_recon) %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  dplyr::filter(species == "beech", id %in% sample_id) %>%
  dplyr::select(-id) %>% 
  rabmp::prepare_data(x = "x", y = "y", species = "species", type = "type", dbh = "dbh")

rm(pattern_1999_recon)

#### Mortality probs ####
mort_prob <- purrr::map_dbl(0:120, function(x) 
  rabmp:::rcpp_calculate_mortality_probs(species = "beech", 
                                         dbh = x, 
                                         int_beech_early = parameters_beech_fitted$mort_int_beech_early,
                                         dbh_beech_early = parameters_beech_fitted$mort_dbh_beech_early,
                                         int_beech_late = parameters_beech_fitted$mort_int_beech_late,
                                         dbh_beech_late = parameters_beech_fitted$mort_dbh_beech_late,
                                         dinc_beech = parameters_beech_fitted$mort_dinc_beech,
                                         int_ash = parameters_beech_fitted$mort_int_ash,
                                         dbh_ash = parameters_beech_fitted$mort_dbh_ash,
                                         int_others = parameters_beech_fitted$mort_int_others,
                                         dbh_others = parameters_beech_fitted$mort_dbh_others))

plot_mort_prob <- ggplot() + 
  geom_line(aes(x = 0:120, y = mort_prob)) + 
  labs(x = "DBH [cm]", y = "Mortality probability") + 
  theme_classic(base_size = 15) 

#### Save plots ####
overwrite <- FALSE

suppoRt::save_ggplot(plot = plot_mort_prob, 
                    path = "Figures/Appendix",
                    filename = "ggplot_structure_mortality.png", 
                    dpi = 300, height = 10, width = 12.5, units = "cm", 
                    overwrite = overwrite)
