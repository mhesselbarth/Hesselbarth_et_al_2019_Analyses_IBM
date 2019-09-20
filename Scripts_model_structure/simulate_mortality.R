#### Import libraries and data ####

# load packages
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(spatstat)
library(tidyverse)

# import parameters
parameters <- rabmp::read_parameters("Data/Input/parameters_beech.txt")

# load data
input_data <- dplyr::filter(rabmp::example_input_data, 
                            spec == "beech", Class == "adult")

# prepare data for rabmp
input_data <- rabmp::prepare_data(data = input_data, 
                                  x = "x_coord", y = "y_coord",
                                  species = "spec", type = "Class", dbh = "bhd")

#### Mortality probs ####
mort_prob <- purrr::map_dbl(0:120, function(x) 
  rabmp:::rcpp_calculate_mortality_probs(species = "Beech", 
                                         dbh = x, 
                                         int_beech_early = parameters$mort_int_beech_early,
                                         dbh_beech_early = parameters$mort_dbh_beech_early,
                                         int_beech_late = parameters$mort_int_beech_late,
                                         dbh_beech_late = parameters$mort_dbh_beech_late,
                                         dinc_beech = parameters$mort_dinc_beech,
                                         int_ash = parameters$mort_int_ash,
                                         dbh_ash = parameters$mort_dbh_ash,
                                         int_others = parameters$mort_int_others,
                                         dbh_others = parameters$mort_dbh_others))

plot_mort_prob <- ggplot() + 
  geom_line(aes(x = 0:120, y = mort_prob)) + 
  labs(x = "DBH [cm]", y = "Mortality probability") + 
  theme_classic(base_size = 15) 

#### Save plots ####
suppoRt::save_ggplot(plot = plot_mort_prob, 
                    path = "Figures/Appendix",
                    filename = "ggplot_structure_mortality.png", 
                    dpi = 300, width = 15, height = 7.5, units = "cm")
