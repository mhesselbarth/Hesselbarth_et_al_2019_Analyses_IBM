###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Model structure Growth ####

#### Import libraries and data ####

# load packages
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(spatstat)
library(tidyverse)

# import parameters
parameters_beech_fitted <- rabmp::read_parameters("Data/Input/parameters_beech_fitted.txt")

parameters_beech_fitted$growth_mod <- 1

pattern_1999_recon <- readr::read_rds("Data/Input/beech_1999_rec.rds")

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

#### Plot potential growth ####
pot_growth <- purrr::map_dbl(0:125, function(x) {
  calculate_growth(dbh = x, parameters = parameters_beech_fitted)
  })

plot_potential <- ggplot() +
  geom_line(aes(x = 0:125, y = pot_growth, col = "Potential growth")) +
  scale_color_manual(values = c("Potential growth" = "black"), name = "") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 125)) +
  labs(x = "DBH [cm]", y = "DBH increment [cm]") +
  theme_classic(base_size = 15) + 
  theme(legend.position = "bottom")

#### Simulate acutal growth ####
model_data <- rabmp::update_i(data = input_data)
model_data <- rabmp::simulate_ci(data = model_data, parameters = parameters_beech_fitted)
model_data <- rabmp::simulate_growth(data = model_data, parameters = parameters_beech_fitted)

model_data <- tibble::tibble(dbh_old = model_data$dbh[model_data$i == 0],
                             dbh_new = model_data$dbh[model_data$i == 1], 
                             ci = model_data$ci[model_data$i == 1]) %>% 
  dplyr::mutate(dbh_inc = dbh_new - dbh_old)

plot_actual <- ggplot() + 
  geom_point(data = model_data, 
             aes(x = dbh_old, y = dbh_inc, col = ci), pch = 1, size = 2) +
  geom_line(aes(x = 0:125, y = pot_growth)) +
  scale_color_viridis_c(name = "CI", option = "A") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 125)) +
  labs(x = "DBH [cm]", y = "DBH increment [cm]") + 
  theme_classic(base_size = 15) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))

#### Save ggplot ####
suppoRt::save_ggplot(plot = plot_potential, 
                    filename = "ggplot_structure_potential.png", 
                    path = "Figures/Appendix", 
                    dpi = 300, width = 15, height = 7.5, units = "cm")

suppoRt::save_ggplot(plot = plot_actual, 
                    filename = "ggplot_structure_actual.png", 
                    path = "Figures/Appendix", 
                    dpi = 300, width = 15, height = 7.5, units = "cm")
