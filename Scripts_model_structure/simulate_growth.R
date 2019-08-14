#### Import libraries and data ####

# load packages
library(helpeR) # devtools::install_github("mhesselbarth/helpeR")
library(rabmp)
library(spatstat)
library(tidyverse)

# import parameters
parameters <- rabmp::read_parameters("Data/parameters.txt")

# load data
input_data <- dplyr::filter(rabmp::example_input_data, 
                            spec == "Beech", Class == "Adult")

# prepare data for rabmp
input_data <- rabmp::prepare_data(data = input_data, 
                                  x = "x_coord", y = "y_coord",
                                  species = "spec", type = "Class", dbh = "bhd")


#### Plot potential growth ####
pot_growth <- purrr::map_dbl(0:80, function(x) calculate_growth(dbh = x, parameters = parameters))

plot_potential <- ggplot() +
  geom_line(aes(x = 0:80, y = pot_growth, col = "black")) +
  scale_color_manual(values = c("black" = "black"), name = "CI") + 
  labs(x = "DBH [cm]", y = "DBH increment [cm]") +
  theme_classic(base_size = 15)

#### Simulate acutal growth ####
parameters$growth_mod <- 1

model_data <- update_i(data = input_data)
model_data <- rabmp::simulate_ci(data = model_data, parameters = parameters)
model_data <- rabmp::simulate_growth(data = model_data, parameters = parameters)

model_data_long <- tidyr::unnest(model_data)

model_data_mod <- tibble::tibble(dbh_old = model_data_long$dbh[model_data_long$i == 0],
                                 dbh_new = model_data_long$dbh[model_data_long$i == 1], 
                                 ci = model_data_long$ci[model_data_long$i == 1]) %>% 
  dplyr::mutate(dbh_inc = dbh_new - dbh_old)

plot_actual <- ggplot() + 
  geom_point(data = model_data_mod, 
             aes(x = dbh_old, y = dbh_inc, col = ci)) +
  geom_line(aes(x = 0:80, y = pot_growth)) +
  scale_color_viridis_c(name = "CI", option = "A") +
  lims(x = c(0, 80)) +
  labs(x = "DBH [cm]", y = "DBH increment [cm]") + 
  theme_classic(base_size = 15)

#### Save ggplot ####
helpeR::save_ggplot(plot = plot_potential, 
                    filename = "plot_potential.png", 
                    path = "Figures/", 
                    dpi = 300, width = 15, height = 7.5, units = "cm")

helpeR::save_ggplot(plot = plot_actual, 
                    filename = "plot_growth.png", 
                    path = "Figures/", 
                    dpi = 300, width = 15, height = 7.5, units = "cm")
