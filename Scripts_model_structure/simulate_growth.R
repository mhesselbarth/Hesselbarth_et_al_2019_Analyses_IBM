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


#### Plot potential growth ####
pot_growth <- purrr::map_dbl(0:80, function(x) {
  calculate_growth(dbh = x, parameters = parameters)
  })

plot_potential <- ggplot() +
  geom_line(aes(x = 0:80, y = pot_growth, col = "black")) +
  scale_color_manual(values = c("black" = "black"), name = "CI") + 
  labs(x = "DBH [cm]", y = "DBH increment [cm]") +
  scale_y_continuous(limits = c(0, 2.5)) + 
  theme_classic(base_size = 15)

#### Simulate acutal growth ####
# parameters$growth_mod <- 1

model_data <- update_i(data = input_data)
model_data <- rabmp::simulate_ci(data = model_data, parameters = parameters)
model_data <- rabmp::simulate_growth(data = input_data, parameters = parameters)

# model_data_long <- tidyr::unnest(model_data)

model_data <- tibble::tibble(dbh_old = model_data$dbh[model_data$i == 0],
                             dbh_new = model_data$dbh[model_data$i == 1], 
                             ci = model_data$ci[model_data$i == 1]) %>% 
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
suppoRt::save_ggplot(plot = plot_potential, 
                    filename = "ggplot_structure_potential.png", 
                    path = "Figures/Appendix", 
                    dpi = 300, width = 15, height = 7.5, units = "cm")

suppoRt::save_ggplot(plot = plot_actual, 
                    filename = "ggplot_structure_actual.png", 
                    path = "Figures/Appendix", 
                    dpi = 300, width = 15, height = 7.5, units = "cm")
