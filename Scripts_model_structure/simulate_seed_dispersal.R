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

#### Plot Probability ####
distance_density <- rabmp:::rcpp_random_distance(number_seeds = 1000000, 
                                                 species = "Beech", max_dist = 120) %>% 
  tibble::tibble(prob = .)

plot_dist_density <- ggplot(data = distance_density) + 
  geom_density(aes(prob)) +   
  labs(x = "Distance [m]", y = "Density") + 
  theme_classic(base_size = 15)



#### Plot point pattern ####

set.seed(42)
id <- sample(1:nrow(input_data), size = 20)

seedlings <- simulate_seed_dispersal(data = input_data[id, ], parameters = parameters, 
                                     plot_area = spatstat::owin(xrange = c(0, 500), 
                                                                yrange = c(0, 500))) %>%
  tidyr::unnest()

plot_seed_pattern <- ggplot(data = seedlings) +  
  geom_rect(aes(xmin = 0, xmax = 500, ymin = 0, ymax = 500), 
            fill = NA, col = "black") + 
  geom_point(aes(x = x, y = y, shape = type, size = dbh)) + 
  scale_size_continuous(name = "DBH [cm]") +
  scale_shape_manual(values = c(19, 1), name = "Life stage") +

  coord_equal() + 
  theme_void(base_size = 15)


#### Save plots #### 
helpeR::save_ggplot(plot = plot_dist_density, 
                    filename = "plot_dist_density.png", 
                    path = "Figures/", 
                    dpi = 300, width = 15, height = 7.5, units = "cm")

helpeR::save_ggplot(plot = plot_seed_pattern, 
                    filename = "plot_seed_pattern.png", 
                    path = "Figures/", 
                    dpi = 300, width = 15, height = 15, units = "cm")


