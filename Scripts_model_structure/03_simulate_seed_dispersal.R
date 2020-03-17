###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Model structure Seed dispersal ####

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

plot_area <- pattern_1999_recon$window

#### Pre-processing of input data ####
set.seed(42)
sample_id <- sample(1:pattern_1999_recon$n, size = 100)

input_data <- tibble::as_tibble(pattern_1999_recon) %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  dplyr::filter(species == "beech", id %in% sample_id) %>%
  dplyr::select(-id) %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh")

rm(pattern_1999_recon)

#### Plot Probability ####
distance_density <- rabmp:::rcpp_random_distance(number_seeds = 1000000, 
                                                 beta = parameters_beech_fitted$seed_beta,
                                                 max_dist = parameters_beech_fitted$seed_max_dist) %>% 
  tibble::tibble(prob = .)

plot_dist_density <- ggplot(data = distance_density) + 
  geom_density(aes(prob)) +   
  labs(x = "Distance [m]", y = "Density") + 
  theme_classic(base_size = 15) +
  theme(text = element_text(family  = "Calibri Light"))

suppoRt::save_ggplot(plot = plot_dist_density, 
                     filename = "plot_dist_density.png", 
                     path = "C:/Users/Maximilian/ownCloud/13_Thesis_defense/Figures/",
                     dpi = dpi, units = units,
                     width = 50, height = 50, 
                     overwrite = FALSE)

#### Plot point pattern ####
seedlings <- rabmp::simulate_seed_dispersal_biotic(data = input_data, 
                                     parameters = parameters_beech_fitted, 
                                     plot_area = plot_area) %>%
  tibble::as_tibble() %>% 
  dplyr::mutate(type = dplyr::case_when(type %in% c("adult", "sapling") ~ "adult", 
                                        type == "seedling" ~ "seedling"))

plot_seed_pattern <- ggplot(data = seedlings) +  
  geom_point(aes(x = x, y = y, shape = type, size = dbh)) + 
  geom_polygon(data = tibble::as_tibble(plot_area), aes(x = x, y = y), 
               col = "black", fill = NA) + 
  scale_size_continuous(name = "DBH [cm]") +
  scale_shape_manual(values = c(1, 3), name = "Life stage") +
  guides(size = FALSE) +
  coord_equal() + 
  theme_void(base_size = 15) + 
  theme(legend.position = "bottom", 
        text = element_text(family  = "Calibri Light"))

#### Save plots #### 
overwrite <- FALSE

ggplot_overall <- plot_dist_density + plot_seed_pattern + 
  patchwork::plot_layout(ncol = 2, nrow = 1,
                         widths = c(0.5, 0.5), heights = c(0.5, 0.5)) + 
  patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")", 
                             theme = theme(text = element_text(family = "Calibri Light")))

suppoRt::save_ggplot(plot = ggplot_overall, 
                     filename = "ggplot_seed_overall.png", 
                     path = "C:/Users/Maximilian/ownCloud/13_Disputation/Figures/",
                     dpi = dpi, units = units,
                     width = 200, height = 125, 
                     overwrite = FALSE)


suppoRt::save_ggplot(plot = plot_dist_density, 
                     filename = "ggplot_structure_seed_density.png", 
                     path = "Figures/Appendix", 
                     dpi = 300, height = 10, width = 12.5, units = "cm", 
                     overwrite = overwrite)

suppoRt::save_ggplot(plot = plot_seed_pattern, 
                     filename = "ggplot_structure_seed_pattern.png", 
                     path = "Figures/Appendix", 
                     dpi = 300, height = 12.5, width = 12.5, units = "cm", 
                     overwrite = overwrite)
