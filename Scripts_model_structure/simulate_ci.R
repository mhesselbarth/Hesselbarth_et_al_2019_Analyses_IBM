#### Import libraries and data ####

# load packages
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(spatstat)
library(tidyverse)

# import parameters
parameters_beech_fitted <- rabmp::read_parameters("Data/Input/parameters_beech_fitted.txt")

pattern_1999_recon <- readr::read_rds("Data/Input/beech_1999_rec.rds")

plot_area <- tibble::as_tibble(pattern_1999_recon$window)

#### Pre-processing of input data ####
set.seed(42)
sample_id <- sample(1:pattern_1999_recon$n, size = 250)

input_data <- tibble::as_tibble(pattern_1999_recon) %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  dplyr::filter(species == "beech", id %in% sample_id) %>%
  dplyr::select(-id) %>% 
  rabmp::prepare_data(x = "x", y = "y", species = "species", type = "type", dbh = "dbh")

rm(pattern_1999_recon)

#### Plot kernel ####
# calculate CI index dbh = 20
ci_index_20 <- purrr::map_dbl(0:20, function(x) 
  (20 ^ parameters_beech_fitted$ci_alpha) * exp(-(x / (20 ^ parameters_beech_fitted$ci_beta))))
# ci_index_20 <- ci_index_20 / (5 ^ parameters$ci_alpha + ci_index_20)

# calculate CI index dbh = 10 
ci_index_10 <- purrr::map_dbl(0:20, function(x) 
  (10 ^ parameters_beech_fitted$ci_alpha) * exp(-(x / (10 ^ parameters_beech_fitted$ci_beta))))
# ci_index_10 <- ci_index_10 / (5 ^ parameters$ci_alpha + ci_index_10)

# create ggplot
plot_kernel <- ggplot() + 
  geom_line(aes(x = -20:20, y = c(rev(ci_index_10), ci_index_10[-1]),
            linetype = "DBH 10")) + 
  geom_line(aes(x = -20:20, y = c(rev(ci_index_20), ci_index_20[-1]), 
            linetype = "DBH 20")) + 
  scale_linetype_manual(name = "DBH [cm]", 
                        values = c("DBH 10" = 2, "DBH 20" = 1)) + 
  scale_x_continuous(breaks = seq(-20, 20, 5),
                     labels = abs(seq(-20, 20, 5))) + 
  labs(x = "Distance", y = "Exponential competition kernel") + 
  theme_classic(base_size = 15)

#### Plot point pattern ####
# calculat ci
data_ci <- rabmp::simulate_ci(input_data, parameters = parameters_beech_fitted) %>% 
  tibble::as_tibble()

# create ggplot without ci
plot_pattern <- ggplot(data = data_ci) + 
  geom_polygon(data = plot_area, aes(x = x, y = y), 
               col = "black", fill = NA) + 
  geom_point(aes(x = x, y = y, size = dbh), pch = 1) + 
  scale_size_continuous(name = "DBH [cm]") + 
  coord_equal() + 
  theme_void(base_size = 15)

# create ggplot with ci
plot_pattern_ci <- ggplot(data = data_ci) + 
  geom_polygon(data = plot_area, aes(x = x, y = y), 
               col = "black", fill = NA) + 
  geom_point(aes(x = x, y = y, col = ci, size = dbh), pch = 1) + 
  scale_size_continuous(name = "DBH [cm]") + 
  scale_color_viridis_c(name = "Competition index", option = "A") + 
  coord_equal() + 
  theme_void(base_size = 15)

#### Save plots #### 
suppoRt::save_ggplot(plot = plot_kernel,
                    filename = "ggplot_structure_kernel.png",
                    path = "Figures/Appendix",
                    dpi = 300, height = 7.5, width = 15, units = "cm")

suppoRt::save_ggplot(plot = plot_pattern,
                    filename = "ggplot_structure_pattern.png",
                    path = "Figures/Appendix",
                    dpi = 300, height = 15, width = 15, units = "cm")

suppoRt::save_ggplot(plot = plot_pattern_ci,
                    filename = "ggplot_structure_pattern_ci.png",
                    path = "Figures/Appendix",
                    dpi = 300, height = 15, width = 15, units = "cm")

