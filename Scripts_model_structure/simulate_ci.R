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

#### Plot kernel ####
# calculate CI index dbh = 20
ci_index_20 <- purrr::map_dbl(0:20, function(x) (20 ^ parameters$ci_alpha) * exp(-(x / (20 ^ parameters$ci_beta))))
# ci_index_20 <- ci_index_20 / (5 ^ parameters$ci_alpha + ci_index_20)

# calculate CI index dbh = 10 
ci_index_10 <- purrr::map_dbl(0:20, function(x) (10 ^ parameters$ci_alpha) * exp(-(x / (10 ^ parameters$ci_beta))))
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
data_ci <- simulate_ci(input_data, parameters = parameters) %>% 
  tidyr::unnest()

# create ggplot without ci
plot_pattern <- ggplot(data = data_ci) + 
  geom_rect(aes(xmin = 0, xmax = 500, ymin = 0, ymax = 500), 
            fill = NA, col = "black") + 
  geom_point(aes(x = x, y = y, size = dbh), pch = 1) + 
  scale_size_continuous(name = "DBH [cm]") + 
  coord_equal() + 
  theme_void(base_size = 15)

# create ggplot with ci
plot_pattern_ci <- ggplot(data = data_ci) + 
  geom_rect(aes(xmin = 0, xmax = 500, ymin = 0, ymax = 500), 
            fill = NA, col = "black") + 
  geom_point(aes(x = x, y = y, col = ci, size = dbh), pch = 1) + 
  scale_size_continuous(name = "DBH [cm]") + 
  scale_color_viridis_c(name = "Competition index", option = "A") + 
  coord_equal() + 
  theme_void(base_size = 15)

#### Save plots #### 
helpeR::save_ggplot(plot = plot_kernel,
                    filename = "plot_kernel.png",
                    path = "Figures/Appendix",
                    dpi = 300, height = 7.5, width = 15, units = "cm")

helpeR::save_ggplot(plot = plot_pattern,
                    filename = "plot_pattern.png",
                    path = "Figures/Appendix",
                    dpi = 300, height = 15, width = 15, units = "cm")

helpeR::save_ggplot(plot = plot_pattern_ci,
                    filename = "plot_pattern_ci.png",
                    path = "Figures/Appendix",
                    dpi = 300, height = 15, width = 15, units = "cm")

