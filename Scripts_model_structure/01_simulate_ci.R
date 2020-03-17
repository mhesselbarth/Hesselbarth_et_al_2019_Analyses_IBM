###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Model structure CI ####

#### Import libraries and data ####

# load packages #
source("Helper_functions/helper_functions_setup.R")

# import parameters
parameters_fitted_biotic <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt",
                                                   sep = ";")

# import data
beech_1999_ppp_rec <- readr::read_rds("Data/Input/beech_1999_ppp_rec.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

#### Pre-processing of input data ####
set.seed(42)
sample_id <- sample(1:beech_1999_ppp_rec$n, size = 250)

input_data <- tibble::as_tibble(beech_1999_ppp_rec) %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  dplyr::filter(species == "beech", id %in% sample_id) %>%
  dplyr::select(-id) %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh")

#### Plot kernel ####
min_max <- 20

# calculate CI index dbh = 20
ci_index_20 <- purrr::map_dbl(0:min_max, function(x) 
  (20 ^ parameters_fitted_biotic$ci_alpha) * exp(-(x / (20 ^ parameters_fitted_biotic$ci_beta))))
# ci_index_20 <- ci_index_20 / (5 ^ parameters$ci_alpha + ci_index_20)

# calculate CI index dbh = 10 
ci_index_10 <- purrr::map_dbl(0:min_max, function(x) 
  (10 ^ parameters_fitted_biotic$ci_alpha) * exp(-(x / (10 ^ parameters_fitted_biotic$ci_beta))))
# ci_index_10 <- ci_index_10 / (5 ^ parameters$ci_alpha + ci_index_10)

# create ggplot
ggplot_kernel <- ggplot() + 
  geom_line(aes(x = -min_max:min_max, 
                y = c(rev(ci_index_10), ci_index_10[-1]),
            linetype = "dbh 10")) + 
  geom_line(aes(x = -min_max:min_max, 
                y = c(rev(ci_index_20), ci_index_20[-1]), 
            linetype = "dbh 20")) + 
  scale_linetype_manual(name = "dbh [cm]", 
                        values = c("dbh 10" = 2, "dbh 20" = 1)) + 
  scale_x_continuous(breaks = seq(-min_max, min_max, 5),
                     labels = abs(seq(-min_max, min_max, 5))) + 
  labs(x = expression(paste(Distance[ij], " [m]")), y = expression(c[i]^raw)) + 
  guides(linetype = FALSE) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        plot.margin = margin(0, 0, 0, 0, "mm"), 
        text = element_text(family = "Calibri Light"))

suppoRt::save_ggplot(plot = ggplot_kernel, 
                     filename = "ggplot_ci_kernel.png", 
                     path = "C:/Users/Maximilian/ownCloud/13_Thesis_defense/Figures/",
                     dpi = dpi, units = units,
                     width = 50, height = 50, 
                     overwrite = FALSE)

#### Plot point pattern ####
# calculat ci
data_ci <- rabmp::simulate_ci(input_data, parameters = parameters_fitted_biotic) %>% 
  tibble::as_tibble()

# create ggplot without ci
# plot_pattern <- ggplot(data = data_ci) + 
#   geom_polygon(data = spatstat::as.data.frame.owin(plot_area), 
#                aes(x = x, y = y),
#                col = "black", fill = NA) +
#   geom_point(aes(x = x, y = y, size = dbh), pch = 1) + 
#   scale_size_continuous(name = "dbh [cm]") + 
#   coord_equal() + 
#   theme_void(base_size = 15) + 
#   theme(legend.position = "bottom")

# create ggplot with ci
ggplot_pattern_ci <- ggplot(data = data_ci) + 
  geom_polygon(data = spatstat::as.data.frame.owin(plot_area),
               aes(x = x, y = y),
               col = "black", fill = "grey", alpha = 0.65) +
  geom_point(aes(x = x, y = y, col = ci, size = dbh), pch = 19) + 
  scale_size_continuous(name = "dbh [cm]") + 
  scale_color_viridis_c(name = expression(c[i]^trans), option = "C") +
  coord_equal() + 
  guides(size = FALSE, col = FALSE) +
  theme_void(base_size = base_size) + 
  theme(legend.position = "bottom", 
        plot.margin = margin(0, 0, 0, 0, "mm"), 
        text = element_text(family = "Calibri Light"))


#### Overall plot ####
ggplot_overall <- ggplot_kernel + ggplot_pattern_ci + 
  patchwork::plot_layout(ncol = 2, nrow = 1,
                         widths = c(0.5, 0.5), heights = c(0.5, 0.5)) + 
  patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")", 
                             theme = theme(text = element_text(family = "Calibri Light")))

#### Save plots #### 
suppoRt::save_ggplot(plot = ggplot_overall, 
                     filename = "ggplot_ci_overall.png", 
                     path = "C:/Users/Maximilian/ownCloud/13_Disputation/Figures/",
                     dpi = dpi, units = units,
                     width = 200, height = 125, 
                     overwrite = FALSE)
