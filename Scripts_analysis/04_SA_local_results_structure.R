###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Results local SA structure ####

#### Import libraries and data ####

# load packages #
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(spatstat)
library(tidyverse)

source("Helper_functions/helper_functions_sa_structure.R")

#### Import data ####

# import model runs default
sa_default <- readr::read_rds("Data/Output/sa_default_y50_e50_r50.rds")

# import increased parameters
sa_increased_5 <- readr::read_rds("Data/Output/sa_increased_5_y50_e50_r50.rds")
sa_increased_10 <- readr::read_rds("Data/Output/sa_increased_10_y50_e50_r50.rds")

# import decreased parameters
sa_decreased_5 <- readr::read_rds("Data/Output/sa_decreased_5_y50_e50_r50.rds")
sa_decreased_10 <- readr::read_rds("Data/Output/sa_decreased_10_y50_e50_r50.rds")

overwrite <- FALSE

base_size <- 10

#### DBH distribution ####
# set parameters #
by <- 10

# Increased parameters #
sa_dbh_dist_inc_5 <- calc_dbh_dist_sa(default = sa_default, 
                                      changed = sa_increased_5, 
                                      by = by) %>% 
  dplyr::mutate(dbh_class = factor(dbh_class, ordered = TRUE), 
                parameter = factor(parameter), 
                direction = "Increased 5%")

sa_dbh_dist_inc_10 <- calc_dbh_dist_sa(default = sa_default, 
                                       changed = sa_increased_10, 
                                       by = by) %>% 
  dplyr::mutate(dbh_class = factor(dbh_class, ordered = TRUE), 
                parameter = factor(parameter), 
                direction = "Increased 10%")

sa_dbh_dist_inc <- dplyr::bind_rows(sa_dbh_dist_inc_5, 
                                    sa_dbh_dist_inc_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_dbh_dist_inc <- ggplot(data = sa_dbh_dist_inc) + 
  geom_bar(aes(x = dbh_class, y = diff_n_rel * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  facet_wrap(~ parameter, scales = "free_y") + 
  scale_x_discrete(name = "DBH class [cm]",
                   breaks = seq(from = as.numeric(min(sa_dbh_dist_inc$dbh_class)),
                                to = as.numeric(max(sa_dbh_dist_inc$dbh_class)),
                                by = 1),
                   labels = paste0("<", seq(from = as.numeric(min(sa_dbh_dist_inc$dbh_class)) * 10,
                                            to = as.numeric(max(sa_dbh_dist_inc$dbh_class)) * 10,
                                            by = by))) +
  scale_y_continuous(name = "Relative difference [%]") +
  scale_color_viridis_d(name = "Parameter change") +
  scale_fill_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_dbh_dist_inc, 
                     filename = "ggplot_sa_dbh_dist_inc.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300, 
                     overwrite = overwrite)

# Decrease parameters #
sa_dbh_dist_dec_5 <- calc_dbh_dist_sa(default = sa_default, 
                                      changed = sa_decreased_5, 
                                      by = by) %>% 
  dplyr::mutate(dbh_class = factor(dbh_class, ordered = TRUE), 
                parameter = factor(parameter), 
                direction = "Decreased 5%")

sa_dbh_dist_dec_10 <- calc_dbh_dist_sa(default = sa_default, 
                                       changed = sa_decreased_10, 
                                       by = by) %>% 
  dplyr::mutate(dbh_class = factor(dbh_class, ordered = TRUE), 
                parameter = factor(parameter), 
                direction = "Decreased 10%")

sa_dbh_dist_dec <- dplyr::bind_rows(sa_dbh_dist_dec_5, 
                                    sa_dbh_dist_dec_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_dbh_dist_dec <- ggplot(data = sa_dbh_dist_dec) + 
  geom_bar(aes(x = dbh_class, y = diff_n_rel * 100, 
               fill = direction), col = "black", 
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  facet_wrap(~ parameter, scales = "free_y") + 
  scale_x_discrete(name = "DBH class [cm]",
                   breaks = seq(from = as.numeric(min(sa_dbh_dist_dec$dbh_class)),
                                to = as.numeric(max(sa_dbh_dist_dec$dbh_class)),
                                by = 1), 
                   labels = paste0("<", seq(from = as.numeric(min(sa_dbh_dist_dec$dbh_class)) * 10,
                                            to = as.numeric(max(sa_dbh_dist_dec$dbh_class)) * 10,
                                            by = by))) +
  scale_y_continuous(name = "Relative difference [%]") +
  scale_color_viridis_d(name = "Parameter change") +
  scale_fill_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_dbh_dist_dec, 
                     filename = "ggplot_sa_dbh_dist_dec.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300, 
                     overwrite = overwrite)

#### DBH growth ####
# Increased parameters #
sa_growth_inc_5 <- calc_growth_sa(default = sa_default, 
                                  changed = sa_increased_5) %>%
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 5%")

sa_growth_inc_10 <- calc_growth_sa(default = sa_default, 
                                   changed = sa_increased_10) %>%
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 10%")

sa_growth_inc <- dplyr::bind_rows(sa_growth_inc_5, 
                                  sa_growth_inc_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_growth_inc <- ggplot(data = sa_growth_inc) + 
  geom_bar(aes(x = parameter, y = diff_inc * 100, 
               fill = direction), col = "black", 
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter" ,) +
  scale_y_continuous(name = "Relative difference DBH increment [%]", 
                     limits = c(-25, 25)) +
  scale_color_viridis_d(name = "Parameter change") +
  scale_fill_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_growth_inc, 
                     filename = "ggplot_sa_growth_inc.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300, 
                     overwrite = overwrite)

# Decreased parameters #
sa_growth_dec_5 <- calc_growth_sa(default = sa_default, 
                                  changed = sa_decreased_5) %>%
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 5%")

sa_growth_dec_10 <- calc_growth_sa(default = sa_default, 
                                   changed = sa_decreased_10) %>%
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 10%")

sa_growth_dec <- dplyr::bind_rows(sa_growth_dec_5, 
                                  sa_growth_dec_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_growth_dec <- ggplot(data = sa_growth_dec) + 
  geom_bar(aes(x = parameter, y = diff_inc * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter" ,) +
  scale_y_continuous(name = "Relative difference DBH increment [%]", 
                     limits = c(-30, 30)) +
  scale_color_viridis_d(name = "Parameter change") +
  scale_fill_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_growth_dec, 
                     filename = "ggplot_sa_growth_dec.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300, 
                     overwrite = overwrite)

#### n died ####
# Increased parameters #
sa_died_inc_5 <- calc_died_sa(default = sa_default,
                              changed = sa_increased_5) %>%
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 5%")

sa_died_inc_10 <- calc_died_sa(default = sa_default,
                               changed = sa_increased_10) %>%
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 10%")

sa_died_inc <- dplyr::bind_rows(sa_died_inc_5, 
                                sa_died_inc_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_died_inc <- ggplot(data = sa_died_inc) + 
  geom_bar(aes(x = parameter, y = diff_n_rel * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter" ,) +
  scale_y_continuous(name = "Relative difference n died [%]",
                     limits = c(-40, 40)) +
  scale_color_viridis_d(name = "Parameter change") +
  scale_fill_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_died_inc, 
                     filename = "ggplot_sa_died_inc.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300, 
                     overwrite = overwrite)

# Decreased parameters #
sa_died_dec_5 <- calc_died_sa(default = sa_default,
                              changed = sa_decreased_5) %>%
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 5%")

sa_died_dec_10 <- calc_died_sa(default = sa_default,
                               changed = sa_decreased_10) %>%
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 10%")

sa_died_dec <- dplyr::bind_rows(sa_died_dec_5, 
                                sa_died_dec_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_died_dec <- ggplot(data = sa_died_dec) + 
  geom_bar(aes(x = parameter, y = diff_n_rel * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter" ,) +
  scale_y_continuous(name = "Relative difference n died [%]",
                     limits = c(-20, 60)) +
  scale_color_viridis_d(name = "Parameter change") +
  scale_fill_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_died_dec, 
                     filename = "ggplot_sa_died_dec.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300, 
                     overwrite = overwrite)
