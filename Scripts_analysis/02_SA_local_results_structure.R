#### Import libraries and data ####

# load packages #
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(sensitivity)
library(spatstat)
library(tidyverse)

source("Scripts_analysis/helper_functions_sa_structure.R")

#### Import data ####

# import model runs default
sa_default <- readr::read_rds("Data/Output/sa_default.rds")

# import increased parameters
sa_increased_5 <- readr::read_rds("Data/Output/sa_increased_5.rds")
sa_increased_10 <- readr::read_rds("Data/Output/sa_increased_10.rds")

# import decreased parameters
sa_decreased_5 <- readr::read_rds("Data/Output/sa_decreased_5.rds")
sa_decreased_10 <- readr::read_rds("Data/Output/sa_decreased_10.rds")

overwrite <- TRUE

#### DBH distribution ####
smaller_small <- 1

bigger_large <- 1
bigger_small <- 0

by_large <- 10
by_small <- 0.1

#### Default parameters ####
sa_default_dbh_dist_large <- calc_dbh_dist(data = sa_default, bigger = bigger_large, 
                                           by = by_large) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel))

sa_default_dbh_dist_small <- calc_dbh_dist(data = sa_default,
                                           bigger = bigger_small, smaller = smaller_small,
                                           by = by_small) %>%
  dplyr::group_by(parameter, dbh_class) %>%
  dplyr::summarise(mean = mean(n_rel),
                   min = min(n_rel),
                   max = max(n_rel),
                   sd = sd(n_rel))

#### Calculate changed parameters ####
# Increased #
sa_inc_5_dbh_dist_large <- calc_dbh_dist(data = sa_increased_5, 
                                         bigger = bigger_large,
                                         by = by_large) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "5%")

sa_inc_5_dbh_dist_small <- calc_dbh_dist(data = sa_increased_5,
                                         bigger = bigger_small, smaller = smaller,
                                         by = by_small) %>%
  dplyr::group_by(parameter, dbh_class) %>%
  dplyr::summarise(mean = mean(n_rel),
                   min = min(n_rel),
                   max = max(n_rel),
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Increased",
                strength = "5%")

sa_inc_10_dbh_dist_large <- calc_dbh_dist(data = sa_increased_10, 
                                          bigger = bigger_large, 
                                          by = by_large) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "10%")

sa_inc_10_dbh_dist_small <- calc_dbh_dist(data = sa_increased_10,
                                          bigger = bigger_small, smaller = smaller,
                                          by = by_small) %>%
  dplyr::group_by(parameter, dbh_class) %>%
  dplyr::summarise(mean = mean(n_rel),
                   min = min(n_rel),
                   max = max(n_rel),
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Increased",
                strength = "10%")

# Decreased #
sa_dec_5_dbh_dist_large <- calc_dbh_dist(data = sa_decreased_5, 
                                         bigger = bigger_large, 
                                         by = by_large) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "5%")

sa_dec_5_dbh_dist_small <- calc_dbh_dist(data = sa_decreased_5,
                                         bigger = bigger_small, smaller = smaller,
                                         by = by_small) %>%
  dplyr::group_by(parameter, dbh_class) %>%
  dplyr::summarise(mean = mean(n_rel),
                   min = min(n_rel),
                   max = max(n_rel),
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Decreased",
                strength = "5%")

sa_dec_10_dbh_dist_large <- calc_dbh_dist(data = sa_decreased_10, 
                                          bigger = bigger_large,
                                          by = by_large) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "10%")

sa_dec_10_dbh_dist_small <- calc_dbh_dist(data = sa_decreased_10,
                                          bigger = bigger_small, smaller = smaller,
                                          by = by_small) %>%
  dplyr::group_by(parameter, dbh_class) %>%
  dplyr::summarise(mean = mean(n_rel),
                   min = min(n_rel),
                   max = max(n_rel),
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Decreased",
                strength = "10%")

#### Plot results ####
sa_overall_dbh_dist_large <- dplyr::bind_rows(sa_inc_5_dbh_dist_large, 
                                              sa_inc_10_dbh_dist_large, 
                                              sa_dec_5_dbh_dist_large,
                                              sa_dec_10_dbh_dist_large) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength), 
                parameter = forcats::as_factor(parameter))

min_class_large <- min(sa_overall_dbh_dist_large$dbh_class)
max_class_large <- max(sa_overall_dbh_dist_large$dbh_class)

sa_overall_dbh_dist_small <- dplyr::bind_rows(sa_inc_5_dbh_dist_small,
                                              sa_inc_10_dbh_dist_small,
                                              sa_dec_5_dbh_dist_small,
                                              sa_dec_10_dbh_dist_small) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength),
                parameter = forcats::as_factor(parameter))

min_class_small <- min(sa_overall_dbh_dist_small$dbh_class)
max_class_small <- max(sa_overall_dbh_dist_small$dbh_class)

sa_ggplot_dbh_dist_large <- ggplot(data = sa_overall_dbh_dist_large) + 
  geom_bar(data = sa_default_dbh_dist_large, 
           aes(x = dbh_class, y = mean), stat = "identity") +
  geom_errorbar(data = sa_default_dbh_dist_large, 
                aes(x = dbh_class, ymin = min, ymax = max), width = 0.5) + 
  geom_point(aes(x = dbh_class, y = mean, col = parameter), pch = 19) +
  facet_wrap(~ combined) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  scale_x_continuous(name = "DBH class [cm]",
                     breaks = seq(from = min_class_large,
                                  to = max_class_large,
                                  by = 1),
                     labels = paste0(">", seq(from = min_class_large,
                                              to = max_class_large * by_large,
                                              by = by_large))) +
  scale_y_continuous(name = "Relative frequency") +
  theme_classic(base_size = 15)

sa_ggplot_dbh_dist_small <- ggplot(data = sa_overall_dbh_dist_small) +
  geom_bar(data = sa_default_dbh_dist_small,
           aes(x = dbh_class, y = mean), stat = "identity") +
  geom_errorbar(data = sa_default_dbh_dist_small,
                aes(x = dbh_class, ymin = min, ymax = max), width = 0.5) +
  geom_point(aes(x = dbh_class, y = mean, col = parameter), pch = 19) +
  facet_wrap(~ combined) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  scale_x_continuous(name = "DBH class [cm]",
                     breaks = seq(from = min_class_small,
                                  to = max_class_small,
                                  by = 1),
                     labels = paste0(">", seq(from = bigger_small + by_small,
                                              to = smaller_small - by_small,
                                              by = by_small))) +
  scale_y_continuous(name = "Relative frequency") +
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_dbh_dist_large, 
                     filename = "sa_ggplot_dbh_dist_large.png", path = "Figures/", 
                     dpi = 300, width = 35, height = 15, units = "cm", 
                     overwrite = overwrite)

suppoRt::save_ggplot(plot = sa_ggplot_dbh_dist_small,
                     filename = "sa_ggplot_dbh_dist_small.png", path = "Figures/",
                     dpi = 300, width = 35, height = 15, units = "cm",
                     overwrite = overwrite)

#### DBH growth ####
#### Default parameters ####
sa_default_growth <- calc_growth(data = sa_default) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(mean = mean(dbh_inc), 
                   min = min(dbh_inc), 
                   max = max(dbh_inc), 
                   sd = sd(dbh_inc))

#### Calculate changed parameters ####
sa_inc_5_growth <- calc_growth(data = sa_increased_5) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(dbh_inc) - sa_default_growth$mean, 
                   diff_rel = diff_abs / sa_default_growth$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "5%")

sa_inc_10_growth <- calc_growth(data = sa_increased_10) %>% 
  dplyr::group_by(parameter)  %>% 
  dplyr::summarise(diff_abs = mean(dbh_inc) - sa_default_growth$mean, 
                   diff_rel = diff_abs / sa_default_growth$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "10%")

sa_dec_5_growth <- calc_growth(data = sa_decreased_5) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(dbh_inc) - sa_default_growth$mean, 
                   diff_rel = diff_abs / sa_default_growth$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "5%")

sa_dec_10_growth <- calc_growth(data = sa_decreased_10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(dbh_inc) - sa_default_growth$mean, 
                   diff_rel = diff_abs / sa_default_growth$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "10%")

#### Plot results ####
sa_overall_growth <- dplyr::bind_rows(sa_inc_5_growth, 
                                        sa_inc_10_growth, 
                                        sa_dec_5_growth,
                                        sa_dec_10_growth) %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength))

sa_ggplot_growth <- ggplot(data = sa_overall_growth) + 
  geom_bar(aes(x = parameter, y = diff_rel, fill = strength),
           position = "dodge", col = "black", stat = "identity") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  facet_wrap(~ direction) +   
  coord_flip() +
  scale_fill_viridis_d(name = "Parameter\nchange", option = "D") + 
  scale_y_continuous(breaks = seq(-10, 10, 5)) + 
  labs(x = "Parameter", y = "Relative difference DBH growth/year [%]") + 
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_growth, 
                     filename = "sa_ggplot_growth.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)

#### N died ####
#### Default parameters ####
sa_default_died <- calc_died(data = sa_default) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(mean = mean(n_died), 
                   min = min(n_died), 
                   max = max(n_died), 
                   sd = sd(n_died))

#### Calculate changed parameters ####
sa_inc_5_died <- calc_died(data = sa_increased_5) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(n_died) - sa_default_died$mean, 
                   diff_rel = diff_abs / sa_default_died$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "5%")

sa_inc_10_died <- calc_died(data = sa_increased_10) %>% 
  dplyr::group_by(parameter)  %>% 
  dplyr::summarise(diff_abs = mean(n_died) - sa_default_died$mean, 
                   diff_rel = diff_abs / sa_default_died$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "10%")

sa_dec_5_died <- calc_died(data = sa_decreased_5) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(n_died) - sa_default_died$mean, 
                   diff_rel = diff_abs / sa_default_died$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "5%")

sa_dec_10_died <- calc_died(data = sa_decreased_10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(n_died) - sa_default_died$mean, 
                   diff_rel = diff_abs / sa_default_died$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "10%")

#### Plot results ####
sa_overall_died <- dplyr::bind_rows(sa_inc_5_died, 
                                    sa_inc_10_died, 
                                    sa_dec_5_died,
                                    sa_dec_10_died) %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength))

sa_ggplot_died <- ggplot(data = sa_overall_died) + 
  geom_bar(aes(x = parameter, y = diff_rel, fill = strength),
           position = "dodge", col = "black", stat = "identity") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  facet_wrap(~ direction) +   
  coord_flip() +
  scale_fill_viridis_d(name = "Parameter\nchange", option = "D") + 
  scale_y_continuous(breaks = seq(-30, 30, 10)) +
  labs(x = "Parameter", y = "Relative difference n died [%]") + 
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_died, 
                     filename = "sa_ggplot_died.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)
