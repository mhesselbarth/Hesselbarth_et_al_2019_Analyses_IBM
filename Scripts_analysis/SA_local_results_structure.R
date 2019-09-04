#### Import libraries and data ####

# load packages #
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(sensitivity)
library(spatstat)
library(tidyverse)

source("Scripts/helper_functions.R")

#### Import data ####

# import model runs default
sa_default <- readr::read_rds("Data/Output/sa_default.rds")

# import increased parameters
sa_increased_5 <- readr::read_rds("Data/Output/sa_increased_5.rds")
sa_increased_10 <- readr::read_rds("Data/Output/sa_increased_10.rds")

# import decreased parameters
sa_decreased_5 <- readr::read_rds("Data/Output/sa_decreased_5.rds")
sa_decreased_10 <- readr::read_rds("Data/Output/sa_decreased_10.rds")

threshold <- 5
overwrite <- TRUE

#### DBH distribution ####
#### Default parameters ####

sa_default_dbh_dist <- calc_dbh_dist(data = sa_default, threshold = threshold) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel))

#### Calculate changed parameters ####
sa_inc_5_dbh_dist <- calc_dbh_dist(data = sa_increased_5, threshold = threshold) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "5%")

sa_inc_10_dbh_dist <- calc_dbh_dist(data = sa_increased_10, threshold = threshold) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "10%")

sa_dec_5_dbh_dist <- calc_dbh_dist(data = sa_decreased_5, threshold = threshold) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "5%")

sa_dec_10_dbh_dist <- calc_dbh_dist(data = sa_decreased_10, threshold = threshold) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(mean = mean(n_rel), 
                   min = min(n_rel), 
                   max = max(n_rel), 
                   sd = sd(n_rel)) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "10%")

#### Plot results ####
sa_overall_dbh_dist <- dplyr::bind_rows(sa_inc_5_dbh_dist, 
                                        sa_inc_10_dbh_dist, 
                                        sa_dec_5_dbh_dist,
                                        sa_dec_10_dbh_dist) %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength))

sa_ggplot_dbh_dist <- ggplot(data = sa_overall_dbh_dist) + 
  geom_bar(data = sa_default_dbh_dist, aes(x = dbh_class, y = mean), stat = "identity") +
  geom_line(aes(x = dbh_class, y = mean, col = parameter), linetype = 1) +
  facet_wrap(~ combined) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  scale_x_continuous(name = "DBH class [cm]",
                     breaks = seq(from = threshold,
                                  to = 120,
                                  by = 10),
                     labels = paste0(">" ,seq(from = threshold,
                                              to = 120,
                                              by = 10))) +
  scale_y_continuous(name = "Relative frequency") +
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_dbh_dist, 
                     filename = "sa_ggplot_dbh_dist.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
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
