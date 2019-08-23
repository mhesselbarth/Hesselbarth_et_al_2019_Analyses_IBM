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
sa_default <- readr::read_rds("Data/Output/Example/sa_default.rds")

# import increased parameters
sa_increased_5 <- readr::read_rds("Data/Output/Example/sa_increased_5.rds")
sa_increased_10 <- readr::read_rds("Data/Output/Example/sa_increased_10.rds")

# import decreased parameters
sa_decreased_5 <- readr::read_rds("Data/Output/Example/sa_decreased_5.rds")
sa_decreased_10 <- readr::read_rds("Data/Output/Example/sa_decreased_10.rds")

# get observation window
window <- spatstat::owin(xrange = c(0, 500), yrange = c(0, 500))

# calulate r 
r <- seq(from = 0,
         to = spatstat::rmax.rule(W = window,
                                  lambda = nrow(dplyr::filter(sa_default[[1]], 
                                                              i == max(i))) / 
                                    spatstat::area(window)),
         length.out = 515)

# set parameters
overwrite <- TRUE

#### Pair-correlation function ####
#### Calculate default result #### 
sa_default_pcf <- calc_pcf(data = sa_default, correction = "Ripley", r = r,
                           window = window) %>% 
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   min = min(pcf), 
                   max = max(pcf), 
                   n = n())

#### Calculate changed parameters ####
sa_inc_5_pcf <- calc_pcf(data = sa_increased_5, correction = "Ripley", r = r,
                         window = window) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(pcf),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "5%")

sa_inc_10_pcf <- calc_pcf(data = sa_increased_10, correction = "Ripley", r = r,
                          window = window) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(pcf),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "10%")

sa_dec_5_pcf <- calc_pcf(data = sa_decreased_5, correction = "Ripley", r = r,
                         window = window) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(pcf),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "5%")

sa_dec_10_pcf <- calc_pcf(data = sa_decreased_10, correction = "Ripley", r = r,
                          window = window) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(pcf),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "10%")

#### Plot results ####
sa_overall_pcf <- dplyr::bind_rows(sa_inc_5_pcf, 
                                   sa_inc_10_pcf, 
                                   sa_dec_5_pcf,
                                   sa_dec_10_pcf) %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength))

# create plot
sa_ggplot_pcf <- ggplot(data = sa_overall_pcf) + 
  geom_ribbon(data = sa_default_pcf,
              aes(x = r, ymin = min, ymax = max), fill = "grey") +
  geom_line(aes(x = r, y = pcf, col = parameter)) +
  facet_wrap(~ combined) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  labs(x = "r [m]", y = "pcf(r)") + 
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_pcf, 
                     filename = "sa_ggplot_pcf.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)

#### Nearest neighbor distribution function ####
#### Calculate default result #### 
sa_default_nnd <- calc_nnd(data = sa_default, correction = "km", r = r, 
                           window = window) %>% 
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   min = min(km), 
                   max = max(km), 
                   n = n())

#### Calculate changed parameters ####
sa_inc_5_nnd <- calc_nnd(data = sa_increased_5, correction = "km", r = r, 
                         window = window) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "5%")

sa_inc_10_nnd <- calc_nnd(data = sa_increased_10, correction = "km", r = r, 
                          window = window) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "10%")

sa_dec_5_nnd <- calc_nnd(data = sa_decreased_5, correction = "km", r = r, 
                         window = window) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "5%")

sa_dec_10_nnd <- calc_nnd(data = sa_decreased_10, correction = "km", r = r, 
                          window = window) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "10%")

#### Plot results ####
sa_overall_nnd <- dplyr::bind_rows(sa_inc_5_nnd, 
                                   sa_inc_10_nnd, 
                                   sa_dec_5_nnd,
                                   sa_dec_10_nnd) %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength))

sa_ggplot_nnd <- ggplot(data = sa_overall_nnd) + 
  geom_ribbon(data = sa_default_nnd, aes(x = r, ymin = min, ymax = max), 
              fill = "grey") +
  geom_line(aes(x = r, y = iso, col = parameter)) +
  facet_wrap(~ direction) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  labs(x = "r [m]", y = "nnd(r)") + 
  coord_cartesian(xlim = c(0, 10)) +
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_nnd, 
                     filename = "sa_ggplot_nnd.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)

# #### Mark-correlation function ####
# 
# # STH WRONG HERE! # 
# 
# #### Calculate default result #### 
# 
# sa_default_kmm <- calc_kmm(data = sa_default, correction = "Ripley",
#                            window = window) %>% 
#   # dplyr::mutate(id = as.numeric(id)) %>%
#   dplyr::group_by(r) %>% 
#   dplyr::summarise(theo = mean(theo), 
#                    min = min(iso), 
#                    max = max(iso), 
#                    n = n())
# 
# #### Calculate changed parameters ####
# 
# sa_inc_5_kmm <- calc_kmm(data = sa_increased_5, correction = "Ripley", 
#                          window = window) %>% 
#   # dplyr::mutate(id = as.numeric(id)) %>% 
#   dplyr::group_by(parameter, r) %>% 
#   dplyr::summarise(theo = mean(theo), 
#                    iso = mean(iso),
#                    n = n()) %>%
#   dplyr::mutate(direction = "Increased 5%")
# 
# sa_inc_10_kmm <- calc_kmm(data = sa_increased_10, correction = "Ripley", 
#                           window = window) %>% 
#   # dplyr::mutate(id = as.numeric(id)) %>% 
#   dplyr::group_by(parameter, r) %>% 
#   dplyr::summarise(theo = mean(theo), 
#                    iso = mean(iso),
#                    n = n()) %>%
#   dplyr::mutate(direction = "Increased 10%")
# 
# sa_dec_5_kmm <- calc_kmm(data = sa_decreased_5, correction = "Ripley", 
#                          window = window) %>% 
#   # dplyr::mutate(id = as.numeric(id)) %>% 
#   dplyr::group_by(parameter, r) %>% 
#   dplyr::summarise(theo = mean(theo), 
#                    iso = mean(iso),
#                    n = n()) %>%
#   dplyr::mutate(direction = "Decreased 5%")
# 
# sa_dec_10_kmm <- calc_kmm(data = sa_decreased_10, correction = "Ripley", 
#                           window = window) %>% 
#   # dplyr::mutate(id = as.numeric(id)) %>% 
#   dplyr::group_by(parameter, r) %>% 
#   dplyr::summarise(theo = mean(theo), 
#                    iso = mean(iso),
#                    n = n()) %>%
#   dplyr::mutate(direction = "Decreased 10%")
# 
# sa_overall_kmm <- dplyr::bind_rows(sa_inc_5_kmm,
#                                    sa_inc_10_kmm,
#                                    sa_dec_5_kmm,
#                                    sa_dec_10_kmm) %>%
#   dplyr::mutate(direction = forcats::as_factor(direction))
# 
# sa_ggplot_kmm <- ggplot(data = sa_overall_kmm) + 
#   geom_ribbon(data = sa_default_kmm, aes(x = r, ymin = min, ymax = max), 
#               fill = "grey") +
#   # geom_line(aes(x = r, y = iso, col = parameter)) +
#   facet_wrap(~ direction) +
#   scale_colour_viridis_d(name = "Parameter", option = "D") +
#   labs(x = "r [m]", y = "kmm(r)") + 
#   coord_cartesian(xlim = c(0, 10)) +
#   theme_classic(base_size = 15)
# 
# suppoRt::save_ggplot(plot = sa_ggplot_kmm, 
#                      filename = "sa_ggplot_kmm.png", path = "Figures/", 
#                      dpi = 300, width = 30, height = 15, units = "cm", 
#                      overwrite = overwrite)


#### Clark and Evans Index ####
#### Calculate default result #### 
sa_default_clark <- calc_clark(data = sa_default, correction = "cdf",
                               window = window) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(mean = mean(clark), 
                   min = min(clark), 
                   max = max(clark),
                   n = n())

#### Calculate changed parameters ####
sa_inc_5_clark <- calc_clark(data = sa_increased_5, correction = "cdf", 
                             window = window) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark$mean, 
                   diff_rel = diff_abs / sa_default_clark$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "5%")

sa_inc_10_clark <- calc_clark(data = sa_increased_10, correction = "cdf", 
                          window = window) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark$mean, 
                   diff_rel = diff_abs / sa_default_clark$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "10%")

sa_dec_5_clark <- calc_clark(data = sa_decreased_5, correction = "cdf", 
                         window = window) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark$mean, 
                   diff_rel = diff_abs / sa_default_clark$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "5%")

sa_dec_10_clark <- calc_clark(data = sa_decreased_10, correction = "cdf", 
                              window = window) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark$mean, 
                   diff_rel = diff_abs / sa_default_clark$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "10%")

#### Plot results ####
sa_overall_clark <- dplyr::bind_rows(sa_inc_5_clark, 
                                     sa_inc_10_clark, 
                                     sa_dec_5_clark,
                                     sa_dec_10_clark) %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength))

sa_ggplot_clark <- ggplot(data = sa_overall_clark) + 
  geom_bar(aes(x = parameter, y = diff_rel, fill = strength), 
           stat = "identity", width = 0.75, position = "dodge",
           col = "black") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  facet_wrap(~ direction) +   
  coord_flip() +
  scale_y_continuous(breaks = seq(-10, 10, 5)) + 
  scale_fill_viridis_d(name = "Parameter\nchange", option = "D") + 
  labs(x = "Parameter", y = "Relative difference CE Index [%]") + 
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_clark, 
                     filename = "sa_ggplot_clark.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)
