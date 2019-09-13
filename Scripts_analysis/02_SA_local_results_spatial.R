#### Import libraries and data ####

# load packages #
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(sensitivity)
library(spatstat)
library(tidyverse)

source("Scripts_analysis/helper_functions_sa_spatial.R")

#### Import data ####
pattern_1999_recon <- readr::read_rds("Data/Input/pattern_1999_reconstructed.rds")

# import model runs default
sa_default <- readr::read_rds("Data/Output/sa_default.rds")

# import increased parameters
sa_increased_5 <- readr::read_rds("Data/Output/sa_increased_5.rds")
sa_increased_10 <- readr::read_rds("Data/Output/sa_increased_10.rds")

# import decreased parameters
sa_decreased_5 <- readr::read_rds("Data/Output/sa_decreased_5.rds")
sa_decreased_10 <- readr::read_rds("Data/Output/sa_decreased_10.rds")

#### Set parameters ####
# get observation window
window <- pattern_1999_recon$window

# rm(pattern_1999_recon)

# calulate r 
r <- seq(from = 0,
         to = spatstat::rmax.rule(W = window,
                                  lambda = nrow(dplyr::filter(sa_default[[1]], 
                                                              i == max(i))) / 
                                    spatstat::area(window)),
         length.out = 515)

# set parameters
overwrite <- FALSE

#### Pair-correlation function ####
#### Calculate default result #### 
sa_default_pcf_large <- calc_pcf(data = sa_default, correction = "Ripley", r = r,
                                 fast = FALSE, divisor = "d",
                                 window = window, bigger = 10) %>% 
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   min = min(iso), 
                   max = max(iso), 
                   n = n())

sa_default_pcf_small <- calc_pcf(data = sa_default, correction = "Ripley", r = r,
                                 fast = FALSE, divisor = "d",
                                 window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   min = min(iso), 
                   max = max(iso), 
                   n = n())

# rm(sa_default)

#### Calculate changed parameters ####
# Increased #
sa_inc_5_pcf_large <- calc_pcf(data = sa_increased_5, correction = "Ripley", r = r,
                            fast = FALSE, divisor = "d",
                            window = window, bigger = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "5%")

sa_inc_5_pcf_small <- calc_pcf(data = sa_increased_5, correction = "Ripley", r = r,
                               fast = FALSE, divisor = "d",
                               window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "5%")

sa_inc_10_pcf_large <- calc_pcf(data = sa_increased_10, correction = "Ripley", r = r,
                             fast = FALSE, divisor = "d",
                             window = window, bigger = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "10%")

sa_inc_10_pcf_small <- calc_pcf(data = sa_increased_10, correction = "Ripley", r = r,
                                fast = FALSE, divisor = "d",
                                window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "10%")

# rm(list = c("sa_increased_5", "sa_increased_10"))

# Decreased #
sa_dec_5_pcf_large <- calc_pcf(data = sa_decreased_5, correction = "Ripley", r = r,
                               fast = FALSE, divisor = "d",
                               window = window, bigger = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "5%")

sa_dec_5_pcf_small <- calc_pcf(data = sa_decreased_5, correction = "Ripley", r = r,
                               fast = FALSE, divisor = "d",
                               window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "5%")

sa_dec_10_pcf_large <- calc_pcf(data = sa_decreased_10, correction = "Ripley", r = r,
                                fast = FALSE, divisor = "d",
                                window = window, bigger = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "10%")

sa_dec_10_pcf_small <- calc_pcf(data = sa_decreased_10, correction = "Ripley", r = r,
                                fast = FALSE, divisor = "d",
                                window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "10%")

# rm(list = c("sa_decreased_5", "sa_decreased_10"))

#### Plot results ####
sa_overall_pcf_large <- dplyr::bind_rows(sa_inc_5_pcf_large, 
                                         sa_inc_10_pcf_large, 
                                         sa_dec_5_pcf_large,
                                         sa_dec_10_pcf_large) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength), 
                parameter = forcats::as_factor(parameter))

sa_overall_pcf_small <- dplyr::bind_rows(sa_inc_5_pcf_small, 
                                         sa_inc_10_pcf_small, 
                                         sa_dec_5_pcf_small,
                                         sa_dec_10_pcf_small) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength), 
                parameter = forcats::as_factor(parameter))

# create plot
sa_ggplot_pcf_large <- ggplot(data = sa_overall_pcf_large) + 
  geom_ribbon(data = sa_default_pcf_large,
              aes(x = r, ymin = min, ymax = max), fill = "grey") +
  geom_line(aes(x = r, y = pcf, col = parameter)) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ combined) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  labs(x = "r [m]", y = "pcf(r)") + 
  theme_classic(base_size = 15)

sa_ggplot_pcf_small <- ggplot(data = sa_overall_pcf_small) + 
  geom_ribbon(data = sa_default_pcf_small,
              aes(x = r, ymin = min, ymax = max), fill = "grey") +
  geom_line(aes(x = r, y = pcf, col = parameter)) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ combined) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  labs(x = "r [m]", y = "pcf(r)") + 
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_pcf_large, 
                     filename = "sa_ggplot_pcf_large.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)

suppoRt::save_ggplot(plot = sa_ggplot_pcf_small, 
                     filename = "sa_ggplot_pcf_small.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)

#### Nearest neighbor distribution function ####
#### Calculate default result #### 
sa_default_nnd_large <- calc_nnd(data = sa_default, correction = "km", r = r, 
                              window = window, bigger = 10) %>% 
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   min = min(km), 
                   max = max(km), 
                   n = n())

sa_default_nnd_small <- calc_nnd(data = sa_default, correction = "km", r = r, 
                                 window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   min = min(km), 
                   max = max(km), 
                   n = n())

#### Calculate changed parameters ####
# Increased # 
sa_inc_5_nnd_large <- calc_nnd(data = sa_increased_5, correction = "km", r = r, 
                               window = window, bigger = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "5%")

sa_inc_5_nnd_small <- calc_nnd(data = sa_increased_5, correction = "km", r = r, 
                               window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "5%")

sa_inc_10_nnd_large <- calc_nnd(data = sa_increased_10, correction = "km", r = r, 
                                window = window, bigger = 5) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "10%")

sa_inc_10_nnd_small <- calc_nnd(data = sa_increased_10, correction = "km", r = r, 
                                window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased", 
                strength = "10%")

# rm(list = c("sa_increased_5", "sa_increased_10"))

# Decreased # 
sa_dec_5_nnd_large <- calc_nnd(data = sa_decreased_5, correction = "km", r = r, 
                               window = window, bigger = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "5%")

sa_dec_5_nnd_small <- calc_nnd(data = sa_decreased_5, correction = "km", r = r, 
                               window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "5%")

sa_dec_10_nnd_large <- calc_nnd(data = sa_decreased_10, correction = "km", r = r, 
                             window = window, bigger = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "10%")

sa_dec_10_nnd_small <- calc_nnd(data = sa_decreased_10, correction = "km", r = r, 
                                window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased", 
                strength = "10%")

# rm(list = c("sa_decreased_5", "sa_decreased_10"))

#### Plot results ####
sa_overall_nnd_large <- dplyr::bind_rows(sa_inc_5_nnd_large, 
                                         sa_inc_10_nnd_large, 
                                         sa_dec_5_nnd_large,
                                         sa_dec_10_nnd_large) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength), 
                parameter = forcats::as_factor(parameter))

sa_overall_nnd_small <- dplyr::bind_rows(sa_inc_5_nnd_small, 
                                         sa_inc_10_nnd_small, 
                                         sa_dec_5_nnd_small,
                                         sa_dec_10_nnd_small) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength), 
                parameter = forcats::as_factor(parameter))

sa_ggplot_nnd_large <- ggplot(data = sa_overall_nnd_large) + 
  geom_ribbon(data = sa_default_nnd_large, aes(x = r, ymin = min, ymax = max), 
              fill = "grey") +
  geom_line(aes(x = r, y = iso, col = parameter)) +
  facet_wrap(~ combined) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  labs(x = "r [m]", y = "nnd(r)") + 
  coord_cartesian(xlim = c(0, 10)) +
  theme_classic(base_size = 15)

sa_ggplot_nnd_small <- ggplot(data = sa_overall_nnd_small) + 
  geom_ribbon(data = sa_default_nnd_small, aes(x = r, ymin = min, ymax = max), 
              fill = "grey") +
  geom_line(aes(x = r, y = iso, col = parameter)) +
  facet_wrap(~ combined) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  labs(x = "r [m]", y = "nnd(r)") + 
  coord_cartesian(xlim = c(0, 10)) +
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_nnd_large, 
                     filename = "sa_ggplot_nnd_large.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)

suppoRt::save_ggplot(plot = sa_ggplot_nnd_small, 
                     filename = "sa_ggplot_nnd_small.png", path = "Figures/", 
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
sa_default_clark_large <- calc_clark(data = sa_default, correction = "cdf",
                                     window = window, bigger = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(mean = mean(clark), 
                   min = min(clark), 
                   max = max(clark),
                   n = n())

sa_default_clark_small <- calc_clark(data = sa_default, correction = "cdf",
                                     window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(mean = mean(clark), 
                   min = min(clark), 
                   max = max(clark),
                   n = n())

# rm(sa_default)

#### Calculate changed parameters ####
# Increased #
sa_inc_5_clark_large <- calc_clark(data = sa_increased_5, correction = "cdf", 
                                   window = window, bigger = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark_large$mean, 
                   diff_rel = diff_abs / sa_default_clark_large$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "5%")

sa_inc_5_clark_small <- calc_clark(data = sa_increased_5, correction = "cdf", 
                                   window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark_small$mean, 
                   diff_rel = diff_abs / sa_default_clark_small$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "5%")

sa_inc_10_clark_large <- calc_clark(data = sa_increased_10, correction = "cdf", 
                                    window = window, bigger = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark_large$mean, 
                   diff_rel = diff_abs / sa_default_clark_large$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "10%")

sa_inc_10_clark_small <- calc_clark(data = sa_increased_10, correction = "cdf", 
                                    window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark_small$mean, 
                   diff_rel = diff_abs / sa_default_clark_small$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Increased", 
                strength = "10%")

# rm(list = c("sa_increased_5", "sa_increased_10"))

# Decreased #
sa_dec_5_clark_large <- calc_clark(data = sa_decreased_5, correction = "cdf", 
                                   window = window, bigger = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark_large$mean, 
                   diff_rel = diff_abs / sa_default_clark_large$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "5%")

sa_dec_5_clark_small <- calc_clark(data = sa_decreased_5, correction = "cdf", 
                                   window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark_small$mean, 
                   diff_rel = diff_abs / sa_default_clark_small$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "5%")

sa_dec_10_clark_large <- calc_clark(data = sa_decreased_10, correction = "cdf", 
                                    window = window, bigger = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark_large$mean, 
                   diff_rel = diff_abs / sa_default_clark_large$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "10%")

sa_dec_10_clark_small <- calc_clark(data = sa_decreased_10, correction = "cdf", 
                                    window = window, bigger = 1, smaller = 10) %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_abs = mean(clark) - sa_default_clark_small$mean, 
                   diff_rel = diff_abs / sa_default_clark_small$mean * 100, 
                   n = n()) %>% 
  dplyr::mutate(diff_scl = diff_rel / max(abs(diff_rel)), 
                direction = "Decreased", 
                strength = "10%")

# rm(list = c("sa_decreased_5", "sa_decreased_10"))

#### Plot results ####
sa_overall_clark_large <- dplyr::bind_rows(sa_inc_5_clark_large, 
                                           sa_inc_10_clark_large, 
                                           sa_dec_5_clark_large,
                                           sa_dec_10_clark_large) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength), 
                parameter = forcats::as_factor(parameter))

sa_overall_clark_small <- dplyr::bind_rows(sa_inc_5_clark_small, 
                                           sa_inc_10_clark_small, 
                                           sa_dec_5_clark_small,
                                           sa_dec_10_clark_small) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(combined = forcats::as_factor(paste(direction, strength)),
                direction = forcats::as_factor(direction),
                strength = forcats::as_factor(strength), 
                parameter = forcats::as_factor(parameter))

sa_ggplot_clark_large <- ggplot(data = sa_overall_clark_large) + 
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

sa_ggplot_clark_small <- ggplot(data = sa_overall_clark_small) + 
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

suppoRt::save_ggplot(plot = sa_ggplot_clark_large, 
                     filename = "sa_ggplot_clark_large.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)

suppoRt::save_ggplot(plot = sa_ggplot_clark_small, 
                     filename = "sa_ggplot_clark_small.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)
