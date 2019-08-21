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

# get observation window
window <- spatstat::owin(xrange = c(0, 500), yrange = c(0, 500))

#### Pair-correlation function ####

#### Calculate default result #### 

sa_default_pcf <- calc_pcf(data = sa_default, correction = "Ripley", 
                           window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>%
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   min = min(iso), 
                   max = max(iso), 
                   n = n())

#### Calculate changed parameters ####

sa_inc_5_pcf <- calc_pcf(data = sa_increased_5, correction = "Ripley", 
                         window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased 5%")

sa_inc_10_pcf <- calc_pcf(data = sa_increased_10, correction = "Ripley", 
                          window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased 10%")

sa_dec_5_pcf <- calc_pcf(data = sa_decreased_5, correction = "Ripley", 
                         window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased 5%")

sa_dec_10_pcf <- calc_pcf(data = sa_decreased_10, correction = "Ripley", 
                          window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased 10%")

sa_overall_pcf <- dplyr::bind_rows(sa_inc_5_pcf, 
                                   sa_inc_10_pcf, 
                                   sa_dec_5_pcf,
                                   sa_dec_10_pcf)

sa_ggplot_pcf <- ggplot(data = sa_overall_pcf) + 
  geom_ribbon(data = sa_default_pcf, 
              aes(x = r, ymin = min, ymax = max), fill = "grey") +
  geom_line(aes(x = r, y = iso, col = parameter)) + 
  facet_wrap(~ direction) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  labs(x = "r [m]", y = "pcf(r)") + 
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_pcf, 
                     filename = "sa_ggplot_pcf.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = TRUE)

#### Nearest neighbor distribution function ####

#### Calculate default result #### 

sa_default_nnd <- calc_nnd(data = sa_default, correction = "km",
                           window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>%
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   min = min(km), 
                   max = max(km), 
                   n = n())

#### Calculate changed parameters ####

sa_inc_5_nnd <- calc_nnd(data = sa_increased_5, correction = "Ripley", 
                         window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased 5%")

sa_inc_10_nnd <- calc_nnd(data = sa_increased_10, correction = "Ripley", 
                          window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased 10%")

sa_dec_5_nnd <- calc_nnd(data = sa_decreased_5, correction = "Ripley", 
                         window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased 5%")

sa_dec_10_nnd <- calc_nnd(data = sa_decreased_10, correction = "Ripley", 
                          window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(km),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased 10%")

sa_overall_nnd <- dplyr::bind_rows(sa_inc_5_nnd, 
                                   sa_inc_10_nnd, 
                                   sa_dec_5_nnd,
                                   sa_dec_10_nnd)

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
                     overwrite = TRUE)

#### Mark-correlation function ####

#### Calculate default result #### 

sa_default_kmm <- calc_kmm(data = sa_default, correction = "Ripley",
                           window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>%
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   min = min(iso), 
                   max = max(iso), 
                   n = n())

#### Calculate changed parameters ####

sa_inc_5_kmm <- calc_kmm(data = sa_increased_5, correction = "Ripley", 
                         window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased 5%")

sa_inc_10_kmm <- calc_kmm(data = sa_increased_10, correction = "Ripley", 
                          window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Increased 10%")

sa_dec_5_kmm <- calc_kmm(data = sa_decreased_5, correction = "Ripley", 
                         window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased 5%")

sa_dec_10_kmm <- calc_kmm(data = sa_decreased_10, correction = "Ripley", 
                          window = window) %>% 
  # dplyr::mutate(id = as.numeric(id)) %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   iso = mean(iso),
                   n = n()) %>%
  dplyr::mutate(direction = "Decreased 10%")

sa_overall_kmm <- dplyr::bind_rows(sa_inc_5_kmm, 
                                   sa_inc_10_kmm, 
                                   sa_dec_5_kmm,
                                   sa_dec_10_kmm)

sa_ggplot_kmm <- ggplot(data = sa_overall_kmm) + 
  geom_ribbon(data = sa_default_kmm, aes(x = r, ymin = min, ymax = max), 
              fill = "grey") +
  geom_line(aes(x = r, y = iso, col = parameter)) +
  facet_wrap(~ direction) +
  scale_colour_viridis_d(name = "Parameter", option = "D") +
  labs(x = "r [m]", y = "kmm(r)") + 
  coord_cartesian(xlim = c(0, 10)) +
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = sa_ggplot_kmm, 
                     filename = "sa_ggplot_kmm.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = TRUE)
