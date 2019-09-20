###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Results local SA spatial ####

#### Import libraries and data ####

# load packages #
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(sensitivity)
library(spatstat)
library(tidyverse)

source("Helper_functions/helper_functions_sa_spatial.R")

#### Import data ####
beech_1999_rec <- readr::read_rds("Data/Input/beech_1999_rec.rds")

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
window <- readr::read_rds("Data/Raw/plot_area_owin.rds")

# rm(pattern_1999_recon)

# set parameters
overwrite <- FALSE
base_size <- 15

#### Pair-correlation function ####

# set parameters #
correction <- "good"
divisor <- "d" 
fast <- FALSE
r <- seq(from = 0, to = 50, length.out = 525)

# r <- seq(from = 0,
#          to = spatstat::rmax.rule(W = window,
#                                   lambda = nrow(dplyr::filter(sa_default[[1]], 
#                                                               i == max(i))) / 
#                                     spatstat::area(window)),
#          length.out = 525)

# increased parameters #
sa_pcf_increased_5 <- calc_pcf(default = sa_default,
                               changed = sa_increased_5,
                               correction = correction,
                               window = window, r = r, 
                               fast = fast) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>%
  dplyr::summarise(theo.default = mean(theo.default, na.rm = TRUE),
                   theo.changed = mean(theo.changed, na.rm = TRUE),
                   pcf_diff = mean(pcf_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(parameter = factor(parameter),
                direction = "Increased 5%")

sa_pcf_increased_10 <- calc_pcf(default = sa_default,
                                changed = sa_increased_10,
                                correction = correction,
                                window = window, r = r, 
                                fast = fast) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>%
  dplyr::summarise(theo.default = mean(theo.default, na.rm = TRUE),
                   theo.changed = mean(theo.changed, na.rm = TRUE),
                   pcf_diff = mean(pcf_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(parameter = factor(parameter),
                direction = "Increased 10%")

sa_pcf_increased <- dplyr::bind_rows(sa_pcf_increased_5, sa_pcf_increased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_pcf_inc <- ggplot(data = sa_pcf_increased) + 
  geom_line(aes(x = r, y = pcf_diff * 100, col = parameter)) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~ direction) +
  scale_x_continuous(name = "r [m]" ,) +
  scale_y_continuous(name = "Relative difference pcf(r) [%]") +
  scale_color_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_pcf_inc, 
                     filename = "ggplot_sa_pcf_inc.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300)

# decreased parameters #
sa_pcf_decreased_5 <- calc_pcf(default = sa_default,
                               changed = sa_decreased_5,
                               correction = correction,
                               window = window, r = r, 
                               fast = fast) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>%
  dplyr::summarise(theo.default = mean(theo.default, na.rm = TRUE),
                   theo.changed = mean(theo.changed, na.rm = TRUE),
                   pcf_diff = mean(pcf_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(parameter = factor(parameter),
                direction = "Decreased 5%")

sa_pcf_decreased_10 <- calc_pcf(default = sa_default,
                                changed = sa_decreased_10,
                                correction = correction,
                                window = window, r = r, 
                                fast = fast) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>%
  dplyr::summarise(theo.default = mean(theo.default, na.rm = TRUE),
                   theo.changed = mean(theo.changed, na.rm = TRUE),
                   pcf_diff = mean(pcf_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(parameter = factor(parameter),
                direction = "Decreased 10%")

sa_pcf_decreased <- dplyr::bind_rows(sa_pcf_decreased_5, sa_pcf_decreased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_pcf_dec <- ggplot(data = sa_pcf_decreased) + 
  geom_line(aes(x = r, y = pcf_diff * 100, col = parameter)) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~ direction) +
  scale_x_continuous(name = "r [m]" ,) +
  scale_y_continuous(name = "Relative difference pcf(r) [%]") +
  scale_color_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_pcf_dec, 
                     filename = "ggplot_sa_pcf_dec.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300)

#### Nearest neighbor distribution function ####

# set parameters #
correction <- "km"
r <- seq(from = 0, to = 10, length.out = 525)

# increased parameters #
sa_nnd_increased_5 <- calc_nnd(default = sa_default,
                               changed = sa_increased_5,
                               correction = correction,
                               window = window, r = r) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo.default = mean(theo.default, na.rm = TRUE), 
                   theo.changed = mean(theo.changed, na.rm = TRUE), 
                   nnd_diff = mean(nnd_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 5%")

sa_nnd_increased_10 <- calc_nnd(default = sa_default,
                                changed = sa_increased_10,
                                correction = correction,
                                window = window, r = r) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo.default = mean(theo.default, na.rm = TRUE), 
                   theo.changed = mean(theo.changed, na.rm = TRUE), 
                   nnd_diff = mean(nnd_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 10%")

sa_nnd_increased <- dplyr::bind_rows(sa_nnd_increased_5, sa_nnd_increased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_nnd_inc <- ggplot(data = sa_nnd_increased) + 
  geom_line(aes(x = r, y = nnd_diff * 100, col = parameter)) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~ direction) +
  scale_x_continuous(name = "r [m]" ,) +
  scale_y_continuous(name = "Relative difference G(r) [%]") +
  scale_color_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_nnd_inc, 
                     filename = "ggplot_sa_nnd_inc.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300)

# decreased values #
sa_nnd_decreased_5 <- calc_nnd(default = sa_default,
                               changed = sa_decreased_5,
                               correction = correction,
                               window = window, r = r) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo.default = mean(theo.default, na.rm = TRUE), 
                   theo.changed = mean(theo.changed, na.rm = TRUE), 
                   nnd_diff = mean(nnd_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 5%")

sa_nnd_decreased_10 <- calc_nnd(default = sa_default,
                                changed = sa_decreased_10,
                                correction = correction,
                                window = window, r = r) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(theo.default = mean(theo.default, na.rm = TRUE), 
                   theo.changed = mean(theo.changed, na.rm = TRUE), 
                   nnd_diff = mean(nnd_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 10%")

sa_nnd_decreased <- dplyr::bind_rows(sa_nnd_decreased_5, sa_nnd_decreased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_nnd_dec <- ggplot(data = sa_nnd_decreased) + 
  geom_line(aes(x = r, y = nnd_diff * 100, col = parameter)) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~ direction) +
  scale_x_continuous(name = "r [m]" ,) +
  scale_y_continuous(name = "Relative difference G(r) [%]") +
  scale_color_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_nnd_dec, 
                     filename = "ggplot_sa_nnd_dec.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300)

### Mark-correlation function ####

# set parameters #
correction <- "Ripley"
r <- seq(from = 0, to = 50, length.out = 525)

# r <- seq(from = 0,
#          to = spatstat::rmax.rule(W = window,
#                                   lambda = nrow(dplyr::filter(sa_default[[1]], 
#                                                               i == max(i))) / 
#                                     spatstat::area(window)),
#          length.out = 525)

# increased values #
sa_kmm_increased_5 <- calc_kmm(default = sa_default,
                               changed = sa_increased_5,
                               window = window, r = r, 
                               correction = correction) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(kmm_diff = mean(kmm_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 5%")

sa_kmm_increased_10 <- calc_kmm(default = sa_default,
                                changed = sa_increased_10,
                                window = window, r = r, 
                                correction = correction) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(kmm_diff = mean(kmm_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 10%")

sa_kmm_increased <- dplyr::bind_rows(sa_kmm_increased_5, sa_kmm_increased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_kmm_inc <- ggplot(data = sa_kmm_increased) + 
  geom_line(aes(x = r, y = kmm_diff * 100, col = parameter)) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~ direction) +
  scale_x_continuous(name = "r [m]" ,) +
  scale_y_continuous(name = "Relative difference kmm(r) [%]") +
  scale_color_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_kmm_inc, 
                     filename = "ggplot_sa_kmm_inc.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300)

# decreased values #
sa_kmm_decreased_5 <- calc_kmm(default = sa_default,
                               changed = sa_decreased_5,
                               window = window, r = r, 
                               correction = correction) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(kmm_diff = mean(kmm_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 5%")

sa_kmm_decreased_10 <- calc_kmm(default = sa_default,
                                changed = sa_decreased_10,
                                window = window, r = r, 
                                correction = correction) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter, r) %>% 
  dplyr::summarise(kmm_diff = mean(kmm_diff, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 10%")

sa_kmm_decreased <- dplyr::bind_rows(sa_kmm_decreased_5, sa_kmm_decreased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_kmm_dec <- ggplot(data = sa_kmm_decreased) + 
  geom_line(aes(x = r, y = kmm_diff * 100, col = parameter)) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~ direction) +
  scale_x_continuous(name = "r [m]" ,) +
  scale_y_continuous(name = "Relative difference kmm(r) [%]") +
  scale_color_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_kmm_dec, 
                     filename = "ggplot_sa_kmm_dec.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300)

#### Clark and Evans Index ####

# increased parameters #
sa_clark_increased_5 <- calc_clark(default = sa_default,
                                   changed = sa_increased_5,
                                   window = window) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_ce = mean(diff_ce, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 5%")

sa_clark_increased_10 <- calc_clark(default = sa_default,
                                    changed = sa_increased_10,
                                    window = window) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_ce = mean(diff_ce, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Increased 10%")

sa_clark_increased <- dplyr::bind_rows(sa_clark_increased_5, sa_clark_increased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_clark_inc <-  ggplot(data = sa_clark_increased) + 
  geom_bar(aes(x = parameter, y = diff_ce * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter" ,) +
  scale_y_continuous(name = "Relative difference CE index [%]", 
                     limits = c(-12.5, 12.5)) +
  scale_color_viridis_d(name = "Parameter change") +
  scale_fill_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_clark_inc, 
                     filename = "ggplot_sa_clark_inc.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300)

# decreased values #
sa_clark_decreased_5 <- calc_clark(default = sa_default,
                                   changed = sa_decreased_5,
                                   window = window) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_ce = mean(diff_ce, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 5%")

sa_clark_decreased_10 <- calc_clark(default = sa_default,
                                    changed = sa_decreased_10,
                                    window = window) %>% 
  dplyr::bind_rows(.id = "parameter") %>% 
  dplyr::group_by(parameter) %>% 
  dplyr::summarise(diff_ce = mean(diff_ce, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = factor(parameter), 
                direction = "Decreased 10%")

sa_clark_decreased <- dplyr::bind_rows(sa_clark_decreased_5, sa_clark_decreased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_clark_dec <-  ggplot(data = sa_clark_decreased) + 
  geom_bar(aes(x = parameter, y = diff_ce * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#FDE725FF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter" ,) +
  scale_y_continuous(name = "Relative difference CE index [%]", 
                     limits = c(-12.5, 12.5)) +
  scale_color_viridis_d(name = "Parameter change") +
  scale_fill_viridis_d(name = "Parameter change") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_clark_dec, 
                     filename = "ggplot_sa_clark_dec.png", 
                     path = "Figures/", 
                     width = 29.7, height = 21.0, units = "cm", dpi = 300)
