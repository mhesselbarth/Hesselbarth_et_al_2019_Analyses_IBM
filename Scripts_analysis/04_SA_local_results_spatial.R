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
library(spatstat)
library(tidyverse)

source("Helper_functions/helper_functions_sa_spatial.R")

# import data #
beech_1999_rec <- readr::read_rds("Data/Input/beech_1999_rec_ppp.rds")

# import model runs default
sa_default <- readr::read_rds("Data/Output/sa_default_y50_e50_r50.rds")

# import increased parameters
sa_increased_5 <- readr::read_rds("Data/Output/sa_increased_5_y50_e50_r50.rds")
sa_increased_10 <- readr::read_rds("Data/Output/sa_increased_10_y50_e50_r50.rds")

# import decreased parameters
sa_decreased_5 <- readr::read_rds("Data/Output/sa_decreased_5_y50_e50_r50.rds")
sa_decreased_10 <- readr::read_rds("Data/Output/sa_decreased_10_y50_e50_r50.rds")

# get observation window
window <- readr::read_rds("Data/Raw/plot_area_owin.rds")

# rm(pattern_1999_recon)

# set parameters
overwrite <- FALSE
base_size <- 12.5

#### Preprocess data ####
#### Preprocess data ####
sa_increased_5 <- sa_increased_5[!names(sa_increased_5) %in% c("ci_max_dist", 
                                                               "seed_max_dist", 
                                                               "growth_mod")]

sa_increased_10 <- sa_increased_10[!names(sa_increased_10) %in% c("ci_max_dist", 
                                                                  "seed_max_dist", 
                                                                  "growth_mod")]

sa_decreased_5 <- sa_decreased_5[!names(sa_decreased_5) %in% c("ci_max_dist", 
                                                               "seed_max_dist", 
                                                               "growth_mod")]

sa_decreased_10 <- sa_decreased_10[!names(sa_decreased_10) %in% c("ci_max_dist", 
                                                                  "seed_max_dist",
                                                                  "growth_mod")]

parameter_levels <- rev(c("ci_alpha", "ci_beta", 
                          "growth_assymp", "growth_rate", "growth_infl", 
                          "seed_str", "seed_empty", "seed_success", "seed_beta", 
                          "mort_dbh_early", "mort_dbh_late", 
                          "mort_int_early", "mort_int_late", 
                          "mort_dinc"))

########################
####                ####
#### Integral value ####
####                ####
########################

#### Pair-correlation function ####

# set parameters #
correction <- "good"
divisor <- "d" 
r <- seq(from = 0, to = 25, length.out = 525)

# increased parameters #
sa_pcf_increased_5 <- calc_pcf_sa_int(default = sa_default,
                                      changed = sa_increased_5,
                                      correction = correction,
                                      divisor = divisor,
                                      window = window, r = r) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Increased 5%")

sa_pcf_increased_10 <- calc_pcf_sa_int(default = sa_default,
                                       changed = sa_increased_10,
                                       correction = correction,
                                       divisor = divisor,
                                       window = window, r = r) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Increased 10%")

sa_pcf_increased <- dplyr::bind_rows(sa_pcf_increased_5, sa_pcf_increased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_pcf_inc <- ggplot(data = sa_pcf_increased) + 
  geom_bar(aes(x = parameter, y = pcf_mean * 100, fill = direction),
           col = "black", stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = expression(paste("Relative difference", integral(g(r), 0, r))),
                     breaks = seq(-12.5, 12.5, 2.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_pcf_inc, 
                     filename = "ggplot_sa_pcf_inc.png", 
                     path = "Figures/", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300)

# decreased parameters #
sa_pcf_decreased_5 <- calc_pcf_sa_int(default = sa_default,
                                      changed = sa_decreased_5,
                                      correction = correction,
                                      divisor = divisor,
                                      window = window, r = r) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Decreased 5%")

sa_pcf_decreased_10 <- calc_pcf_sa_int(default = sa_default,
                                       changed = sa_decreased_10,
                                       correction = correction,
                                       divisor = divisor,
                                       window = window, r = r) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Decreased 10%")

sa_pcf_decreased <- dplyr::bind_rows(sa_pcf_decreased_5, sa_pcf_decreased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_pcf_dec <- ggplot(data = sa_pcf_decreased) + 
  geom_bar(aes(x = parameter, y = pcf_mean * 100, fill = direction),
           col = "black", stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = expression(paste("Relative difference", integral(g(r), 0, r))),
                     breaks = seq(-12.5, 12.5, 2.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_pcf_dec, 
                     filename = "ggplot_sa_pcf_dec.png", 
                     path = "Figures/", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300)

#### Nearest neighbor distribution function ####

# set parameters #
correction <- "km"
r <- seq(from = 0, to = 10, length.out = 525)

# increased parameters #
sa_nnd_increased_5 <- calc_nnd_sa_int(default = sa_default,
                                      changed = sa_increased_5,
                                      correction = correction,
                                      window = window, r = r) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Increased 5%")

sa_nnd_increased_10 <- calc_nnd_sa_int(default = sa_default,
                                       changed = sa_increased_10,
                                       correction = correction,
                                       window = window, r = r) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels), 
                direction = "Increased 10%")

sa_nnd_increased <- dplyr::bind_rows(sa_nnd_increased_5, sa_nnd_increased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_nnd_inc <- ggplot(data = sa_nnd_increased) + 
  geom_bar(aes(x = parameter, y = nnd_mean * 100, fill = direction),
           col = "black", stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = expression(paste("Relative difference", integral(G(r), 0, r))),
                     breaks = seq(-12.5, 12.5, 2.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_nnd_inc, 
                     filename = "ggplot_sa_nnd_inc.png", 
                     path = "Figures/", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300)

# decreased values #
sa_nnd_decreased_5 <- calc_nnd_sa_int(default = sa_default,
                                      changed = sa_decreased_5,
                                      correction = correction,
                                      window = window, r = r) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Decreased 5%")

sa_nnd_decreased_10 <- calc_nnd_sa_int(default = sa_default,
                                       changed = sa_decreased_10,
                                       correction = correction,
                                       window = window, r = r) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Decreased 10%")

sa_nnd_decreased <- dplyr::bind_rows(sa_nnd_decreased_5, sa_nnd_decreased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_nnd_dec <- ggplot(data = sa_nnd_decreased) + 
  geom_bar(aes(x = parameter, y = nnd_mean * 100, fill = direction),
           col = "black", stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = expression(paste("Relative difference", integral(G(r), 0, r))),
                     breaks = seq(-12.5, 12.5, 2.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_nnd_dec, 
                     filename = "ggplot_sa_nnd_dec.png", 
                     path = "Figures/", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300)

### Mark-correlation function ####

# set parameters #
correction <- "Ripley"
r <- seq(from = 0, to = 50, length.out = 525)

# increased values #
sa_kmm_increased_5 <- calc_kmm_sa_int(default = sa_default,
                                      changed = sa_increased_5,
                                      window = window, r = r, 
                                      correction = correction) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Increased 5%")

sa_kmm_increased_10 <- calc_kmm_sa_int(default = sa_default,
                                       changed = sa_increased_10,
                                       window = window, r = r, 
                                       correction = correction) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Increased 10%")

sa_kmm_increased <- dplyr::bind_rows(sa_kmm_increased_5, sa_kmm_increased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_kmm_inc <- ggplot(data = sa_kmm_increased) + 
  geom_bar(aes(x = parameter, y = kmm_mean * 100, fill = direction),
           col = "black", stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = expression(paste("Relative difference", integral(kmm(r), 0, r))),
                     breaks = seq(-12.5, 12.5, 2.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_kmm_inc, 
                     filename = "ggplot_sa_kmm_inc.png", 
                     path = "Figures/", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300)

# decreased values #
sa_kmm_decreased_5 <- calc_kmm_sa_int(default = sa_default,
                                      changed = sa_decreased_5,
                                      window = window, r = r, 
                                      correction = correction) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Decreased 5%")

sa_kmm_decreased_10 <- calc_kmm_sa_int(default = sa_default,
                                       changed = sa_decreased_10,
                                       window = window, r = r, 
                                       correction = correction) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Decreased 10%")

sa_kmm_decreased <- dplyr::bind_rows(sa_kmm_decreased_5, sa_kmm_decreased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_kmm_dec <- ggplot(data = sa_kmm_decreased) + 
  geom_bar(aes(x = parameter, y = kmm_mean * 100, fill = direction),
           col = "black", stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = expression(paste("Relative difference", integral(kmm(r), 0, r))),
                     breaks = seq(-12.5, 12.5, 2.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_kmm_dec, 
                     filename = "ggplot_sa_kmm_dec.png", 
                     path = "Figures/", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300)

#### Clark and Evans Index ####

# increased parameters #
sa_clark_increased_5 <- calc_clark_sa(default = sa_default,
                                      changed = sa_increased_5,
                                      window = window) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Increased 5%")

sa_clark_increased_10 <- calc_clark_sa(default = sa_default,
                                       changed = sa_increased_10,
                                       window = window) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Increased 10%")

sa_clark_increased <- dplyr::bind_rows(sa_clark_increased_5, sa_clark_increased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Increased 5%", 
                                              "Increased 10%")))

ggplot_sa_clark_inc <- ggplot(data = sa_clark_increased) + 
  geom_bar(aes(x = parameter, y = ce_mean * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter" ) +
  scale_y_continuous(name = "Relative difference CE index [%]", 
                     limits = c(-12.5, 12.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_clark_inc, 
                     filename = "ggplot_sa_clark_inc.png", 
                     path = "Figures/", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300)

# decreased values #
sa_clark_decreased_5 <- calc_clark_sa(default = sa_default,
                                      changed = sa_decreased_5,
                                      window = window) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels), 
                direction = "Decreased 5%")

sa_clark_decreased_10 <- calc_clark_sa(default = sa_default,
                                       changed = sa_decreased_10,
                                       window = window) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = "Decreased 10%")

sa_clark_decreased <- dplyr::bind_rows(sa_clark_decreased_5, sa_clark_decreased_10) %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("Decreased 5%", 
                                              "Decreased 10%")))

ggplot_sa_clark_dec <- ggplot(data = sa_clark_decreased) + 
  geom_bar(aes(x = parameter, y = ce_mean * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter" ) +
  scale_y_continuous(name = "Relative difference CE index [%]", 
                     limits = c(-12.5, 12.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_clark_dec, 
                     filename = "ggplot_sa_clark_dec.png", 
                     path = "Figures/", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300)

########################
####                ####
#### Function value ####
####                ####
########################

# #### Pair-correlation function ####
# 
# # set parameters #
# correction <- "good"
# divisor <- "d"
# fast <- FALSE
# r <- seq(from = 0, to = 25, length.out = 525)
# 
# # increased parameters #
# sa_pcf_increased_5 <- calc_pcf_sa_fun(default = sa_default,
#                                       changed = sa_increased_5,
#                                       correction = correction,
#                                       divisor = divisor,
#                                       window = window, r = r) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Increased 5%")
# 
# sa_pcf_increased_10 <- calc_pcf_sa_fun(default = sa_default,
#                                        changed = sa_increased_10,
#                                        correction = correction,
#                                        divisor = divisor,
#                                        window = window, r = r) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Increased 10%")
# 
# sa_pcf_increased <- dplyr::bind_rows(sa_pcf_increased_5, sa_pcf_increased_10) %>%
#   dplyr::mutate(direction = factor(direction,
#                                    levels = c("Increased 5%",
#                                               "Increased 10%")))
# 
# ggplot_sa_pcf_inc <- ggplot(data = sa_pcf_increased) +
#   geom_line(aes(x = r, y = pcf_diff * 100, col = parameter)) +
#   geom_hline(yintercept = 0) +
#   facet_wrap(~ direction) +
#   scale_x_continuous(name = "r [m]" ,) +
#   scale_y_continuous(name = "Relative difference pcf(r) [%]") +
#   scale_color_viridis_d(name = "Parameter change") +
#   theme_classic(base_size = base_size) +
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_pcf_inc,
#                      filename = "ggplot_sa_pcf_inc.png",
#                      path = "Figures/",
#                      width = 29.7, height = 21.0, units = "cm", dpi = 300)
# 
# # decreased parameters #
# sa_pcf_decreased_5 <- calc_pcf_sa_fun(default = sa_default,
#                                       changed = sa_decreased_5,
#                                       correction = correction,
#                                       divisor = divisor,
#                                       window = window, r = r) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Decreased 5%")
# 
# sa_pcf_decreased_10 <- calc_pcf_sa_fun(default = sa_default,
#                                        changed = sa_decreased_10,
#                                        correction = correction,
#                                        divisor = divisor,
#                                        window = window, r = r) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Decreased 10%")
# 
# sa_pcf_decreased <- dplyr::bind_rows(sa_pcf_decreased_5, sa_pcf_decreased_10) %>%
#   dplyr::mutate(direction = factor(direction,
#                                    levels = c("Decreased 5%",
#                                               "Decreased 10%")))
# 
# ggplot_sa_pcf_dec <- ggplot(data = sa_pcf_decreased) +
#   geom_line(aes(x = r, y = pcf_diff * 100, col = parameter)) +
#   geom_hline(yintercept = 0) +
#   facet_wrap(~ direction) +
#   scale_x_continuous(name = "r [m]" ,) +
#   scale_y_continuous(name = "Relative difference pcf(r) [%]") +
#   scale_color_viridis_d(name = "Parameter change") +
#   theme_classic(base_size = base_size) +
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_pcf_dec,
#                      filename = "ggplot_sa_pcf_dec.png",
#                      path = "Figures/",
#                      width = 29.7, height = 21.0, units = "cm", dpi = 300)
# 
# #### Nearest neighbor distribution function ####
# 
# # set parameters #
# correction <- "km"
# r <- seq(from = 0, to = 10, length.out = 525)
# 
# # increased parameters #
# sa_nnd_increased_5 <- calc_nnd_sa_fun(default = sa_default,
#                                       changed = sa_increased_5,
#                                       correction = correction,
#                                       window = window, r = r) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Increased 5%")
# 
# sa_nnd_increased_10 <- calc_nnd_sa_fun(default = sa_default,
#                                        changed = sa_increased_10,
#                                        correction = correction,
#                                        window = window, r = r) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Increased 10%")
# 
# sa_nnd_increased <- dplyr::bind_rows(sa_nnd_increased_5, sa_nnd_increased_10) %>%
#   dplyr::mutate(direction = factor(direction,
#                                    levels = c("Increased 5%",
#                                               "Increased 10%")))
# 
# ggplot_sa_nnd_inc <- ggplot(data = sa_nnd_increased) +
#   geom_line(aes(x = r, y = nnd_diff * 100, col = parameter)) +
#   geom_hline(yintercept = 0) +
#   facet_wrap(~ direction) +
#   scale_x_continuous(name = "r [m]" ,) +
#   scale_y_continuous(name = "Relative difference G(r) [%]") +
#   scale_color_viridis_d(name = "Parameter change") +
#   theme_classic(base_size = base_size) +
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_nnd_inc,
#                      filename = "ggplot_sa_nnd_inc.png",
#                      path = "Figures/",
#                      width = 29.7, height = 21.0, units = "cm", dpi = 300)
# 
# # decreased values #
# sa_nnd_decreased_5 <- calc_nnd_sa_fun(default = sa_default,
#                                       changed = sa_decreased_5,
#                                       correction = correction,
#                                       window = window, r = r) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Decreased 5%")
# 
# sa_nnd_decreased_10 <- calc_nnd_sa_fun(default = sa_default,
#                                        changed = sa_decreased_10,
#                                        correction = correction,
#                                        window = window, r = r) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Decreased 10%")
# 
# sa_nnd_decreased <- dplyr::bind_rows(sa_nnd_decreased_5, sa_nnd_decreased_10) %>%
#   dplyr::mutate(direction = factor(direction,
#                                    levels = c("Decreased 5%",
#                                               "Decreased 10%")))
# 
# ggplot_sa_nnd_dec <- ggplot(data = sa_nnd_decreased) +
#   geom_line(aes(x = r, y = nnd_diff * 100, col = parameter)) +
#   geom_hline(yintercept = 0) +
#   facet_wrap(~ direction) +
#   scale_x_continuous(name = "r [m]" ,) +
#   scale_y_continuous(name = "Relative difference G(r) [%]") +
#   scale_color_viridis_d(name = "Parameter change") +
#   theme_classic(base_size = base_size) +
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_nnd_dec,
#                      filename = "ggplot_sa_nnd_dec.png",
#                      path = "Figures/",
#                      width = 29.7, height = 21.0, units = "cm", dpi = 300)
# 
# ### Mark-correlation function ####
# 
# # set parameters #
# correction <- "Ripley"
# r <- seq(from = 0, to = 50, length.out = 525)
# 
# # increased values #
# sa_kmm_increased_5 <- calc_kmm_sa_fun(default = sa_default,
#                                       changed = sa_increased_5,
#                                       window = window, r = r,
#                                       correction = correction) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Increased 5%")
# 
# sa_kmm_increased_10 <- calc_kmm_sa_fun(default = sa_default,
#                                        changed = sa_increased_10,
#                                        window = window, r = r,
#                                        correction = correction) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Increased 10%")
# 
# sa_kmm_increased <- dplyr::bind_rows(sa_kmm_increased_5, sa_kmm_increased_10) %>%
#   dplyr::mutate(direction = factor(direction,
#                                    levels = c("Increased 5%",
#                                               "Increased 10%")))
# 
# ggplot_sa_kmm_inc <- ggplot(data = sa_kmm_increased) +
#   geom_line(aes(x = r, y = kmm_diff * 100, col = parameter)) +
#   geom_hline(yintercept = 0) +
#   facet_wrap(~ direction) +
#   scale_x_continuous(name = "r [m]" ,) +
#   scale_y_continuous(name = "Relative difference kmm(r) [%]") +
#   scale_color_viridis_d(name = "Parameter change") +
#   theme_classic(base_size = base_size) +
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_kmm_inc,
#                      filename = "ggplot_sa_kmm_inc.png",
#                      path = "Figures/",
#                      width = 29.7, height = 21.0, units = "cm", dpi = 300)
# 
# # decreased values #
# sa_kmm_decreased_5 <- calc_kmm_sa_fun(default = sa_default,
#                                       changed = sa_decreased_5,
#                                       window = window, r = r,
#                                       correction = correction) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Decreased 5%")
# 
# sa_kmm_decreased_10 <- calc_kmm_sa_fun(default = sa_default,
#                                        changed = sa_decreased_10,
#                                        window = window, r = r,
#                                        correction = correction) %>%
#   dplyr::mutate(parameter = factor(parameter),
#                 direction = "Decreased 10%")
# 
# sa_kmm_decreased <- dplyr::bind_rows(sa_kmm_decreased_5, sa_kmm_decreased_10) %>%
#   dplyr::mutate(direction = factor(direction,
#                                    levels = c("Decreased 5%",
#                                               "Decreased 10%")))
# 
# ggplot_sa_kmm_dec <- ggplot(data = sa_kmm_decreased) +
#   geom_line(aes(x = r, y = kmm_diff * 100, col = parameter)) +
#   geom_hline(yintercept = 0) +
#   facet_wrap(~ direction) +
#   scale_x_continuous(name = "r [m]" ,) +
#   scale_y_continuous(name = "Relative difference kmm(r) [%]") +
#   scale_color_viridis_d(name = "Parameter change") +
#   theme_classic(base_size = base_size) +
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_kmm_dec,
#                      filename = "ggplot_sa_kmm_dec.png",
#                      path = "Figures/",
#                      width = 29.7, height = 21.0, units = "cm", dpi = 300)
