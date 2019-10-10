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

base_size <- 12.5

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

# reverse parameter levels because of coord_flip()
parameter_levels <- rev(c("ci_alpha", "ci_beta", 
                          "growth_assymp", "growth_rate", "growth_infl", 
                          "seed_str", "seed_empty", "seed_success", "seed_beta", 
                          "mort_dbh_early", "mort_dbh_late", 
                          "mort_int_early", "mort_int_late", 
                          "mort_dinc"))

#### DBH distribution ####
# set parameters #
by <- 10

# Increased parameters #
sa_dbh_dist_inc_5 <- calc_dbh_dist_sa(default = sa_default, 
                                      changed = sa_increased_5, 
                                      by = by) %>% 
  dplyr::mutate(dbh_class = factor(dbh_class, ordered = TRUE), 
                parameter = factor(parameter, 
                                   levels = parameter_levels), 
                direction = "Increased 5%")

sa_dbh_dist_inc_10 <- calc_dbh_dist_sa(default = sa_default, 
                                       changed = sa_increased_10, 
                                       by = by) %>% 
  dplyr::mutate(dbh_class = factor(dbh_class, ordered = TRUE), 
                parameter = factor(parameter, 
                                   levels = parameter_levels), 
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
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) + 
  scale_x_discrete(name = "DBH class [cm]",
                   breaks = seq(from = as.numeric(min(sa_dbh_dist_inc$dbh_class)),
                                to = as.numeric(max(sa_dbh_dist_inc$dbh_class)),
                                by = 1),
                   labels = paste0("<", seq(from = as.numeric(min(sa_dbh_dist_inc$dbh_class)) * 10,
                                            to = as.numeric(max(sa_dbh_dist_inc$dbh_class)) * 10,
                                            by = by))) +
  scale_y_continuous(name = "Relative difference [%]") +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5))

suppoRt::save_ggplot(plot = ggplot_sa_dbh_dist_inc, 
                     filename = "ggplot_sa_dbh_dist_inc.png", 
                     path = "Figures/SA", 
                     width = 21.0, height = 29.7, units = "cm", dpi = 300, 
                     overwrite = overwrite)

# Decrease parameters #
sa_dbh_dist_dec_5 <- calc_dbh_dist_sa(default = sa_default, 
                                      changed = sa_decreased_5, 
                                      by = by) %>% 
  dplyr::mutate(dbh_class = factor(dbh_class, ordered = TRUE), 
                parameter = factor(parameter, 
                                   levels = parameter_levels), 
                direction = "Decreased 5%")

sa_dbh_dist_dec_10 <- calc_dbh_dist_sa(default = sa_default, 
                                       changed = sa_decreased_10, 
                                       by = by) %>% 
  dplyr::mutate(dbh_class = factor(dbh_class, ordered = TRUE), 
                parameter = factor(parameter, 
                                   levels = parameter_levels), 
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
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) + 
  scale_x_discrete(name = "DBH class [cm]",
                   breaks = seq(from = as.numeric(min(sa_dbh_dist_dec$dbh_class)),
                                to = as.numeric(max(sa_dbh_dist_dec$dbh_class)),
                                by = 1), 
                   labels = paste0("<", seq(from = as.numeric(min(sa_dbh_dist_dec$dbh_class)) * 10,
                                            to = as.numeric(max(sa_dbh_dist_dec$dbh_class)) * 10,
                                            by = by))) +
  scale_y_continuous(name = "Relative difference [%]") +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5))

suppoRt::save_ggplot(plot = ggplot_sa_dbh_dist_dec, 
                     filename = "ggplot_sa_dbh_dist_dec.png", 
                     path = "Figures/SA", 
                     width = 21.0, height = 29.7, units = "cm", dpi = 300, 
                     overwrite = overwrite)

#### DBH growth ####
# Increased parameters #
sa_growth_inc_5 <- calc_growth_sa(default = sa_default, 
                                  changed = sa_increased_5) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels), 
                direction = "Increased 5%")

sa_growth_inc_10 <- calc_growth_sa(default = sa_default, 
                                   changed = sa_increased_10) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels), 
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
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = "Difference mean annual DBH increment [%]", 
                     breaks = seq(-20, 20, 5), 
                     limits = c(-22.5, 22.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_growth_inc, 
                     filename = "ggplot_sa_growth_inc.png", 
                     path = "Figures/SA", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300, 
                     overwrite = overwrite)

# Decreased parameters #
sa_growth_dec_5 <- calc_growth_sa(default = sa_default, 
                                  changed = sa_decreased_5) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels), 
                direction = "Decreased 5%")

sa_growth_dec_10 <- calc_growth_sa(default = sa_default, 
                                   changed = sa_decreased_10) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels), 
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
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter" ,) +
  scale_y_continuous(name = "Difference mean annual DBH increment [%]", 
                     breaks = seq(-20, 20, 5), 
                     limits = c(-22.5, 22.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_growth_dec, 
                     filename = "ggplot_sa_growth_dec.png", 
                     path = "Figures/SA", 
                     width = 17.5, height = 12.5 , units = "cm", dpi = 300, 
                     overwrite = overwrite)

#### n died ####
# Increased parameters #
sa_died_inc_5 <- calc_died_sa(default = sa_default,
                              changed = sa_increased_5) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),  
                direction = "Increased 5%")

sa_died_inc_10 <- calc_died_sa(default = sa_default,
                               changed = sa_increased_10) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels), 
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
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = "Difference n died [%]",
                     breaks = seq(-35, 35, 5),
                     limits = c(-37.5, 37.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_died_inc, 
                     filename = "ggplot_sa_died_inc.png", 
                     path = "Figures/SA", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300, 
                     overwrite = overwrite)

# Decreased parameters #
sa_died_dec_5 <- calc_died_sa(default = sa_default,
                              changed = sa_decreased_5) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),  
                direction = "Decreased 5%")

sa_died_dec_10 <- calc_died_sa(default = sa_default,
                               changed = sa_decreased_10) %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels), 
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
  geom_hline(yintercept = -10, linetype = 2, col = "#21908CFF") + 
  geom_hline(yintercept = -5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 0, linetype = 1) + 
  geom_hline(yintercept = 5, linetype = 2, col = "#440154FF") + 
  geom_hline(yintercept = 10, linetype = 2, col = "#21908CFF") + 
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = "Difference n died [%]",
                     breaks = seq(-60, 60, 5),
                     limits = c(-62.5, 62.5)) +
  scale_fill_manual(name = "Parameter change", 
                    values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_died_dec, 
                     filename = "ggplot_sa_died_dec.png", 
                     path = "Figures/SA", 
                     width = 17.5, height = 12.5, units = "cm", dpi = 300, 
                     overwrite = overwrite)
