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
source("Helper_functions/helper_functions_setup.R")
source("Helper_functions/helper_functions_sa_structure.R")

#### Import data ####

# import model runs default
sa_default <- readr::read_rds("Data/Output/SA/sa_default_y50_e5_r50.rds")

# import increased parameters
sa_increased_5 <- readr::read_rds("Data/Output/SA/sa_increased_5_y50_e5_r50.rds")
sa_increased_10 <- readr::read_rds("Data/Output/SA/sa_increased_10_y50_e5_r50.rds")

# import decreased parameters
sa_decreased_5 <- readr::read_rds("Data/Output/SA/sa_decreased_5_y50_e5_r50.rds")
sa_decreased_10 <- readr::read_rds("Data/Output/SA/sa_decreased_10_y50_e5_r50.rds")

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
  dplyr::mutate(direction = "Increased +5%")

sa_dbh_dist_inc_10 <- calc_dbh_dist_sa(default = sa_default, 
                                       changed = sa_increased_10, 
                                       by = by) %>% 
  dplyr::mutate(direction = "Increased +10%")

# Decrease parameters #
sa_dbh_dist_dec_5 <- calc_dbh_dist_sa(default = sa_default, 
                                      changed = sa_decreased_5, 
                                      by = by) %>% 
  dplyr::mutate(direction = "Decreased -5%")

sa_dbh_dist_dec_10 <- calc_dbh_dist_sa(default = sa_default, 
                                       changed = sa_decreased_10, 
                                       by = by) %>% 
  dplyr::mutate(direction = "Decreased -10%")

sa_dbh_dist <- dplyr::bind_rows(sa_dbh_dist_inc_5, 
                                sa_dbh_dist_inc_10,
                                sa_dbh_dist_dec_5, 
                                sa_dbh_dist_dec_10) %>% 
  dplyr::mutate(dbh_class = factor(dbh_class, ordered = TRUE), 
                parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = factor(direction, 
                                   levels = c("Decreased -10%", 
                                              "Decreased -5%", 
                                              "Increased +5%", 
                                              "Increased +10%")))

ggplot_sa_dbh_dist <- ggplot(data = sa_dbh_dist) + 
  geom_bar(aes(x = dbh_class, y = diff_n_rel * 100, 
               fill = direction), col = "black", 
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#440154FF") +
  geom_hline(yintercept = -5, linetype = 2, col = "#31688EFF") +
  geom_hline(yintercept = 0, linetype = 1) +
  geom_hline(yintercept = 5, linetype = 2, col = "#35B779FF") +
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") +
  facet_wrap(~ parameter, scales = "free_y", ncol = 4) + 
  scale_x_discrete(name = "dbh class [cm]",
                   breaks = seq(from = as.numeric(min(sa_dbh_dist_dec$dbh_class)),
                                to = as.numeric(max(sa_dbh_dist_dec$dbh_class)),
                                by = 1), 
                   labels = paste0("<", seq(from = as.numeric(min(sa_dbh_dist$dbh_class)) * 10,
                                            to = as.numeric(max(sa_dbh_dist$dbh_class)) * 10,
                                            by = by))) +
  scale_y_continuous(name = "Relative difference [%]") +
  scale_fill_manual(name = "Parameter change",
                    values = c("#440154FF", "#31688EFF", 
                               "#35B779FF","#FDE725FF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5))

suppoRt::save_ggplot(plot = ggplot_sa_dbh_dist, 
                     filename = "ggplot_sa_dbh_dist.png", 
                     path = "Figures/Appendix", 
                     dpi = dpi,
                     width = height_full, height = width_full, units = units,
                     overwrite = overwrite)

#### DBH growth ####

# Increased parameters #
sa_growth_inc_5 <- calc_growth_sa(default = sa_default, 
                                  changed = sa_increased_5) %>%
  dplyr::mutate(direction = "Increased +5%")

sa_growth_inc_10 <- calc_growth_sa(default = sa_default, 
                                   changed = sa_increased_10) %>%
  dplyr::mutate(direction = "Increased +10%")

# Decreased parameters #
sa_growth_dec_5 <- calc_growth_sa(default = sa_default, 
                                  changed = sa_decreased_5) %>%
  dplyr::mutate(direction = "Decreased -5%")

sa_growth_dec_10 <- calc_growth_sa(default = sa_default, 
                                   changed = sa_decreased_10) %>%
  dplyr::mutate(direction = "Decreased -10%")

sa_growth <- dplyr::bind_rows(sa_growth_inc_5,
                              sa_growth_inc_10,
                              sa_growth_dec_5, 
                              sa_growth_dec_10) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = factor(direction, 
                                   levels = c("Decreased -10%", 
                                              "Decreased -5%", 
                                              "Increased +5%", 
                                              "Increased +10%")))

ggplot_sa_growth <- ggplot(data = sa_growth) + 
  geom_bar(aes(x = parameter, y = diff_inc * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#440154FF") +
  geom_hline(yintercept = -5, linetype = 2, col = "#31688EFF") +
  geom_hline(yintercept = 0, linetype = 1) +
  geom_hline(yintercept = 5, linetype = 2, col = "#35B779FF") +
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") +
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = "Difference mean annual DBH increment [%]", 
                     breaks = seq(-20, 20, 5), 
                     limits = c(-22.5, 22.5)) +
  scale_fill_manual(name = "Parameter change",
                    values = c("#440154FF", "#31688EFF", 
                               "#35B779FF","#FDE725FF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_growth, 
                     filename = "ggplot_sa_growth.png", 
                     path = "Figures/Appendix/", 
                     dpi = dpi,
                     width = width_full, height = height_small, units = units, 
                     overwrite = overwrite)

#### n died ####

# Increased parameters #
sa_died_inc_5 <- calc_died_sa(default = sa_default,
                              changed = sa_increased_5) %>%
  dplyr::mutate(direction = "Increased +5%")

sa_died_inc_10 <- calc_died_sa(default = sa_default,
                               changed = sa_increased_10) %>%
  dplyr::mutate(direction = "Increased +10%")

# Decreased parameters #
sa_died_dec_5 <- calc_died_sa(default = sa_default,
                              changed = sa_decreased_5) %>%
  dplyr::mutate(direction = "Decreased -5%")

sa_died_dec_10 <- calc_died_sa(default = sa_default,
                               changed = sa_decreased_10) %>%
  dplyr::mutate(direction = "Decreased -10%")

sa_died <- dplyr::bind_rows(sa_died_inc_5, 
                            sa_died_inc_10, 
                            sa_died_dec_5, 
                            sa_died_dec_10) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                direction = factor(direction, 
                                   levels = c("Decreased -10%", 
                                              "Decreased -5%", 
                                              "Increased +5%", 
                                              "Increased +10%")))

ggplot_sa_died <- ggplot(data = sa_died) + 
  geom_bar(aes(x = parameter, y = diff_n_rel * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#440154FF") +
  geom_hline(yintercept = -5, linetype = 2, col = "#31688EFF") +
  geom_hline(yintercept = 0, linetype = 1) +
  geom_hline(yintercept = 5, linetype = 2, col = "#35B779FF") +
  geom_hline(yintercept = 10, linetype = 2, col = "#FDE725FF") +
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = "Difference n died [%]",
                     breaks = seq(-45, 45, 5),
                     limits = c(-45, 45)) +
  scale_fill_manual(name = "Parameter change",
                    values = c("#440154FF", "#31688EFF", 
                               "#35B779FF","#FDE725FF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_died, 
                     filename = "ggplot_sa_died.png", 
                     path = "Figures/Appendix/",     
                     dpi = dpi,
                     width = width_full, height = height_small, units = units, 
                     overwrite = overwrite)
