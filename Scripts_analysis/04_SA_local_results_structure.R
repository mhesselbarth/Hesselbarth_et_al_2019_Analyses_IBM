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
                          "seed_str", "seed_empty", "seed_success", "seed_eta", 
                          "mort_dbh_early", "mort_dbh_late", 
                          "mort_int_early", "mort_int_late", 
                          "mort_dinc"))

#### Preprocess data ####
# filter size classes # 

sa_default_sapling <- purrr::map(sa_default, function(x)
  dplyr::filter(x, type == "sapling"))

sa_default_adult <- purrr::map(sa_default, function(x) 
  dplyr::filter(x, type == "adult"))

rm(sa_default)

sa_increased_5_sapling <- purrr::map(sa_increased_5, function(x) 
  dplyr::filter(x, type == "sapling"))

sa_increased_5_adult <- purrr::map(sa_increased_5, function(x) 
  dplyr::filter(x, type == "adult"))

rm(sa_increased_5)

sa_increased_10_sapling <- purrr::map(sa_increased_10, function(x) 
  dplyr::filter(x, type == "sapling"))

sa_increased_10_adult <- purrr::map(sa_increased_10, function(x) 
  dplyr::filter(x, type == "adult"))

rm(sa_increased_10)

sa_decreased_5_sapling <- purrr::map(sa_decreased_5, function(x) 
  dplyr::filter(x, type == "sapling"))

sa_decreased_5_adult <- purrr::map(sa_decreased_5, function(x) 
  dplyr::filter(x, type == "adult"))

rm(sa_decreased_5)

sa_decreased_10_sapling <- purrr::map(sa_decreased_10, function(x) 
  dplyr::filter(x, type == "sapling"))

sa_decreased_10_adult <- purrr::map(sa_decreased_10, function(x) 
  dplyr::filter(x, type == "adult"))

rm(sa_decreased_10)

#### n individuals ####

# Increased parameters #
sa_individuals_inc_5_sapling <- calc_n_sa(default = sa_default_sapling,
                                          changed = sa_increased_5_sapling) %>%
  dplyr::mutate(size = "sapling", 
                direction = "Increased +5%")

sa_individuals_inc_5_adult <- calc_n_sa(default = sa_default_adult,
                                        changed = sa_increased_5_adult) %>%
  dplyr::mutate(size = "adult", 
                direction = "Increased +5%")

sa_individuals_inc_10_sapling <- calc_n_sa(default = sa_default_sapling,
                                           changed = sa_increased_10_sapling) %>%
  dplyr::mutate(size = "sapling", 
                direction = "Increased +10%")

sa_individuals_inc_10_adult <- calc_n_sa(default = sa_default_adult,
                                         changed = sa_increased_10_adult) %>%
  dplyr::mutate(size = "adult", 
                direction = "Increased +10%")

# Decreased parameters #
sa_individuals_dec_5_sapling <- calc_n_sa(default = sa_default_sapling,
                                          changed = sa_decreased_5_sapling) %>%
  dplyr::mutate(size = "sapling", 
                direction = "Decreased -5%")

sa_individuals_dec_5_adult <- calc_n_sa(default = sa_default_adult,
                                        changed = sa_decreased_5_adult) %>%
  dplyr::mutate(size = "adult", 
                direction = "Decreased -5%")

sa_individuals_dec_10_sapling <- calc_n_sa(default = sa_default_sapling,
                                           changed = sa_decreased_10_sapling) %>%
  dplyr::mutate(size = "sapling", 
                direction = "Decreased -10%")

sa_individuals_dec_10_adult <- calc_n_sa(default = sa_default_adult,
                                         changed = sa_decreased_10_adult) %>%
  dplyr::mutate(size = "adult", 
                direction = "Decreased -10%")

# combine to one dataframe # 
sa_individuals <- dplyr::bind_rows(sa_individuals_inc_5_sapling,
                                   sa_individuals_inc_5_adult,
                                   sa_individuals_inc_10_sapling,
                                   sa_individuals_inc_10_adult,
                                   sa_individuals_dec_5_sapling,
                                   sa_individuals_dec_5_adult,
                                   sa_individuals_dec_10_sapling,
                                   sa_individuals_dec_10_adult) %>% 
  dplyr::mutate(parameter = dplyr::case_when(parameter == "seed_beta" ~  "seed_eta", 
                                             TRUE ~ parameter)) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                size = factor(size, levels = c("sapling", "adult"), 
                              labels = c("Sapling", "Adult")),
                direction = factor(direction, 
                                   levels = c("Decreased -10%", 
                                              "Decreased -5%",
                                              "Increased +5%", 
                                              "Increased +10%")))

ggplot_sa_individuals <- ggplot(data = sa_individuals) + 
  geom_bar(aes(x = parameter, y = diff_n * 100, 
               fill = direction), col = "black",
           stat = "identity", position = "dodge") + 
  geom_hline(yintercept = -10, linetype = 2, col = "#000004FF") +
  geom_hline(yintercept = -5, linetype = 2, col = "#781C6DFF") +
  geom_hline(yintercept = 0, linetype = 1) +
  geom_hline(yintercept = 5, linetype = 2, col = "#ED7953FF") +
  geom_hline(yintercept = 10, linetype = 2, col = "#FCFFA4FF") +
  coord_flip() +
  facet_wrap(~ size) +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = "Difference individuals [%]",
                     breaks = seq(-100, 100, 25),
                     limits = c(-100, 100)) +
  scale_fill_manual(name = "Parameter change",
                    values = c("#000004FF", "#781C6DFF" ,
                               "#ED6925FF", "#FCFFA4FF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_individuals, 
                     filename = "ggplot_sa_individuals.png", 
                     path = "Figures/Appendix/",     
                     dpi = dpi,
                     width = width_full, height = height_full * (2/3), units = units, 
                     overwrite = overwrite)

#### DBH growth ####
# 
# # Increased parameters #
# sa_growth_inc_5 <- calc_growth_sa(default = sa_default, 
#                                   changed = sa_increased_5) %>%
#   dplyr::mutate(direction = "Increased +5%")
# 
# sa_growth_inc_10 <- calc_growth_sa(default = sa_default, 
#                                    changed = sa_increased_10) %>%
#   dplyr::mutate(direction = "Increased +10%")
# 
# # Decreased parameters #
# sa_growth_dec_5 <- calc_growth_sa(default = sa_default, 
#                                   changed = sa_decreased_5) %>%
#   dplyr::mutate(direction = "Decreased -5%")
# 
# sa_growth_dec_10 <- calc_growth_sa(default = sa_default, 
#                                    changed = sa_decreased_10) %>%
#   dplyr::mutate(direction = "Decreased -10%")
# 
# sa_growth <- dplyr::bind_rows(sa_growth_inc_5,
#                               sa_growth_inc_10,
#                               sa_growth_dec_5, 
#                               sa_growth_dec_10) %>% 
#   dplyr::mutate(parameter = factor(parameter, 
#                                    levels = parameter_levels),
#                 direction = factor(direction, 
#                                    levels = c("Decreased -10%", 
#                                               "Decreased -5%", 
#                                               "Increased +5%", 
#                                               "Increased +10%")))
# 
# ggplot_sa_growth <- ggplot(data = sa_growth) + 
#   geom_bar(aes(x = parameter, y = diff_inc * 100, 
#                fill = direction), col = "black",
#            stat = "identity", position = "dodge") + 
#   geom_hline(yintercept = -10, linetype = 2, col = "#0D0887FF") +
#   geom_hline(yintercept = -5, linetype = 2, col = "#9C179EFF") +
#   geom_hline(yintercept = 0, linetype = 1) +
#   geom_hline(yintercept = 5, linetype = 2, col = "#ED7953FF") +
#   geom_hline(yintercept = 10, linetype = 2, col = "#F0F921FF") +
#   coord_flip() +
#   scale_x_discrete(name = "Parameter") +
#   scale_y_continuous(name = "Difference mean annual DBH increment [%]", 
#                      breaks = seq(-20, 20, 5), 
#                      limits = c(-22.5, 22.5)) +
#   scale_fill_manual(name = "Parameter change",
#                     values = c("#0D0887FF", "#9C179EFF" ,
#                                "#ED7953FF", "#F0F921FF")) +
#   theme_classic(base_size = base_size) + 
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_growth, 
#                      filename = "ggplot_sa_growth.png", 
#                      path = "Figures/Appendix/", 
#                      dpi = dpi,
#                      width = width_full, height = height_full, units = units, 
#                      overwrite = overwrite)
# 
# #### n died ####
# 
# # Increased parameters #
# sa_died_inc_5 <- calc_died_sa(default = sa_default,
#                               changed = sa_increased_5) %>%
#   dplyr::mutate(direction = "Increased +5%")
# 
# sa_died_inc_10 <- calc_died_sa(default = sa_default,
#                                changed = sa_increased_10) %>%
#   dplyr::mutate(direction = "Increased +10%")
# 
# # Decreased parameters #
# sa_died_dec_5 <- calc_died_sa(default = sa_default,
#                               changed = sa_decreased_5) %>%
#   dplyr::mutate(direction = "Decreased -5%")
# 
# sa_died_dec_10 <- calc_died_sa(default = sa_default,
#                                changed = sa_decreased_10) %>%
#   dplyr::mutate(direction = "Decreased -10%")
# 
# sa_died <- dplyr::bind_rows(sa_died_inc_5, 
#                             sa_died_inc_10, 
#                             sa_died_dec_5, 
#                             sa_died_dec_10) %>% 
#   dplyr::mutate(parameter = factor(parameter, 
#                                    levels = parameter_levels),
#                 direction = factor(direction, 
#                                    levels = c("Decreased -10%", 
#                                               "Decreased -5%", 
#                                               "Increased +5%", 
#                                               "Increased +10%")))
# 
# ggplot_sa_died <- ggplot(data = sa_died) + 
#   geom_bar(aes(x = parameter, y = diff_n_rel * 100, 
#                fill = direction), col = "black",
#            stat = "identity", position = "dodge") + 
#   geom_hline(yintercept = -10, linetype = 2, col = "#0D0887FF") +
#   geom_hline(yintercept = -5, linetype = 2, col = "#9C179EFF") +
#   geom_hline(yintercept = 0, linetype = 1) +
#   geom_hline(yintercept = 5, linetype = 2, col = "#ED7953FF") +
#   geom_hline(yintercept = 10, linetype = 2, col = "#F0F921FF") +
#   coord_flip() +
#   scale_x_discrete(name = "Parameter") +
#   scale_y_continuous(name = "Difference n died [%]",
#                      breaks = seq(-45, 45, 5),
#                      limits = c(-45, 45)) +
#   scale_fill_manual(name = "Parameter change",
#                     values = c("#0D0887FF", "#9C179EFF" ,
#                                "#ED7953FF", "#F0F921FF")) +
#   theme_classic(base_size = base_size) + 
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_died, 
#                      filename = "ggplot_sa_died.png", 
#                      path = "Figures/Appendix/",     
#                      dpi = dpi,
#                      width = width_full, height = height_small, units = units, 
#                      overwrite = overwrite)
