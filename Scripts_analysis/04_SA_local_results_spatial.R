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
source("Helper_functions/helper_functions_setup.R")
source("Helper_functions/helper_functions_sa_spatial.R")

# import data #
beech_1999_rec <- readr::read_rds("Data/Input/beech_1999_rec_ppp.rds")

# import model runs default
sa_default <- readr::read_rds("Data/Output/SA/sa_default_y50_e5_r50.rds")

# import increased parameters
sa_increased_5 <- readr::read_rds("Data/Output/SA/sa_increased_5_y50_e5_r50.rds")
sa_increased_10 <- readr::read_rds("Data/Output/SA/sa_increased_10_y50_e5_r50.rds")

# import decreased parameters
sa_decreased_5 <- readr::read_rds("Data/Output/SA/sa_decreased_5_y50_e5_r50.rds")
sa_decreased_10 <- readr::read_rds("Data/Output/SA/sa_decreased_10_y50_e5_r50.rds")

# get observation window
window <- readr::read_rds("Data/Raw/plot_area_owin.rds")

# rm(pattern_1999_recon)

#### Preprocess data ####
parameter_levels <- rev(c("ci_alpha", "ci_beta", 
                          "growth_assymp", "growth_rate", "growth_infl", 
                          "seed_str", "seed_empty", "seed_success", "seed_beta", 
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

########################
####                ####
#### Integral value ####
####                ####
########################

# #### Nearest neighbor distribution function ####
# 
# # set parameters #
# correction_nnd <- "km"
# r_nnd <- seq(from = 0, to = 10, length.out = 525)
# 
# # increased parameters #
# sa_nnd_increased_5_sapling <- calc_nnd_sa_int(default = sa_default_sapling,
#                                               changed = sa_increased_5_sapling,
#                                               correction = correction_nnd,
#                                               r = r_nnd,
#                                               window = window) %>% 
#   dplyr::mutate(direction = "Increased +5%", 
#                 size = "sapling")
# 
# sa_nnd_increased_5_adult <- calc_nnd_sa_int(default = sa_default_adult,
#                                       changed = sa_increased_5_adult,
#                                       correction = correction_nnd,
#                                       r = r_nnd,
#                                       window = window) %>% 
#   dplyr::mutate(direction = "Increased +5%", 
#                 size = "adult")
# 
# sa_nnd_increased_10_sapling <- calc_nnd_sa_int(default = sa_default_sapling,
#                                                changed = sa_increased_10_sapling,
#                                                correction = correction_nnd,
#                                                r = r_nnd,
#                                                window = window) %>% 
#   dplyr::mutate(size = "sapling",
#                 direction = "Increased +10%")
# 
# sa_nnd_increased_10_adult <- calc_nnd_sa_int(default = sa_default_adult,
#                                              changed = sa_increased_10_adult,
#                                              correction = correction_nnd,
#                                              r = r_nnd,
#                                              window = window) %>% 
#   dplyr::mutate(size = "adult",
#                 direction = "Increased +10%")
# 
# # decreased values #
# sa_nnd_decreased_5_sapling <- calc_nnd_sa_int(default = sa_default_sapling,
#                                               changed = sa_decreased_5_sapling,
#                                               correction = correction_nnd,
#                                               r = r_nnd,
#                                               window = window) %>% 
#   dplyr::mutate(direction = "Decreased -5%", 
#                 size = "sapling")
# 
# sa_nnd_decreased_5_adult <- calc_nnd_sa_int(default = sa_default_adult,
#                                             changed = sa_decreased_5_adult,
#                                             correction = correction_nnd,
#                                             r = r_nnd,
#                                             window = window) %>% 
#   dplyr::mutate(direction = "Decreased -5%", 
#                 size = "adult")
# 
# sa_nnd_decreased_10_sapling <- calc_nnd_sa_int(default = sa_default_sapling,
#                                                changed = sa_decreased_10_sapling,
#                                                correction = correction_nnd,
#                                                r = r_nnd,
#                                                window = window) %>% 
#   dplyr::mutate(size = "sapling",
#                 direction = "Decreased -10%")
# 
# sa_nnd_decreased_10_adult <- calc_nnd_sa_int(default = sa_default_adult,
#                                              changed = sa_decreased_10_adult,
#                                              correction = correction_nnd,
#                                              r = r_nnd,
#                                              window = window) %>% 
#   dplyr::mutate(size = "adult",
#                 direction = "Decreased -10%")
# 
# sa_nnd <- dplyr::bind_rows(sa_nnd_increased_5_sapling,
#                            sa_nnd_increased_5_adult,
#                            sa_nnd_increased_10_sapling,
#                            sa_nnd_increased_10_adult,
#                            sa_nnd_decreased_5_sapling,
#                            sa_nnd_decreased_5_adult,
#                            sa_nnd_decreased_10_sapling, 
#                            sa_nnd_increased_10_adult) %>% 
#   dplyr::mutate(parameter = factor(parameter, 
#                                    levels = parameter_levels),
#                 size = factor(size, levels = c("sapling", "adult"), 
#                               labels = c("Sapling", "Adult")),
#                 direction = factor(direction, 
#                                    levels = c("Decreased -10%", 
#                                               "Decreased -5%",
#                                               "Increased +5%", 
#                                               "Increased +10%")))
# 
# ggplot_sa_nnd <- ggplot(data = sa_nnd) + 
#   geom_bar(aes(x = parameter, y = nnd_mean * 100, fill = direction),
#            col = "black", stat = "identity", position = position_dodge()) +
#   facet_wrap(~ size) +
#   geom_hline(yintercept = -10, linetype = 2, col = "#0D0887FF") +
#   geom_hline(yintercept = -5, linetype = 2, col = "#9C179EFF") +
#   geom_hline(yintercept = 0, linetype = 1) +
#   geom_hline(yintercept = 5, linetype = 2, col = "#ED7953FF") +
#   geom_hline(yintercept = 10, linetype = 2, col = "#F0F921FF") +
#   coord_flip() +
#   scale_x_discrete(name = "Parameter") +
#   scale_y_continuous(name = expression(paste("Relative difference", integral(G(r), 0, r))),
#                      breaks = seq(-20, 20, 5)) +
#   scale_fill_manual(name = "Parameter change",
#                     values = c("#0D0887FF", "#9C179EFF" ,
#                                "#ED7953FF", "#F0F921FF")) +
#   theme_classic(base_size = base_size) + 
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_nnd, 
#                      filename = "ggplot_sa_nnd.png", 
#                      path = "Figures/Appendix/", 
#                      dpi = 300,
#                      width = width_full, height = height_small, units = units)

#### Pair-correlation function ####

# set parameters #
r_pcf <- seq(from = 0, to = 50, length.out = 525)
correction_pcf <- "Ripley"
stoyan_pcf <- 0.25
divisor_pcf <- "d"

# increased parameters #
sa_pcf_increased_5_sapling <- calc_pcf_sa_int(default = sa_default_sapling,
                                              changed = sa_increased_5_sapling,
                                              correction = correction_pcf,
                                              divisor = divisor_pcf,
                                              stoyan = stoyan_pcf,
                                              r = r_pcf,
                                              window = window) %>% 
  dplyr::mutate(size = "sapling", 
                direction = "Increased +5%")

sa_pcf_increased_5_adult <- calc_pcf_sa_int(default = sa_default_adult,
                                            changed = sa_increased_5_adult,
                                            correction = correction_pcf,
                                            divisor = divisor_pcf,
                                            stoyan = stoyan_pcf,
                                            r = r_pcf,
                                            window = window) %>% 
  dplyr::mutate(size = "adult", 
                direction = "Increased +5%")

sa_pcf_increased_10_sapling <- calc_pcf_sa_int(default = sa_default_sapling,
                                               changed = sa_increased_10_sapling,
                                               correction = correction_pcf,
                                               divisor = divisor_pcf,
                                               stoyan = stoyan_pcf,
                                               r = r_pcf,
                                               window = window) %>% 
  dplyr::mutate(size = "sapling", 
                direction = "Increased +10%")

sa_pcf_increased_10_adult <- calc_pcf_sa_int(default = sa_default_adult,
                                             changed = sa_increased_10_adult,
                                             correction = correction_pcf,
                                             divisor = divisor_pcf,
                                             stoyan = stoyan_pcf,
                                             r = r_pcf,
                                             window = window) %>% 
  dplyr::mutate(size = "adult", 
                direction = "Increased +10%")

# decreased parameters #
sa_pcf_decreased_5_sapling <- calc_pcf_sa_int(default = sa_default_sapling,
                                              changed = sa_decreased_5_sapling,
                                              correction = correction_pcf,
                                              divisor = divisor_pcf,
                                              stoyan = stoyan_pcf,
                                              r = r_pcf,
                                              window = window) %>% 
  dplyr::mutate(size = "sapling", 
                direction = "Decreased -5%")

sa_pcf_decreased_5_adult <- calc_pcf_sa_int(default = sa_default_adult,
                                            changed = sa_decreased_5_adult,
                                            correction = correction_pcf,
                                            divisor = divisor_pcf,
                                            stoyan = stoyan_pcf,
                                            r = r_pcf,
                                            window = window) %>% 
  dplyr::mutate(size = "adult", 
                direction = "Decreased -5%")

sa_pcf_decreased_10_sapling <- calc_pcf_sa_int(default = sa_default_sapling,
                                               changed = sa_decreased_10_sapling,
                                               correction = correction_pcf,
                                               divisor = divisor_pcf,
                                               stoyan = stoyan_pcf,
                                               r = r_pcf,
                                               window = window) %>% 
  dplyr::mutate(size = "sapling", 
                direction = "Decreased -10%")

sa_pcf_decreased_10_adult <- calc_pcf_sa_int(default = sa_default_adult,
                                             changed = sa_decreased_10_adult,
                                             correction = correction_pcf,
                                             divisor = divisor_pcf,
                                             stoyan = stoyan_pcf,
                                             r = r_pcf,
                                             window = window) %>% 
  dplyr::mutate(size = "adult", 
                direction = "Decreased -10%")


sa_pcf <- dplyr::bind_rows(sa_pcf_increased_5_sapling,
                           sa_pcf_increased_5_adult,
                           sa_pcf_increased_10_sapling,
                           sa_pcf_increased_10_adult,
                           sa_pcf_decreased_5_sapling,
                           sa_pcf_decreased_5_adult,
                           sa_pcf_decreased_10_sapling,
                           sa_pcf_decreased_10_adult,) %>% 
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = parameter_levels),
                size = factor(size, levels = c("sapling", "adult"), 
                              labels = c("Sapling", "Adult")),
                direction = factor(direction, 
                                   levels = c("Decreased -10%", 
                                              "Decreased -5%",
                                              "Increased +5%", 
                                              "Increased +10%")))

ggplot_sa_pcf <- ggplot(data = sa_pcf) + 
  geom_bar(aes(x = parameter, y = pcf_mean * 100, fill = direction),
           col = "black", stat = "identity", position = position_dodge()) +
  facet_wrap(~ size) +
  geom_hline(yintercept = -10, linetype = 2, col = "#0D0887FF") +
  geom_hline(yintercept = -5, linetype = 2, col = "#9C179EFF") +
  geom_hline(yintercept = 0, linetype = 1) +
  geom_hline(yintercept = 5, linetype = 2, col = "#ED7953FF") +
  geom_hline(yintercept = 10, linetype = 2, col = "#F0F921FF") +
  coord_flip() +
  scale_x_discrete(name = "Parameter") +
  scale_y_continuous(name = expression(paste("Relative difference", integral(g(r), 0, r))),
                     breaks = seq(-12.5, 12.5, 2.5)) +
  scale_fill_manual(name = "Parameter change",
                    values = c("#0D0887FF", "#9C179EFF" ,
                               "#ED7953FF", "#F0F921FF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sa_pcf, 
                     filename = "ggplot_sa_pcf.png", 
                     path = "Figures/Appendix/",
                     dpi = dpi, 
                     width = width_full, height = height_full * (2/3),
                     units = units, 
                     overwrite = overwrite)

### Mark-correlation function ####
# 
# # set parameters #
# correction_kmm <- "Ripley"
# r_kmm <- seq(from = 0, to = 50, length.out = 525)
# 
# # increased values #
# sa_kmm_increased_5_sapling <- calc_kmm_sa_int(default = sa_default_sapling,
#                                               changed = sa_increased_5_sapling,
#                                               r = r_kmm, 
#                                               correction = correction_kmm, 
#                                               window = window) %>% 
#   dplyr::mutate(size = "sapling", 
#                 direction = "Increased +5%")
# 
# sa_kmm_increased_5_adult <- calc_kmm_sa_int(default = sa_default_adult,
#                                             changed = sa_increased_5_adult,
#                                             r = r_kmm, 
#                                             correction = correction_kmm, 
#                                             window = window) %>% 
#   dplyr::mutate(size = "adult", 
#                 direction = "Increased +5%")
# 
# sa_kmm_increased_10_sapling <- calc_kmm_sa_int(default = sa_default_sapling,
#                                                changed = sa_increased_10_sapling,
#                                                r = r_kmm, 
#                                                correction = correction_kmm, 
#                                                window = window) %>% 
#   dplyr::mutate(size = "sapling", 
#                 direction = "Increased +10%")
# 
# sa_kmm_increased_10_adult <- calc_kmm_sa_int(default = sa_default_adult,
#                                              changed = sa_increased_10_adult,
#                                              r = r_kmm, 
#                                              correction = correction_kmm, 
#                                              window = window) %>% 
#   dplyr::mutate(size = "adult", 
#                 direction = "Increased +10%")
# 
# # decreased parameters 
# sa_kmm_decreased_5_sapling <- calc_kmm_sa_int(default = sa_default_sapling,
#                                               changed = sa_decreased_5_sapling,
#                                               r = r_kmm, 
#                                               correction = correction_kmm, 
#                                               window = window) %>% 
#   dplyr::mutate(size = "sapling", 
#                 direction = "Decreased -5%")
# 
# sa_kmm_decreased_5_adult <- calc_kmm_sa_int(default = sa_default_adult,
#                                             changed = sa_decreased_5_adult,
#                                             r = r_kmm, 
#                                             correction = correction_kmm, 
#                                             window = window) %>% 
#   dplyr::mutate(size = "adult", 
#                 direction = "Decreased -5%")
# 
# sa_kmm_decreased_10_sapling <- calc_kmm_sa_int(default = sa_default_sapling,
#                                                changed = sa_decreased_10_sapling,
#                                                r = r_kmm, 
#                                                correction = correction_kmm, 
#                                                window = window) %>% 
#   dplyr::mutate(size = "sapling", 
#                 direction = "Decreased -10%")
# 
# sa_kmm_decreased_10_adult <- calc_kmm_sa_int(default = sa_default_adult,
#                                              changed = sa_decreased_10_adult,
#                                              r = r_kmm, 
#                                              correction = correction_kmm, 
#                                              window = window) %>% 
#   dplyr::mutate(size = "adult", 
#                 direction = "Decreased -10%")
# 
# sa_kmm <- dplyr::bind_rows(sa_kmm_increased_5_sapling,
#                            sa_kmm_increased_5_adult,
#                            sa_kmm_increased_10_sapling,
#                            sa_kmm_increased_10_adult,
#                            sa_kmm_decreased_5_sapling,
#                            sa_kmm_decreased_5_adult,
#                            sa_kmm_decreased_10_sapling,
#                            sa_kmm_decreased_10_adult,) %>% 
#   dplyr::mutate(parameter = factor(parameter, 
#                                    levels = parameter_levels),
#                 size = factor(size, levels = c("sapling", "adult"), 
#                               labels = c("Sapling", "Adult")),
#                 direction = factor(direction, 
#                                    levels = c("Decreased -10%", 
#                                               "Decreased -5%",
#                                               "Increased +5%", 
#                                               "Increased +10%")))
# 
# ggplot_sa_kmm <- ggplot(data = sa_kmm) + 
#   geom_bar(aes(x = parameter, y = kmm_mean * 100, fill = direction),
#            col = "black", stat = "identity", position = position_dodge()) +
#   facet_wrap(~ size) + 
#   geom_hline(yintercept = -10, linetype = 2, col = "#0D0887FF") +
#   geom_hline(yintercept = -5, linetype = 2, col = "#9C179EFF") +
#   geom_hline(yintercept = 0, linetype = 1) +
#   geom_hline(yintercept = 5, linetype = 2, col = "#ED7953FF") +
#   geom_hline(yintercept = 10, linetype = 2, col = "#F0F921FF") +
#   coord_flip() +
#   scale_x_discrete(name = "Parameter") +
#   scale_y_continuous(name = expression(paste("Relative difference", integral(kmm(r), 0, r))),
#                      breaks = seq(-12.5, 12.5, 2.5)) +
#   scale_fill_manual(name = "Parameter change", 
#                     values = c("#0D0887FF", "#9C179EFF" ,
#                                "#ED7953FF", "#F0F921FF")) +
#   theme_classic(base_size = base_size) + 
#   theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = ggplot_sa_kmm_dec, 
#                      filename = "ggplot_sa_kmm_dec.png", 
#                      path = "Figures/SA", 
#                      dpi = dpi,
#                      width = width_full, height = height_small, units = units)
