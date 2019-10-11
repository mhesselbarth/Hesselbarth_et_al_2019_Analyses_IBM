###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Comparison spatial model spatial ####

#### Import libraries and data ####

# load packages #
library(rabmp)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(spatstat)
library(tidyverse)

source("Helper_functions/helper_functions_comparison_spatial.R")

pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds") %>% 
  spatstat::subset.ppp(species == "beech" & dbh_99 > 1)

pattern_2007 <- readr::read_rds("Data/Raw/pattern_2007_ppp.rds") %>% 
  spatstat::subset.ppp(species == "beech" & dbh_07 > 1)

pattern_2013 <- readr::read_rds("Data/Raw/pattern_2013_ppp.rds") %>% 
  spatstat::subset.ppp(species == "beech" & dbh_13 > 1)

model_run_y50_e5_r50 <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_real_b.rds")

names(model_run_y50_e5_r50) <- rep("Biotic model", 
                                    times = length(model_run_y50_e5_r50))

window <- readr::read_rds("Data/Raw/plot_area_owin.rds")

parameters_fitted_biotic <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt", 
                                                   sep = ";")

overwrite <- FALSE

base_size <- 10

#### Calculate Nearest-neighbor distribution function ####
# set pa.rameters #
r <- seq(from = 0, to = 10, length.out = 525)
correction <- "km"

# calculate NND #
nnd_model <- calc_nnd_comp(data = model_run_y50_e5_r50, 
                           window = window, r = r, correction = correction) %>% 
  dplyr::mutate(data_type = "Biotic model")

nnd_1999 <- pattern_1999 %>% 
  spatstat::Gest(r = r, correction = correction) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction) %>% 
  dplyr::mutate(data_type = "Observed data 1999")

nnd_2007 <- spatstat::subset.ppp(pattern_2007, 
                                 inside_fence == 0 & type != "dead") %>% 
  spatstat::Gest(r = r, correction = correction) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction) %>% 
  dplyr::mutate(data_type = "Observed data 2007")

nnd_2013 <- spatstat::subset.ppp(pattern_2013, 
                                 inside_fence == 0 & type != "dead") %>% 
  spatstat::Gest(r = r, correction = correction) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction) %>% 
  dplyr::mutate(data_type = "Observed data 2013")

# combine to one df #
nnd_overall <- dplyr::bind_rows(nnd_1999, 
                                nnd_2007, 
                                nnd_2013) %>% 
  dplyr::mutate(data_type = factor(data_type, 
                                   levels = c("Observed data 1999",
                                              "Observed data 2007",
                                              "Observed data 2013")))

# create plot #
ggplot_biotic_nnd <- ggplot(data = nnd_overall) + 
  geom_ribbon(data = nnd_model, aes(x = r, ymin = fun_lo, ymax = fun_hi), 
              fill = "grey") +
  geom_line(data = nnd_model, aes(x = r, y = fun_mean)) +
  geom_line(aes(x = r, y = km, col = data_type), size = 0.75) +
  scale_color_viridis_d(name = "") +
  labs(x = "r [m]", y = "G(r)") +
  theme_classic(base_size = 10) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_biotic_nnd, 
                     filename = "ggplot_biotic_nnd.png",
                     path = "Figures/Biotic_model/", 
                     dpi = 300, width = 15, height = 12.5, units = "cm", 
                     overwrite = overwrite)

#### Pair-correlation function ####
# set parameters #
r <- seq(from = 0, to = 20, length.out = 525)
correction <- "Ripley"

# calculate pcf #
pcf_model <- calc_pcf_comp(data = model_run_y50_e5_r50, 
                           window = window, r = r,
                           correction = correction) %>% 
  dplyr::mutate(data_type = "Biotic model")

pcf_1999 <- pattern_1999 %>% 
  spatstat::pcf(r = r, correction = correction, divisor = "d") %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type = "Observed data 1999")

pcf_2007 <- spatstat::subset.ppp(pattern_2007, 
                                 inside_fence == 0 & type != "dead") %>% 
  spatstat::pcf(r = r, correction = correction, divisor = "d") %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type = "Observed data 2007")

pcf_2013 <- spatstat::subset.ppp(pattern_2013, 
                                 inside_fence == 0 & type != "dead") %>% 
  spatstat::pcf(r = r, correction = correction, divisor = "d") %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type = "Observed data 2013")

# combine to one df #
pcf_overall <- dplyr::bind_rows(pcf_1999, 
                                pcf_2007, 
                                pcf_2013) %>% 
  dplyr::mutate(data_type = factor(data_type, 
                                   levels = c("Observed data 1999",
                                              "Observed data 2007",
                                              "Observed data 2013")))

# create plot #
ggplot_biotic_pcf <- ggplot(data = pcf_overall) + 
  geom_ribbon(data = pcf_model, aes(x = r, ymin = fun_lo, ymax = fun_hi),
              fill = "grey") +
  geom_line(data = pcf_model, aes(x = r, y = fun_mean)) +
  geom_line(aes(x = r, y = pcf, col = data_type), size = 0.75) +
  geom_hline(yintercept = 1, linetype = 2, size = 0.25) +
  scale_color_viridis_d(name = "") +
  labs(x = "r [m]", y = "g(r)") +
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_biotic_pcf, 
                     filename = "ggplot_biotic_pcf.png", 
                     path = "Figures/Biotic_model/", 
                     dpi = 300, width = 15, height = 12.5, units = "cm", 
                     overwrite = overwrite)

#### Mark correlation function ####
# set parameters #
r <- seq(from = 0, to = 20, length.out = 525)
correction <- "Ripley"

# calculate kmm #
kmm_model <- calc_kmm_comp(data = model_run_y50_e5_r50, 
                           window = window, r = r, 
                           correction = correction) %>% 
  dplyr::mutate(data_type = "Biotic model")

kmm_1999 <- spatstat::subset.ppp(pattern_1999, select = dbh_99) %>% 
  spatstat::markcorr(r = r, correction = correction) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::mutate(data_type = "Observed data 1999")

kmm_2007 <- spatstat::subset.ppp(pattern_2007, 
                                 inside_fence == 0 & type != "dead") %>%
  spatstat::subset.ppp(select = dbh_07) %>% 
  spatstat::markcorr(r = r, correction = correction) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::mutate(data_type = "Observed data 2007")

kmm_2013 <- spatstat::subset.ppp(pattern_2013, 
                                 inside_fence == 0 & type != "dead" & !is.na(dbh_13)) %>% 
  spatstat::subset.ppp(select = dbh_13) %>% 
  spatstat::markcorr(r = r, correction = correction) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::mutate(data_type = "Observed data 2013")

# combine to one df #
kmm_overall <- dplyr::bind_rows(kmm_1999, 
                                kmm_2007, 
                                kmm_2013) %>% 
  dplyr::mutate(data_type = factor(data_type, 
                                   levels = c("Observed data 1999",
                                              "Observed data 2007",
                                              "Observed data 2013")))

# create plot #
ggplot_biotic_kmm <- ggplot(data = kmm_overall) + 
  geom_ribbon(data = kmm_model, aes(x = r, ymin = fun_lo, ymax = fun_hi), 
              fill = "grey") +
  geom_line(data = kmm_model, aes(x = r, y = fun_mean)) +
  geom_line(aes(x = r, y = kmm, col = data_type), size = 0.75) +
  geom_hline(yintercept = 1, linetype = 2, size = 0.25) +
  scale_color_viridis_d(name = "") +
  labs(x = "r [m]", y = "kmm(r)") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_biotic_kmm, 
                     filename = "ggplot_biotic_kmm.png", 
                     path = "Figures/Biotic_model/", 
                     dpi = 300, width = 15, height = 12.5, units = "cm", 
                     overwrite = overwrite)

#### CI index ####
from <- 0.25
to <- 1

ci_model <- calc_ci_comp(data = model_run_y50_e5_r50, 
                         parameters = parameters_fitted_biotic, 
                         from = from, to = to) %>% 
  dplyr::mutate(data_type = "Biotic model")


ci_1999 <- spatstat::subset.ppp(pattern_1999, select = dbh_99) %>% 
  data.table::as.data.table() %>% 
  as.matrix() %>% 
  rabmp:::rcpp_calculate_ci(alpha = parameters_fitted_biotic$ci_alpha,
                            beta = parameters_fitted_biotic$ci_beta,
                            max_dist = parameters_fitted_biotic$ci_max_dist) %>% 
  density(from = from, to = to, n = 512)

ci_1999 <- tibble::tibble(x = ci_1999$x, y = ci_1999$y) %>% 
  dplyr::mutate(data_type = "Observed data 1999")
  

ci_2007 <- spatstat::subset.ppp(pattern_2007, 
                                inside_fence == 0 & type != "dead", 
                                select = dbh_07) %>%
  data.table::as.data.table() %>% 
  as.matrix() %>% 
  rabmp:::rcpp_calculate_ci(alpha = parameters_fitted_biotic$ci_alpha,
                            beta = parameters_fitted_biotic$ci_beta,
                            max_dist = parameters_fitted_biotic$ci_max_dist) %>% 
  density(from = from, to = to, n = 512)

ci_2007 <- tibble::tibble(x = ci_2007$x, y = ci_2007$y) %>% 
  dplyr::mutate(data_type = "Observed data 2007")


ci_2013 <- spatstat::subset.ppp(pattern_2013, 
                                inside_fence == 0 & type != "dead" & !is.na(dbh_13), 
                                select = dbh_13) %>% 
  data.table::as.data.table() %>% 
  as.matrix() %>% 
  rabmp:::rcpp_calculate_ci(alpha = parameters_fitted_biotic$ci_alpha,
                            beta = parameters_fitted_biotic$ci_beta,
                            max_dist = parameters_fitted_biotic$ci_max_dist) %>% 
  density(from = from, to = to, n = 512)

ci_2013 <- tibble::tibble(x = ci_2013$x, y = ci_2013$y) %>% 
  dplyr::mutate(data_type = "Observed data 2013")

# combine to one df #
ci_overall <- dplyr::bind_rows(ci_1999, 
                               ci_2007, 
                               ci_2013) %>% 
  dplyr::mutate(data_type = factor(data_type, 
                                   levels = c("Observed data 1999",
                                              "Observed data 2007",
                                              "Observed data 2013")))

ggplot_biotic_ci <- ggplot(data = ci_overall) + 
  geom_density(aes(x = x, y = y, fill = data_type),
               stat = "identity", col = "black", alpha = 0.5) + 
  geom_density(data = ci_model,
               aes(x = x, y = fun_mean, fill = data_type),
               stat = "identity", col = "black", alpha = 0.5) +
  scale_fill_viridis_d(name = "") +
  labs(x = "Competition value", y = "Density") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_biotic_ci, 
                     filename = "ggplot_biotic_ci.png", 
                     path = "Figures/Biotic_model/", 
                     dpi = 300, width = 15, height = 12.5, units = "cm", 
                     overwrite = overwrite)
