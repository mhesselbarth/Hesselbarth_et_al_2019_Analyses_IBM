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

window <- readr::read_rds("Data/Raw/plot_area_owin.rds")

model_run_reco <- readr::read_rds("Data/Output/model_run_y50_e10_r50_reco_a.rds")

names(model_run_reco) <- rep("Abiotic model", times = length(model_run_reco))

overwrite <- FALSE

#### Calculate Nearest-neighbor distribution function ####
# set parameters #
r <- seq(from = 0, to = 10, length.out = 525)
correction <- "km"

# calculate NND #
nnd_model <- calc_nnd(data = model_run_reco, 
                      window = window, r = r, correction = correction) %>% 
  dplyr::bind_rows() %>% 
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo),
                   km_mean = mean(km),
                   km_lo = min(km),
                   km_hi = max(km)) %>%
  dplyr::mutate(data_type = "Abiotic model")

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
ggplot_abiotic_nnd <- ggplot(data = nnd_overall) + 
  geom_ribbon(data = nnd_model, aes(x = r, ymin = km_lo, ymax = km_hi), 
              fill = "grey") +
  # geom_line(data = nnd_model, aes(x = r, y = km_mean)) +
  geom_line(aes(x = r, y = km, col = data_type), size = 0.75) +
  scale_color_viridis_d(name = "") +
  labs(x = "r [m]", y = "G(r)") +
  theme_classic(base_size = 10) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_abiotic_nnd, 
                     filename = "ggplot_abiotic_nnd.png", path = "Figures/", 
                     dpi = 300, width = 15, height = 12.5, units = "cm", 
                     overwrite = overwrite)

#### Pair-correlation function ####
# set parameters #
r <- seq(from = 0, to = 20, length.out = 525)
correction <- "Ripley"
fast <- FALSE
divisor <- "d"

# calculate pcf #
pcf_model <- calc_pcf(data = model_run_reco, 
                      window = window, r = r, fast = fast, 
                      correction = correction, divisor = divisor) %>% 
  dplyr::bind_rows() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   pcf_mean = mean(pcf), 
                   pcf_lo = min(pcf),
                   pcf_hi = max(pcf)) %>% 
  dplyr::mutate(data_type = "Abiotic model")

pcf_1999 <- pattern_1999 %>% 
  spatstat::pcf(r = r, correction = correction, divisor = divisor) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type = "Observed data 1999")

pcf_2007 <- spatstat::subset.ppp(pattern_2007, 
                                 inside_fence == 0 & type != "dead") %>% 
  spatstat::pcf(r = r, correction = correction, divisor = divisor) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type = "Observed data 2007")

pcf_2013 <- spatstat::subset.ppp(pattern_2013, 
                                 inside_fence == 0 & type != "dead") %>% 
  spatstat::pcf(r = r, correction = correction, divisor = divisor) %>% 
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
ggplot_abiotic_pcf <- ggplot(data = pcf_overall) + 
  geom_ribbon(data = pcf_model, aes(x = r, ymin = pcf_lo, ymax = pcf_hi),
              fill = "grey") +
  geom_line(aes(x = r, y = pcf, col = data_type), size = 0.75) +
  geom_hline(yintercept = 1, linetype = 2, size = 0.25) +
  scale_color_viridis_d(name = "") +
  labs(x = "r [m]", y = "g(r)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_abiotic_pcf, 
                     filename = "ggplot_abiotic_pcf.png", path = "Figures/", 
                     dpi = 300, width = 15, height = 12.5, units = "cm", 
                     overwrite = overwrite)

#### Mark correlation function ####
# set parameters #
r <- seq(from = 0, to = 20, length.out = 525)
correction <- "Ripley"
fast <- FALSE
divisor <- "d"

# calculate kmm #
kmm_model <- calc_kmm(data = model_run_reco, 
                      window = window, r = r, 
                      correction = correction) %>% 
  dplyr::bind_rows() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::group_by(r) %>% 
  dplyr::summarise(theo = mean(theo), 
                   kmm_mean = mean(kmm),
                   kmm_lo = min(kmm),
                   kmm_hi = max(kmm)) %>% 
  dplyr::mutate(data_type = "Abiotic model")

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
ggplot_abiotic_kmm <- ggplot(data = kmm_overall) + 
  geom_ribbon(data = kmm_model, aes(x = r, ymin = kmm_lo, ymax = kmm_hi), 
              fill = "grey") +
  geom_line(aes(x = r, y = kmm, col = data_type), size = 0.75) +
  geom_hline(yintercept = 1, linetype = 2, size = 0.25) +
  scale_color_viridis_d(name = "") +
  labs(x = "r [m]", y = "kmm(r)") +
  theme_classic(base_size = 10) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_abiotic_kmm, 
                     filename = "ggplot_abiotic_kmm.png", path = "Figures/", 
                     dpi = 300, width = 15, height = 12.5, units = "cm", 
                     overwrite = overwrite)
