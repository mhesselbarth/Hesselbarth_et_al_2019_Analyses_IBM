###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Comparison model spatial ####

#### Import libraries and data ####
source("Helper_functions/helper_functions_setup.R")
source("Helper_functions/helper_functions_comparison_spatial.R")

pattern_2007 <- readr::read_rds("Data/Raw/pattern_2007_ppp.rds")
pattern_2013 <- readr::read_rds("Data/Raw/pattern_2013_ppp.rds")

model_run_y50_e5_r50_biotic <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_real_b.rds")[1:10]
model_run_y50_e5_r50_abiotic <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_real_a.rds")[1:10]

window <- readr::read_rds("Data/Raw/plot_area_owin.rds")

parameters_fitted_biotic <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt",
                                                   sep = ";")

parameters_fitted_abiotic <- rabmp::read_parameters("Data/Input/parameters_fitted_abiotic.txt",
                                                    sep = ";")

#### Preprocess data ####

# set names # 

names(model_run_y50_e5_r50_biotic) <- rep("Biotic model", 
                                          times = length(model_run_y50_e5_r50_biotic))

names(model_run_y50_e5_r50_abiotic) <- rep("Abiotic model", 
                                           times = length(model_run_y50_e5_r50_abiotic))

# filter size classes # 
model_run_y50_e5_r50_biotic_sapling <- purrr::map(model_run_y50_e5_r50_biotic, 
                                                  function(x) 
                                                    dplyr::filter(x, type == "sapling"))

model_run_y50_e5_r50_biotic_adult <- purrr::map(model_run_y50_e5_r50_biotic, 
                                                function(x) 
                                                  dplyr::filter(x, type == "adult"))

model_run_y50_e5_r50_abiotic_sapling <- purrr::map(model_run_y50_e5_r50_abiotic, 
                                                   function(x) 
                                                     dplyr::filter(x, type == "sapling"))

model_run_y50_e5_r50_abiotic_adult <- purrr::map(model_run_y50_e5_r50_abiotic, 
                                                 function(x) 
                                                   dplyr::filter(x, type == "adult"))

pattern_2007_sapling <- spatstat::subset.ppp(pattern_2007, 
                                             dbh_07 > 1 & dbh_07 <= 10 & 
                                               inside_fence == 0)

pattern_2007_adult <- spatstat::subset.ppp(pattern_2007, dbh_07 > 10 & 
                                             inside_fence == 0)

pattern_2013_sapling <- spatstat::subset.ppp(pattern_2013, 
                                             dbh_13 > 1 & dbh_13 <= 10 & 
                                               inside_fence == 0)

pattern_2013_adult <- spatstat::subset.ppp(pattern_2013, dbh_13 > 10 & 
                                             inside_fence == 0)

#### Nearest-neighbor distribution function ####
r_nnd <- seq(from = 0, to = 10, length.out = 525)
correction_nnd <- "km"

# calculate NND #
nnd_model_biotic_sapling <- calc_nnd_comp(data = model_run_y50_e5_r50_biotic_sapling, 
                                          r = r_nnd, correction = correction_nnd, 
                                          window = window) %>% 
  dplyr::mutate(data_type_model = "Biotic model", 
                size_model = "Sapling")

nnd_model_biotic_adult <- calc_nnd_comp(data = model_run_y50_e5_r50_biotic_adult, 
                                        r = r_nnd, correction = correction_nnd, 
                                        window = window) %>% 
  dplyr::mutate(data_type_model = "Biotic model", 
                size_model = "Adult")

nnd_model_abiotic_sapling <- calc_nnd_comp(data = model_run_y50_e5_r50_abiotic_sapling, 
                                           r = r_nnd, correction = correction_nnd, 
                                           window = window) %>% 
  dplyr::mutate(data_type_model = "Abiotic model", 
                size_model = "Sapling")

nnd_model_abiotic_adult <- calc_nnd_comp(data = model_run_y50_e5_r50_abiotic_adult, 
                                         r = r_nnd, correction = correction_nnd, 
                                         window = window) %>% 
  dplyr::mutate(data_type_model = "Abiotic model", 
                size_model = "Adult")

nnd_2007_sapling <- spatstat::subset.ppp(pattern_2007_sapling, 
                                         inside_fence == 0 & type != "dead") %>% 
  spatstat::Gest(r = r_nnd, correction = correction_nnd) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction_nnd) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size_field = "Sapling")

nnd_2007_adult <- spatstat::subset.ppp(pattern_2007_adult, 
                                       inside_fence == 0 & type != "dead") %>% 
  spatstat::Gest(r = r_nnd, correction = correction_nnd) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction_nnd) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size_field = "Adult")

nnd_2013_sapling <- spatstat::subset.ppp(pattern_2013_sapling, 
                                         inside_fence == 0 & type != "dead") %>% 
  spatstat::Gest(r = r_nnd, correction = correction_nnd) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction_nnd) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size_field = "Sapling")

nnd_2013_adult <- spatstat::subset.ppp(pattern_2013_adult, 
                                       inside_fence == 0 & type != "dead") %>% 
  spatstat::Gest(r = r_nnd, correction = correction_nnd) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction_nnd) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size_field = "Adult")

# combine to one df #
nnd_overall_model <- dplyr::bind_rows(nnd_model_biotic_sapling,
                                      nnd_model_biotic_adult,
                                      nnd_model_abiotic_sapling,
                                      nnd_model_abiotic_adult) %>% 
  dplyr::mutate(data_type_model = factor(data_type_model, 
                                         levels = c("Biotic model",
                                                    "Abiotic model")), 
                size_model = factor(size_model, levels = c("Sapling", "Adult")))

nnd_overall_field <- dplyr::bind_rows(nnd_2007_sapling,
                                      nnd_2007_adult,
                                      nnd_2013_sapling, 
                                      nnd_2013_adult) %>% 
  dplyr::mutate(data_type_field = factor(data_type_field, 
                                         levels = c("Field data 2007",
                                                    "Field data 2013")), 
                size_field = factor(size_field, levels = c("Sapling", "Adult")))

# create plot #
ggplot_nnd <- ggplot(data = nnd_overall_model) + 
  geom_ribbon(aes(x = r, ymin = fun_lo, ymax = fun_hi, fill = size_model), 
              alpha = 0.5) +
  geom_line(data = nnd_overall_field, 
            aes(x = r, y = km, col = size_field, linetype = data_type_field)) +
  facet_wrap(~ data_type_model) + 
  scale_color_viridis_d(name = "") +
  scale_fill_viridis_d(name = "") +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  labs(x = "r [m]", y = "G(r)") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_nnd,
                     filename = "ggplot_nnd.png",
                     path = "Figures/",
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units,
                     overwrite = overwrite)

#### Pair-correlation function #### 
r_pcf <- seq(from = 0, to = 50, length.out = 525)
correction_pcf <- "Ripley"
stoyan_pcf <- 0.25
divisor_pcf <- "d"

# calculate pcf #
pcf_model_biotic_sapling <- calc_pcf_comp(data = model_run_y50_e5_r50_biotic_sapling,
                                          r = r_pcf, correction = correction_pcf,
                                          window = window, stoyan = stoyan_pcf, 
                                          divisor = divisor_pcf) %>% 
  dplyr::mutate(data_type_model = "Biotic model", 
                size_model = "Sapling")

pcf_model_biotic_adult <- calc_pcf_comp(data = model_run_y50_e5_r50_biotic_adult,
                                        r = r_pcf, correction = correction_pcf,
                                        window = window, stoyan = stoyan_pcf, 
                                        divisor = divisor_pcf) %>% 
  dplyr::mutate(data_type_model = "Biotic model", 
                size_model = "Adult")

pcf_model_abiotic_sapling <- calc_pcf_comp(data = model_run_y50_e5_r50_abiotic_sapling,
                                           r = r_pcf, correction = correction_pcf,
                                           window = window, stoyan = stoyan_pcf, 
                                           divisor = divisor_pcf) %>% 
  dplyr::mutate(data_type_model = "Abiotic model", 
                size_model = "Sapling")

pcf_model_abiotic_adult <- calc_pcf_comp(data = model_run_y50_e5_r50_abiotic_adult,
                                         r = r_pcf, correction = correction_pcf,
                                         window = window, stoyan = stoyan_pcf, 
                                         divisor = divisor_pcf) %>% 
  dplyr::mutate(data_type_model = "Abiotic model", 
                size_model = "Adult")

pcf_2007_sapling <- spatstat::pcf(pattern_2007_sapling, 
                                  r = r_pcf, correction = correction_pcf, 
                                  divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size_field = "Sapling")

pcf_2007_adult <- spatstat::pcf(pattern_2007_adult, 
                                r = r_pcf, correction = correction_pcf, 
                                divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size_field = "Adult")

pcf_2013_sapling <- spatstat::pcf(pattern_2013_sapling, 
                                  r = r_pcf, correction = correction_pcf, 
                                  divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size_field = "Sapling")

pcf_2013_adult <- spatstat::pcf(pattern_2013_adult, 
                                r = r_pcf, correction = correction_pcf, 
                                divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size_field = "Adult")

pcf_overall_model <- dplyr::bind_rows(pcf_model_biotic_sapling,
                                      pcf_model_biotic_adult,
                                      pcf_model_abiotic_sapling,
                                      pcf_model_abiotic_adult) %>% 
  dplyr::mutate(data_type_model = factor(data_type_model, 
                                         levels = c("Biotic model",
                                                    "Abiotic model")),
                size_model = factor(size_model, levels = c("Sapling", "Adult")))

pcf_overall_field <- dplyr::bind_rows(pcf_2007_sapling,
                                      pcf_2007_adult,
                                      pcf_2013_sapling,
                                      pcf_2013_adult) %>% 
  dplyr::mutate(data_type_field = factor(data_type_field, 
                                         levels = c("Field data 2007",
                                                    "Field data 2013")),
                size_field = factor(size_field, levels = c("Sapling", "Adult")))

# create plot #
ggplot_pcf <- ggplot(data = pcf_overall_model) + 
  geom_ribbon(aes(x = r, ymin = fun_lo, ymax = fun_hi, fill = size_model),
              alpha = 0.5) +
  geom_line(data = pcf_overall_field,
            aes(x = r, y = pcf, col = size_field, linetype = data_type_field)) +
  geom_hline(yintercept = 1, linetype = 3, size = 0.25) +
  facet_wrap(~ data_type_model, scales = "free_y") +
  scale_color_viridis_d(name = "") +
  scale_fill_viridis_d(name = "") +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  labs(x = "r [m]", y = "g(r)") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_pcf,
                     filename = "ggplot_pcf.png",
                     path = "Figures/",
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units,
                     overwrite = overwrite)

#### Mark-correlation function #### 
r_kmm <- seq(from = 0, to = 50, length.out = 525)
correction_kmm <- "Ripley"

# calculate kmm #
kmm_model_biotic_sapling <- calc_kmm_comp(data = model_run_y50_e5_r50_biotic_sapling,
                                          r = r_kmm, correction = correction_kmm,
                                          window = window, stoyan = stoyan_kmm) %>% 
  dplyr::mutate(data_type_model = "Biotic model", 
                size_model = "Sapling")

kmm_model_biotic_adult <- calc_kmm_comp(data = model_run_y50_e5_r50_biotic_adult,
                                        r = r_kmm, correction = correction_kmm,
                                        window = window) %>% 
  dplyr::mutate(data_type_model = "Biotic model", 
                size_model = "Adult")

kmm_model_abiotic_sapling <- calc_kmm_comp(data = model_run_y50_e5_r50_abiotic_sapling,
                                           r = r_kmm, correction = correction_kmm,
                                           window = window) %>% 
  dplyr::mutate(data_type_model = "Abiotic model", 
                size_model = "Sapling")

kmm_model_abiotic_adult <- calc_kmm_comp(data = model_run_y50_e5_r50_abiotic_adult,
                                         r = r_kmm, correction = correction_kmm,
                                         window = window) %>% 
  dplyr::mutate(data_type_model = "Abiotic model", 
                size_model = "Adult")

kmm_2007_sapling <- spatstat::markcorr(spatstat::subset.ppp(pattern_2007_sapling, 
                                                            select = dbh_07), 
                                       r = r_kmm, correction = correction_kmm) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size_field = "Sapling")

kmm_2007_adult <- spatstat::markcorr(spatstat::subset.ppp(pattern_2007_adult, 
                                                          select = dbh_07), 
                                     r = r_kmm, correction = correction_kmm) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size_field = "Adult")

kmm_2013_sapling <- spatstat::markcorr(spatstat::subset.ppp(pattern_2013_sapling, 
                                                            select = dbh_13), 
                                       r = r_kmm, correction = correction_kmm) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size_field = "Sapling")

kmm_2013_adult <- spatstat::markcorr(spatstat::subset.ppp(pattern_2013_adult, 
                                                          select = dbh_13), 
                                     r = r_kmm, correction = correction_kmm) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size_field = "Adult")

kmm_overall_model <- dplyr::bind_rows(kmm_model_biotic_sapling,
                                      kmm_model_biotic_adult,
                                      kmm_model_abiotic_sapling,
                                      kmm_model_abiotic_adult) %>% 
  dplyr::mutate(data_type_model = factor(data_type_model, 
                                         levels = c("Biotic model",
                                                    "Abiotic model")),
                size_model = factor(size_model, levels = c("Sapling", "Adult")))

kmm_overall_field <- dplyr::bind_rows(kmm_2007_sapling,
                                      kmm_2007_adult,
                                      kmm_2013_sapling,
                                      kmm_2013_adult) %>% 
  dplyr::mutate(data_type_field = factor(data_type_field, 
                                         levels = c("Field data 2007",
                                                    "Field data 2013")),
                size_field = factor(size_field, levels = c("Sapling", "Adult")))

# create plot #
ggplot_kmm <- ggplot(data = kmm_overall_model) + 
  geom_ribbon(aes(x = r, ymin = fun_lo, ymax = fun_hi, fill = size_model),
              alpha = 0.5) +
  geom_line(data = kmm_overall_field,
            aes(x = r, y = kmm, col = size_field, linetype = data_type_field)) +
  geom_hline(yintercept = 1, linetype = 3, size = 0.25) +
  facet_wrap(~ data_type_model, scales = "free_y") + 
  scale_color_viridis_d(name = "") +
  scale_fill_viridis_d(name = "") +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  coord_cartesian(ylim = c(0.5, 1.75)) +
  labs(x = "r [m]", y = "kmm(r)") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_kmm,
                     filename = "ggplot_kmm.png",
                     path = "Figures/",
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units,
                     overwrite = overwrite)


#### CI index ####
from_ci <- 0.25
to_ci <- 1

ci_model_biotic <- calc_ci_comp(data = model_run_y50_e5_r50_biotic, 
                                parameters = parameters_fitted_biotic, 
                                from = from_ci, to = to_ci) %>% 
  dplyr::mutate(data_type_model = "Biotic model")

ci_model_abiotic <- calc_ci_comp(data = model_run_y50_e5_r50_abiotic, 
                                 parameters = parameters_fitted_abiotic, 
                                 from = from_ci, to = to_ci) %>% 
  dplyr::mutate(data_type_model = "Abiotic model")

ci_2007 <- spatstat::subset.ppp(pattern_2007, 
                                inside_fence == 0 & type != "dead", 
                                select = dbh_07) %>%
  data.table::as.data.table() %>% 
  as.matrix() %>% 
  rabmp:::rcpp_calculate_ci(alpha = parameters_fitted_biotic$ci_alpha,
                            beta = parameters_fitted_biotic$ci_beta,
                            max_dist = parameters_fitted_biotic$ci_max_dist) %>% 
  density(from = from_ci, to = to_ci, n = 512)

ci_2007 <- tibble::tibble(x = ci_2007$x, y = ci_2007$y) %>% 
  dplyr::mutate(data_type_field = "Field data 2007")

ci_2013 <- spatstat::subset.ppp(pattern_2013, 
                                inside_fence == 0 & type != "dead" & !is.na(dbh_13), 
                                select = dbh_13) %>% 
  data.table::as.data.table() %>% 
  as.matrix() %>% 
  rabmp:::rcpp_calculate_ci(alpha = parameters_fitted_biotic$ci_alpha,
                            beta = parameters_fitted_biotic$ci_beta,
                            max_dist = parameters_fitted_biotic$ci_max_dist) %>% 
  density(from = from_ci, to = to_ci, n = 512)

ci_2013 <- tibble::tibble(x = ci_2013$x, y = ci_2013$y) %>% 
  dplyr::mutate(data_type_field = "Field data 2013")

# combine to one df #
ci_overall_model <- dplyr::bind_rows(ci_model_biotic,
                                     ci_model_abiotic) %>% 
  dplyr::mutate(data_type_model = factor(data_type_model, 
                                         levels = c("Biotic model",
                                                    "Abiotic model")))

ci_overall_field <- dplyr::bind_rows(ci_2007, 
                                     ci_2013) %>% 
  dplyr::mutate(data_type_field = factor(data_type_field, 
                                         levels = c("Field data 2007",
                                                    "Field data 2013")))

ggplot_ci <- ggplot(data = ci_overall_model) + 
  geom_density(data = ci_overall_field,
               aes(x = x, y = y, fill = data_type_field),
               stat = "identity", alpha = 0.5, col = NA) +
  geom_density(aes(x = x, y = fun_mean, col = "Model"),
               stat = "identity", alpha = 0.5) +
  scale_fill_viridis_d(name = "") +
  scale_color_manual(name = "", values = "black") +
  facet_wrap(~ data_type_model) + 
  labs(x = "Competition value", y = "Density") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_ci, 
                     filename = "ggplot_ci.png", 
                     path = "Figures/", 
                     dpi = 300, 
                     width = width_full, height = height_small, units = units, 
                     overwrite = overwrite)