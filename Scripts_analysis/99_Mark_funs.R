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

beech_2007_ppp <- readr::read_rds("Data/Input/beech_2007_ppp.rds")

beech_2013_ppp <- readr::read_rds("Data/Input/beech_2013_ppp.rds")

model_run_y50_e5_r50_biotic <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_real_b.rds")
model_run_y50_e5_r50_abiotic <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_real_a.rds")

window <- readr::read_rds("Data/Raw/plot_area_owin.rds")

parameters_fitted_biotic <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt",
                                                   sep = ";")

parameters_fitted_abiotic <- rabmp::read_parameters("Data/Input/parameters_fitted_abiotic_real.txt",
                                                    sep = ";")

#### Preprocess data ####
# set names # 

names(model_run_y50_e5_r50_biotic) <- rep("Biotic model version", 
                                          times = length(model_run_y50_e5_r50_biotic))

names(model_run_y50_e5_r50_abiotic) <- rep("Combined model version", 
                                           times = length(model_run_y50_e5_r50_abiotic))

sim_i <- 50

#### Mark-correlation function ####

r_kmm <- seq(from = 0, to = 50, length.out = 525)
correction_kmm <- "Ripley"

# calculate kmm #
kmm_model_biotic <- calc_kmm_comp(data = model_run_y50_e5_r50_biotic,
                                          sim_i = sim_i,
                                          r = r_kmm, correction = correction_kmm,
                                          window = window) %>% 
  dplyr::mutate(data_type_model = "Biotic model version")

kmm_model_abiotic <- calc_kmm_comp(data = model_run_y50_e5_r50_abiotic,
                                           sim_i = sim_i,
                                           r = r_kmm, correction = correction_kmm,
                                           window = window) %>% 
  dplyr::mutate(data_type_model = "Combined model version")


kmm_2007 <- spatstat::markcorr(spatstat::subset.ppp(beech_2007_ppp, 
                                                            select = dbh_07), 
                                       r = r_kmm, correction = correction_kmm) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007")

kmm_2013 <- spatstat::markcorr(spatstat::subset.ppp(beech_2013_ppp, 
                                                            select = dbh_13), 
                                       r = r_kmm, correction = correction_kmm) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "kmm")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013")

kmm_overall_model <- dplyr::bind_rows(kmm_model_biotic,
                                      kmm_model_abiotic) %>% 
  dplyr::mutate(data_type_model = factor(data_type_model, 
                                         levels = c("Biotic model version",
                                                    "Combined model version")))

kmm_overall_field <- dplyr::bind_rows(kmm_2007,
                                      kmm_2013) %>% 
  dplyr::mutate(data_type_field = factor(data_type_field, 
                                         levels = c("Field data 2007",
                                                    "Field data 2013")))

# create plot #
ggplot_kmm <- ggplot(data = kmm_overall_model) + 
  geom_ribbon(aes(x = r, ymin = fun_lo, ymax = fun_hi, fill = "Simulation model"),
              alpha = 0.5) +
  geom_line(data = kmm_overall_field,
            aes(x = r, y = kmm, linetype = data_type_field)) +
  geom_hline(yintercept = 1, linetype = 3, size = 0.25) +
  facet_wrap(~ data_type_model, scales = "free_y") + 
  scale_fill_manual(name = "", values = "grey25") +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  coord_cartesian(ylim = c(0.5, 1.75)) +
  labs(x = "r [m]", y = expression(italic(k[mm](r)))) +
  guides(colour = FALSE) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))




#### Mark-correlation function #### 
r_gmm <- seq(from = 0, to = 50, length.out = 525)
correction_gmm <- "Ripley"

# calculate gmm #
gmm_model_biotic <- calc_vario_comp(data = model_run_y50_e5_r50_biotic,
                                    sim_i = sim_i,
                                    r = r_gmm, correction = correction_gmm,
                                    window = window) %>% 
  dplyr::mutate(data_type_model = "Biotic model version")

gmm_model_abiotic <- calc_vario_comp(data = model_run_y50_e5_r50_abiotic,
                                     sim_i = sim_i,
                                     r = r_gmm, correction = correction_gmm,
                                     window = window) %>% 
  dplyr::mutate(data_type_model = "Combined model version")

gmm_2007 <- spatstat::markvario(spatstat::subset.ppp(beech_2007_ppp, 
                                                     select = dbh_07), 
                                r = r_gmm, correction = correction_gmm, 
                                normalise = TRUE) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "gmm")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007")

gmm_2013 <- spatstat::markvario(spatstat::subset.ppp(beech_2013_ppp, 
                                                     select = dbh_13), 
                                r = r_gmm, correction = correction_gmm, 
                                normalise = TRUE) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "gmm")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013")

gmm_overall_model <- dplyr::bind_rows(gmm_model_biotic,
                                      gmm_model_abiotic) %>% 
  dplyr::mutate(data_type_model = factor(data_type_model, 
                                         levels = c("Biotic model version",
                                                    "Combined model version")))

gmm_overall_field <- dplyr::bind_rows(gmm_2007,
                                      gmm_2013) %>% 
  dplyr::mutate(data_type_field = factor(data_type_field, 
                                         levels = c("Field data 2007",
                                                    "Field data 2013")))

# create plot #
ggplot_gmm <- ggplot(data = gmm_overall_model) + 
  geom_ribbon(aes(x = r, ymin = fun_lo, ymax = fun_hi, fill = "Simulation model"),
              alpha = 0.5) +
  geom_line(data = gmm_overall_field,
            aes(x = r, y = gmm, linetype = data_type_field)) +
  geom_hline(yintercept = 1, linetype = 3, size = 0.25) +
  facet_wrap(~ data_type_model, scales = "free_y") + 
  scale_fill_manual(name = "", values = "grey25") +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  # coord_cartesian(ylim = c(0.5, 1.75)) +
  labs(x = "r [m]", y = expression(italic(gamma[mm](r)))) +
  guides(colour = FALSE) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_gmm,
                     filename = "ggplot_gmm_y50.png",
                     path = "Figures/",
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units,
                     overwrite = overwrite)

