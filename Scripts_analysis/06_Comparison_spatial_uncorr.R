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

model_run_y50_e5_r50_corr <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_real_a.rds")
model_run_y50_e5_r50_uncorr <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_uncorr_a.rds")

window <- readr::read_rds("Data/Raw/plot_area_owin.rds")

parameters_fitted_abiotic_real <- rabmp::read_parameters("Data/Input/parameters_fitted_abiotic_real.txt",
                                                         sep = ";")

parameters_fitted_abiotic_reco <- rabmp::read_parameters("Data/Input/parameters_fitted_abiotic_reco.txt",
                                                         sep = ";")

#### Preprocess data ####

# set names # 

names(model_run_y50_e5_r50_corr) <- rep("Correlated data", 
                                        times = length(model_run_y50_e5_r50_corr))

names(model_run_y50_e5_r50_uncorr) <- rep("Uncorrelated data", 
                                          times = length(model_run_y50_e5_r50_uncorr))

# filter size classes # 
model_run_y50_e5_r50_corr_sapling <- purrr::map(model_run_y50_e5_r50_corr, 
                                                function(x) 
                                                  dplyr::filter(x, type == "sapling"))

model_run_y50_e5_r50_corr_adult <- purrr::map(model_run_y50_e5_r50_corr, 
                                                function(x) 
                                                  dplyr::filter(x, type == "adult"))

model_run_y50_e5_r50_uncorr_sapling <- purrr::map(model_run_y50_e5_r50_uncorr, 
                                                  function(x) 
                                                    dplyr::filter(x, type == "sapling"))

model_run_y50_e5_r50_uncorr_adult <- purrr::map(model_run_y50_e5_r50_uncorr, 
                                                function(x) 
                                                  dplyr::filter(x, type == "adult"))

#### Pair-correlation function #### 
r_pcf <- seq(from = 0, to = 50, length.out = 525)
correction_pcf <- "Ripley"
stoyan_pcf <- 0.25
divisor_pcf <- "d"

# calculate pcf #
pcf_model_corr_sapling <- calc_pcf_comp(data = model_run_y50_e5_r50_corr_sapling,
                                        sim_i = 50,
                                        r = r_pcf, correction = correction_pcf,
                                        window = window, stoyan = stoyan_pcf, 
                                        divisor = divisor_pcf) %>% 
  dplyr::mutate(data_type_model = "Correlated data", 
                size_model = "Sapling")

pcf_model_corr_adult <- calc_pcf_comp(data = model_run_y50_e5_r50_corr_adult,
                                      sim_i = 50,
                                      r = r_pcf, correction = correction_pcf,
                                      window = window, stoyan = stoyan_pcf, 
                                      divisor = divisor_pcf) %>% 
  dplyr::mutate(data_type_model = "Correlated data", 
                size_model = "Adult")

pcf_model_uncorr_sapling <- calc_pcf_comp(data = model_run_y50_e5_r50_uncorr_sapling,
                                          sim_i = 50,
                                          r = r_pcf, correction = correction_pcf,
                                          window = window, stoyan = stoyan_pcf, 
                                          divisor = divisor_pcf) %>% 
  dplyr::mutate(data_type_model = "Uncorrelated data", 
                size_model = "Sapling")

pcf_model_uncorr_adult <- calc_pcf_comp(data = model_run_y50_e5_r50_uncorr_adult,
                                      sim_i = 50,
                                      r = r_pcf, correction = correction_pcf,
                                      window = window, stoyan = stoyan_pcf, 
                                      divisor = divisor_pcf) %>% 
  dplyr::mutate(data_type_model = "Uncorrelated data", 
                size_model = "Adult")

pcf_overall <- dplyr::bind_rows(pcf_model_corr_sapling,
                                pcf_model_corr_adult,
                                pcf_model_uncorr_sapling,
                                pcf_model_uncorr_adult) %>% 
  dplyr::mutate(data_type_model = factor(data_type_model, 
                                         levels = c("Correlated data",
                                                    "Uncorrelated data")),
                size_model = factor(size_model, levels = c("Sapling", "Adult")))

# create plot #
ggplot_pcf <- ggplot(data = pcf_overall) + 
  geom_ribbon(aes(x = r, ymin = fun_lo, ymax = fun_hi, fill = size_model),
              alpha = 0.5) +
  facet_wrap(~ data_type_model, scales = "free_y") +
  scale_fill_manual(name = "", values = c("Sapling" = "#0D0887FF", 
                                          "Adult" = "#CC4678FF")) +
  scale_color_manual(name = "", values = c("Sapling" = "#0D0887FF", 
                                           "Adult" = "#CC4678FF")) +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  labs(x = "r [m]", y = expression(italic(g(r)))) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_pcf,
                     filename = "ggplot_pcf_uncorr.png",
                     path = "Figures/Appendix/",
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units,
                     overwrite = overwrite)

#### Mark-correlation function #### 
r_kmm <- seq(from = 0, to = 50, length.out = 525)
correction_kmm <- "Ripley"

kmm_model_corr_sapling <- calc_kmm_comp(data = model_run_y50_e5_r50_corr_sapling,
                                        sim_i = 50,
                                        r = r_kmm, correction = correction_kmm,
                                        window = window) %>% 
  dplyr::mutate(data_type_model = "Correlated data", 
                size_model = "Sapling")

kmm_model_corr_adult <- calc_kmm_comp(data = model_run_y50_e5_r50_corr_adult,
                                      sim_i = 50,
                                      r = r_kmm, correction = correction_kmm,
                                      window = window) %>% 
  dplyr::mutate(data_type_model = "Correlated data", 
                size_model = "Adult")

kmm_model_uncorr_sapling <- calc_kmm_comp(data = model_run_y50_e5_r50_uncorr_sapling,
                                          sim_i = 50,
                                          r = r_kmm, correction = correction_kmm,
                                          window = window) %>% 
  dplyr::mutate(data_type_model = "Uncorrelated data", 
                size_model = "Sapling")

kmm_model_uncorr_adult <- calc_kmm_comp(data = model_run_y50_e5_r50_uncorr_adult,
                                        sim_i = 50,
                                        r = r_kmm, correction = correction_kmm,
                                        window = window) %>% 
  dplyr::mutate(data_type_model = "Uncorrelated data", 
                size_model = "Adult")

kmm_overall <- dplyr::bind_rows(kmm_model_corr_sapling,
                                kmm_model_corr_adult,
                                kmm_model_uncorr_sapling,
                                kmm_model_uncorr_adult) %>% 
  dplyr::mutate(data_type_model = factor(data_type_model, 
                                         levels = c("Correlated data",
                                                    "Uncorrelated data")),
                size_model = factor(size_model, levels = c("Sapling", "Adult")))

# create plot #
ggplot_kmm <- ggplot(data = kmm_overall) + 
  geom_ribbon(aes(x = r, ymin = fun_lo, ymax = fun_hi, fill = size_model),
              alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = 3, size = 0.25) +
  facet_wrap(~ data_type_model, scales = "free_y") + 
  scale_fill_manual(name = "", values = c("Sapling" = "#0D0887FF", 
                                          "Adult" = "#CC4678FF")) +
  scale_color_manual(name = "", values = c("Sapling" = "#0D0887FF", 
                                           "Adult" = "#CC4678FF")) +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  coord_cartesian(ylim = c(0.5, 1.75)) +
  labs(x = "r [m]", y = expression(italic(k[mm](r)))) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_kmm,
                     filename = "ggplot_kmm_uncorr.png",
                     path = "Figures/Appendix/",
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units,
                     overwrite = overwrite)

#### CI index ####
from_ci <- 0.25
to_ci <- 1

ci_model_corr_sapling <- calc_ci_comp(data = model_run_y50_e5_r50_corr_sapling,
                                      sim_i = 50,
                                      parameters = parameters_fitted_biotic_real, 
                                      from = from_ci, to = to_ci) %>% 
  dplyr::mutate(data_type_model = "Correlated data", 
                size_model = "Sapling")

ci_model_corr_adult <- calc_ci_comp(data = model_run_y50_e5_r50_corr_adult,
                                    sim_i = 50,
                                    parameters = parameters_fitted_biotic_real, 
                                    from = from_ci, to = to_ci) %>% 
  dplyr::mutate(data_type_model = "Correlated data", 
                size_model = "Adult")

ci_model_uncorr_sapling <- calc_ci_comp(data = model_run_y50_e5_r50_uncorr_sapling,
                                        sim_i = 50,
                                        parameters = parameters_fitted_biotic_reco, 
                                        from = from_ci, to = to_ci) %>% 
  dplyr::mutate(data_type_model = "Uncorrelated data", 
                size_model = "Sapling")

ci_model_uncorr_adult <- calc_ci_comp(data = model_run_y50_e5_r50_uncorr_adult,
                                      sim_i = 50,
                                      parameters = parameters_fitted_biotic_reco, 
                                      from = from_ci, to = to_ci) %>% 
  dplyr::mutate(data_type_model = "Uncorrelated data", 
                size_model = "Adult")

# combine to one df #
ci_overall <- dplyr::bind_rows(ci_model_corr_sapling,
                               ci_model_corr_adult,
                               ci_model_uncorr_sapling,
                               ci_model_uncorr_adult) %>% 
  dplyr::mutate(data_type_model = factor(data_type_model, 
                                         levels = c("Correlated datal",
                                                    "Uncorrelated data")), 
                size_model = factor(size_model, levels = c("Sapling", "Adult")))

ggplot_ci <- ggplot(data = ci_overall) + 
  geom_density(aes(x = x, y = fun_mean, fill = size_model),
               col = NA, stat = "identity", alpha = 0.5) +
  scale_fill_manual(name = "", values = c("Sapling" = "#0D0887FF", 
                                          "Adult" = "#CC4678FF")) +
  scale_color_manual(name = "", values = c("Sapling" = "#0D0887FF", 
                                           "Adult" = "#CC4678FF")) +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  facet_wrap(~ data_type_model) + 
  labs(x = "Competition value", y = "Density") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_ci, 
                     filename = "ggplot_ci_uncorr.png", 
                     path = "Figures/Appendix/", 
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units, 
                     overwrite = overwrite)
