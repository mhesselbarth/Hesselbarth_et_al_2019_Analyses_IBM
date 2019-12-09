###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Sobol indice abiotic ####

#### load packages and data ####
source("Helper_functions/helper_functions_setup.R")
source("Helper_functions/helper_functions_sa_sobol.R")

parameters_fitted_abiotic <- rabmp::read_parameters("Data/Input/parameters_fitted_abiotic_real.txt",
                                                    sep = ";")

pattern_1999 <- readr::read_rds("Data/Input/beech_1999_ppp.rds")

abiotic_habitats_real <- readr::read_rds("Data/Input/abiotic_cond_real_model.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

years <- 50
save_each <- years

n <- 100 # 250

set.seed(42)

#### Pre-process data ####
pattern_1999_dt <- tibble::as_tibble(pattern_1999) %>%
  dplyr::select(x, y, species, dbh_99, type) %>% 
  dplyr::filter(species == "beech") %>%
  dplyr::mutate(type = "adult") %>% 
  dplyr::select(-species) %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh_99")

#### Number of individuals ####

set.seed(42)

# sample parameters #
param_set_1_indiv <- tgp::lhs(n = n, rect = matrix(data = c(-0.0824235, -0.0261765,  # growth_abiotic +- SE * 1.96
                                                            -2.4, -2.1, # mort_dbh_early_low
                                                            -2.1, -1.9, # mort_dbh_early_high
                                                            0.033 ,0.052, # mort_dbh_late_low
                                                            0.052, 0.07, # mort_dbh_late_high
                                                            -2.4, -1.4, # mort_dinc_low, 
                                                            -1.4, -0.45, # mort_dinc_high
                                                            1.2, 1.8, # mort_int_early_low
                                                            1.8, 2.5, # mort_int_early_high
                                                            -10, -8.9, # mort_int_late_low
                                                            -8.9, -7.8, # mort_int_late_high
                                                            0.003838327, 0.005457279, # seed_success_high
                                                            0.005457279, 0.009846355), # seed_success_low
                                                   ncol = 2, byrow = TRUE))

param_set_2_indiv <- tgp::lhs(n = n, rect = matrix(data = c(-0.0824235, -0.0261765,  # growth_abiotic +- SE * 1.96
                                                            -2.4, -2.1, # mort_dbh_early_low
                                                            -2.1, -1.9, # mort_dbh_early_high
                                                            0.033 ,0.052, # mort_dbh_late_low
                                                            0.052, 0.07, # mort_dbh_late_high
                                                            -2.4, -1.4, # mort_dinc_low, 
                                                            -1.4, -0.45, # mort_dinc_high
                                                            1.2, 1.8, # mort_int_early_low
                                                            1.8, 2.5, # mort_int_early_high
                                                            -10, -8.9, # mort_int_late_low
                                                            -8.9, -7.8, # mort_int_late_high
                                                            0.003838327, 0.005457279, # seed_success_high
                                                            0.005457279, 0.009846355), # seed_success_low
                                                   ncol = 2, byrow = TRUE))

# create an instance of the class sobol #
sobol_model_indiv <- sensitivity::sobol2007(model = NULL, 
                                            X1 =  data.frame(param_set_1_indiv), 
                                            X2 =  data.frame(param_set_2_indiv), 
                                            nboot = 10000) 

# get parameter combinitations from sobol model
param_sampled_indiv <- purrr::map(seq_len(nrow(sobol_model_indiv$X)), 
                                  function(x) as.numeric(sobol_model_indiv$X[x, ]))

# # get the simulated model response #
# simulation_results_indiv_abiotic <-
#   suppoRt::submit_to_cluster(calc_sobol_indiv_abiotic,
#                              x = param_sampled_indiv,
#                              const = list(data = pattern_1999_dt,
#                                           parameters = parameters_fitted_abiotic,
#                                           abiotic = abiotic_habitats_real$scaled,
#                                           plot_area = plot_area,
#                                           years = years,
#                                           save_each = save_each),
#                              n_jobs = length(param_sampled_indiv),
#                              log_worker = TRUE,
#                              template = list(job_name = "sobol_indiv",
#                                              walltime = "04:00:00",
#                                              queue = "medium",
#                                              mem_cpu = "4096",
#                                              log_file = "sobol_indiv.log"))
# 
# # flatten to vector #
# simulation_results_indiv_abiotic <- purrr::flatten_dbl(simulation_results_indiv_abiotic)
# 
# # save results #
# suppoRt::save_rds(object = simulation_results_indiv_abiotic,
#                   filename = "sa_simulation_results_indiv_abiotic.rds",
#                   path = "Data/Output/SA/",
#                   overwrite = overwrite)

simulation_results_indiv_abiotic <- readr::read_rds("Data/Output/SA/sa_simulation_results_indiv_abiotic.rds")
simulation_results_indiv_centered <- simulation_results_indiv_abiotic - mean(simulation_results_indiv_abiotic)

# add the simulation results to the sobol instance #   
sensitivity::tell(sobol_model_indiv, simulation_results_indiv_centered)

# convert to df # 
sobol_model_indiv_df <- tibble::as_tibble(sobol_model_indiv$S) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("growth_abiotic", 
                              "mort_dbh_early_low", "mort_dbh_early_high", 
                              "mort_dbh_late_low", "mort_dbh_late_high", 
                              "mort_dinc_low", "mort_dinc_high", 
                              "mort_int_early_low", "mort_int_early_high", 
                              "mort_int_late_low", "mort_int_late_high", 
                              "seed_success_high", "seed_success_low"),
                effect = "Main effect", 
                output = "Number of individuals")

sobol_model_indiv_df <- tibble::as_tibble(sobol_model_indiv$T) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("growth_abiotic", 
                              "mort_dbh_early_low", "mort_dbh_early_high", 
                              "mort_dbh_late_low", "mort_dbh_late_high", 
                              "mort_dinc_low", "mort_dinc_high", 
                              "mort_int_early_low", "mort_int_early_high", 
                              "mort_int_late_low", "mort_int_late_high", 
                              "seed_success_high", "seed_success_low"),
                effect = "Total effect", 
                output = "Number of individuals") %>% 
  dplyr::bind_rows(sobol_model_indiv_df, .) %>% 
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0, 
                                         value > 1 ~ 1, 
                                         value > 0 & value < 1 ~ value), 
                min_ci = dplyr::case_when(min_ci < 0 ~ 0, 
                                          min_ci > 0 ~ min_ci), 
                max_ci = dplyr::case_when(max_ci > 1 ~ 1, 
                                          max_ci < 1 ~ max_ci))

#### Pair-correlation function ####

set.seed(42)

# sample parameters #
param_set_1_pcf <- tgp::lhs(n = n, rect = matrix(data = c(-0.0824235, -0.0261765,  # growth_abiotic +- SE * 1.96
                                                          -2.4, -2.1, # mort_dbh_early_low
                                                          -2.1, -1.9, # mort_dbh_early_high
                                                          0.033 ,0.052, # mort_dbh_late_low
                                                          0.052, 0.07, # mort_dbh_late_high
                                                          -2.4, -1.4, # mort_dinc_low, 
                                                          -1.4, -0.45, # mort_dinc_high
                                                          1.2, 1.8, # mort_int_early_low
                                                          1.8, 2.5, # mort_int_early_high
                                                          -10, -8.9, # mort_int_late_low
                                                          -8.9, -7.8, # mort_int_late_high
                                                          0.003838327, 0.005457279, # seed_success_high
                                                          0.005457279, 0.009846355), # seed_success_low
                                                 ncol = 2, byrow = TRUE))

param_set_2_pcf <- tgp::lhs(n = n, rect = matrix(data = c(-0.0824235, -0.0261765,  # growth_abiotic +- SE * 1.96
                                                          -2.4, -2.1, # mort_dbh_early_low
                                                          -2.1, -1.9, # mort_dbh_early_high
                                                          0.033 ,0.052, # mort_dbh_late_low
                                                          0.052, 0.07, # mort_dbh_late_high
                                                          -2.4, -1.4, # mort_dinc_low, 
                                                          -1.4, -0.45, # mort_dinc_high
                                                          1.2, 1.8, # mort_int_early_low
                                                          1.8, 2.5, # mort_int_early_high
                                                          -10, -8.9, # mort_int_late_low
                                                          -8.9, -7.8, # mort_int_late_high
                                                          0.003838327, 0.005457279, # seed_success_high
                                                          0.005457279, 0.009846355), # seed_success_low
                                                 ncol = 2, byrow = TRUE))

# create an instance of the class sobol #
sobol_model_pcf <- sensitivity::sobol2007(model = NULL, 
                                          X1 =  data.frame(param_set_1_pcf), 
                                          X2 =  data.frame(param_set_1_pcf), 
                                          nboot = 10000) 

# get parameter combinitations from sobol model
param_sampled_pcf <- purrr::map(seq_len(nrow(sobol_model_pcf$X)), 
                                function(x) as.numeric(sobol_model_pcf$X[x, ]))

# # get the simulated model response #
# simulation_results_pcf_abiotic <-
#   suppoRt::submit_to_cluster(calc_sobol_pcf_abiotic,
#                              x = param_sampled_pcf,
#                              const = list(data = pattern_1999_dt,
#                                           parameters = parameters_fitted_abiotic,
#                                           abiotic = abiotic_habitats_real$scaled,
#                                           plot_area = plot_area,
#                                           years = years,
#                                           save_each = save_each),
#                              n_jobs = length(param_sampled_pcf),
#                              log_worker = TRUE,
#                              template = list(job_name = "sobol_pcf",
#                                              walltime = "24:00:00",
#                                              queue = "medium",
#                                              mem_cpu = "8192",
#                                              log_file = "sobol_pcf.log"))
# 
# # flatten to vector #
# simulation_results_pcf_abiotic <- purrr::flatten_dbl(simulation_results_pcf_abiotic)
# 
# # save results #
# suppoRt::save_rds(object = simulation_results_pcf_abiotic,
#                   filename = "sa_simulation_results_pcf_abiotic.rds",
#                   path = "Data/Output/SA/",
#                   overwrite = overwrite)

simulation_results_pcf_abiotic <- readr::read_rds("Data/Output/SA/sa_simulation_results_pcf_abiotic.rds")
simulation_results_pcf_centered <- simulation_results_pcf_abiotic - mean(simulation_results_pcf_abiotic)

# add the simulation results to the sobol instance #   
sensitivity::tell(sobol_model_pcf, simulation_results_pcf_centered)

# convert to df # 
sobol_model_pcf_df <- tibble::as_tibble(sobol_model_pcf$S) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("growth_abiotic", 
                              "mort_dbh_early_low", "mort_dbh_early_high", 
                              "mort_dbh_late_low", "mort_dbh_late_high", 
                              "mort_dinc_low", "mort_dinc_high", 
                              "mort_int_early_low", "mort_int_early_high", 
                              "mort_int_late_low", "mort_int_late_high", 
                              "seed_success_high", "seed_success_low"),
                effect = "Main effect", 
                output = "Integral pair-correlation function")

sobol_model_pcf_df <- tibble::as_tibble(sobol_model_pcf$T) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("growth_abiotic", 
                              "mort_dbh_early_low", "mort_dbh_early_high", 
                              "mort_dbh_late_low", "mort_dbh_late_high", 
                              "mort_dinc_low", "mort_dinc_high", 
                              "mort_int_early_low", "mort_int_early_high", 
                              "mort_int_late_low", "mort_int_late_high", 
                              "seed_success_high", "seed_success_low"),
                effect = "Total effect", 
                output = "Integral pair-correlation function") %>% 
  dplyr::bind_rows(sobol_model_pcf_df, .) %>% 
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0, 
                                         value > 1 ~ 1, 
                                         value > 0 & value < 1 ~ value), 
                min_ci = dplyr::case_when(min_ci < 0 ~ 0, 
                                          min_ci > 0 ~ min_ci), 
                max_ci = dplyr::case_when(max_ci > 1 ~ 1, 
                                          max_ci < 1 ~ max_ci))

#### Overall results ####
sobol_model_overall_df <- dplyr::bind_rows(sobol_model_indiv_df, 
                                           sobol_model_pcf_df) %>% 
  dplyr::mutate(effect = factor(effect, levels = c("Main effect", 
                                                   "Total effect")), 
                output = factor(output, levels = c("Number of individuals",
                                                   "Integral pair-correlation function")))

dplyr::group_by(sobol_model_overall_df, effect, output) %>% 
  dplyr::summarise(value = round(sum(value), digits = 2))

dplyr::filter(sobol_model_overall_df, output == "Number of individuals")
dplyr::filter(sobol_model_overall_df, output == "Integral pair-correlation function")


ggplot_sobol <- ggplot(data = sobol_model_overall_df) + 
  geom_point(aes(x = parameter, y = value, col = effect),
             size = 3.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x  = parameter, ymin = min_ci, ymax = max_ci,col = effect),
                width = 0.1, position = position_dodge(width = 0.5),
                size = 0.25) +
  facet_wrap(~ output, scales = "free_x") +
  # scale_color_viridis_d(name = "", option = "C") +
  scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF", 
                                           "Total effect" = "#CC4678FF")) +
  scale_y_continuous(name = "Effect strength", limits = c(0, 1)) +
  scale_x_discrete(name = "Parameter") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1))

suppoRt::save_ggplot(plot = ggplot_sobol, 
                     filename = "ggplot_sobol_abiotic.png", 
                     path = "Figures/", 
                     units = units, dpi = dpi, 
                     width = width_full, height = height_small, 
                     overwrite = overwrite)
