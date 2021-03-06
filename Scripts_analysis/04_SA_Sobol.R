###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Sobol indice ####

#### load packages and data ####
source("Helper_functions/helper_functions_setup.R")
source("Helper_functions/helper_functions_sa_sobol.R")

parameters_fitted_biotic <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt",
                                                   sep = ";")

pattern_1999 <- readr::read_rds("Data/Input/beech_1999_ppp.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

years <- 50
save_each <- years

n <- 250

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
param_set_1_indiv <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                            0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                            1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                            -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                            1.2, 2.5,  # mort_int_early CI from Holzwarth et al.
                                                            -10.0, -7.8), # mort_int_late
                                                   ncol = 2, byrow = TRUE))

param_set_2_indiv <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                            0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                            1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                            -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                            1.2, 2.5,  # mort_int_early CI from Holzwarth et al.
                                                            -10.0, -7.8), # mort_int_late
                                                   ncol = 2, byrow = TRUE))

# create an instance of the class sobol #
sobol_model_indiv_adult <- sensitivity::sobol2007(model = NULL, 
                                                  X1 = data.frame(param_set_1_indiv), 
                                                  X2 = data.frame(param_set_2_indiv), 
                                                  nboot = 10000) 

sobol_model_indiv_sapling <- sensitivity::sobol2007(model = NULL, 
                                                    X1 = data.frame(param_set_1_indiv), 
                                                    X2 = data.frame(param_set_2_indiv), 
                                                    nboot = 10000) 

# get parameter combinitations from sobol model
param_sampled_indiv <- purrr::map(seq_len(nrow(sobol_model_indiv_adult$X)), 
                                  function(x) as.numeric(sobol_model_indiv_adult$X[x, ]))

# # get the simulated model response #
# simulation_results_indiv <-
#   suppoRt::submit_to_cluster(calc_sobol_indiv,
#                              x = param_sampled_indiv,
#                              const = list(data = pattern_1999_dt,
#                                           parameters = parameters_fitted_biotic,
#                                           plot_area = plot_area,
#                                           years = years,
#                                           save_each = save_each),
#                              n_jobs = length(param_sampled_indiv),
#                              template = list(job_name = "sobol_indiv",
#                                              walltime = "06:00:00",
#                                              queue = "medium",
#                                              service = "normal",
#                                              mem_cpu = "8192",
#                                              log_file = "sobol_indiv.log"))
# 
# # save results #
# suppoRt::save_rds(object = simulation_results_indiv,
#                   filename = "sa_simulation_results_indiv.rds",
#                   path = "Data/Output/SA/",
#                   overwrite = overwrite)

# read data 
simulation_results_indiv <- readr::read_rds("Data/Output/SA/sa_simulation_results_indiv.rds")

# get only number of individuals #
simulation_results_indiv_adult <- purrr::map_dbl(simulation_results_indiv, 
                                                 function(x) x$adult)

simulation_results_indiv_sapling <- purrr::map_dbl(simulation_results_indiv, 
                                                 function(x) x$sapling)

# center results by mean #
simulation_results_indiv_adult_centered <- simulation_results_indiv_adult - 
  mean(simulation_results_indiv_adult)

simulation_results_indiv_sapling_centered <- simulation_results_indiv_sapling - 
  mean(simulation_results_indiv_sapling)

# add the simulation results to the sobol instance #   
sensitivity::tell(sobol_model_indiv_adult, simulation_results_indiv_adult_centered)

sensitivity::tell(sobol_model_indiv_sapling, simulation_results_indiv_sapling_centered)

# convert to df # 
sobol_model_indiv_adult_df <- tibble::as_tibble(sobol_model_indiv_adult$S) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early", "mort_int_late"), 
                effect = "Main effect", 
                output = "Number of individuals", 
                group = "Adults")

sobol_model_indiv_adult_df <- tibble::as_tibble(sobol_model_indiv_adult$T) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early", "mort_int_late"), 
                effect = "Total effect", 
                output = "Number of individuals", 
                group = "Adults") %>% 
  dplyr::bind_rows(sobol_model_indiv_adult_df, .) %>% 
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0, 
                                          value > 1 ~ 1, 
                                          value > 0 & value < 1 ~ value), 
                min_ci = dplyr::case_when(min_ci < 0 ~ 0, 
                                          min_ci > 0 ~ min_ci), 
                max_ci = dplyr::case_when(max_ci > 1 ~ 1, 
                                          max_ci < 1 ~ max_ci))

sobol_model_indiv_sapling_df <- tibble::as_tibble(sobol_model_indiv_sapling$S) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early", "mort_int_late"), 
                effect = "Main effect", 
                output = "Number of individuals", 
                group = "Saplings")

sobol_model_indiv_sapling_df <- tibble::as_tibble(sobol_model_indiv_sapling$T) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early", "mort_int_late"), 
                effect = "Total effect", 
                output = "Number of individuals", 
                group = "Saplings") %>% 
  dplyr::bind_rows(sobol_model_indiv_sapling_df, .) %>% 
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
# sample parameters #
param_set_1_pcf <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                          0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                          1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                          -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                          1.2, 2.5,  # mort_int_early CI from Holzwarth et al.
                                                          -10.0, -7.8), # mort_int_late
                                                   ncol = 2, byrow = TRUE))

param_set_2_pcf <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                          0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                          1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                          -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                          1.2, 2.5,  # mort_int_early CI from Holzwarth et al.
                                                          -10.0, -7.8), # mort_int_late
                                                   ncol = 2, byrow = TRUE))

# create an instance of the class sobol #
sobol_model_pcf_adult <- sensitivity::sobol2007(model = NULL, 
                                                X1 =  data.frame(param_set_1_pcf), 
                                                X2 =  data.frame(param_set_1_pcf), 
                                                nboot = 10000) 

sobol_model_pcf_sapling <- sensitivity::sobol2007(model = NULL, 
                                                  X1 =  data.frame(param_set_1_pcf), 
                                                  X2 =  data.frame(param_set_1_pcf), 
                                                  nboot = 10000) 

# get parameter combinitations from sobol model
param_sampled_pcf <- purrr::map(seq_len(nrow(sobol_model_pcf_adult$X)), 
                                function(x) as.numeric(sobol_model_pcf_adult$X[x, ]))

# # get the simulated model response #
# simulation_results_pcf <-
#   suppoRt::submit_to_cluster(calc_sobol_pcf,
#                              x = param_sampled_pcf,
#                              const = list(data = pattern_1999_dt,
#                                           parameters = parameters_fitted_biotic,
#                                           plot_area = plot_area,
#                                           years = years,
#                                           save_each = save_each),
#                              n_jobs = length(param_sampled_pcf),
#                              log_worker = TRUE,
#                              template = list(job_name = "sobol_pcf",
#                                              walltime = "48:00:00",
#                                              queue = "medium",
#                                              service = "normal",
#                                              mem_cpu = "24576",
#                                              log_file = "sobol_pcf.log"))
# 
# # save results #
# suppoRt::save_rds(object = simulation_results_pcf,
#                   filename = "sa_simulation_results_pcf.rds",
#                   path = "Data/Output/SA/",
#                   overwrite = overwrite)

# read data
simulation_results_pcf <- readr::read_rds("Data/Output/SA/sa_simulation_results_pcf.rds")

# get only number of individuals #
simulation_results_pcf_adult <- purrr::map_dbl(simulation_results_pcf, 
                                                 function(x) x$adult)

simulation_results_pcf_sapling <- purrr::map_dbl(simulation_results_pcf, 
                                                   function(x) x$sapling)

# center results by mean
simulation_results_pcf_adult_centered <- simulation_results_pcf_adult -
  mean(simulation_results_pcf_adult)

simulation_results_pcf_sapling_centered <- simulation_results_pcf_sapling - 
  mean(simulation_results_pcf_sapling)

# add the simulation results to the sobol instance #   
sensitivity::tell(sobol_model_pcf_adult, simulation_results_pcf_adult_centered)

sensitivity::tell(sobol_model_pcf_sapling, simulation_results_pcf_sapling_centered)

# convert to df # 
sobol_model_pcf_adult_df <- tibble::as_tibble(sobol_model_pcf_adult$S) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early", "mort_int_late"), 
                effect = "Main effect", 
                output = "Summarised pair-correlation function", 
                group = "Adults")

sobol_model_pcf_adult_df <- tibble::as_tibble(sobol_model_pcf_adult$T) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early", "mort_int_late"), 
                effect = "Total effect", 
                output = "Summarised pair-correlation function", 
                group = "Adults") %>% 
  dplyr::bind_rows(sobol_model_pcf_adult_df, .) %>% 
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0, 
                                         value > 1 ~ 1, 
                                         value > 0 & value < 1 ~ value), 
                min_ci = dplyr::case_when(min_ci < 0 ~ 0, 
                                          min_ci > 0 ~ min_ci), 
                max_ci = dplyr::case_when(max_ci > 1 ~ 1, 
                                          max_ci < 1 ~ max_ci))

# convert to df # 
sobol_model_pcf_sapling_df <- tibble::as_tibble(sobol_model_pcf_sapling$S) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early", "mort_int_late"), 
                effect = "Main effect", 
                output = "Summarised pair-correlation function", 
                group = "Saplings")

sobol_model_pcf_sapling_df <- tibble::as_tibble(sobol_model_pcf_sapling$T) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early", "mort_int_late"), 
                effect = "Total effect", 
                output = "Summarised pair-correlation function", 
                group = "Saplings") %>% 
  dplyr::bind_rows(sobol_model_pcf_sapling_df, .) %>% 
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0, 
                                         value > 1 ~ 1, 
                                         value > 0 & value < 1 ~ value), 
                min_ci = dplyr::case_when(min_ci < 0 ~ 0, 
                                          min_ci > 0 ~ min_ci), 
                max_ci = dplyr::case_when(max_ci > 1 ~ 1, 
                                          max_ci < 1 ~ max_ci))

#### Overall results ####
sobol_model_overall_df <- dplyr::bind_rows(sobol_model_indiv_adult_df, 
                                           sobol_model_indiv_sapling_df,
                                           sobol_model_pcf_adult_df, 
                                           sobol_model_pcf_sapling_df) %>% 
  dplyr::mutate(effect = factor(effect, levels = c("Main effect", 
                                                   "Total effect")), 
                output = factor(output, levels = c("Number of individuals",
                                                   "Summarised pair-correlation function")), 
                group = factor(group, levels = c("Saplings", "Adults")), 
                parameter = factor(parameter, levels = c("ci_alpha", "ci_beta", "growth_infl", 
                                                         "mort_dbh_early", "mort_int_early", "mort_int_late")))

dplyr::group_by(sobol_model_overall_df, effect, output, group) %>% 
  dplyr::summarise(value = round(sum(value), digits = 2))

dplyr::filter(sobol_model_overall_df, output == "Number of individuals")
dplyr::filter(sobol_model_overall_df, output == "Summarised pair-correlation function")

ggplot_sobol <- ggplot(data = sobol_model_overall_df) + 
  geom_point(aes(x = parameter, y = value, col = effect),
             size = 3.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x  = parameter, ymin = min_ci, ymax = max_ci,col = effect),
                width = 0.1, position = position_dodge(width = 0.5),
                size = 0.25) +
  facet_wrap(~ output + group, scales = "free_x") +
  scale_color_manual(name = "", values = c("Main effect" = "#0D0887FF",
                                           "Total effect" = "#ED7953FF")) +
  scale_y_continuous(name = "Effect strength", limits = c(0, 0.75)) +
  scale_x_discrete(name = "Parameter") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1))

suppoRt::save_ggplot(plot = ggplot_sobol, 
                     filename = "ggplot_sobol.png", 
                     path = "Figures/", 
                     units = units, dpi = dpi, 
                     width = width_full, height = height_full * 2/3, 
                     overwrite = overwrite)
