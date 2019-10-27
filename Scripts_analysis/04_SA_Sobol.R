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

n <- 250 # 100

set.seed(42)

#### Pre-process data ####
pattern_1999_dt <- tibble::as_tibble(pattern_1999) %>%
  dplyr::select(x, y, species, dbh_99, type) %>% 
  dplyr::filter(species == "beech") %>%
  dplyr::mutate(type = "adult") %>% 
  dplyr::select(-species) %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh_99")

#### Number of individuals ####

# sample parameters #
param_set_1_indiv <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                            0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                            1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                            -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                            1.2, 2.5),  # mort_int_early CI from Holzwarth et al.
                                                   ncol = 2, byrow = TRUE))

param_set_2_indiv <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                            0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                            1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                            -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                            1.2, 2.5), # mort_int_early CI from Holzwarth et al.
                                                   ncol = 2, byrow = TRUE))

# create an instance of the class sobol #
sobol_model_indiv <- sensitivity::sobol2007(model = NULL, 
                                            X1 =  data.frame(param_set_1_indiv), 
                                            X2 =  data.frame(param_set_2_indiv), 
                                            conf = 0.95,
                                            nboot = 10000) 

# get parameter combinitations from sobol model
param_sampled_indiv <- purrr::map(seq_len(nrow(sobol_model_indiv$X)), 
                                  function(x) as.numeric(sobol_model_indiv$X[x, ]))

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
#                                              walltime = "03:00:00",
#                                              queue = "medium", 
#                                              service = "normal",
#                                              mem_cpu = "2048", 
#                                              log_file = "sobol_indiv.log"))
# 
# # flatten to vector #
# simulation_results_indiv <- purrr::flatten_dbl(simulation_results_indiv)
# 
# # save results #
# suppoRt::save_rds(object = simulation_results_indiv,
#                   filename = "sa_simulation_results_indiv.rds",
#                   path = "Data/Output/SA/",
#                   overwrite = overwrite)

simulation_results_indiv <- readr::read_rds("Data/Output/SA/sa_simulation_results_indiv.rds")

# add the simulation results to the sobol instance #   
sensitivity::tell(sobol_model_indiv, simulation_results_indiv)

# convert to df # 
sobol_model_indiv_df <- tibble::as_tibble(sobol_model_indiv$S) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early"), 
                effect = "First order indices", 
                output = "Number of individuals")

sobol_model_indiv_df <- tibble::as_tibble(sobol_model_indiv$T) %>% 
  purrr::set_names(c("value", "bias", "std_error", "min_ci", "max_ci")) %>% 
  dplyr::mutate(parameter = c("ci_alpha", "ci_beta", "growth_infl", 
                              "mort_dbh_early", "mort_int_early"), 
                effect = "Total indices", 
                output = "Number of individuals") %>% 
  dplyr::bind_rows(sobol_model_indiv_df, .) %>% 
  dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0, 
                                          value > 1 ~ 1, 
                                          value > 0 & value < 1 ~ value), 
                min_ci = dplyr::case_when(min_ci < 0 ~ 0, 
                                          min_ci > 0 ~ min_ci), 
                max_ci = dplyr::case_when(max_ci > 1 ~ 1, 
                                          max_ci < 1 ~ max_ci))

ggplot_sobol_individual <- ggplot(data = sobol_model_indiv_df) + 
  geom_point(aes(x = parameter, y = value, col = effect), 
             size = 5, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(x  = parameter, ymin = min_ci, ymax = max_ci,col = effect), 
                width = 0.1, position = position_dodge(width = 0.5)) +
  scale_color_viridis_d(name = "", option = "C") + 
  scale_y_continuous(name = "Effect strength", limits = c(0, 1)) + 
  scale_x_discrete(name = "Parameter") +
  theme(legend.position = "bottom") +
  theme_classic(base_size = base_size)

suppoRt::save_ggplot(plot = ggplot_sobol_individual, 
                     filename = "ggplot_sobol_individual.png", 
                     path = "Figures/", 
                     dpi = dpi,
                     width = width_full, height = height_small, units = units, 
                     overwrite = overwrite)

#### Pair-correlation function ####

set.seed(42)

# sample parameters #
param_set_1_pcf <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                          0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                          1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                          -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                          -10, -7.8),  # mort_int_late CI from Holzwarth et al.
                                                 ncol = 2, byrow = TRUE))

param_set_2_pcf <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                          0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                          1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                          -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                          -10, -7.8),  # mort_int_late CI from Holzwarth et al.
                                                 ncol = 2, byrow = TRUE))

# create an instance of the class sobol #
sobol_model_pcf <- sensitivity::sobol2007(model = NULL, 
                                            X1 =  data.frame(param_set_1_pcf), 
                                            X2 =  data.frame(param_set_1_pcf), 
                                            conf = 0.95,
                                            nboot = 10000) 

# get parameter combinitations from sobol model
param_sampled_pcf <- purrr::map(seq_len(nrow(sobol_model_pcf$X)), 
                                  function(x) as.numeric(sobol_model_pcf$X[x, ]))

# get the simulated model response #
simulation_results_pcf <- 
  suppoRt::submit_to_cluster(calc_sobol_pcf, 
                             x = param_sampled_pcf, 
                             const = list(data = pattern_1999_dt,
                                          parameters = parameters_fitted_biotic,
                                          plot_area = plot_area,
                                          years = years,
                                          save_each = save_each),
                             n_jobs = length(param_sampled_pcf),
                             template = list(job_name = "sobol_pcf",
                                             walltime = "02:00:00",
                                             queue = "medium", 
                                             service = "short",
                                             mem_cpu = "2048", 
                                             log_file = "sobol_pcf.log"))

# flatten to vector #
simulation_results_pcf <- purrr::flatten_dbl(simulation_results_pcf)

suppoRt::save_rds(object = simulation_results_pcf,
                  filename = "simulation_results_pcf.rds",
                  path = "Data/Output/SA/",
                  overwrite = overwrite)

simulation_results_indiv <- readr::read_rds("Data/Output/SA/sa_simulation_results_indiv.rds")

# add the simulation results to the sobol instance #   
sensitivity::tell(sobol_model_pcf, simulation_results_pcf)

# save results #
suppoRt::save_rds(object = sobol_model_pcf, 
                  filename = "sobol_model_pcf.rds", 
                  path = "Data/Output/", 
                  overwrite = overwrite)
