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

pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

years <- 50
save_each <- years

n <- 50

set.seed(42)

#### Pre-process data ####
pattern_1999_dt <- tibble::as_tibble(pattern_1999) %>%
  dplyr::select(x, y, species, dbh_99, type) %>% 
  dplyr::filter(species == "beech") %>%
  dplyr::mutate(type = "adult") %>% 
  dplyr::select(-species) %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh_99")

param_set_1 <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                      0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                      1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                      -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                      1.2, 2.5), 
                                             ncol = 2, byrow = TRUE)) # mort_int_early CI from Holzwarth et al.

param_set_2 <- tgp::lhs(n = n, rect = matrix(data = c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                      0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                      1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                      -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                      1.2, 2.5), 
                                             ncol = 2, byrow = TRUE)) # mort_int_early CI from Holzwarth et al.

#### Number of individuals ####

# create an instance of the class sobol #
sobol_model_indiv <- sensitivity::soboljansen(model = NULL, 
                                              X1 =  data.frame(param_set_1), 
                                              X2 =  data.frame(param_set_2), 
                                              conf = 0.95,
                                              nboot = 5000) 

# get parameter combinitations from sobol model
param_sampled_indiv <- purrr::map(seq_len(nrow(sobol_model_indiv$X)), 
                                  function(x) as.numeric(sobol_model_indiv$X[x, ]))

# get the simulated model response #
simulation_results_indiv <- 
  suppoRt::submit_to_cluster(calc_sobol_indiv, 
                             x = param_sampled_indiv, 
                             const = list(data = pattern_1999_dt,
                                          parameters = parameters_fitted_biotic,
                                          plot_area = plot_area,
                                          years = years,
                                          save_each = save_each),
                             n_jobs = length(param_sampled),
                             template = list(job_name = "sobol_indiv",
                                             walltime = "02:00:00",
                                             queue = "medium", 
                                             service = "short",
                                             mem_cpu = "1024", 
                                             log_file = "sobol_indiv.log"))

# flatten to vector #
simulation_results_indiv <- purrr::flatten_dbl(simulation_results_indiv)

# add the simulation results to the sobol instance #   
sensitivity::tell(sobol_model_indiv, simulation_results_indiv)

# save results #
suppoRt::save_rds(object = sobol_model_indiv, 
                  filename = "sobol_model_indiv.rds", 
                  path = "Data/Output/", 
                  overwrite = overwrite)

#### Pair-correlation function ####
# create an instance of the class sobol #
sobol_model_pcf <- sensitivity::soboljansen(model = NULL, 
                                            X1 =  data.frame(param_set_1), 
                                            X2 =  data.frame(param_set_2), 
                                            conf = 0.95,
                                            nboot = 5000) 

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
                             n_jobs = length(param_sampled),
                             template = list(job_name = "sobol_pcf",
                                             walltime = "02:00:00",
                                             queue = "medium", 
                                             service = "short",
                                             mem_cpu = "1024", 
                                             log_file = "sobol_pcf.log"))

# flatten to vector #
simulation_results_pcf <- purrr::flatten_dbl(simulation_results_pcf)

# add the simulation results to the sobol instance #   
sensitivity::tell(sobol_model_pcf, simulation_results_pcf)

# save results #
suppoRt::save_rds(object = sobol_model_pcf, 
                  filename = "sobol_model_pcf.rds", 
                  path = "Data/Output/", 
                  overwrite = overwrite)
