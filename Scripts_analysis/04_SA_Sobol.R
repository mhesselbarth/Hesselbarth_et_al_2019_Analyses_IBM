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

parameters_fitted_abiotic <- rabmp::read_parameters("Data/Input/parameters_fitted_abiotic.txt",
                                                    sep = ";")

pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

years <- 5
save_each <- years

#### Pre-process data ####
pattern_1999_dt <- tibble::as_tibble(pattern_1999) %>%
  dplyr::select(x, y, species, dbh_99, type) %>% 
  dplyr::filter(species == "beech") %>%
  dplyr::mutate(type = "adult") %>% 
  dplyr::select(-species) %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh_99")

#### Number of individuals ####
param_set_1 <- tgp::lhs(n = 5, rect = matrix(c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                1.2, 2.5), ncol = 2)) # mort_int_early CI from Holzwarth et al.

param_set_2 <- tgp::lhs(n = 5, rect = matrix(c(0.9879718, 1.125665, # ci_alpha +- SE * 1.96
                                                0.4242094, 0.4458642, # ci_beta +- SE * 1.96
                                                1.251156, 1.44794, # growth_infl +- SE * 1.96
                                                -2.4, -1.9, # mort_dbh_early CI from Holzwarth et al.
                                                1.2, 2.5), ncol = 2)) # mort_int_early CI from Holzwarth et al.

#### Run Sobol analysis ####

# create an instance of the class sobol #
sobol_model <- sensitivity::soboljansen(model = NULL, 
                                        X1 = param_set_1, 
                                        X2 = param_set_2, 
                                        conf = 0.95,
                                        nboot = 1000) 

# get parameter combinitations from sobol model
param_sampled <- purrr::map(seq_len(nrow(sobol_model$X)), 
                            function(x) sobol_model$X[x, ])

# get the simulated model response #
simulation_resuls <- suppoRt::submit_to_cluster(calc_sobol_indiv, 
                                                param = param_sampled, 
                                                const = list(data = pattern_1999_dt,
                                                             plot_area = plot_area,
                                                             years = years,
                                                             save_each = save_each),
                                                n_jobs = length(param_sampled),
                                                template = list(job_name = "sobol_indiv",
                                                                walltime = "00:30:00",
                                                                queue = "medium", 
                                                                service = "short",
                                                                mem_cpu = "1024", 
                                                                log_file = "sobol_indiv.log"))

# add the simulation results to the sobol instance #   
sensitivity::tell(sobol_model, simulation_resuls)

# show result #
print(sobol_model)
