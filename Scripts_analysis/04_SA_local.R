###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Local sensitivity analysis ####

#### Import libraries and data ####

# load packages #
source("Helper_functions/helper_functions_setup.R")
source("Helper_functions/helper_functions_sa.R")

# import data #
parameters_beech_fitted <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt", 
                                                  sep = ";")

pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

#### Set SA parameters ####
repetitions <- 50 # 25

plot_area <- pattern_1999$window
years <- 50 # 50
save_each <- 5
return_seedlings <- FALSE
return_tibble <- TRUE
return_nested <- FALSE
verbose <- FALSE

#### Create parameters ####

# increase parameters by 5% and 10%
parameters_beech_inc_5 <- change_parameters(x = parameters_beech_fitted, 
                                            change = 0.05) %>% 
  rep(each = repetitions)

parameters_beech_inc_10 <- change_parameters(x = parameters_beech_fitted, 
                                             change = 0.1) %>% 
  rep(each = repetitions)

# decrease parameters by 5% and 10%
parameters_beech_dec_5 <- change_parameters(x = parameters_beech_fitted, 
                                            change = -0.05) %>% 
  rep(each = repetitions)

parameters_beech_dec_10 <- change_parameters(x = parameters_beech_fitted, 
                                             change = -0.1) %>% 
  rep(each = repetitions)

parameters_beech_default <- list(parameters_beech_fitted) %>% 
  rep(each = repetitions)

# remove parameters that don't have an influence
parameters_beech_inc_5 <- parameters_beech_inc_5[!names(parameters_beech_inc_5) 
                                                 %in% c("ci_max_dist",
                                                        "seed_max_dist", 
                                                        "growth_mod")]

parameters_beech_inc_10 <- parameters_beech_inc_10[!names(parameters_beech_inc_10) 
                                                   %in% c("ci_max_dist", 
                                                          "seed_max_dist", 
                                                          "growth_mod")]

parameters_beech_dec_5 <- parameters_beech_dec_5[!names(parameters_beech_dec_5) 
                                                 %in% c("ci_max_dist", 
                                                        "seed_max_dist", 
                                                        "growth_mod")]

parameters_beech_dec_10 <- parameters_beech_dec_10[!names(parameters_beech_dec_10) 
                                                   %in% c("ci_max_dist", 
                                                          "seed_max_dist",
                                                          "growth_mod")]

#### Pre-processing of input data ####
data <- tibble::as_tibble(pattern_1999) %>%
  dplyr::select(x, y, species, dbh_99, type) %>% 
  dplyr::filter(species == "beech", type != "dead") %>%
  dplyr::mutate(type = "adult") %>% 
  dplyr::select(-species) %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh_99")

# rm(pattern_1999)

#### Default parameters ####
# sa_default <- purrr::map(parameters_beech_default, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    return_seedlings = return_seedlings,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_default_y50_e5_r50 <- suppoRt::submit_to_cluster(rabmp::run_model_biotic,
                                                    parameters = parameters_beech_default,
                                                    const = list(data = data,
                                                                 plot_area = plot_area,
                                                                 years = years,
                                                                 save_each = save_each,
                                                                 return_seedlings = return_seedlings,
                                                                 return_nested = return_nested,
                                                                 return_tibble = return_tibble,
                                                                 verbose = verbose),
                                                    n_jobs = length(parameters_beech_default),
                                                    template = list(job_name = "sa_default",
                                                                    walltime = "01:00:00",
                                                                    queue = "medium", 
                                                                    service = "short",
                                                                    mem_cpu = "1024", 
                                                                    log_file = "sa_default_y50_e5_r50.log"))

names(sa_default_y50_e5_r50) <- rep("default", times = repetitions)

suppoRt::save_rds(object = sa_default_y50_e5_r50,
                  filename = "sa_default_y50_e5_r50.rds",
                  path = "Data/Output/SA/")

# rm(sa_default_y50_e5_r50)

#### Increased parameters #### 
# sa_increased_5 <- purrr::map(parameters_beech_inc_5, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    return_seedlings = return_seedlings,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_increased_5_y50_e5_r50 <- suppoRt::submit_to_cluster(rabmp::run_model_biotic,
                                                        parameters = parameters_beech_inc_5,
                                                        const = list(data = data,
                                                                     plot_area = plot_area,
                                                                     years = years,
                                                                     save_each = save_each,
                                                                     return_seedlings = return_seedlings,
                                                                     return_nested = return_nested,
                                                                     return_tibble = return_tibble,
                                                                     verbose = verbose),
                                                        n_jobs = length(parameters_beech_inc_5),
                                                        template = list(job_name = "sa_inc_5",
                                                                        walltime = "01:00:00",
                                                                        queue = "medium", 
                                                                        service = "short",
                                                                        mem_cpu = "1024",  
                                                                        log_file = "sa_inc_5_y50_e5_r50.log"))

names(sa_increased_5_y50_e5_r50) <- names(parameters_beech_inc_5)

suppoRt::save_rds(object = sa_increased_5_y50_e5_r50,
                  filename = "sa_increased_5_y50_e5_r50.rds",
                  path = "Data/Output/SA/")

# rm(sa_increased_5_y50_e5_r50)

# sa_increased_10 <- purrr::map(parameters_beech_inc_10, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    return_seedlings = return_seedlings,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_increased_10_y50_e5_r50 <- suppoRt::submit_to_cluster(rabmp::run_model_biotic,
                                                         parameters = parameters_beech_inc_10,
                                                         const = list(data = data,
                                                                      plot_area = plot_area,
                                                                      years = years,
                                                                      save_each = save_each,
                                                                      return_seedlings = return_seedlings,
                                                                      return_nested = return_nested,
                                                                      return_tibble = return_tibble,
                                                                      verbose = verbose),
                                                         n_jobs = length(parameters_beech_inc_10),
                                                         template = list(job_name = "sa_inc_10",
                                                                         walltime = "01:00:00",
                                                                         queue = "medium", 
                                                                         service = "short",
                                                                         mem_cpu = "1024", 
                                                                         log_file = "sa_inc_10_y50_e5_r50.log"))

names(sa_increased_10_y50_e5_r50) <- names(parameters_beech_inc_10)

suppoRt::save_rds(object = sa_increased_10_y50_e5_r50,
                  filename = "sa_increased_10_y50_e5_r50.rds",
                  path = "Data/Output/SA/")

# rm(sa_increased_10_y50_e5_r50)

#### Decreased parameters ####
# sa_decreased_5 <- purrr::map(parameters_beech_dec_5, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    return_seedlings = return_seedlings,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_decreased_5_y50_e5_r50 <- suppoRt::submit_to_cluster(rabmp::run_model_biotic,
                                                        parameters = parameters_beech_dec_5,
                                                        const = list(data = data,
                                                                     plot_area = plot_area,
                                                                     years = years,
                                                                     save_each = save_each,
                                                                     return_seedlings = return_seedlings,
                                                                     return_nested = return_nested,
                                                                     return_tibble = return_tibble,
                                                                     verbose = verbose),
                                                        n_jobs = length(parameters_beech_dec_5),
                                                        template = list(job_name = "sa_dec_5",
                                                                        walltime = "01:00:00",
                                                                        queue = "medium", 
                                                                        service = "short",
                                                                        mem_cpu = "1024", 
                                                                        log_file = "sa_dec_5_y50_e5_r50.log"))

names(sa_decreased_5_y50_e5_r50) <- names(parameters_beech_dec_5)

suppoRt::save_rds(object = sa_decreased_5_y50_e5_r50,
                  filename = "sa_decreased_5_y50_e5_r50.rds",
                  path = "Data/Output/SA/")

# rm(sa_decreased_5_y50_e5_r50)

# sa_decreased_10 <- purrr::map(parameters_beech_dec_10, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    return_seedlings = return_seedlings,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_decreased_10_y50_e5_r50 <- suppoRt::submit_to_cluster(rabmp::run_model_biotic,
                                                         parameters = parameters_beech_dec_10,
                                                         const = list(data = data,
                                                                      plot_area = plot_area,
                                                                      years = years,
                                                                      save_each = save_each,
                                                                      return_seedlings = return_seedlings,
                                                                      return_nested = return_nested,
                                                                      return_tibble = return_tibble,
                                                                      verbose = verbose),
                                                         n_jobs = length(parameters_beech_dec_10),
                                                         template = list(job_name = "sa_dec_10",
                                                                         walltime = "01:00:00",
                                                                         queue = "medium", 
                                                                         service = "short",
                                                                         mem_cpu = "1024", 
                                                                         log_file = "sa_dec_10_y50_e5_r50.log"))

names(sa_decreased_10_y50_e5_r50) <- names(parameters_beech_dec_10)

suppoRt::save_rds(object = sa_decreased_10_y50_e5_r50,
                  filename = "sa_decreased_10_y50_e5_r50.rds",
                  path = "Data/Output/SA/")

# rm(sa_decreased_10_y50_e5_r50)
