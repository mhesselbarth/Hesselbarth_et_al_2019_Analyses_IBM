###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Systematic find sigma ####

#### Import libraries ####

# load packages #
source("Helper_functions/helper_functions_setup.R")
source("Helper_functions/helper_function_explore_sigma.R")

#### Import data ####
beech_1999_ppp <- readr::read_rds("Data/Input/beech_1999_ppp.rds")

beech_2007_ppp <- readr::read_rds("Data/Input/beech_2007_ppp.rds")
beech_2007_sapling_ppp <- readr::read_rds("Data/Input/beech_2007_sapling_ppp.rds")
beech_2007_adult_ppp <- readr::read_rds("Data/Input/beech_2007_adult_ppp.rds")

beech_2013_ppp <- readr::read_rds("Data/Input/beech_2013_ppp.rds")
beech_2013_sapling_ppp <- readr::read_rds("Data/Input/beech_2013_sapling_ppp.rds")
beech_2013_adult_ppp <- readr::read_rds("Data/Input/beech_2013_adult_ppp.rds")

parameters_fitted_abiotic <- rabmp::read_parameters("Data/Input/parameters_fitted_abiotic.txt",
                                                    sep = ";")

#### Pre-process data ####
beech_2007_df <- tibble::as_tibble(beech_2007_ppp)
beech_2013_df <- tibble::as_tibble(beech_2013_ppp)

plot_area <- beech_1999_ppp$window

probs <- list(c(0.1, 0.9), c(0.25, 0.75))
probs_id <- c(1, 2)

sigma <- seq(from = 5, to = 75, by = 5)

sigma_probs <- suppoRt::expand_grid_unique(x = sigma, 
                                           y = probs_id)


years <- 50
save_each <- 50

#### Run systematic sigma exploratation ####
model_runs_sigma <- suppoRt::submit_to_cluster(explore_sigma,
                                               sigma = sigma_probs[, 1],
                                               probs_id = sigma_probs[, 2],
                                               const = list(data = beech_1999_ppp,
                                                            parameters = parameters_fitted_abiotic,
                                                            probs = probs,
                                                            plot_area = plot_area,
                                                            years = years,
                                                            save_each = save_each),
                                               n_jobs = nrow(sigma_probs),
                                               template = list(job_name = "systematic",
                                                               walltime = "02:00:00",
                                                               queue = "medium",
                                                               service = "short",
                                                               mem_cpu = "2048",
                                                               log_file = "systematic.log"))

suppoRt::save_rds(object = model_runs_sigma,
                  filename = "model_runs_sigma.rds",
                  path = "Data/Output/model_runs/", 
                  overwrite = TRUE)
