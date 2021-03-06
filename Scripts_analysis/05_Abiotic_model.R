###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Run abiotic model ####

#### Import libraries and data ####

# load packages #
source("Helper_functions/helper_functions_setup.R")

parameters_fitted_abiotic_real <- rabmp::read_parameters("Data/Input/parameters_fitted_abiotic_real.txt", 
                                                    sep = ";")

# beech_1999_rec_ppp <- readr::read_rds("Data/Input/beech_1999_rec_ppp.rds")

beech_1999_ppp <- readr::read_rds("Data/Input/beech_1999_ppp.rds")

abiotic_habitats_real <- readr::read_rds("Data/Input/abiotic_cond_real_model.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

#### Set SA parameters ####
repetitions <- 50 # 50

years <- rep(x = 50, times = repetitions) # 50
probs <- c(0.25, 0.75)
save_each <- 5
verbose <- FALSE

#### Pre-processing of input data ####

# data_reconstruction <- tibble::as_tibble(beech_1999_rec_ppp) %>%
#   dplyr::select(-species) %>% 
#   rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh")

data_real <- tibble::as_tibble(beech_1999_ppp) %>%
  dplyr::select(x, y, species, dbh_99, type) %>% 
  dplyr::mutate(type = "adult") %>% 
  dplyr::select(-species) %>% 
  rabmp::prepare_data(x = "x", y = "y", type = "type", dbh = "dbh_99")

# rm(beech_1999_rec_ppp)

rm(beech_1999_ppp)

#### Run model ####

# # reconstructed data #
# model_run_y50_e5_r50_reco_a <- suppoRt::submit_to_cluster(rabmp::run_model_abiotic,
#                                                           years = years,
#                                                           const = list(data = data_reconstruction,
#                                                                        parameters = parameters_fitted_abiotic,
#                                                                        abiotic = abiotic_habitats_reco$scaled,
#                                                                        probs = probs,
#                                                                        plot_area = plot_area,
#                                                                        save_each = save_each,
#                                                                        verbose = verbose),
#                                                           n_jobs = length(years),
#                                                           template = list(job_name = "abiotic_reco",
#                                                                           walltime = "02:00:00",
#                                                                           queue = "medium", 
#                                                                           service = "short",
#                                                                           mem_cpu = "2048", 
#                                                                           log_file = "abiotic_reco.log"))
# 
# suppoRt::save_rds(object = model_run_y50_e5_r50_reco_a,
#                   filename = "model_run_y50_e5_r50_reco_a.rds",
#                   path = "Data/Output/model_runs")

# rm(model_run_y50_e5_r50_reco_a)

# real world data #
model_run_y50_e5_r50_real_a <- suppoRt::submit_to_cluster(rabmp::run_model_abiotic,
                                                          years = years,
                                                          const = list(data = data_real,
                                                                       parameters = parameters_fitted_abiotic_real,
                                                                       abiotic = abiotic_habitats_real$scaled,
                                                                       probs = probs,
                                                                       plot_area = plot_area,
                                                                       save_each = save_each,
                                                                       verbose = verbose),
                                                          n_jobs = length(years),
                                                          template = list(job_name = "abiotic_real",
                                                                          walltime = "02:00:00",
                                                                          queue = "medium", 
                                                                          service = "short",
                                                                          mem_cpu = "4096", 
                                                                          log_file = "abiotic_real.log"))

suppoRt::save_rds(object = model_run_y50_e5_r50_real_a,
                  filename = "model_run_y50_e5_r50_real_a.rds",
                  path = "Data/Output/model_runs")

# rm(model_run_y50_e5_r50_real_a)
