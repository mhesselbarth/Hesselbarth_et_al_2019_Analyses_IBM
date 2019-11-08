###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Run biotic model ####

#### Import libraries and data ####

# load packages #
source("Helper_functions/helper_functions_setup.R")

parameters_fitted_biotic <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt",
                                                   sep = ";")

# beech_1999_rec_ppp <- readr::read_rds("Data/Input/beech_1999_rec_ppp.rds")

beech_1999_ppp <- readr::read_rds("Data/Input/beech_1999_ppp.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

#### Set SA parameters ####
repetitions <- 50 # 50

years <- rep(x = 50, times = repetitions) # 50
save_each <- 5
return_nested <- FALSE
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
# model_run_y50_e5_r50_reco_b <- suppoRt::submit_to_cluster(rabmp::run_model_biotic,
#                                                            years = years,
#                                                            const = list(data = data_reconstruction,
#                                                                         parameters = parameters_fitted_biotic,
#                                                                         plot_area = plot_area,
#                                                                         save_each = save_each,
#                                                                         return_nested = return_nested,
#                                                                         verbose = verbose),
#                                                            n_jobs = length(years),
#                                                            template = list(job_name = "y50_e5_r50_reco",
#                                                                            walltime = "02:00:00",
#                                                                            queue = "medium", 
#                                                                            service = "short",
#                                                                            mem_cpu = "1024", 
#                                                                            log_file = "y50_e5_r50_reco.log"))
# 
# suppoRt::save_rds(object = model_run_y50_e5_r50_reco_b,
#                   filename = "model_run_y50_e5_r50_reco_b.rds",
#                   path = "Data/Output/model_runs")

# rm(model_run_y50_e5_r50_reco_b)

# real world data #
model_run_y50_e5_r50_real_b <- suppoRt::submit_to_cluster(rabmp::run_model_biotic,
                                                           years = years,
                                                           const = list(data = data_real,
                                                                        parameters = parameters_fitted_biotic,
                                                                        plot_area = plot_area,
                                                                        save_each = save_each,
                                                                        return_nested = return_nested,
                                                                        verbose = verbose),
                                                           n_jobs = length(years),
                                                           template = list(job_name = "y50_e5_r50_real",
                                                                           walltime = "02:00:00",
                                                                           queue = "medium", 
                                                                           service = "short",
                                                                           mem_cpu = "1024", 
                                                                           log_file = "y50_e5_r50_real.log"))

suppoRt::save_rds(object = model_run_y50_e5_r50_real_b,
                  filename = "model_run_y50_e5_r50_real_b.rds",
                  path = "Data/Output/model_runs")

# rm(model_run_y50_e5_r50_real_b)
