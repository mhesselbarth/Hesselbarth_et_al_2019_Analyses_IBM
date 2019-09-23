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
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(spatstat)
library(tidyverse)

parameters_beech_default <- rabmp::read_parameters("Data/Input/parameters_beech.txt", return_list = TRUE)

pattern_1999_recon <- readr::read_rds("Data/Input/beech_1999_rec.rds")

#### Set SA parameters ####
repetitions <- 50 # 50

plot_area <- pattern_1999_recon$window
years <- rep(x = 100, times = repetitions) # 50
save_each <- 100
return_nested <- FALSE
verbose <- FALSE

#### Pre-processing of input data ####

data <- tibble::as_tibble(pattern_1999_recon) %>%
  dplyr::filter(species == "beech") %>%
  rabmp::prepare_data(x = "x", y = "y", species = "species", type = "type", dbh = "dbh")

rm(pattern_1999_recon)

#### Run model ####
# model_run_y100_r50_e5 <- purrr::map(years, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = parameters_beech_default,
#                    plot_area = plot_area,
#                    years = x,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

model_run_y100_e100_r50 <- suppoRt::submit_to_cluster(rabmp::run_model,
                                                    years = years,
                                                    const = list(data = data,
                                                                 parameters = parameters_beech_default,
                                                                 plot_area = plot_area,
                                                                 save_each = save_each,
                                                                 return_nested = return_nested,
                                                                 verbose = verbose),
                                                    n_jobs = length(years),
                                                    template = list(job_name = "y100_e100_r50",
                                                                    walltime = "12:00:00",
                                                                    queue = "medium", 
                                                                    mem_cpu = "4096", 
                                                                    log_file = "y100_e100_r50.log"))

suppoRt::save_rds(object = model_run_y100_e100_r50,
                  filename = "model_run_y100_e100_r50.rds",
                  path = "Data/Output/")

# rm(model_run_y100_e100_r50)
