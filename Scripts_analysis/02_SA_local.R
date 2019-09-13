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
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(sensitivity)
library(spatstat)
library(tidyverse)

parameters_beech_default <- rabmp::read_parameters("Data/Input/parameters_beech.txt", return_list = TRUE)

beech_1999_rec <- readr::read_rds("Data/Input/beech_1999_rec.rds")

source("Helper_functions/helper_functions_sa.R")

#### Set SA parameters ####
repetitions <- 10 # 50

plot_area <- beech_1999_rec$window
years <- 50 # 50
save_each <- 50
return_nested <- FALSE
verbose <- FALSE

#### Create parameters ####
# increase parameters by 5% and 10%
parameters_beech_inc_5 <- change_parameters(x = parameters_beech_default, 
                                            change = 0.05) %>% 
  rep(each = repetitions)

parameters_beech_inc_10 <- change_parameters(x = parameters_beech_default, 
                                             change = 0.1) %>% 
  rep(each = repetitions)
# decrease parameters by 5% and 10%
parameters_beech_dec_5 <- change_parameters(x = parameters_beech_default, 
                                            change = -0.05) %>% 
  rep(each = repetitions)

parameters_beech_dec_10 <- change_parameters(x = parameters_beech_default, 
                                             change = -0.1) %>% 
  rep(each = repetitions)

parameters_beech_default <- list(parameters_beech_default) %>% 
  rep(each = repetitions)

#### Pre-processing of input data ####

data <- tibble::as_tibble(beech_1999_rec) %>%
  rabmp::prepare_data(x = "x", y = "y", species = "species", type = "type", dbh = "dbh")

rm(pattern_1999_recon)

#### Default parameters ####
# sa_default <- purrr::map(parameters_beech_default, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_default <- suppoRt::submit_to_cluster(rabmp::run_model,
                                         parameters = parameters_beech_default,
                                         const = list(data = data,
                                                      plot_area = plot_area,
                                                      years = years,
                                                      save_each = save_each,
                                                      return_nested = return_nested,
                                                      verbose = verbose),
                                         n_jobs = length(parameters_beech_default),
                                         template = list(job_name = "sa_default",
                                                         walltime = "06:00:00",
                                                         queue = "medium", 
                                                         mem_cpu = "3072", 
                                                         log_file = "sa_default.log"))

names(sa_default) <- rep("default", times = repetitions)

suppoRt::save_rds(object = sa_default,
                  filename = "sa_default.rds",
                  path = "Data/Output/")

rm(sa_default)

#### Increased parameters #### 
# sa_increased_5 <- purrr::map(parameters_beech_inc_5, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_increased_5 <- suppoRt::submit_to_cluster(rabmp::run_model,
                                             parameters = parameters_beech_inc_5,
                                             const = list(data = data,
                                                          plot_area = plot_area,
                                                          years = years,
                                                          save_each = save_each,
                                                          return_nested = return_nested,
                                                          verbose = verbose),
                                             n_jobs = length(parameters_beech_inc_5),
                                             template = list(job_name = "sa_inc_5",
                                                             walltime = "06:00:00",
                                                             queue = "medium", 
                                                             mem_cpu = "3072", 
                                                             log_file = "sa_inc_5.log"))

names(sa_increased_5) <- names(parameters_beech_inc_5)

suppoRt::save_rds(object = sa_increased_5,
                  filename = "sa_increased_5.rds",
                  path = "Data/Output/")

rm(sa_increased_5)

# sa_increased_10 <- purrr::map(parameters_beech_inc_10, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_increased_10 <- suppoRt::submit_to_cluster(rabmp::run_model,
                                              parameters = parameters_beech_inc_10,
                                              const = list(data = data,
                                                           plot_area = plot_area,
                                                           years = years,
                                                           save_each = save_each,
                                                           return_nested = return_nested,
                                                           verbose = verbose),
                                              n_jobs = length(parameters_beech_inc_10),
                                              template = list(job_name = "sa_inc_10",
                                                              walltime = "06:00:00",
                                                              queue = "medium", 
                                                              mem_cpu = "3072", 
                                                              log_file = "sa_inc_10.log"))

names(sa_increased_10) <- names(parameters_beech_inc_10)

suppoRt::save_rds(object = sa_increased_10,
                  filename = "sa_increased_10.rds",
                  path = "Data/Output/")

rm(sa_increased_10)

#### Decreased parameters ####
# sa_decreased_5 <- purrr::map(parameters_beech_dec_5, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_decreased_5 <- suppoRt::submit_to_cluster(rabmp::run_model,
                                             parameters = parameters_beech_dec_5,
                                             const = list(data = data,
                                                          plot_area = plot_area,
                                                          years = years,
                                                          save_each = save_each,
                                                          return_nested = return_nested,
                                                          verbose = verbose),
                                             n_jobs = length(parameters_beech_dec_5),
                                             template = list(job_name = "sa_dec_5",
                                                             walltime = "06:00:00",
                                                             queue = "medium", 
                                                             mem_cpu = "3072", 
                                                             log_file = "sa_dec_5.log"))

names(sa_decreased_5) <- names(parameters_beech_dec_5)

suppoRt::save_rds(object = sa_decreased_5,
                  filename = "sa_decreased_5.rds",
                  path = "Data/Output/")

rm(sa_decreased_5)

# sa_decreased_10 <- purrr::map(parameters_beech_dec_10, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = x,
#                    plot_area = plot_area,
#                    years = years,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

sa_decreased_10 <- suppoRt::submit_to_cluster(rabmp::run_model,
                                              parameters = parameters_beech_dec_10,
                                              const = list(data = data,
                                                           plot_area = plot_area,
                                                           years = years,
                                                           save_each = save_each,
                                                           return_nested = return_nested,
                                                           verbose = verbose),
                                              n_jobs = length(parameters_beech_dec_10),
                                              template = list(job_name = "sa_dec_10",
                                                              walltime = "06:00:00",
                                                              queue = "medium", 
                                                              mem_cpu = "3072", 
                                                              log_file = "sa_dec_10.log"))

names(sa_decreased_10) <- names(parameters_beech_dec_10)

suppoRt::save_rds(object = sa_decreased_10,
                  filename = "sa_decreased_10.rds",
                  path = "Data/Output/")

rm(sa_decreased_10)
