###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Helper functions model comparison structure #### 

calc_n_comp <- function(data, verbose = TRUE) {
  
  n_data <- length(data)
  
  result <- purrr::map(seq_along(data), function(x) {
    
    # print progress
    if (verbose) {
      message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
    }
    
    # get living trees of last time step
    temp_data <- dplyr::filter(data[[x]], type != "dead")
    
    dplyr::group_by(temp_data, i, type) %>% 
      dplyr::summarise(n = dplyr::n())
  })
  
  if (verbose) { 
    message("")
  }
  
  result <- dplyr::bind_rows(result) %>% 
    dplyr::group_by(i, type) %>% 
    dplyr::summarise(n_mean = mean(n), 
                     n_min = min(n), 
                     n_max = max(n), 
                     n_sd = sd(n))
  
  return(result)
}

# calc_n_ingrowth_comp <- function(data, verbose = TRUE) {
#   
#   n_data <- length(data)
#   
#   result <- purrr::map(seq_along(data), function(x) {
#     
#     # print progress
#     if (verbose) {
#       message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
#     }
#     
#     temp_data_start <- dplyr::filter(data[[x]], i == 0, type != "dead")
#     
#     temp_data_end <- dplyr::filter(data[[x]], i != 0, type != "dead", 
#                                    !id %in% temp_data_start$id)
#     
#     dplyr::group_by(temp_data_end, i) %>% 
#       dplyr::summarise(n = dplyr::n())
#   })
#   
#   if (verbose) { 
#     message("")
#   }
#   
#   result <- dplyr::bind_rows(result) %>% 
#     dplyr::group_by(i) %>% 
#     dplyr::summarise(n_mean = mean(n), 
#                      n_min = min(n), 
#                      n_max = max(n), 
#                      n_sd = sd(n))
#   
#   return(result)
# }
# 
# calc_n_dead_comp <- function(data, by, verbose = TRUE) {
#   
#   n_data <- length(data)
#   
#   # loop through all data
#   result <- purrr::map(seq_along(data), function(x) {
#     
#     # print progress
#     if (verbose) {
#       message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
#     }
#     
#     # get living trees of last time step
#     temp_data <- dplyr::filter(data[[x]], type == "dead")
#     
#     dplyr::mutate(temp_data, i = cut(i, breaks = seq(from = 0, to = 50, by = 10))) %>% 
#       dplyr::group_by(i) %>% 
#       dplyr::summarise(n = dplyr::n())
#   })
#   
#   if (verbose) {
#     message("")
#   }
#   
#   return(result)
# }

calc_dbh_dist_comp <- function(data, sim_i, by, verbose = TRUE) {
  
  n_data <- length(data)
  
  # loop through all data
  result <- purrr::map(seq_along(data), function(x) {
 
    # print progress
    if (verbose) {
      message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
    }
     
    # get living trees of last time step
    temp_data <- dplyr::filter(data[[x]], i == sim_i, type != "dead")
    
    # count number of trees for relative n
    n_points <- nrow(temp_data)
      
    # classify to dbh class and count n
    dplyr::mutate(temp_data, dbh_class = cut(dbh, breaks = seq(from = 0, 
                                                               to = max(dbh) + by, 
                                                               by = by), 
                                  labels = FALSE)) %>% 
      dplyr::group_by(dbh_class) %>% 
      dplyr::summarise(n = dplyr::n()) %>% 
      dplyr::mutate(n_rel = n / n_points)
  })
  
  if (verbose) { 
    message("")
  }
  
  return(result)
}

#### helper function model comparison #### 
calc_growth_comp <- function(data, sim_i, by, verbose = TRUE) {
  
  n_data <- length(data)
  
  # loop through all data
  result <- purrr::map(seq_along(data), function(x) {
    
    # print progress
    if (verbose) {
      message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
    }

    # get current data
    temp_data <- dplyr::filter(data[[x]], i <= sim_i)
    
    # how many years did each tree survive
    temp_years_living <- dplyr::group_by(temp_data, id) %>% 
      dplyr::summarise(min_i = min(i),
                       max_i = max(i), 
                       years_lived = max_i - min_i)
    
    # get only living trees at start of simulation
    temp_start <- dplyr::filter(temp_data,
                                i == min(i), type != "dead") %>% 
      dplyr::select(id, dbh)

    # get only last year of each tree
    temp_end <- dplyr::group_by(temp_data, id) %>% 
      dplyr::top_n(n = 1, wt = i) %>% 
      dplyr::select(id, dbh)
    
    # combine data to one df
    dplyr::left_join(x = temp_start, 
                     y = temp_end, 
                     by = "id",
                     suffix = c(".start", ".end")) %>% 
      dplyr::left_join(y = temp_years_living, by = "id") %>% 
      dplyr::filter(years_lived > 0) %>% 
      # tidyr::replace_na(replace = list(dbh.start = 0)) %>% 
      dplyr::mutate(dbh_inc = (dbh.end - dbh.start) / years_lived) %>% 
      dplyr::mutate(dbh_class = cut(dbh.start, breaks = seq(from = 0, 
                                                            to = max(dbh.start) + by, 
                                                            by = by), 
                                    labels = FALSE))
  })
  
  if (verbose) {
    message("")
  }
  
  return(result)
}
