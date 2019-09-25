#### helper function model comparison #### 
calc_dbh_dist <- function(data, by, verbose = TRUE) {
  
  n_data <- length(data)
  
  # loop through all data
  result <- purrr::map(seq_along(data), function(x) {
 
    # print progress
    if (verbose) {
      message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
    }
     
    # get living trees of last time step
    temp_data <- dplyr::filter(data[[x]], i == max(i), type != "dead")
    
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
calc_growth <- function(data, by, verbose = TRUE) {
  
  n_data <- length(data)
  
  # loop through all data
  result <- purrr::map(seq_along(data), function(x) {
    
    # print progress
    if (verbose) {
      message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
    }

    # get current data
    temp_data <- data[[x]]
    
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

calc_died <- function(data, by, verbose = TRUE) {
  
  n_data <- length(data)
  
  # loop through all data
  result <- purrr::map(seq_along(data), function(x) {
    
    # print progress
    if (verbose) {
      message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
    }
    
    # get number of living trees
    temp_living <- dplyr::filter(data[[x]], i == min(i), type != "dead") %>%
      dplyr::mutate(dbh_class = cut(dbh, breaks = seq(from = 0, 
                                                      to = max(dbh) + by, 
                                                      by = by), 
                                    labels = FALSE)) %>% 
      dplyr::group_by(dbh_class) %>% 
      dplyr::summarise(n_living = dplyr::n())
    
    # get number of trees that died absolute and relative to starting living trees
    dplyr::filter(data[[x]], i != min(i), type == "dead") %>% 
      dplyr::mutate(dbh_class = cut(dbh, breaks = seq(from = 0, 
                                                            to = max(dbh) + by, 
                                                            by = by), 
                                    labels = FALSE)) %>% 
      dplyr::group_by(dbh_class) %>% 
      dplyr::summarise(n_died = n()) %>% 
      dplyr::left_join(y = temp_living, by = "dbh_class") %>% 
      tidyr::replace_na(replace = list(n_died = 0, n_living = 0)) %>% 
      dplyr::mutate(n_died_rel = n_died / n_living, 
                    n_died_rel_total = n_died / sum(temp_living$n_living))
  })
  
  if (verbose) {
    message("")
  }
  
  return(result)
}
