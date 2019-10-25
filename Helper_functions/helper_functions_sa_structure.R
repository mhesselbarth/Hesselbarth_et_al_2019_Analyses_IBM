###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Helper functions SA structure #### 

calc_n_sa <- function(default, changed, 
                      verbose = TRUE) {
  
  # get name of changed to get changed parameters
  names_parameters <- names(changed)
  
  # get length of input for printig
  n_default <- length(default)
  n_changed <- length(changed)
  
  # repeat counter 1:repetitions of default data for each changed parameter
  counter_default <- rep(x = 1:n_default, times = length(unique(names(changed))))
  
  # loop through all input data
  individuals_default <- purrr::map(seq_along(default), function(x) {
    
    if (verbose) {
      message("\r> Progress (default): ", x, "/", n_default, appendLF = FALSE)
    }
    
    # get number of living trees
    temp_default <- dplyr::filter(default[[x]], i == max(i), type != "dead") %>% 
      nrow()
  })
  
  if (verbose) {
    message("")
  }
  
  # loop through all input data
  individuals_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      message("\r> Progress (default): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    # get default data
    temp_default <- individuals_default[[counter_default[x]]]
    
    # get number of living trees
    temp_changed <- dplyr::filter(changed[[x]], i == max(i), type != "dead") %>%
      nrow()
    
    tibble::tibble(n_changed = temp_changed, n_default = temp_default) %>% 
      dplyr::mutate(diff_n = (n_changed - n_default) / n_default)
  })
  
  # add changed parameters as names
  names(individuals_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  individuals_changed <- dplyr::bind_rows(individuals_changed, 
                                          .id = "parameter") %>% 
    dplyr::group_by(parameter) %>% 
    dplyr::summarise(diff_n = mean(suppoRt::replace_infinite(diff_n), 
                                   na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  return(individuals_changed)
}

calc_growth_sa <- function(default, changed, 
                           verbose = TRUE) {
  
  # get name of changed to get changed parameters
  names_parameters <- names(changed)
  
  # get length of input for printig
  n_default <- length(default)
  n_changed <- length(changed)
  
  # repeat counter 1:repetitions of default data for each changed parameter
  counter_default <- rep(x = 1:n_default, times = length(unique(names(changed))))
  
  # loop through all input data
  growth_default <- purrr::map(seq_along(default), function(x) {
    
    # print progress
    if (verbose) {
      message("\r> Progress (default): ", x, "/", n_default, appendLF = FALSE)
    }
    
    # get current data
    temp_data <- default[[x]]
    
    # how many years did each tree survive
    temp_years_living <- dplyr::group_by(temp_data, id) %>% 
      dplyr::summarise(min_i = min(i),
                       max_i = max(i), 
                       years_lived = max_i - min_i)
    
    # get id of trees that lived more than 1 year
    temp_id_lived <- dplyr::filter(temp_years_living, years_lived > 0) %>% 
      dplyr::pull(id)
    
    # filter data 
    temp_data <- dplyr::filter(temp_data, id %in% temp_id_lived)
    
    # get only living trees at start of simulation
    temp_start <- dplyr::filter(temp_data, i == min(i), type != "dead") %>% 
      dplyr::select(id, dbh)
    
    # get only last year of each tree
    temp_end <- dplyr::group_by(temp_data, id) %>% 
      dplyr::top_n(n = 1, wt = i) %>% 
      dplyr::select(id, dbh)
    
    # combine data to one df
    dplyr::full_join(x = temp_start, 
                     y = temp_end, 
                     by = "id", suffix = c(".start", ".end")) %>% 
      dplyr::left_join(y = temp_years_living, by = "id") %>% 
      tidyr::replace_na(replace = list(dbh.start = 0)) %>% 
      dplyr::mutate(dbh_inc = (dbh.end - dbh.start) / years_lived) %>% 
      dplyr::pull(dbh_inc) %>% 
      mean(na.rm = TRUE)
  })
  
  if (verbose) {
    message("")
  }
  
  growth_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    # get default data
    temp_default <- growth_default[[counter_default[x]]]
    
    # get data of last time step
    temp_data <- changed[[x]]
    
    # how many years did each tree survive
    temp_years_living <- dplyr::group_by(temp_data, id) %>% 
      dplyr::summarise(min_i = min(i),
                       max_i = max(i), 
                       years_lived = max_i - min_i)
    
    # get id of trees that lived more than 1 year
    temp_id_lived <- dplyr::filter(temp_years_living, years_lived > 0) %>% 
      dplyr::pull(id)
    
    # filter data 
    temp_data <- dplyr::filter(temp_data, id %in% temp_id_lived)
    
    # get only living trees at start of simulation
    temp_start <- dplyr::filter(temp_data, i == min(i), type != "dead") %>% 
      dplyr::select(id, dbh)
    
    # get only last year of each tree
    temp_end <- dplyr::group_by(temp_data, id) %>% 
      dplyr::top_n(n = 1, wt = i) %>% 
      dplyr::select(id, dbh)
    
    # combine data to one df and calculate difference to default
    dplyr::full_join(x = temp_start, y = temp_end, 
                     by = "id", suffix = c(".start", ".end")) %>% 
      dplyr::left_join(y = temp_years_living, by = "id") %>% 
      tidyr::replace_na(replace = list(dbh.start = 0)) %>% 
      dplyr::mutate(dbh_inc = (dbh.end - dbh.start) / years_lived) %>% 
      dplyr::pull(dbh_inc) %>% 
      mean(na.rm = TRUE) %>% 
      tibble::tibble(diff_inc = .) %>% 
      dplyr::mutate(diff_inc = (diff_inc - temp_default) / temp_default)
  })
  
  # add changed parameters as names
  names(growth_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  growth_changed <- dplyr::bind_rows(growth_changed, .id = "parameter") %>% 
    dplyr::group_by(parameter) %>% 
    dplyr::summarise(diff_inc = mean(suppoRt::replace_infinite(diff_inc), 
                                     na.rm = TRUE)) %>%
    dplyr::ungroup() 
  
  return(growth_changed)
}

calc_died_sa <- function(default, changed, 
                         verbose = TRUE) {
  
  # get name of changed to get changed parameters
  names_parameters <- names(changed)
  
  # get length of input for printig
  n_default <- length(default)
  n_changed <- length(changed)
  
  # repeat counter 1:repetitions of default data for each changed parameter
  counter_default <- rep(x = 1:n_default, times = length(unique(names(changed))))
  
  # loop through all input data
  died_default <- purrr::map(seq_along(default), function(x) {
    
    if (verbose) {
      message("\r> Progress (default): ", x, "/", n_default, appendLF = FALSE)
    }
  
    # get number of living trees
    temp_living <- dplyr::filter(default[[x]], i == min(i), type != "dead") %>% 
      nrow()
    
    # get number of trees that died absolute and relative to starting living trees
    dplyr::filter(default[[x]], i != min(i), type == "dead") %>% 
      dplyr::summarise(n = n()) %>% 
      dplyr::mutate(n_rel = n / temp_living)
  })
  
  if (verbose) {
    message("")
  }
  
  # loop through all input data
  died_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      message("\r> Progress (default): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    # get default data
    temp_default <- died_default[[counter_default[x]]]
    
    # get number of living trees
    temp_living <- dplyr::filter(changed[[x]], i == min(i), type != "dead") %>% 
      nrow()
    
    # get number of trees that died absolute and relative to starting living trees
    temp_died <- dplyr::filter(changed[[x]], i != min(i), type == "dead") %>% 
      dplyr::summarise(n = n()) %>% 
      dplyr::mutate(n_rel = n / temp_living)
    
    dplyr::bind_cols(temp_died, temp_default) %>%
      purrr::set_names(c("n.changed", "n_rel.changed", 
                         "n.default", "n_rel.default" )) %>% 
      dplyr::mutate(diff_n = (n.changed - n.default) / n.default, 
                    diff_n_rel = (n_rel.changed - n_rel.default) / n_rel.default)
  })
  
  # add changed parameters as names
  names(died_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  
  died_changed <- dplyr::bind_rows(died_changed, .id = "parameter") %>% 
    dplyr::group_by(parameter) %>% 
    dplyr::summarise(diff_n = mean(suppoRt::replace_infinite(diff_n), 
                                   na.rm = TRUE),
                     diff_n_rel = mean(suppoRt::replace_infinite(diff_n_rel), 
                                       na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  return(died_changed)
}
