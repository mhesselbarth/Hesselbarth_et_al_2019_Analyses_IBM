#### Helper functions #### 
#### SA_local_results_spatial.R ####

calc_pcf <- function(default, changed, 
                     window, r, fast = FALSE,
                     verbose = TRUE, ...) {

  # check correction
  if (length(correction) != 1) {
    stop("Please select one correction.")
  }
    
  # get name of changed to get changed parameters
  names_parameters <- names(changed)
  
  # get length of input for printig
  n_default <- length(default)
  n_changed <- length(changed)

  # repeat counter 1:repetitions of default data for each changed parameter
  counter_default <- rep(x = 1:n_default, times = length(unique(names(changed))))
  
  pcf_default <- purrr::map(seq_along(default), function(x) {

    if (verbose) {
      
      message("\r> Progress (default): ", x, "/", n_default, appendLF = FALSE)
    }
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = default[[x]]$x, y = default[[x]]$y,
                              window = window)

    # calculate pcf 
    if (!fast) {
      temp_sf <- spatstat::pcf(temp_ppp, r = r, ...)
    } 
    
    else {
      temp_sf <- onpoint::estimate_pcf_fast(temp_ppp, r = r, ...)
    }

    # convert to data frame
    tibble::as_tibble(temp_sf) %>% 
      purrr::set_names(c("r", "theo", "pcf"))
  })
  
  if (verbose) {
    
    message("")
  }
  
  # get pcf of changed data
  pcf_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = changed[[x]]$x, y = changed[[x]]$y,
                              window = window)
    
    # calculate pcf
    if (!fast) {
      temp_sf <- spatstat::pcf(temp_ppp, r = r, ...)
    } 
    
    else {
      temp_sf <- onpoint::estimate_pcf_fast(temp_ppp, r = r, ...)
    }
    
    # get default pcf values
    temp_sf_default <- pcf_default[[counter_default[x]]]
    
    # calculate difference to default
    tibble::as_tibble(temp_sf) %>% 
      purrr::set_names(c("r", "theo", "pcf")) %>% 
      dplyr::left_join(y = temp_sf_default, by = "r", 
                       suffix = c(".default", ".changed")) %>% 
      dplyr::mutate(pcf_diff = (pcf.changed - pcf.default) / pcf.default)
  })
  
  # add changed parameters as names
  names(pcf_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  return(pcf_changed)
}

calc_nnd <- function(default, changed, 
                     window, r, verbose = TRUE, ...) {

  # check correction
  if (length(correction) != 1) {
    stop("Please select one correction.")
  }
  
  # get name of changed to get changed parameters
  names_parameters <- names(changed)
  
  # get length of input for printig
  n_default <- length(default)
  n_changed <- length(changed)
  
  # repeat counter 1:repetitions of default data for each changed parameter
  counter_default <- rep(x = 1:n_default, times = length(unique(names(changed))))
  
  nnd_default <- purrr::map(seq_along(default), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (default): ", x, "/", n_default, appendLF = FALSE)
    }
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = default[[x]]$x, y = default[[x]]$y,
                              window = window)
    
    # calculate nnd 
    temp_sf <- spatstat::Gest(temp_ppp, r = r, ...)
    
    # convert to data frame
    tibble::as_tibble(temp_sf) %>% 
      dplyr::select(r, theo, correction) %>% 
      purrr::set_names(c("r", "theo", "nnd"))
  })
  
  if (verbose) {
    message("")
  }
  
  # get pcf of changed data
  nnd_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = changed[[x]]$x, y = changed[[x]]$y,
                              window = window)
    
    # calculate nnd 
    temp_sf <- spatstat::Gest(temp_ppp, r = r, ...)
    
    # get default pcf values
    temp_sf_default <- nnd_default[[counter_default[x]]]

    # calculate difference to default
    tibble::as_tibble(temp_sf) %>% 
      dplyr::select(r, theo, correction) %>% 
      purrr::set_names(c("r", "theo", "nnd")) %>%  
      dplyr::left_join(y = temp_sf_default, by = "r", 
                       suffix = c(".default", ".changed")) %>% 
      dplyr::mutate(nnd_diff = (nnd.changed - nnd.default) / nnd.default)
  })
  
  # add changed parameters as names
  names(nnd_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  return(nnd_changed)
}

calc_kmm <-  function(default, changed, 
                      window, r, fast = FALSE,
                      verbose = TRUE, ...) {
  
  # check correction
  if (length(correction) != 1) {
    stop("Please select one correction.")
  }
  
  # get name of changed to get changed parameters
  names_parameters <- names(changed)
  
  # get length of input for printig
  n_default <- length(default)
  n_changed <- length(changed)
  
  # repeat counter 1:repetitions of default data for each changed parameter
  counter_default <- rep(x = 1:n_default, times = length(unique(names(changed))))
  
  kmm_default <- purrr::map(seq_along(default), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (default): ", x, "/", n_default, appendLF = FALSE)
    }
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = default[[x]]$x, y = default[[x]]$y,
                              window = window)

    # check which points are inside
    temp_inside <- spatstat::inside.owin(x = default[[x]]$x, y = default[[x]]$y,
                                         w = window)

    temp_default <- default[[x]][temp_inside, ]
    
    spatstat::marks(temp_ppp) <- temp_default$dbh
    
    temp_sf <- spatstat::markcorr(temp_ppp, r = r, ...)
    
    # convert to data frame
    tibble::as_tibble(temp_sf) %>% 
      purrr::set_names(c("r", "theo", "kmm"))
  })
  
  if (verbose) {
    
    message("")
  }
  
  # get pcf of changed data
  kmm_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = changed[[x]]$x, y = changed[[x]]$y,
                              window = window)
    
    # check which points are inside
    temp_inside <- spatstat::inside.owin(x = changed[[x]]$x, y = changed[[x]]$y,
                                         w = window)
    
    temp_default <- changed[[x]][temp_inside, ]
    
    spatstat::marks(temp_ppp) <- temp_default$dbh
    
    temp_sf <- spatstat::markcorr(temp_ppp, r = r, ...)
    
    # get default pcf values
    temp_sf_default <- kmm_default[[counter_default[x]]]
    
    # calculate difference to default
    tibble::as_tibble(temp_sf) %>%
      purrr::set_names(c("r", "theo", "kmm")) %>% 
      dplyr::left_join(y = temp_sf_default, by = "r", 
                       suffix = c(".default", ".changed")) %>% 
      dplyr::mutate(kmm_diff = (kmm.changed - kmm.default) / kmm.default)
  })
  
  # add changed parameters as names
  names(kmm_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  return(kmm_changed)
}

calc_clark <- function(default, changed, 
                       window, verbose = TRUE) {

  # get name of changed to get changed parameters
  names_parameters <- names(changed)
  
  # get length of input for printig
  n_default <- length(default)
  n_changed <- length(changed)
  
  # repeat counter 1:repetitions of default data for each changed parameter
  counter_default <- rep(x = 1:n_default, times = length(unique(names(changed))))
  
  ce_default <- purrr::map(seq_along(default), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (default): ", x, "/", n_default, appendLF = FALSE)
    }
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = default[[x]]$x, y = default[[x]]$y,
                              window = window)
    
    # calculate CE 
    spatstat::clarkevans(temp_ppp, correction = "cdf")
  })
  
  if (verbose) {

    message("")
  }
  
  # get pcf of changed data
  ce_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = changed[[x]]$x, y = changed[[x]]$y,
                              window = window)
    
    # calculate CE 
    temp_ce <- spatstat::clarkevans(temp_ppp, correction = "cdf")
    
    # get default pcf values
    temp_ce_default <- ce_default[[counter_default[x]]]
    
    temp_diff <- tibble::tibble(diff_ce = (temp_ce - temp_ce_default) / temp_ce_default)
    
    return(temp_diff)
  })
  
  # add changed parameters as names
  names(ce_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  return(ce_changed)
}

# calc_pcf <- function(data, window, r, 
#                      fast = TRUE,
#                      smaller = NULL, bigger = NULL, ...) {
#   
#   names_input <- names(data)
#   
#   length_input <- length(data)
#   
#   purrr::map_dfr(seq_along(data), function(x) {
#     
#     if (is.null(smaller) && is.null(bigger)) {
#       
#       data <- dplyr::filter(data[[x]], i == max(i))
#     }
#     
#     else{
#       
#       if (!is.null(smaller) && is.null(bigger)) { 
#         
#         data <- dplyr::filter(data[[x]], i == max(i), dbh <= smaller)
#       }
#       
#       else if (is.null(smaller) && !is.null(bigger)) {
#         
#         data <- dplyr::filter(data[[x]], i == max(i), dbh >= bigger)
#       }
#       
#       else if (!is.null(smaller) && !is.null(bigger)) {
#         
#         if (smaller <= bigger) {
#         
#           stop("Not possible to provide 'smaller' and 'bigger'.")
#         }
#         
#         data <- dplyr::filter(data[[x]], i == max(i), 
#                               dbh >= bigger, dbh <= smaller)
#       }
#     }
#     
#     message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
#     
#     data <- spatstat::ppp(x = data$x, y = data$y, 
#                           window = window)
#     
#     if (!fast) {
#       sf <- spatstat::pcf(data, r = r, ...)
#     }
#     
#     else {
#       sf <- shar::estimate_pcf_fast(data, r = r, ...)
#     }
#     
#     tibble::add_column(tibble::as_tibble(sf), 
#                        parameter = names_input[[x]], .before = 1)
#     
#   }, .id = "id") 
# }
# 
# calc_nnd <- function(data, window, r, smaller = NULL, bigger = NULL, ...) {
#   
#   names_input <- names(data)
#   
#   length_input <- length(data)
#   
#   purrr::map_dfr(seq_along(data), function(x) {
#     
#     if (is.null(smaller) && is.null(bigger)) {
#       
#       data <- dplyr::filter(data[[x]], i == max(i))
#     }
#     
#     else{
#       
#       if (!is.null(smaller) && is.null(bigger)) { 
#         
#         data <- dplyr::filter(data[[x]], i == max(i), dbh <= smaller)
#       }
#       
#       else if (is.null(smaller) && !is.null(bigger)) {
#         
#         data <- dplyr::filter(data[[x]], i == max(i), dbh >= bigger)
#       }
#       
#       else if (!is.null(smaller) && !is.null(bigger)) {
#         
#         if (smaller <= bigger) {
#           
#           stop("Not possible to provide 'smaller' and 'bigger'.")
#         }
#         
#         data <- dplyr::filter(data[[x]], i == max(i), 
#                               dbh >= bigger, dbh <= smaller)
#       }
#     }
#     
#     message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
#     
#     data <- spatstat::ppp(x = data$x, y = data$y, 
#                           window = window)
#     
#     sf <- spatstat::Gest(data, r = r, ...)
#     
#     tibble::add_column(tibble::as_tibble(sf), 
#                        parameter = names_input[[x]], .before = 1)
#     
#   }, .id = "id") 
# }
# 
# calc_kmm <- function(data, window, r, smaller = NULL, bigger = NULL, ...) {
# 
#   names_input <- names(data)
# 
#   length_input <- length(data)
# 
#   purrr::map_dfr(seq_along(data), function(x) {
# 
#     if (is.null(smaller) && is.null(bigger)) {
#       
#       data <- dplyr::filter(data[[x]], i == max(i))
#     }
#     
#     else{
#       
#       if (!is.null(smaller) && is.null(bigger)) { 
#         
#         data <- dplyr::filter(data[[x]], i == max(i), dbh <= smaller)
#       }
#       
#       else if (is.null(smaller) && !is.null(bigger)) {
#         
#         data <- dplyr::filter(data[[x]], i == max(i), dbh >= bigger)
#       }
#       
#       else if (!is.null(smaller) && !is.null(bigger)) {
#         
#         if (smaller <= bigger) {
#           
#           stop("Not possible to provide 'smaller' and 'bigger'.")
#         }
#         
#         data <- dplyr::filter(data[[x]], i == max(i), 
#                               dbh >= bigger, dbh <= smaller)
#       }
#     }
#     
#     message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
#   
#     data_ppp <- spatstat::ppp(x = data$x, y = data$y,
#                               window = window)
#     
#     data_inside <- spatstat::inside.owin(x = data$x, y = data$y,
#                                          w = window)
#     
#     data <- data[data_inside, ]
# 
#     spatstat::marks(data_ppp) <- data$dbh
# 
#     sf <- spatstat::markcorr(data_ppp, r = r, ...)
# 
#     tibble::add_column(tibble::as_tibble(sf),
#                        parameter = names_input[[x]], .before = 1)
# 
#   }, .id = "id")
# }
# 
# 
# calc_clark <- function(data, window, smaller = NULL, bigger = NULL, ...) {
#   
#   names_input <- names(data)
#   
#   length_input <- length(data)
#   
#   purrr::map_dfr(seq_along(data), function(x) {
#     
#     if (is.null(smaller) && is.null(bigger)) {
#       
#       data <- dplyr::filter(data[[x]], i == max(i))
#     }
#     
#     else{
#       
#       if (!is.null(smaller) && is.null(bigger)) { 
#         
#         data <- dplyr::filter(data[[x]], i == max(i), dbh <= smaller)
#       }
#       
#       else if (is.null(smaller) && !is.null(bigger)) {
#         
#         data <- dplyr::filter(data[[x]], i == max(i), dbh >= bigger)
#       }
#       
#       else if (!is.null(smaller) && !is.null(bigger)) {
#         
#         if (smaller <= bigger) {
#           
#           stop("Not possible to provide 'smaller' and 'bigger'.")
#         }
#         
#         data <- dplyr::filter(data[[x]], i == max(i), 
#                               dbh >= bigger, dbh <= smaller)
#       }
#     }
#     
#     message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
#     
#     data <- spatstat::ppp(x = data$x, y = data$y, 
#                           window = window)
#     
#     index <- spatstat::clarkevans(data, ...)
#     
#     tibble::add_column(tibble::tibble(clark = index), 
#                        parameter = names_input[[x]], .before = 1)
#     
#   }, .id = "id") 
# }
