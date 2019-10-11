#### helper function model comparison spatial ### 

calc_nnd_comp <- function(data, window, r, correction, verbose = TRUE, ...) {
  
  # get length of input for printig
  n_data <- length(data)
  
  nnd_result <- purrr::map(seq_along(data), function(x) {
    
    if (verbose) {
      message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
    }
    
    temp_data <- dplyr::filter(data[[x]], i == max(i), type != "dead")
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_data$x, y = temp_data$y,
                              window = window)
    
    # calculate nnd 
    temp_sf <- spatstat::Gest(temp_ppp, r = r, correction = correction, ...)
    
    # convert to data frame
    tibble::as_tibble(temp_sf)
  })
  
  if (verbose) {
    message("")
  }
  
  nnd_result <- dplyr::bind_rows(nnd_result) %>% 
    dplyr::group_by(r) %>% 
    dplyr::summarise(theo = mean(theo),
                     fun_mean = mean(!!dplyr::sym(correction)),
                     fun_lo = min(!!dplyr::sym(correction)),
                     fun_hi = max(!!dplyr::sym(correction)))
  
  return(nnd_result)
}

calc_pcf_comp <- function(data, window, verbose = TRUE, ...) {
  
  # get length of input for printig
  n_data <- length(data)
  
  pcf_result <- purrr::map(seq_along(data), function(x) {
    
    if (verbose) {
      message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
    }
    
    # get data of last timestep
    temp_data <- dplyr::filter(data[[x]], i == max(i), type != "dead")
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_data$x, y = temp_data$y,
                              window = window)
    
    temp_sf <- spatstat::pcf(temp_ppp, ...)
    
    # convert to data frame
    tibble::as_tibble(temp_sf) %>% 
      purrr::set_names(c("r", "theo", "pcf"))
  })

  if (verbose) {
    message("")
  }
  
  pcf_result <- dplyr::bind_rows(pcf_result) %>% 
    dplyr::group_by(r) %>% 
    dplyr::summarise(theo = mean(theo),
                     fun_mean = mean(pcf),
                     fun_lo = min(pcf),
                     fun_hi = max(pcf))
  
  return(pcf_result)
}

calc_kmm_comp <- function(data, window, r, verbose = TRUE, ...) {
  
  # get length of input for printig
  n_data <- length(data)

  kmm_data <- purrr::map(seq_along(data), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (default): ", x, "/", n_data, appendLF = FALSE)
    }
    
    # get data of last timestep
    temp_data <- dplyr::filter(data[[x]], i == max(i), type != "dead")
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_data$x, y = temp_data$y,
                              window = window)
    
    # check which points are inside
    temp_inside <- spatstat::inside.owin(x = temp_data$x, y = temp_data$y,
                                         w = window)
    
    temp_data <- temp_data[temp_inside, ]
    
    spatstat::marks(temp_ppp) <- temp_data$dbh
    
    temp_sf <- spatstat::markcorr(temp_ppp, r = r, ...)
    
    # convert to data frame
    tibble::as_tibble(temp_sf) %>% 
      purrr::set_names(c("r", "theo", "kmm"))
  })
  
  if (verbose) {
    message("")
  }
  
  kmm_data <- dplyr::bind_rows(kmm_data) %>% 
    dplyr::group_by(r) %>% 
    dplyr::summarise(theo = mean(theo), 
                     fun_mean = mean(kmm),
                     fun_lo = min(kmm),
                     fun_hi = max(kmm))
  
  return(kmm_data)
}

calc_ci_comp <- function(data, parameters, from = 0, to = 1, verbose = TRUE) {
  
  # get length of input for printig
  n_data <- length(data)
  
  ci_data <- purrr::map(seq_along(data), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (default): ", x, "/", n_data, appendLF = FALSE)
    }
    
    # get data of last timestep
    temp_data <- dplyr::filter(data[[x]], i == max(i), type != "dead")
    
    # check which points are inside
    temp_inside <- spatstat::inside.owin(x = temp_data$x, y = temp_data$y,
                                         w = window)
    
    temp_data <- temp_data[temp_inside, ]
    
    # convert to data.table
    temp_dt <- data.table::as.data.table(temp_data)
    
    # calculate competition
    ci <- rabmp:::rcpp_calculate_ci(matrix = as.matrix(temp_dt[, list(x, y, dbh)]),
                                    alpha = parameters$ci_alpha,
                                    beta = parameters$ci_beta,
                                    max_dist = parameters$ci_max_dist)
    
    ci_density <- density(ci, from = from, to = to, n = 512)
    
    tibble::tibble(x = ci_density$x, y = ci_density$y)
  })
  
  if (verbose) {
    message("")
  }
  
  ci_data <- dplyr::bind_rows(ci_data) %>% 
    dplyr::group_by(x) %>% 
    dplyr::summarise(fun_mean = mean(y),
                     fun_lo = min(y),
                     fun_hi = max(y))
  
  return(ci_data)
}
