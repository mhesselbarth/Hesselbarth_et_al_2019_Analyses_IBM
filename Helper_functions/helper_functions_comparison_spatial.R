#### helper function model comparison spatial 3### 

calc_nnd <- function(data, window, r, verbose = TRUE, ...) {
  
  # get length of input for printig
  n_data <- length(data)
  
  nnd_result <- purrr::map(seq_along(data), function(x, ...) {
    
    if (verbose) {
      message("\r> Progress: ", x, "/", n_data, appendLF = FALSE)
    }
    
    temp_data <- dplyr::filter(data[[x]], i == max(i), type != "dead")
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_data$x, y = temp_data$y,
                              window = window)
    
    # calculate nnd 
    temp_sf <- spatstat::Gest(temp_ppp, r = r, ...)
    
    # convert to data frame
    tibble::as_tibble(temp_sf)
  })
  
  if (verbose) {
    message("")
  }
  
  return(nnd_result)
}

calc_pcf <- function(data, 
                     window, r, fast = FALSE,
                     verbose = TRUE, ...) {
  
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
    
    # calculate pcf 
    if (!fast) {
      temp_sf <- spatstat::pcf(temp_ppp, r = r, ...)
    } 
    
    else {
      temp_sf <- onpoint::estimate_pcf_fast(temp_ppp, r = r, ...)
    }
    
    # convert to data frame
    tibble::as_tibble(temp_sf)
  })
  
  if (verbose) {
    message("")
  }
  
  return(pcf_result)
}

calc_kmm <-  function(data, window, r, verbose = TRUE, ...) {
  
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
  
  return(kmm_data)
}
