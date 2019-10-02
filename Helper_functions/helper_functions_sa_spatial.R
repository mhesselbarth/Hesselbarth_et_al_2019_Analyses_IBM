###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Helper functions SA spatial #### 

#### PCF #### 
calc_pcf_sa_fun <- function(default, changed, window, r, verbose = TRUE, ...) {

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
    
    # get only latest time step of data
    temp_default <- dplyr::filter(default[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_default$x, y = temp_default$y,
                              window = window)

    # calculate pcf 
   temp_sf <- spatstat::pcf(temp_ppp, r = r, ...)

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
    
    # get only latest time step of data
    temp_changed <- dplyr::filter(changed[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
                              window = window)
    
    # calculate pcf
    temp_sf <- spatstat::pcf(temp_ppp, r = r, ...)

    # get default pcf values
    temp_sf_default <- pcf_default[[counter_default[x]]]
    
    # calculate difference to default
    tibble::as_tibble(temp_sf) %>% 
      purrr::set_names(c("r", "theo", "pcf")) %>% 
      dplyr::left_join(y = temp_sf_default, by = "r", 
                       suffix = c("_default", "_changed")) %>% 
      dplyr::mutate(pcf_diff = (pcf_changed - pcf_default) / pcf_default)
  })
  
  # add changed parameters as names
  names(pcf_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  pcf_changed <- dplyr::bind_rows(pcf_changed, .id = "parameter") %>% 
    dplyr::group_by(parameter, r) %>%
    dplyr::summarise(theo = 1,
                     pcf_diff = mean(pcf_diff, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(pcf_changed)
}

calc_pcf_sa_int <- function(default, changed, window, r, verbose = TRUE, ...) {
  
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
    
    # get only latest time step of data
    temp_default <- dplyr::filter(default[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_default$x, y = temp_default$y,
                              window = window)
    
    # calculate pcf
    temp_sf <- spatstat::pcf(temp_ppp, r = r, ...)
    
    # get function of pcf
    temp_fun <- spatstat::as.function.fv(temp_sf)
    
    # get integral
    integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
  })
  
  if (verbose) {
    
    message("")
  }
  
  # get pcf of changed data
  pcf_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    # get only latest time step of data
    temp_changed <- dplyr::filter(changed[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
                              window = window)
    
    # calculate pcf
    temp_sf <- spatstat::pcf(temp_ppp, r = r, ...)

    # get function of pcf
    temp_fun <- spatstat::as.function.fv(temp_sf)
    
    # calculate integral
    temp_pcf <- integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
    
    # get default pcf values
    temp_pcf_default <- pcf_default[[counter_default[x]]]
    
    # get relative difference between temp and default
    tibble::tibble(pcf_diff = (temp_pcf - temp_pcf_default) / temp_pcf_default)
  })
  
  # add changed parameters as names
  names(pcf_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  pcf_changed <- dplyr::bind_rows(pcf_changed, .id = "parameter") %>% 
    dplyr::group_by(parameter) %>%
    dplyr::summarise(pcf_mean = mean(pcf_diff, na.rm = TRUE), 
                     pcf_min = min(pcf_diff, na.rm = TRUE), 
                     pcf_max = max(pcf_diff, na.rm = TRUE), 
                     pcf_sd = sd(pcf_diff, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(pcf_changed)
}

#### Nearest neigbor ####
calc_nnd_sa_fun <- function(default, changed, 
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
    
    temp_default <- dplyr::filter(default[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_default$x, y = temp_default$y,
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
    
    temp_changed <- dplyr::filter(changed[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
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
                       suffix = c("_default", "_changed")) %>% 
      dplyr::mutate(nnd_diff = (nnd_changed - nnd_default) / nnd_default)
  })
  
  # add changed parameters as names
  names(nnd_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  nnd_changed <- dplyr::bind_rows(nnd_changed, .id = "parameter") %>% 
    dplyr::group_by(parameter, r) %>% 
    dplyr::summarise(theo_default = mean(theo_default, na.rm = TRUE), 
                     theo_changed = mean(theo_changed, na.rm = TRUE), 
                     nnd_diff = mean(nnd_diff, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(nnd_changed)
}

calc_nnd_sa_int <- function(default, changed, 
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
    
    temp_default <- dplyr::filter(default[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_default$x, y = temp_default$y,
                              window = window)
    
    # calculate nnd 
    temp_sf <- spatstat::Gest(temp_ppp, r = r, ...)
    
    # get function of nnd
    temp_fun <- spatstat::as.function.fv(temp_sf)
    
    # get integral
    integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
  })
  
  if (verbose) {
    message("")
  }
  
  # get pcf of changed data
  nnd_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    temp_changed <- dplyr::filter(changed[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
                              window = window)
    
    # calculate nnd 
    temp_sf <- spatstat::Gest(temp_ppp, r = r, ...)
    
    # convert to nnd function
    temp_fun <- spatstat::as.function.fv(temp_sf)
    
    # calculate integral
    temp_nnd <- integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
    
    # get default nnd values
    temp_nnd_default <- nnd_default[[counter_default[x]]]
    
    # get relative difference between temp and default
    tibble::tibble(nnd_diff = (temp_nnd - temp_nnd_default) / temp_nnd_default)
  })
  
  # add changed parameters as names
  names(nnd_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  nnd_changed <- dplyr::bind_rows(nnd_changed, .id = "parameter") %>% 
    dplyr::group_by(parameter) %>% 
    dplyr::summarise(nnd_mean = mean(nnd_diff, na.rm = TRUE), 
                     nnd_min = min(nnd_diff, na.rm = TRUE), 
                     nnd_max = max(nnd_diff, na.rm = TRUE), 
                     nnd_sd = sd(nnd_diff, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(nnd_changed)
}

#### Mark-correlation function ####
calc_kmm_sa_fun <- function(default, changed, 
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
    
    temp_default <- dplyr::filter(default[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = default[[x]]$x, y = default[[x]]$y,
                              window = window)

    # check which points are inside
    temp_inside <- spatstat::inside.owin(x = temp_default$x, y = temp_default$y,
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
    
    temp_changed <- dplyr::filter(changed[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
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
                       suffix = c("_default", "_changed")) %>% 
      dplyr::mutate(kmm_diff = (kmm_changed - kmm_default) / kmm_default)
  })
  
  # add changed parameters as names
  names(kmm_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  kmm_changed <- dplyr::bind_rows(kmm_changed, .id = "parameter") %>% 
    dplyr::group_by(parameter, r) %>% 
    dplyr::summarise(theo = 1, 
                     kmm_diff = mean(kmm_diff, na.rm = TRUE)) %>% 
    dplyr::ungroup() 
  
  return(kmm_changed)
}

calc_kmm_sa_int <- function(default, changed, 
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
    
    temp_default <- dplyr::filter(default[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = default[[x]]$x, y = default[[x]]$y,
                              window = window)
    
    # check which points are inside
    temp_inside <- spatstat::inside.owin(x = temp_default$x, y = temp_default$y,
                                         w = window)
    
    temp_default <- default[[x]][temp_inside, ]
    
    spatstat::marks(temp_ppp) <- temp_default$dbh
    
    temp_sf <- spatstat::markcorr(temp_ppp, r = r, ...)
    
    # get function of kmmr
    temp_fun <- spatstat::as.function.fv(temp_sf)
    
    # get integral
    integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
  })
  
  if (verbose) {
    
    message("")
  }
  
  # get pcf of changed data
  kmm_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    temp_changed <- dplyr::filter(changed[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
                              window = window)
    
    # check which points are inside
    temp_inside <- spatstat::inside.owin(x = changed[[x]]$x, y = changed[[x]]$y,
                                         w = window)
    
    temp_default <- changed[[x]][temp_inside, ]
    
    spatstat::marks(temp_ppp) <- temp_default$dbh
    
    temp_sf <- spatstat::markcorr(temp_ppp, r = r, ...)
    
    # calculate integral
    temp_kmm <- integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
    
    # get default kmm values
    temp_kmm_default <- kmm_default[[counter_default[x]]]
    
    # get relative difference between temp and default
    tibble::tibble(kmm_diff = (temp_kmm - temp_kmm_default) / temp_kmm_default)
  })
  
  # add changed parameters as names
  names(kmm_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  kmm_changed <- dplyr::bind_rows(kmm_changed, .id = "parameter") %>% 
    dplyr::group_by(parameter) %>% 
    dplyr::summarise(kmm_mean = mean(kmm_diff, na.rm = TRUE), 
                     kmm_min = min(kmm_diff, na.rm = TRUE), 
                     kmm_max = max(kmm_diff, na.rm = TRUE), 
                     kmm_sd = sd(kmm_diff, na.rm = TRUE)) %>% 
    dplyr::ungroup() 
  
  return(kmm_changed)
}

#### Clark and Evans Index ####
calc_clark_sa <- function(default, changed, 
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
    
    temp_default <- dplyr::filter(default[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_default$x, y = temp_default$y,
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
    
    temp_changed <- dplyr::filter(changed[[x]], i == max(i))
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
                              window = window)
    
    # calculate CE 
    temp_ce <- spatstat::clarkevans(temp_ppp, correction = "cdf")
    
    # get default pcf values
    temp_ce_default <- ce_default[[counter_default[x]]]
    
    tibble::tibble(ce_diff = (temp_ce - temp_ce_default) / 
                     temp_ce_default)
  })
  
  # add changed parameters as names
  names(ce_changed) <- names_parameters
  
  if (verbose) {
    message("")
  }
  
  ce_changed <- dplyr::bind_rows(ce_changed, .id = "parameter") %>% 
    dplyr::group_by(parameter) %>% 
    dplyr::summarise(ce_mean = mean(ce_diff, na.rm = TRUE), 
                     ce_min = min(ce_diff, na.rm = TRUE), 
                     ce_max = max(ce_diff, na.rm = TRUE), 
                     ce_sd = sd(ce_diff, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(ce_changed)
}
