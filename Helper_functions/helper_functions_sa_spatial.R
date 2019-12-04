###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Helper functions SA spatial #### 

#### PCF #### 
calc_pcf_sa_int <- function(default, changed, window, nsim, verbose = TRUE, ...) {
  
  arguments <- list(...)
  
  if (!"correction" %in% names(arguments)) {
    stop("Please provide correction argument")
  }
  
  else {
    if (length(arguments$correction) > 1) {
      stop("Please provide only one correction argument")
    }
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
    temp_default <- dplyr::filter(default[[x]], i == max(i), type != "dead")
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_default$x, y = temp_default$y,
                              window = window)
    
    # calculate pcf
    temp_env <- spatstat::envelope(Y = temp_ppp, fun = "pcf", nsim = nsim, 
                                   verbose = FALSE,
                                   funargs = arguments)
    
    onpoint::summarise_envelope(x = temp_env)
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
    temp_changed <- dplyr::filter(changed[[x]], i == max(i), type != "dead")
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
                              window = window)
    
    # calculate pcf
    temp_env <- spatstat::envelope(Y = temp_ppp, fun = "pcf", nsim = nsim, 
                                   verbose = FALSE,
                                   funargs = arguments)
    
    temp_pcf <- onpoint::summarise_envelope(x = temp_env)
    
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

calc_nnd_sa_int <- function(default, changed, 
                            window, verbose = TRUE, ...) {
  
  arguments <- list(...)
  
  if (!"correction" %in% names(arguments)) {
    stop("Please provide correction argument")
  }
  
  else {
    if (length(arguments$correction) > 1) {
      stop("Please provide only one correction argument")
    }
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
    
    temp_default <- dplyr::filter(default[[x]], i == max(i), type != "dead")
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_default$x, y = temp_default$y,
                              window = window)
    
    # # calculate nnd 
    # temp_sf <- spatstat::Gest(temp_ppp, r = r, ...)
    temp_sf <- spatstat::Gest(temp_ppp, ...) %>% 
      tibble::as_tibble() %>%
      dplyr::select(r, theo, arguments$correction) %>% 
      purrr::set_names(c("r", "theo", "nnd")) 
    
    # # get function of nnd
    # temp_fun <- spatstat::as.function.fv(temp_sf)
  
    # # get integral
    # integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
    MESS::auc(temp_sf[, 1, drop = TRUE], temp_sf[, 3, drop = TRUE], 
              from = min(temp_sf$r), to = max(temp_sf$r))
  })
  
  if (verbose) {
    message("")
  }
  
  # get pcf of changed data
  nnd_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    temp_changed <- dplyr::filter(changed[[x]], i == max(i), type != "dead")
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
                              window = window)
    
    # # calculate nnd 
    # temp_sf <- spatstat::Gest(temp_ppp, r = r, ...)
    temp_sf <- spatstat::Gest(temp_ppp, ...) %>% 
      tibble::as_tibble() %>%
      dplyr::select(r, theo, arguments$correction) %>% 
      purrr::set_names(c("r", "theo", "nnd")) 
    
    # # get function of nnd
    # temp_fun <- spatstat::as.function.fv(temp_sf)
    
    # # get integral
    # integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
    temp_nnd <- MESS::auc(temp_sf[, 1, drop = TRUE], temp_sf[, 3, drop = TRUE], 
                          from = min(temp_sf$r), to = max(temp_sf$r))
    
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
calc_kmm_sa_int <- function(default, changed, 
                            window, verbose = TRUE, ...) {
  
  arguments <- list(...)
  
  if (!"correction" %in% names(arguments)) {
    stop("Please provide correction argument")
  }
  
  else {
    if (length(arguments$correction) > 1) {
      stop("Please provide only one correction argument")
    }
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
    
    temp_default <- dplyr::filter(default[[x]], i == max(i), type != "dead")
    
    # check which points are inside
    temp_inside <- spatstat::inside.owin(x = temp_default$x, y = temp_default$y,
                                         w = window)
    
    temp_default <- default[[x]][temp_inside, ]
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = default[[x]]$x, y = default[[x]]$y,
                              window = window)
    
    spatstat::marks(temp_ppp) <- temp_default$dbh
    
    # temp_sf <- spatstat::markcorr(temp_ppp, r = r, ...)
    temp_sf <- spatstat::markcorr(temp_ppp, ...) %>% 
      tibble::as_tibble() %>% 
      purrr::set_names(c("r", "theo", "kmm"))
    
    # # get function of kmmr
    # temp_fun <- spatstat::as.function.fv(temp_sf)
    
    # # get integral
    # integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
    MESS::auc(temp_sf$r, temp_sf$kmm, 
              from = min(temp_sf$r), to = max(temp_sf$r))
  })
  
  if (verbose) {
    
    message("")
  }
  
  # get pcf of changed data
  kmm_changed <- purrr::map(seq_along(changed), function(x) {
    
    if (verbose) {
      
      message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
    }
    
    temp_changed <- dplyr::filter(changed[[x]], i == max(i), type != "dead")
    
    # check which points are inside
    temp_inside <- spatstat::inside.owin(x = changed[[x]]$x, y = changed[[x]]$y,
                                         w = window)
    
    temp_changed <- changed[[x]][temp_inside, ]
    
    # convert to ppp
    temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
                              window = window)
    
    spatstat::marks(temp_ppp) <- temp_changed$dbh
    
    # temp_sf <- spatstat::markcorr(temp_ppp, r = r, ...)
    temp_sf <- spatstat::markcorr(temp_ppp, ...) %>% 
      tibble::as_tibble() %>% 
      purrr::set_names(c("r", "theo", "kmm"))
    
    # # calculate integral
    # temp_kmm <- integrate(f = temp_fun, lower = min(r), upper = max(r))[[1]]
    temp_kmm <- MESS::auc(temp_sf$r, temp_sf$kmm, 
                          from = min(temp_sf$r), to = max(temp_sf$r))
    
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
# calc_clark_sa <- function(default, changed, 
#                           window, verbose = TRUE) {
# 
#   # get name of changed to get changed parameters
#   names_parameters <- names(changed)
#   
#   # get length of input for printig
#   n_default <- length(default)
#   n_changed <- length(changed)
#   
#   # repeat counter 1:repetitions of default data for each changed parameter
#   counter_default <- rep(x = 1:n_default, times = length(unique(names(changed))))
#   
#   ce_default <- purrr::map(seq_along(default), function(x) {
#     
#     if (verbose) {
#       
#       message("\r> Progress (default): ", x, "/", n_default, appendLF = FALSE)
#     }
#     
#     temp_default <- dplyr::filter(default[[x]], i == max(i), type != "dead")
#     
#     # convert to ppp
#     temp_ppp <- spatstat::ppp(x = temp_default$x, y = temp_default$y,
#                               window = window)
#     
#     # calculate CE 
#     spatstat::clarkevans(temp_ppp, correction = "cdf")
#   })
#   
#   if (verbose) {
# 
#     message("")
#   }
#   
#   # get pcf of changed data
#   ce_changed <- purrr::map(seq_along(changed), function(x) {
#     
#     if (verbose) {
#       
#       message("\r> Progress (changed): ", x, "/", n_changed, appendLF = FALSE)
#     }
#     
#     temp_changed <- dplyr::filter(changed[[x]], i == max(i), type != "dead")
#     
#     # convert to ppp
#     temp_ppp <- spatstat::ppp(x = temp_changed$x, y = temp_changed$y,
#                               window = window)
#     
#     # calculate CE 
#     temp_ce <- spatstat::clarkevans(temp_ppp, correction = "cdf")
#     
#     # get default pcf values
#     temp_ce_default <- ce_default[[counter_default[x]]]
#     
#     tibble::tibble(ce_diff = (temp_ce - temp_ce_default) / 
#                      temp_ce_default)
#   })
#   
#   # add changed parameters as names
#   names(ce_changed) <- names_parameters
#   
#   if (verbose) {
#     message("")
#   }
#   
#   ce_changed <- dplyr::bind_rows(ce_changed, .id = "parameter") %>% 
#     dplyr::group_by(parameter) %>% 
#     dplyr::summarise(ce_mean = mean(ce_diff, na.rm = TRUE), 
#                      ce_min = min(ce_diff, na.rm = TRUE), 
#                      ce_max = max(ce_diff, na.rm = TRUE), 
#                      ce_sd = sd(ce_diff, na.rm = TRUE)) %>% 
#     dplyr::ungroup()
#   
#   return(ce_changed)
# }
