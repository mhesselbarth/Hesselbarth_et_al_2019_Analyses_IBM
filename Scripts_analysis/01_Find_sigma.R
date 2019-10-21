###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Systematic find sigma ####

#### Import libraries ####

# load packages #
source("Helper_functions/helper_functions_setup.R")
source("Helper_functions/helper_functions_abiotic_conditions.R")

# initialse function potential growth #
fun_potential <- function(dbh, assymp, rate, infl) {     
  
  growth <- assymp * rate * infl * exp(-rate * dbh) * 
    (1 - exp(-rate * dbh)) ^ (infl - 1)
  
  return(growth)
}

# initialse function #
fun_actual <- function(df, par) { 
  
  data_matrix <- as.matrix(df[, c("x", "y", "dbh_99", "growth_pot", "abiotic")])
  
  growth_modelled <- rabmp:::rcpp_calculate_actual_abiotic(matrix = data_matrix, 
                                                           alpha = par[1], 
                                                           beta = par[2],
                                                           mod = 1,
                                                           gamma = par[3],
                                                           max_dist = 30)
  
  difference <- sum(abs(df$growth_full - growth_modelled))
  
  return(difference)
}

explore_sigma <- function(data_1999,
                          data_2013,
                          id_top,
                          sigma,
                          parameters,
                          probs_id,
                          probs,
                          plot_area,
                          years,
                          save_each) {

  # classify abiotic #
  
  # filter data using threshold sapling/adult #
  data_abiotic <- spatstat::subset.ppp(data_1999, dbh_99 > 10 & species == "beech")

  # get intensity
  habitat <- spatstat::density.ppp(data_abiotic, at = "pixel",
                                   weights = data_abiotic$marks$dbh_99,
                                   dimyx = c(645, 609),
                                   kernel = "epanechnikov", sigma = sigma)

  # number of rows added at edges #
  n_rows <- 3

  # convert to raster and add padding #
  habitat <- tibble::as_tibble(habitat)

  habitat <- raster::rasterFromXYZ(habitat)

  habitat <- landscapemetrics::pad_raster(habitat,
                                          pad_raster_value = NA,
                                          pad_raster_cells = n_rows)[[1]]

  # add subsequently number of rows at edge #
  for (i in 1:n_rows) {

    # get all NA cells #
    cells_na <- raster::Which(is.na(habitat),
                              cells = TRUE)

    # get neigbors of NA cells #
    neighbours <- raster::adjacent(x = habitat,
                                   cells = cells_na,
                                   directions = 8)

    # get mean of all neighboring cells
    neighbours_value <- tibble::tibble(from = neighbours[, 1],
                                       to = neighbours[, 2],
                                       x_from = habitat[neighbours[, 1]],
                                       x_to = habitat[neighbours[, 2]])

    neighbours_value <- dplyr::group_by(neighbours_value, from)

    neighbours_value <- dplyr::summarise(neighbours_value,
                                         x = mean(x_to, na.rm = TRUE))

    neighbours_value <- dplyr::filter(neighbours_value, !is.na(x))

    # add values
    habitat[neighbours_value$from] <- neighbours_value$x
  }

  # scale value to -1 to 1 #
  habitat$scaled <- scales::rescale(raster::values(habitat),
                                    to = c(-1, 1), na.rm = TRUE)

  # set names #
  names(habitat) <- c("absolute", "scaled")
  
  # fit parameters #
  data_2013$abiotic <- rabmp::extract_abiotic(data = data.table::data.table(x = data_2013$x, 
                                                                            y = data_2013$y), 
                                              abiotic = habitat$scaled)

  # set starting functions adapted from Pommerening, A., Maleki, K., 2014. #
  # Differences between competition kernels and traditional size-ratio based #
  # competition indices used in forest ecology. For. Ecol. Manage. 331, 135-143. #
  start_values_actual <- c(parameters$ci_alpha, 
                           parameters$ci_beta, 
                           parameters$growth_abiotic)
  
  # fit fun #
  fitted_fun_actual <- optim(par = start_values_actual,
                             fn = fun_actual, 
                             df = dplyr::filter(data_2013, 
                                                !id %in% id_top), 
                             method = "BFGS",
                             control = list(trace = FALSE))
  
  parameters$ci_alpha <- fitted_fun_actual$par[[1]]
  parameters$ci_beta <- fitted_fun_actual$par[[2]]
  parameters$growth_abiotic <- fitted_fun_actual$par[[3]]
  
  # run model #
  
  data_1999 <- tibble::as_tibble(data_1999)
  data_1999 <- dplyr::select(data_1999, 
                        x, y, species, dbh_99, type)
  data_1999 <- dplyr::filter(data_1999, species == "beech")
  data_1999 <- dplyr::mutate(data_1999, type = "adult")
  data_1999 <- dplyr::select(data_1999, -species)
  data_1999 <- rabmp::prepare_data(data_1999,
                                   x = "x", y = "y", type = "type", dbh = "dbh_99")

  result <-  rabmp::run_model_abiotic(data = data_1999,
                                      parameters = parameters,
                                      abiotic = habitat$scaled,
                                      probs = probs[[probs_id]],
                                      plot_area = plot_area,
                                      years = years,
                                      save_each = save_each,
                                      verbose = FALSE)
  
  result$sigma <- sigma
  result$probs <- paste(probs[[probs_id]], collapse = "/")
  
  result$ci_alpha <- parameters$ci_alpha
  result$ci_beta <- parameters$ci_beta
  result$growth_abiotic <- parameters$growth_abiotic

  return(result)
}

#### Import data ####
beech_1999_ppp <- readr::read_rds("Data/Input/beech_1999_ppp.rds")

beech_2007_ppp <- readr::read_rds("Data/Input/beech_2007_ppp.rds")
beech_2007_sapling_ppp <- readr::read_rds("Data/Input/beech_2007_sapling_ppp.rds")
beech_2007_adult_ppp <- readr::read_rds("Data/Input/beech_2007_adult_ppp.rds")

beech_2013_ppp <- readr::read_rds("Data/Input/beech_2013_ppp.rds")
beech_2013_sapling_ppp <- readr::read_rds("Data/Input/beech_2013_sapling_ppp.rds")
beech_2013_adult_ppp <- readr::read_rds("Data/Input/beech_2013_adult_ppp.rds")

parameters_fitted_abiotic <- rabmp::read_parameters("Data/Input/parameters_fitted_abiotic.txt",
                                                    sep = ";")

model_runs_mort <- readr::read_rds("Data/Output/model_runs/model_runs_sigma.rds")

#### Pre-process data ####
beech_2007_df <- tibble::as_tibble(beech_2007_ppp)
beech_2013_df <- tibble::as_tibble(beech_2013_ppp)

plot_area <- beech_1999_ppp$window

parameters_fitted_abiotic$growth_abiotic <- 0

probs <- list(c(0.1, 0.9), c(0.25, 0.75))

combined_ps <- suppoRt::expand_grid_unique(x = 1:2,
                                           y = seq(from = 5, to = 75, by = 5))

probs_id <- combined_ps[, 1]
sigma <- combined_ps[, 2]

years <- 50
save_each <- 50

breaks <- seq(from = 0, to = max(beech_2013_df$dbh_13) + 10, by = 10)


# filter data and calculate mean dbh growth #
beech_2013_df <- dplyr::filter(beech_2013_df, 
                               !is.na(dbh_99), 
                               !is.na(dbh_13),
                               type == "living", 
                               inside_fence == 0) %>% 
  dplyr::mutate(growth_mean = (growth_07 + growth_13) / 2,
                growth_full = (dbh_13 - dbh_99) / 14) %>% 
  dplyr::filter(growth_mean >= 0,
                growth_full >= 0)

# classify into dbh classes and get only trees with highest growth #
breaks <- seq(from = 0, to = max(beech_2013_df$dbh_13) + 10, by = 10)

id_top <- dplyr::mutate(beech_2013_df,
                        dbh_class = cut(dbh_13, breaks = breaks)) %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::top_n(n = n() * 0.05, wt = growth_full) %>%
  dplyr::pull(id)

# calculate potential growth #
beech_2013_df$growth_pot <- fun_potential(dbh = beech_2013_df$dbh_99, 
                                          assymp = parameters_fitted_abiotic$growth_assymp, 
                                          rate = parameters_fitted_abiotic$growth_rate, 
                                          infl = parameters_fitted_abiotic$growth_infl)

#### Run systematic sigma exploratation ####
model_runs_mort <- suppoRt::submit_to_cluster(explore_sigma,
                                              sigma = sigma,
                                              probs_id = probs_id,
                                              const = list(data_1999 = beech_1999_ppp,
                                                           data_2013 = beech_2013_df,
                                                           id_top = id_top,
                                                           parameters = parameters_fitted_abiotic,
                                                           probs = probs,
                                                           plot_area = plot_area,
                                                           years = years,
                                                           save_each = save_each),
                                              export = list(fun_actual = fun_actual),
                                              n_jobs = nrow(combined_ps),
                                              template = list(job_name = "systematic",
                                                              walltime = "02:00:00",
                                                              queue = "medium",
                                                              service = "short",
                                                              mem_cpu = "2048",
                                                              log_file = "systematic.log"))

suppoRt::save_rds(object = model_runs_mort,
                  filename = "model_runs_sigma.rds",
                  path = "Data/Output/model_runs/", 
                  overwrite = overwrite)

##### DBH dist ####
by <- 10

model_runs_dbh <- purrr::map(seq_along(model_runs_mort), function(x) {
  
  message("\r> Progress: ", x, "/", length(model_runs_mort), appendLF = FALSE)
  
  
  # get living trees of last time step
  temp_data <- dplyr::filter(model_runs_mort[[x]], 
                             i == max(i), type != "dead")
  
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

model_runs_dbh <- dplyr::bind_rows(model_runs_dbh, .id = "id") %>% 
  dplyr::mutate(id = as.integer(id))

dbh_dist_2007 <- dplyr::filter(beech_2007_df, 
                               type != "dead", !is.na(dbh_07), inside_fence == 0) %>% 
  dplyr::mutate(dbh_class = cut(dbh_07, breaks = seq(from = 0, 
                                                     to = max(dbh_07) + by, 
                                                     by = by), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n = dplyr::n(), 
                   n_rel = n / nrow(.)) %>% 
  dplyr::mutate(data_type = "Field data 2007")

dbh_dist_2013 <- dplyr::filter(beech_2013_df, 
                               type != "dead", !is.na(dbh_13), inside_fence == 0) %>% 
  dplyr::mutate(dbh_class = cut(dbh_13, breaks = seq(from = 0, 
                                                     to = max(dbh_13) + by, 
                                                     by = by), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n = dplyr::n(), 
                   n_rel = n / nrow(.)) %>% 
  dplyr::mutate(data_type = "Field data 2013")

dbh_dist_field <- dplyr::bind_rows(dbh_dist_2007,
                                   dbh_dist_2013)

ggplot(data = model_runs_dbh) +
  geom_bar(data = dbh_dist_field ,
           aes(x = dbh_class, y = n_rel * 100, fill = data_type),
           position = position_dodge(), stat = "identity") +
  geom_point(aes(x = dbh_class, y = n_rel * 100), pch = "-", 
             size = 5, col = "red") +
  scale_fill_viridis_d(name = "", option = "D") +
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 10)) +
  facet_wrap( ~ id) + 
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

model_runs_dbh_filtered <- dplyr::filter(model_runs_dbh, !id %in% c(1, 2, 5, 6, 
                                                                    7, 8, 9, 10, 
                                                                    11, 12, 16, 
                                                                    17, 18, 19, 
                                                                    20, 21, 22, 
                                                                    23, 24, 25, 
                                                                    26, 27, 28, 
                                                                    29, 30))

ggplot(data = model_runs_dbh_filtered) +
  geom_bar(data = dbh_dist_field ,
           aes(x = dbh_class, y = n_rel * 100, fill = data_type),
           position = position_dodge(), stat = "identity") +
  geom_point(aes(x = dbh_class, y = n_rel * 100), pch = "-", 
             size = 5, col = "red") +
  scale_fill_viridis_d(name = "", option = "D") +
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 10)) +
  facet_wrap( ~ id) + 
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

id_dbh <- unique(model_runs_dbh_filtered$id)

#### Nearest-neighbor distribution function ####
r_nnd <- seq(from = 0, to = 10, length.out = 525)
correction_nnd <- "km"

model_runs_nnd <- purrr::map(seq_along(model_runs_mort), function(x) {
  
  message("\r> Progress: ", x, "/", length(model_runs_mort), appendLF = FALSE)
  
  # get data of last timestep
  temp_data <- dplyr::filter(model_runs_mort[[x]], i == max(i), type != "dead")
  
  saplings <- dplyr::filter(temp_data, type == "sapling")
  
  adults <- dplyr::filter(temp_data, type == "adult")
  
  # convert to ppp
  saplings_ppp <- spatstat::ppp(x = saplings$x, y = saplings$y,
                                window = plot_area)
  
  adults_ppp <- spatstat::ppp(x = adults$x, y = adults$y,
                              window = plot_area)
  
  saplings_sf <- spatstat::Gest(saplings_ppp, correction = correction_nnd, 
                               r = r_nnd) %>% 
    tibble::as_tibble() %>% 
    dplyr::select(r, theo, correction_nnd) %>% 
    purrr::set_names(c("r", "theo", "nnd")) %>% 
    dplyr::mutate(size = "saplings")
  
  adults_sf <- spatstat::Gest(adults_ppp, correction = correction_nnd, 
                              r = r_nnd) %>% 
    tibble::as_tibble() %>% 
    dplyr::select(r, theo, correction_nnd) %>% 
    purrr::set_names(c("r", "theo", "nnd")) %>% 
    dplyr::mutate(size = "adults")
  
  dplyr::bind_rows(saplings_sf,
                   adults_sf)
})

model_runs_nnd <- dplyr::bind_rows(model_runs_nnd, .id = "id") %>% 
  dplyr::mutate(id = as.integer(id))

nnd_2007_sapling <- spatstat::Gest(beech_2007_sapling_ppp, 
                                   r = r_nnd, correction = correction_nnd) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction_nnd) %>%
  purrr::set_names(c("r", "theo", "nnd")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size = "saplings")

nnd_2007_adult <- spatstat::Gest(beech_2007_adult_ppp, 
                                 r = r_nnd, correction = correction_nnd) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction_nnd) %>% 
  purrr::set_names(c("r", "theo", "nnd")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size = "adults")

nnd_2013_sapling <- spatstat::Gest(beech_2013_sapling_ppp, 
                                   r = r_nnd, correction = correction_nnd) %>% 
  tibble::as_tibble() %>%
  dplyr::select(r, theo, correction_nnd) %>%
  purrr::set_names(c("r", "theo", "nnd")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size = "saplings")

nnd_2013_adult <- spatstat::Gest(beech_2013_adult_ppp, 
                                r = r_nnd, correction = correction_nnd) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction_nnd) %>%
  purrr::set_names(c("r", "theo", "nnd")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size = "adults")

nnd_overall_field <- dplyr::bind_rows(nnd_2007_sapling,
                                      nnd_2007_adult,
                                      nnd_2013_sapling,
                                      nnd_2013_adult) %>% 
  dplyr::mutate(data_type_field = factor(data_type_field, 
                                         levels = c("Field data 2007",
                                                    "Field data 2013")),
                size_field = factor(size, levels = c("saplings", "adults")))

ggplot(data = model_runs_nnd) + 
  geom_line(aes(x = r, y = nnd, col = factor(id))) + 
  geom_line(data = nnd_overall_field, 
            aes(x = r, y = nnd, linetype = data_type_field)) +
  facet_wrap(~ size) +
  scale_color_viridis_d() +
  theme_classic()

model_runs_nnd_filtered <- dplyr::filter(model_runs_nnd, !id %in% c(1, 2, 16,
                                                                    6, 7, 8,11,
                                                                    18, 19, 20,
                                                                    21, 22, 23, 
                                                                    24, 25, 26, 
                                                                    27,28, 29, 
                                                                    30))

ggplot(data = model_runs_nnd_filtered) + 
  geom_line(aes(x = r, y = nnd, col = factor(id))) + 
  geom_line(data = nnd_overall_field, 
            aes(x = r, y = nnd, linetype = data_type_field)) +
  facet_wrap(~ size) +
  scale_color_viridis_d() +
  theme_classic()

id_nnd <- unique(model_runs_nnd_filtered$id)

#### Pair correlation function #####
r_pcf <- seq(from = 0, to = 50, length.out = 525)
correction_pcf <- "Ripley"
stoyan_pcf <- 0.25
divisor_pcf <- "d"

model_runs_pcf <- purrr::map(seq_along(model_runs_mort), function(x) {
  
  message("\r> Progress: ", x, "/", length(model_runs_mort), appendLF = FALSE)
  
  # get data of last timestep
  temp_data <- dplyr::filter(model_runs_mort[[x]], i == max(i), type != "dead")
  
  saplings <- dplyr::filter(temp_data, type == "sapling")
  
  adults <- dplyr::filter(temp_data, type == "adult")
  
  # convert to ppp
  saplings_ppp <- spatstat::ppp(x = saplings$x, y = saplings$y,
                                window = plot_area)
  
  adults_ppp <- spatstat::ppp(x = adults$x, y = adults$y,
                              window = plot_area)
  
  saplings_sf <- spatstat::pcf(saplings_ppp, 
                               r = r_pcf, correction = correction_pcf,
                               divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
    tibble::as_tibble() %>% 
    purrr::set_names(c("r", "theo", "pcf")) %>% 
    dplyr::mutate(size = "saplings")
  
  adults_sf <- spatstat::pcf(adults_ppp, 
                             r = r_pcf, correction = correction_pcf,
                             divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
    tibble::as_tibble() %>% 
    purrr::set_names(c("r", "theo", "pcf")) %>% 
    dplyr::mutate(size = "adults")
  
  dplyr::bind_rows(saplings_sf,
                   adults_sf)
})

model_runs_pcf <- dplyr::bind_rows(model_runs_pcf, .id = "id") %>% 
  dplyr::mutate(id = as.integer(id))

pcf_2007_sapling <- spatstat::pcf(beech_2007_sapling_ppp, 
                                  r = r_pcf, correction = correction_pcf, 
                                  divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size = "saplings")

pcf_2007_adult <- spatstat::pcf(beech_2007_adult_ppp, 
                                r = r_pcf, correction = correction_pcf, 
                                divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type_field = "Field data 2007", 
                size = "adults")

pcf_2013_sapling <- spatstat::pcf(beech_2013_sapling_ppp, 
                                  r = r_pcf, correction = correction_pcf, 
                                  divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size = "saplings")

pcf_2013_adult <- spatstat::pcf(beech_2013_adult_ppp, 
                                r = r_pcf, correction = correction_pcf, 
                                divisor = divisor_pcf, stoyan = stoyan_pcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "pcf")) %>% 
  dplyr::mutate(data_type_field = "Field data 2013", 
                size = "adults")

pcf_overall_field <- dplyr::bind_rows(pcf_2007_sapling,
                                      pcf_2007_adult,
                                      pcf_2013_sapling,
                                      pcf_2013_adult) %>% 
  dplyr::mutate(data_type_field = factor(data_type_field, 
                                         levels = c("Field data 2007",
                                                    "Field data 2013")),
                size_field = factor(size, levels = c("saplings", "adults")))

ggplot(model_runs_pcf) + 
  geom_line(aes(x = r, y = pcf, col = factor(id))) + 
  geom_line(data = pcf_overall_field, 
            aes(x = r, y = pcf, linetype = data_type_field)) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ size) +
  scale_color_viridis_d() +
  theme_classic()

model_runs_pcf_filtered <- dplyr::filter(model_runs_pcf, id %in% 
                                           unique(model_runs_dbh_filtered$id))

ggplot(model_runs_pcf_filtered) + 
  geom_line(aes(x = r, y = pcf, col = factor(id))) + 
  geom_line(data = pcf_overall_field, 
            aes(x = r, y = pcf, linetype = data_type_field)) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ size) +
  scale_color_viridis_d() +
  theme_classic()

