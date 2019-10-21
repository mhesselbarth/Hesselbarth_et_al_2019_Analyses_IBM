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

foo <- function(data,
                sigma,
                parameters,
                probs_id,
                probs,
                plot_area,
                years,
                save_each) {

  # filter data using threshold sapling/adult #
  data_abiotic <- spatstat::subset.ppp(data, dbh_99 > 10 & species == "beech")

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

  data <- tibble::as_tibble(data)
  data <- dplyr::select(data, 
                        x, y, species, dbh_99, type)
  data <- dplyr::filter(data, species == "beech")
  data <- dplyr::mutate(data, type = "adult")
  data <- dplyr::select(data, -species)
  data <- rabmp::prepare_data(data,
                              x = "x", y = "y", type = "type", dbh = "dbh_99")

  result <-  rabmp::run_model_abiotic(data = data,
                                      parameters = parameters,
                                      abiotic = habitat$scaled,
                                      probs = probs[[probs_id]],
                                      plot_area = plot_area,
                                      years = years,
                                      save_each = save_each,
                                      verbose = FALSE)
  
  result$sigma <- sigma
  result$probs <- paste(probs[[probs_id]], collapse = "/")

  return(result)
}

#### Import data ####

# import data  #
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

#### Run systematic sigma exploratation ####
model_runs_mort <- suppoRt::submit_to_cluster(foo,
                                              sigma = sigma,
                                              probs_id = probs_id,
                                              const = list(data = pattern_1999,
                                                           parameters = parameters_fitted_abiotic,
                                                           probs = probs,
                                                           plot_area = plot_area,
                                                           years = years,
                                                           save_each = save_each),
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

