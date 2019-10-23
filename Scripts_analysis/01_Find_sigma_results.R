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

#### Load data ####
model_runs_sigma <- readr::read_rds("Data/Output/model_runs/model_runs_sigma.rds")

beech_1999_ppp <- readr::read_rds("Data/Input/beech_1999_ppp.rds")

beech_2007_ppp <- readr::read_rds("Data/Input/beech_2007_ppp.rds")
beech_2007_sapling_ppp <- readr::read_rds("Data/Input/beech_2007_sapling_ppp.rds")
beech_2007_adult_ppp <- readr::read_rds("Data/Input/beech_2007_adult_ppp.rds")

beech_2013_ppp <- readr::read_rds("Data/Input/beech_2013_ppp.rds")
beech_2013_sapling_ppp <- readr::read_rds("Data/Input/beech_2013_sapling_ppp.rds")
beech_2013_adult_ppp <- readr::read_rds("Data/Input/beech_2013_adult_ppp.rds")

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

#### Preprocess data ####
beech_2007_df <- tibble::as_tibble(beech_2007_ppp)
beech_2013_df <- tibble::as_tibble(beech_2013_ppp)

#### number of living trees ####
model_runs_n <- purrr::map(seq_along(model_runs_sigma), function(x) {
  
  # get living trees of last time step
  temp_data <- dplyr::filter(model_runs_sigma[[x]], 
                             i == max(i), type != "dead") %>% 
    dplyr::group_by(type) %>% 
    dplyr::summarise(n = dplyr::n())
})

model_runs_n <- dplyr::bind_rows(model_runs_n, .id = "id") %>% 
  dplyr::mutate(id = as.integer(id))

n_2007 <- dplyr::filter(beech_2007_df, 
                        type != "dead", !is.na(dbh_07), inside_fence == 0) %>% 
  dplyr::mutate(type = dplyr::case_when(dbh_07 > 1 & dbh_07 <= 10 ~ "sapling", 
                                        dbh_07 > 10 ~  "adult")) %>% 
  dplyr::group_by(type) %>% 
  dplyr::summarise(n = dplyr::n()) %>% 
  dplyr::mutate(data_type = "Field data 2007")

n_2013 <- dplyr::filter(beech_2013_df, 
                        type != "dead", !is.na(dbh_13), inside_fence == 0) %>% 
  dplyr::mutate(type = dplyr::case_when(dbh_13 > 1 & dbh_13 <= 10 ~ "sapling", 
                                        dbh_13 > 10 ~  "adult")) %>% 
  dplyr::group_by(type) %>% 
  dplyr::summarise(n = dplyr::n()) %>% 
  dplyr::mutate(data_type = "Field data 2013")

n_field <- dplyr::bind_rows(n_2007,
                            n_2013)

ggplot(data = model_runs_n) +
  geom_bar(data = n_field ,
           aes(x = type, y = n, fill = data_type),
           position = position_dodge(), stat = "identity") +
  geom_point(aes(x = type, y = n), pch = "-", 
             size = 5, col = "red") +
  scale_fill_viridis_d(name = "", option = "D") +
  facet_wrap( ~ id) + 
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

model_runs_n_filtered <- dplyr::filter(model_runs_n, !id %in% c(1, 2, 3, 8, 10, 
                                                                12, 14, 16, 18, 
                                                                20, 22, 24, 26, 
                                                                28, 30))

ggplot(data = model_runs_n_filtered) +
  geom_bar(data = n_field ,
           aes(x = type, y = n, fill = data_type),
           position = position_dodge(), stat = "identity") +
  geom_point(aes(x = type, y = n), pch = "-", 
             size = 5, col = "red") +
  scale_fill_viridis_d(name = "", option = "D") +
  facet_wrap( ~ id) + 
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

id_n <- unique(model_runs_dbh_filtered$id)

##### DBH dist ####
by <- 10

model_runs_dbh <- purrr::map(seq_along(model_runs_sigma), function(x) {
  
  message("\r> Progress: ", x, "/", length(model_runs_sigma), appendLF = FALSE)
  
  
  # get living trees of last time step
  temp_data <- dplyr::filter(model_runs_sigma[[x]], 
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

model_runs_dbh_filtered <- dplyr::filter(model_runs_dbh, !id %in% c(1, 2, 3, 
                                                                    8, 10, 12, 
                                                                    14, 16, 18, 
                                                                    20, 22, 24,
                                                                    26, 28, 30))

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

model_runs_nnd <- purrr::map(seq_along(model_runs_sigma), function(x) {
  
  message("\r> Progress: ", x, "/", length(model_runs_sigma), appendLF = FALSE)
  
  # get data of last timestep
  temp_data <- dplyr::filter(model_runs_sigma[[x]], i == max(i), type != "dead")
  
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

ggplot(data = dplyr::filter(model_runs_nnd, size == "saplings")) + 
  geom_line(aes(x = r, y = nnd)) + 
  geom_line(data = dplyr::filter(nnd_overall_field, size == "saplings"),
            aes(x = r, y = nnd, linetype = data_type_field)) +
  facet_wrap(~ id) +
  scale_color_viridis_d() +
  theme_classic()

model_runs_nnd_filtered <- dplyr::filter(model_runs_nnd, !id %in% c(1, 2, 8, 
                                                                    10, 12, 14, 
                                                                    16, 18, 20,
                                                                    22, 24, 26,
                                                                    28, 30))

ggplot(data = dplyr::filter(model_runs_nnd_filtered, size == "saplings")) + 
  geom_line(aes(x = r, y = nnd)) + 
  geom_line(data = dplyr::filter(nnd_overall_field, size == "saplings"),
            aes(x = r, y = nnd, linetype = data_type_field)) +
  facet_wrap(~ id) +
  scale_color_viridis_d() +
  theme_classic()

id_nnd <- unique(model_runs_nnd_filtered$id)

#### Pair correlation function #####
r_pcf <- seq(from = 0, to = 50, length.out = 525)
correction_pcf <- "Ripley"
stoyan_pcf <- 0.25
divisor_pcf <- "d"

model_runs_pcf <- purrr::map(seq_along(model_runs_sigma), function(x) {
  
  message("\r> Progress: ", x, "/", length(model_runs_sigma), appendLF = FALSE)
  
  # get data of last timestep
  temp_data <- dplyr::filter(model_runs_sigma[[x]], i == max(i), type != "dead")
  
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

ggplot(dplyr::filter(model_runs_pcf, size == "saplings")) + 
  geom_line(aes(x = r, y = pcf), col = "#21908CFF") + 
  geom_line(data = dplyr::filter(pcf_overall_field, size == "saplings"),
            aes(x = r, y = pcf, linetype = data_type_field)) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ id) +
  scale_color_viridis_d() +
  theme_classic()

model_runs_pcf_filtered <- dplyr::filter(model_runs_pcf, !id %in% c(1, 2, 3, 
                                                                    4, 5, 7, 9, 
                                                                    11, 13, 15, 
                                                                    17, 19, 21, 
                                                                    23, 25, 27, 
                                                                    29))

ggplot(dplyr::filter(model_runs_pcf_filtered, size == "saplings")) + 
  geom_line(aes(x = r, y = pcf), col = "#21908CFF") + 
  geom_line(data = dplyr::filter(pcf_overall_field, size == "saplings"),
            aes(x = r, y = pcf, linetype = data_type_field)) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ id) +
  scale_color_viridis_d() +
  theme_classic()

id_pcf <- unique(model_runs_pcf_filtered$id)


#### Find common x
x <- dplyr::intersect(x = dplyr::intersect(x = dplyr::intersect(x = id_pcf,
                                                                y = id_nnd), 
                                           y = id_dbh), y = id_n)

dplyr::select(model_runs_sigma[[x]], sigma, probs)
# sigma = 15, 
# probs = 0.25/0.75
