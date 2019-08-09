
# load packages
library(spatstat)
library(tidyverse)

# import data
pattern_1999 <- readr::read_rds("C:/Users/Lena/Desktop/projektarbeit/Daten/pattern_1999.rds")

pattern_1999_reconstructed <- readr::read_rds("C:/Users/Lena/Desktop/projektarbeit/Daten/pattern_1999_reconstructed_20000.rds")

# remove not needed marks
pattern_1999 <- spatstat::subset.ppp(pattern_1999, select = c("Species", "DBH_07"))

#### 1. compare spatial structure ####

# 1.1 Pair-correlation function #

# calculate pair-correlation function observed pattern
pcf_1999 <- spatstat::pcf(pattern_1999,
                          correction = "best", divisor = "d") %>%
  tibble::as.tibble() %>%
  dplyr::mutate(pattern = "original pattern",
                summary_function = "Pair correlation function")

# calculate pair-correlation function reconstructed pattern
pcf_1999_reconstructed <- spatstat::pcf(pattern_1999_reconstructed,
                                        correction = "best", divisor = "d") %>%
  tibble::as.tibble() %>%
  dplyr::mutate(pattern = "reconstructed pattern",
                summary_function = "Pair correlation function")

# combine to one data frame
pcf_complete <- dplyr::bind_rows(pcf_1999, pcf_1999_reconstructed)

names(pcf_complete)[[3]] <- "value"

# 1.2 Nearest neighbourh distance function #

# calculate nn-function observed pattern
gest_1999 <- spatstat::Gest(pattern_1999,
                            correction = "han") %>%
  tibble::as.tibble() %>%
  dplyr::mutate(pattern = "original pattern",
                summary_function = "Nearest neighbour distribution function")

# calculate pair-correlation function reconstructed pattern
gest_1999_reconstructed <- spatstat::Gest(pattern_1999_reconstructed,
                                          correction = "han") %>%
  tibble::as.tibble() %>%
  dplyr::mutate(pattern = "reconstructed pattern",
                summary_function = "Nearest neighbour distribution function")

# combine to one data frame
gest_complete <- dplyr::bind_rows(gest_1999, gest_1999_reconstructed)

names(gest_complete)[[3]] <- "value"

# 1.3 Plot results #

# create one complete data frame
spatial_complete <- dplyr::bind_rows(pcf_complete, gest_complete)

# plot results
ggplot(spatial_complete) +
  geom_line(aes(x = r, y = theo, col = "CSR"), linetype = 2) +
  geom_line(aes(x = r, y = value, col = pattern)) +
  facet_wrap(~ summary_function, scales = "free") +
  scale_color_manual(values = c("CSR" = "grey",
                                "original pattern" = "black",
                                "reconstructed pattern" = "red"), name = "") +
  labs(x = "r [m]", y = "G(r)") +
  theme_bw()

#### 2. Marks structure ####

# 2.1 DBH structure #

# calculate mark correlation function of DBH original pattern
mark_correlation_1999 <- spatstat::subset.ppp(pattern_1999, select = "DBH_99") %>%
  spatstat::markcorr(correction = "Ripley") %>%
  tibble::as.tibble() %>%
  dplyr::mutate(pattern = "original pattern",
                summary_function = "Mark correlation function")

# calculate mark correlation function of DBH reconstructed pattern
mark_correlation_1999_reconstructed <- spatstat::subset.ppp(pattern_1999_reconstructed, select = "DBH") %>%
  spatstat::markcorr(correction = "Ripley") %>%
  tibble::as.tibble() %>%
  dplyr::mutate(pattern = "reconstructed pattern",
                summary_function = "Mark correlation function")

# combinte to one data frame
mark_correlation_complete <- dplyr::bind_rows(mark_correlation_1999,
                                              mark_correlation_1999_reconstructed)

names(mark_correlation_complete)[[3]] <- "value"

# 2.2 Species structure #

# calculate mark connection function of DBH original pattern
mark_connection_1999 <- spatstat::subset.ppp(pattern_1999, select = "Species") %>%
  spatstat::markconnect(correction = "Ripley") %>%
  tibble::as.tibble() %>%
  dplyr::mutate(pattern = "original pattern",
                summary_function = "Mark connection function")

# calculate mark connection function of DBH reconstructed pattern
mark_connection_1999_reconstructed <- spatstat::subset.ppp(pattern_1999_reconstructed, select = "Species") %>%
  spatstat::markconnect(correction = "Ripley") %>%
  tibble::as.tibble() %>%
  dplyr::mutate(pattern = "reconstructed pattern",
                summary_function = "Mark connection function")

# combinte to one data frame
mark_connection_complete <- dplyr::bind_rows(mark_connection_1999,
                                             mark_connection_1999_reconstructed)

names(mark_connection_complete)[[3]] <- "value"

# 2.3 Plot results #

# create one complete data frame
marks_complete <- dplyr::bind_rows(mark_correlation_complete,
                                   mark_connection_complete)

# plot results
ggplot(marks_complete) +
  geom_line(aes(x = r, y = theo, col = "CSR"), linetype = 2) +
  geom_line(aes(x = r, y = value, col = pattern)) +
  facet_wrap(~ summary_function, scales = "free") +
  scale_color_manual(values = c("CSR" = "grey",
                                "original pattern" = "black",
                                "reconstructed pattern" = "red"), name = "") +
  labs(x = "r [m]", y = "G(r)") +
  theme_bw()
