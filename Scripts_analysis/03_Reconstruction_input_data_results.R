###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Results reconstruct input pattern ####

#### Import libraries and data ####
source("Helper_functions/helper_functions_setup.R")

#### import data ####
pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

# filter for all beech trees
beech_1999 <- spatstat::subset.ppp(pattern_1999, species == "beech")

# inport reconstruted data
beech_1999_rec <- readr::read_rds("Data/Input/beech_1999_rec.rds")

#### Compare non-spatial structure #### 

# make sure number of points is identical
beech_1999$n
beech_1999_rec$n

# make sure only beech
beech_1999$marks$species %>% table()
beech_1999_rec$marks$species %>% table()

# make sure all all living/adult
beech_1999$marks$type %>% table()
beech_1999_rec$marks$type %>% table()

#### DBH distribution ####
dbh_dist <- beech_1999$marks %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(dbh_class = cut(dbh_99, breaks = seq(from = 0, to = 130, by = 5)), 
                data_type = "Input data") %>% 
  dplyr::select(dbh_class, data_type) %>% 
  dplyr::group_by(data_type, dbh_class) %>% 
  dplyr::summarise(n = n())

dbh_dist_rec <- beech_1999_rec$marks %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(dbh_class = cut(dbh, breaks = seq(from = 0, to = 130, by = 5)), 
                data_type = "Reconstructed data") %>% 
  dplyr::select(dbh_class, data_type) %>% 
  dplyr::group_by(data_type, dbh_class) %>% 
  dplyr::summarise(n = n())

dbh_dist_full <- dplyr::bind_rows(dbh_dist, dbh_dist_rec)

ggplot(data = dbh_dist_full) + 
  geom_bar(aes(x = dbh_class, y = n, fill = data_type), 
           stat  = "identity", position = "dodge") + 
  scale_fill_viridis_d(name = "Data type") + 
  theme_classic()

#### Compare spatial structure ####
correction_pcf <- "Ripley"
divisor <- "d"

correction_nnd <- "km"

correction_mcf <- "Ripley"

# calculate pcf # 
pcf_beech <- spatstat::pcf(X = beech_1999, 
                           correction = correction_pcf, 
                           divisor = divisor) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "fun")) %>% 
  dplyr::mutate(data_type = "Input data", 
                sf = "Pair-correlation function")

pcf_beech_rec <- spatstat::pcf(X = beech_1999_rec, 
                               correction = correction_pcf, 
                               divisor = divisor) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "fun")) %>% 
  dplyr::mutate(data_type = "Reconstructed data", 
                sf = "Pair-correlation function")

pcf_full <- dplyr::bind_rows(pcf_beech, pcf_beech_rec)

# calculate NND #
nnd_beech <- spatstat::Gest(X = beech_1999, 
                            correction = correction_nnd) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, correction_nnd) %>% 
  purrr::set_names(c("r", "theo", "fun")) %>% 
  dplyr::mutate(data_type = "Input data", 
                sf = "Neareast neighbor distribution function")

nnd_beech_rec <- spatstat::Gest(X = beech_1999_rec, 
                                correction = correction_nnd) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(r, theo, km) %>% 
  purrr::set_names(c("r", "theo", "fun")) %>% 
  dplyr::mutate(data_type = "Reconstructed data", 
                sf = "Neareast neighbor distribution function")

nnd_full <- dplyr::bind_rows(nnd_beech, nnd_beech_rec)

# calculate mark-corrlation function # 
mc_beech <- spatstat::subset.ppp(x = beech_1999, select = dbh_99) %>% 
  spatstat::markcorr(correction = correction_mcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "fun")) %>% 
  dplyr::mutate(data_type = "Input data", 
                sf = "Mark-correlation function")

mc_beech_rec <- spatstat::subset.ppp(x = beech_1999_rec, select = dbh) %>% 
  spatstat::markcorr(correction = correction_mcf) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("r", "theo", "fun")) %>% 
  dplyr::mutate(data_type = "Reconstructed data", 
                sf = "Mark-correlation function")

mc_full <- dplyr::bind_rows(mc_beech, mc_beech_rec)

# combine all sf #
sf_full <- dplyr::bind_rows(pcf_full, nnd_full, mc_full) %>% 
  dplyr::mutate(sf = factor(sf, levels = c("Neareast neighbor distribution function",
                                           "Pair-correlation function", 
                                           "Mark-correlation function")))

ggplot_sf_intput_recon <- ggplot(data = sf_full) + 
  geom_line(aes(x = r, y = fun, col = data_type)) + 
  geom_line(aes(x = r, y = theo), col = "black", linetype = 2) +
  scale_color_manual(name = "Data type", values = c("#440154FF", "#21908CFF")) + 
  facet_wrap(~ sf, scales = "free", 
             nrow = 3, ncol = 1) + 
  labs(x = "r [m]", y = "f(r)") +
  theme_classic(base_size = 15) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_sf_intput_recon, 
                     filename = "ggplot_sf_intput_recon.png", 
                     path = "Figures/Appendix/", 
                     width = 21.0, height = 29.7, units = "cm", dpi = 300)


