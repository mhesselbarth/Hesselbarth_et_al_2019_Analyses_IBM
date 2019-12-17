###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Appendix ####

# load packages #
source("Helper_functions/helper_functions_setup.R")

# load ppp # 
beech_1999_ppp <- readr::read_rds("Data/Input/beech_1999_ppp.rds")

# filter data #
saplings_1999_ppp <- spatstat::subset.ppp(beech_1999_ppp, dbh_99 > 1 & dbh_99 <= 10)
adults_1999_ppp <- spatstat::subset.ppp(beech_1999_ppp, dbh_99 > 10)

# calculate NND #
nnd_saplings <- spatstat::Gest(X = saplings_1999_ppp, 
                               correction = "km", 
                               r = seq(from  = 0, to = 12.5, length.out = 525)) 

nnd_adults <- spatstat::Gest(X = adults_1999_ppp, 
                             correction = "km", 
                             r = seq(from  = 0, to = 12.5, length.out = 525)) 

# convert to df #
nnd_saplings_df <- tibble::as_tibble(nnd_saplings) %>% 
  dplyr::mutate(type = "Sapling")

nnd_adults_df <- tibble::as_tibble(nnd_adults) %>% 
  dplyr::mutate(type = "Adult")

# combine to one df #
nnd_overall_df <- dplyr::bind_rows(nnd_saplings_df, nnd_adults_df) %>% 
  dplyr::mutate(type = factor(type, 
                              levels = c("Sapling", "Adult")))

# plot result #
ggplot_nnd <- ggplot(data = nnd_overall_df) + 
  geom_line(aes(x = r, y = km, col = type)) + 
  geom_hline(yintercept = 1, linetype = 3, size = 0.25) +
  geom_hline(yintercept = 0.5, linetype = 3, size = 0.25) +
  scale_color_manual(name = "", values = c("Sapling" = "#0D0887FF", 
                                           "Adult" = "#CC4678FF")) + 
  scale_x_continuous(breaks = seq(from = 0, to = 12.5, by = 2.5)) +
  labs(x = "r [m]", y = expression(italic(G(r)))) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_nnd,
                     filename = "nnd_overall_df.png",
                     path = "Figures/Appendix",
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units,
                     overwrite = overwrite)
