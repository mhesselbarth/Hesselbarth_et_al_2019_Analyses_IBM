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

# distance to k neighbours
nn_saplings <- purrr::map_df(1:10, function(x) {
  
  nn_dist <- nndist(saplings_1999_ppp, k = x)
  
  tibble::tibble(nn = x, mean = mean(nn_dist), sd = sd(nn_dist), type = "Sapling")
})

nn_adults <- purrr::map_df(1:10, function(x) {
  
  nn_dist <- nndist(adults_1999_ppp, k = x)
  
  tibble::tibble(nn = x, mean = mean(nn_dist), sd = sd(nn_dist), type = "Adult")
})

# combine to one df #
nn_overall <- dplyr::bind_rows(nn_saplings, nn_adults) %>% 
  dplyr::mutate(nn = factor(nn, ordered = TRUE), 
                type = factor(type, levels = c("Sapling", "Adult")))

# plot result #
ggplot_nn <- ggplot(data = nn_overall) + 
  geom_bar(aes(x = nn, y = mean, fill = type), 
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(x = nn, ymin = mean - sd, ymax = mean + sd, group = type), 
                width = 0.25, position = position_dodge(0.9)) +
  scale_fill_manual(name = "", values = c("Sapling" = "#0D0887FF", 
                                           "Adult" = "#CC4678FF")) + 
  labs(x = expression(paste(italic(k),  " nearest neighbour")), 
       y = "Distance [m]") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

ggplot_nnd_nn <- ggplot_nnd + ggplot_nn + 
  patchwork::plot_layout(ncol = 2, nrow = 1,
                         widths = c(0.5, 0.5), heights = c(0.5, 0.5)) + 
  patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")")

suppoRt::save_ggplot(plot = ggplot_nnd_nn,
                     filename = "ggplot_nnd_nn.png",
                     path = "Figures/Appendix",
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units,
                     overwrite = overwrite)

#### Summarise g(r) ####
set.seed(42)

pattern_complex <- expand.grid(x = seq(from = -10, to = 110, by = 20),
                               y = seq(from = -10, to = 110, by = 20)) %>%
  dplyr::mutate(x = x + runif(n = nrow(.), min = 0, max = 5),
                y = y + runif(n = nrow(.), min = 0, max = 5)) %>%
  spatstat::ppp(x = .$x, y = .$y,
                window = owin(c(0, 100), c(0, 100)))

# conver to df
pattern_complex_df <- tibble::as_tibble(pattern_complex)


pattern_complex_df <- purrr::map_dfr(1:nrow(pattern_complex_df), function(i) {
  
  purrr::map_dfr(1:10, function(j) {
    
    tibble::tibble(x = pattern_complex_df[[i, 1]] + runif(n = 1, min = -5, max = 5), 
                   y = pattern_complex_df[[i, 2]] + runif(n = 1, min = -5, max = 5))
    
  })
  
})

pattern_complex <- spatstat::ppp(x = pattern_complex_df$x, 
                                 y = pattern_complex_df$y, 
                                 window = owin(c(0, 100), c(0, 100)))

complex_env <- spatstat::envelope(pattern_complex, fun = "pcf", 
                                  nsim = 199, nrank = 5,
                                  funargs = list(divisor = "d", 
                                                 correction = "Ripley", 
                                                 stoyan = 0.25))

complex_env_sum <- onpoint::summarise_envelope(complex_env)

complex_env_sum_gg <- plot(complex_env_sum, x_lab = "r [m]", y = expression(italic(g(r))), 
                           base_size = base_size)

suppoRt::save_ggplot(plot = complex_env_sum_gg,
                     filename = "ggplot_summarised_env.png",
                     path = "Figures/Appendix",
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units,
                     overwrite = overwrite)
