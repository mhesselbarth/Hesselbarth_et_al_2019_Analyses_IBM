###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Model structure pattern ####

#### Import libraries and data ####

# load packages
library(extrafont)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(spatstat)
library(tidyverse)

plot_area <- readr::read_rds("Data/Raw/plot_area_owin.rds")

#### Create pattern ####
n <- 250

type <- sample(x = c("living", "dead"), size = n, 
               prob = c(0.7, 0.3), replace = TRUE)

ppp_random <- spatstat::runifpoint(n = n, win = plot_area) %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(type = type)

plot_ppp_total <- ggplot(data = ppp_random) + 
  geom_point(aes(x = x, y = y, shape = type), size = 1.5) +
  geom_polygon(data = spatstat::as.data.frame.owin(plot_area),
               aes(x = x, y = y), col = "black", fill = NA) + 
  scale_shape_manual(values = c(1, 2)) +
  guides(pch = FALSE) +
  coord_equal() +
  theme_void()

#### Save plots ####
suppoRt::save_ggplot(plot = plot_ppp_total, 
                     filename = "plot_ppp_total.png", 
                     path = "C:/Users/Maximilian/ownCloud/13_Thesis_defense/Figures/",
                     dpi = 300, units = "mm",
                     width = 50, height = 50, # width = 200, height = 125, 
                     overwrite = T)
