###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

# load packages #
library(rabmp)
library(raster)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(spatstat)
library(tidyverse)

# set parapemters plotting #
overwrite <- FALSE

base_size <- 12.5

width_full <- 21.0
width_small <- 17.5

height_full <- 29.7
height_small <- 12.5

units <- "cm"

dpi <- 300
