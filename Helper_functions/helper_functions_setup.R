###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

# load packages #
library(data.table)
library(magrittr)
library(landscapemetrics)
library(MESS)
library(onpoint)
library(patchwork)
library(quantreg)
library(rabmp)
library(raster)
library(readr)
library(Rcpp)
library(sensitivity)
library(shar)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(stringr)
library(tgp)
library(spatstat)
library(tidyverse)


# set parapemters plotting #
overwrite <- FALSE

base_size <- 12.5

width_full <- 210
width_small <- 175

height_full <- 297
height_small <- 125

units <- "mm"

dpi <- 300
