###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Preprocess data ####

#### Import libraries ####
source("Helper_functions/helper_functions_setup.R")

#### Import data #### 
pattern_1999_ppp <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")
pattern_2007_ppp <- readr::read_rds("Data/Raw/pattern_2007_ppp.rds")
pattern_2013_ppp <- readr::read_rds("Data/Raw/pattern_2013_ppp.rds")

#### Filter data ####
beech_1999_ppp <- spatstat::subset.ppp(pattern_1999_ppp, 
                                       species == "beech" & type == "living" &
                                       dbh_99 > 1)

beech_2007_ppp <- spatstat::subset.ppp(pattern_2007_ppp, 
                                       species == "beech" & type == "living" &
                                         dbh_07 > 1 & inside_fence == 0)

beech_2013_ppp <- spatstat::subset.ppp(pattern_2013_ppp, 
                                       species == "beech" & type == "living" &
                                         dbh_13 > 1 & inside_fence == 0)

#### Create size classes #### 
beech_2007_sapling_ppp <- spatstat::subset.ppp(beech_2007_ppp, 
                                               dbh_07 > 1 & dbh_07 <= 10 & 
                                             inside_fence == 0)

beech_2007_adult_ppp <- spatstat::subset.ppp(beech_2007_ppp, dbh_07 > 10 & 
                                           inside_fence == 0)

beech_2013_sapling_ppp <- spatstat::subset.ppp(beech_2013_ppp, 
                                               dbh_13 > 1 & dbh_13 <= 10 & 
                                             inside_fence == 0)

beech_2013_adult_ppp <- spatstat::subset.ppp(beech_2013_ppp, dbh_13 > 10 & 
                                           inside_fence == 0)

#### Save data ####
suppoRt::save_rds(object = beech_1999_ppp, 
                  filename = "beech_1999_ppp.rds", 
                  path = "Data/Input/", 
                  overwrite = overwrite)

suppoRt::save_rds(object = beech_2007_ppp, 
                  filename = "beech_2007_ppp.rds", 
                  path = "Data/Input/", 
                  overwrite = overwrite)

suppoRt::save_rds(object = beech_2013_ppp, 
                  filename = "beech_2013_ppp.rds", 
                  path = "Data/Input/", 
                  overwrite = overwrite)

suppoRt::save_rds(object = beech_2007_sapling_ppp, 
                  filename = "beech_2007_sapling_ppp.rds", 
                  path = "Data/Input/", 
                  overwrite = overwrite)

suppoRt::save_rds(object = beech_2007_adult_ppp, 
                  filename = "beech_2007_adult_ppp.rds", 
                  path = "Data/Input/", 
                  overwrite = overwrite)

suppoRt::save_rds(object = beech_2013_sapling_ppp, 
                  filename = "beech_2013_sapling_ppp.rds", 
                  path = "Data/Input/", 
                  overwrite = overwrite)

suppoRt::save_rds(object = beech_2013_adult_ppp, 
                  filename = "beech_2013_adult_ppp.rds", 
                  path = "Data/Input/", 
                  overwrite = overwrite)
