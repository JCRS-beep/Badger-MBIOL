# Modelling our different scenarios, producing outputs for plotting
# 23/02/25
# Can be run entirely in order

# loading required packages
library(ggplot2)
library(gridExtra)
library(tidyr)
library(tidyverse)
library(dplyr)
library(readr)
library(here)


# sourcing required functions
source(here("Code/Functions/01_all_functions.R"))

# data extraction
source(here("Code/02_data_extraction.R"))  # select options 2 (existing data) and 1 (all data)

# setting up parameters and vital rates for demographic model
stages <- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")

Umat <- matrix(0, nrow=4, ncol=4)
rownames(Umat) <- stages
Umat[2,1]<- bright_survival_vec[1]  # yearling f survival
Umat[2,2]<- bright_survival_vec[2]  # adult f survival - could use macdonald 2002 paper values
Umat[4,3]<-bright_survival_vec[3]   # yearling m survival
Umat[4,4]<- bright_survival_vec[4]  # adult m survival

# extract straight from data?
params<- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
                    b= beta,       # calculated from mcdonald 2016
                    rep_K= rogers_k,          # max litter size (K), 
                    h= 6   # harem size per male
)

# projections (model 1) -----
n0 <- c(25, 10, 25, 10)  # vec structure = yf, af, ym, am 

# baseline = 20 year projection
proj0 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 20, 
                  DDapply="Fmat", 
                  intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                  remyear = NULL, 
                  rem_strat =NULL ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL ,
                  return.vec= TRUE, 
                  return.remvec = FALSE) 



# Scenario 1 = 70% removal trial at year 10
proj1 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 20, 
                  DDapply="Fmat", 
                  intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                  remyear = 5, 
                  rem_strat = "random" ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL ,
                  return.vec= TRUE, 
                  return.remvec = TRUE) 



# scenario 2 - biased male removals
proj2 <- rem.proj(Umat,      
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  DDapply="Fmat", 
                  intensity= 70,  # percentage you want REMOVED from pop at time=ry
                  remyear = 10, 
                  rem_strat = "males" ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = 0.15 ,
                  return.vec= TRUE, 
                  return.remvec = TRUE) 


# Scenario 3 - biased female adult removals
proj3 <- rem.proj(Umat,      
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  DDapply="Fmat", 
                  intensity= 70,  # percentage you want REMOVED from pop at time=ry
                  remyear = 10, 
                  rem_strat = 2 ,  #2nd in list = adult fems
                  bias = 0.2,  # bias too strong?
                  return.vec= TRUE, 
                  return.remvec = TRUE) 
# EXTINCTON! but why do individuals that develop into adults not breed at following timestep? issue with setting to 0?


# Repeated projections, do not return long outputs !
# Basline projection analysis
rep_proj0 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                    params, 
                    stagenames = stages,
                    time = 20, 
                    DDapply="Fmat", 
                    intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                    remyear = NULL, 
                    rem_strat =NULL ,  # if specified removals, "adults, females, yearling females... 
                    bias = NULL ,
                    return.vec= TRUE, 
                    return.remvec = FALSE, 
                    reps = 10) # looks better !



# first removal scenario = 70% random
rep_proj1 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                     params, 
                     stagenames = stages,
                     time = 20, 
                     DDapply="Fmat", 
                     intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                     remyear = 5, 
                     rem_strat = "random" ,  # if specified removals, "adults, females, yearling females... 
                     bias = 70 ,
                     return.vec= TRUE, 
                     return.remvec = TRUE, 
                     reps = 100) # looks better !


# Second scenario = 70% male biased 
rep_proj2 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                     params, 
                     stagenames = stages,
                     time = 20, 
                     DDapply="Fmat", 
                     intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                     remyear = 5, 
                     rem_strat = "male" ,  # if specified removals, "adults, females, yearling females... 
                     bias = 70 ,
                     return.vec= TRUE, 
                     return.remvec = TRUE, 
                     reps = 100) 


# Third scenario = 90% female biased 
rep_proj3 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear = 5, 
                         rem_strat = "female" ,  # if specified removals, "adults, females, yearling females... 
                         bias = 70 ,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) 





                           




