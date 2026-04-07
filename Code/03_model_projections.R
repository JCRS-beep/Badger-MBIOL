# Modelling our different scenarios, producing outputs for plotting
# 23/02/25
# Can be run entirely in order

# UPDATES - replace rem.proj with multi.rem 


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
                    h= 6)   # harem size per male


# projections (model 1) -----
n0 <- c(12, 41, 12, 34) # vec structure = yf, af, ym, am. if pop size = 100, ssd gives this vec

# baseline = 20 year projection
proj0 <- multi.rem(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  stagedist = c(0.1, 0.4, 0.1, 0.4),   # approx stage dist 
                  params = params, 
                  stagenames = stages,
                  time = 20, 
                  DDapply="Fmat", 
                  intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                  remyear = NULL, 
                  rem_strat = NULL ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL ,
                  return.vec= TRUE, 
                  return.remvec = FALSE) 

ssd <- ssd(proj0, vis = FALSE)
stagedist <- ssd$stageMat[20,] # final dist used as ssd

# Repeated projections, do not return long outputs !

# need to generate initial vecs in repeatable way
set.seed(123)  # setting our number 
reps <- 100
pop.sizes <- runif(reps, min=25, max=240) 

initials <- matrix(0, nrow = reps, ncol = 4)

for (t in 1:reps){   # loop to fill rows of matrix with vector
  initials[t,] <- floor(stagedist*pop.sizes[t])
}


# Basline projection analysis
rep_proj0 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         initial.vecs = initials,
                         stagedist = stagedist,
                         params = params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                         remyear = NULL, 
                         rem_strat =NULL ,  # if specified removals, "adults, females, yearling females... 
                         bias = NULL ,
                         return.vec= TRUE, 
                         return.remvec = FALSE, 
                         reps = 100) 


# first removal scenario = 70% random
rep_proj1 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         initial.vecs = initials,
                         stagedist = stagedist,
                         params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear = 5, 
                         rem_strat = "random" ,  # if specified removals, "adults, females, yearling females... 
                         bias = NULL ,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) # looks better !
                    


# Second scenario = 70% ADULT male biased 
rep_proj2 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         initial.vecs = initials,
                         stagedist = stagedist,
                         params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear = 5, 
                         rem_strat = 4,  # adult males
                         bias = 0.1,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) 


# Third scenario = 70% ADULT female biased 
rep_proj3 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         initial.vecs = initials,
                         stagedist = stagedist,
                         params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear = 5, 
                         rem_strat = 2 ,  # 2 = Af 
                         bias = 0.1,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) 




# Frequency - if projections are repeated every year for 2 years -----
# Scenario 1 = 70% removal trial at year 10
du_proj1 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                        initial.vecs = initials,
                        stagedist = stagedist,
                        params, 
                        stagenames = stages,
                        time = 20, 
                        DDapply="Fmat", 
                        intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                        remyear = c(5,6),  
                        rem_strat = "random" ,  # if specified removals, "adults, females, yearling females... 
                        bias = NULL ,
                        return.vec= TRUE, 
                        return.remvec = TRUE, 
                        reps = 100) # looks better !



# scenario 2 - biased male removals
du_proj2 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                        initial.vecs = initials,
                        stagedist = stagedist,
                        params = params, 
                        stagenames = stages,
                        time = 20, 
                        DDapply="Fmat", 
                        intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                        remyear = c(5,6),  
                        rem_strat = 4 ,  # if specified removals, "adults, females, yearling females... 
                        bias = 0.1 ,
                        return.vec= TRUE, 
                        return.remvec = TRUE, 
                        reps = 100) # looks better !


# Scenario 3 - biased female adult removals
du_proj3 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                        initial.vecs = initials,
                        stagedist = stagedist,
                        params = params, 
                        stagenames = stages,
                        time = 20, 
                        DDapply="Fmat", 
                        intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                        remyear = c(5,6), 
                        rem_strat = 2 ,  # 2 = Af 
                        bias = 0.1,
                        return.vec= TRUE, 
                        return.remvec = TRUE, 
                        reps = 100) 





# if removals are repeated for 5 years --------
# first removal scenario = 50% random
multi_proj1 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                           initial.vecs = initials,
                           stagedist = stagedist,
                           params = params,  
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear = c(5:9), 
                         rem_strat = "random" ,  # if specified removals, "adults, females, yearling females... 
                         bias = NULL ,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) # looks better !


# Second scenario = 70% ADULT male biased 
multi_proj2 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                           initial.vecs = initials,
                           stagedist = stagedist,
                           params = params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear =  c(5:9), 
                         rem_strat = 4 ,  # if specified removals, "adults, females, yearling females... 
                         bias = 0.1,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) 


# Third scenario = 70% ADULT female biased 
multi_proj3 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                           initial.vecs = initials,
                           stagedist = stagedist,
                         params = params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear =  c(5:9), 
                         rem_strat = 2 ,  # 2 = Af 
                         bias = 0.1,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) 
                           




