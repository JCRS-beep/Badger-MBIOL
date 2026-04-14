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


# parameter definition and df set ups
source(here("Code/03_parameter_rate_setup.R"))   # enter 2 then 1

# everything above this in a seperate script to source?

# projections (model 1) -----
n0 <- c(12, 41, 12, 34) # vec structure = yf, af, ym, am. if pop size = 100, ssd gives this vec

# baseline = 20 year projection
proj0 <- multi.rem(Umat,     
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
stagedist <- round(ssd$stageMat[20,], 3) # final dist used as ssd, rounded to 4 sf

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
rep_proj0 <- repeat.proj(func = "goal", 
                         Umat,    
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
rep_proj1 <- repeat.proj(func = "goal", 
                         Umat,    
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
rep_proj2 <- repeat.proj(func = "goal", 
                         Umat,   
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
rep_proj3 <- repeat.proj(func = "goal", 
                         Umat,    
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
# Scenario 1 = 70% removal trial at year 5 AND 6
du_proj1 <- repeat.proj(func = "goal", 
                        Umat,    
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



# scenario 2 - biased male removals at year 5 AND 6
du_proj2 <- repeat.proj(func = "goal", 
                        Umat,     
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


# Scenario 3 - biased female adult removals at year 5 AND 6
du_proj3 <- repeat.proj(func = "goal", 
                        Umat,     
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
# first removal scenario = 70% random 
multi_proj1 <- repeat.proj(func = "goal", 
                           Umat,     
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
multi_proj2 <- repeat.proj(func = "goal", 
                           Umat,     
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
multi_proj3 <- repeat.proj(func = "goal",
                           Umat,     
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
                           

# continuous removals all years of projection at low levels
# removing 10% every year = needs old rem.proj function OR 15 years of 10 % = intensity of 150%?
constant.proj <- repeat.proj(func = "prop",
                             Umat,      
                             initial.vecs = initials,
                             stagedist = NULL,
                             params = params, 
                             stagenames = stages,
                             time = 20, 
                             DDapply="Fmat", 
                             intensity= 10,  # percentage you want REMOVED from pop at time T=ry
                             remyear =  c(5:19), 
                             rem_strat = "random" ,  # 2 = Af 
                             bias = NULL,
                             return.vec= TRUE, 
                             return.remvec = TRUE, 
                             reps = 100) 


