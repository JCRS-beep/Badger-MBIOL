# Modelling our different scenarios, producing outputs for plotting
# 23/02/25


# loading required packages
library(here)


# sourcing required functions

# everything above this in a seperate script to source?

# projections (model 1) -----
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

ssd0 <- ssd(proj0, vis = FALSE)
stagedist <- round(ssd0$stageMat[20,], 3) # final dist used as ssd, rounded to 3 sf

# Repeated projections, do not return long outputs !


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
cont.proj <- repeat.proj(func = "prop",
                         Umat,      
                         initial.vecs = initials,
                         stagedist,
                         params, 
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


# meta proj scenarios
# baseline - 20 years in patches