# Script for Chrissy
# 15.03.26

# loading data (not needed)
# Custom functions - go through all_functions folder
# Look at projection outputs

# loading required packages ---------------
library(ggplot2)
library(gridExtra)
library(tidyr)
library(tidyverse)
library(dplyr)
library(readr)
library(here)


# sourcing required functions from all functions script
source(here("Code/Functions/01_all_functions.R"))

# data extraction
source(here("Code/02_data_extraction.R"))  # select options 2 (existing data) and 1 (all data)

# setting up parameters and vital rates for demographic model ------
stages <- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")

Umat <- matrix(0, nrow=4, ncol=4)
rownames(Umat) <- stages
Umat[2,1]<- bright_survival_vec[1]  # yearling f survival
Umat[2,2]<- bright_survival_vec[2]  # adult f survival - could use macdonald 2002 paper values
Umat[4,3]<-bright_survival_vec[3]   # yearling m survival
Umat[4,4]<- bright_survival_vec[4]  # adult m survival

# assigning values from data extraction script to parameter dataframe
params<- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
                    b= beta,       # calculated from mcdonald 2016
                    rep_K= rogers_k,          # max litter size (K), 
                    h= 6)  # harem size per male




# projections (model 1 - non-spatial) -----
n0 <- c(5, 10, 5, 10)  # initial population structure = yf, af, ym, am 

# baseline = 20 year projection, no removals
proj0 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 20, 
                  DDapply="Fmat", 
                  intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                  remyear = NULL, 
                  rem_strat = NULL ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL ,
                  return.vec= TRUE, 
                  return.remvec = FALSE) 

# plots for basic visualisations = the trajectory of each scenario investigated
# setting colours for male and female plots
col_vec <- c("#FF6A6A", "#87CEEB")

# stage abundance over time
(proj0_plot <- dd_plot(proj0, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = NULL,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16))

# population size over time
(proj0_Nplot <- dd_plot(proj0, 
                        y_val= "N", 
                        ylab = "Abundance", 
                        xlab = "Time (t)",
                        rem_year = NULL,
                        mytheme = theme_classic(), 
                        cols= col_vec,    # can be vector of cols
                        legend.pos = "top",
                        base_size = 16))


# Scenario 1 = 70% removal trial at year 10 -----
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


# plots
# Stage Abundance
proj1_plot <- dd_plot(proj1, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 5,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16)
# Population size
proj1_Nplot <- dd_plot(proj1, 
                        y_val= "N", 
                        ylab = "Abundance", 
                        xlab = "Time (t)",
                        rem_year = 5,
                        mytheme = theme_classic(), 
                        cols= col_vec,    # can be vector of cols
                        legend.pos = "top",
                        base_size = 16)
grid.arrange(proj1_plot, proj1_Nplot)



# Scenario 2 - biased male removals ----
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

# plots 
# Stage abundance
proj2_plot <- dd_plot(proj2, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16)
# Population size 
proj2_Nplot <- dd_plot(proj2, 
                       y_val= "N", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16)
grid.arrange(proj2_plot, proj2_Nplot)



# Scenario 3 - biased female adult removals --------
proj3 <- rem.proj(Umat,      
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  DDapply="Fmat", 
                  intensity= 80,  # percentage you want REMOVED from pop at time=ry
                  remyear = 10, 
                  rem_strat = 2 ,  #2nd in list = adult fems
                  bias = 0.2,  # bias too strong?
                  return.vec= TRUE, 
                  return.remvec = TRUE)

# plots
# Stage abundance
proj3_plot <- dd_plot(proj3, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16)

# Population size
proj3_Nplot <- dd_plot(proj3, 
                       y_val= "N", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16)
grid.arrange(proj3_plot, proj3_Nplot)


# Output metrics ------
# Lambda 
lambda0 <- growth.rate(proj0, vis = TRUE)
summary(lambda0$lambda)

lambda1 <- growth.rate(proj1, vis = TRUE, rem_year = 5)
summary(lambda1$lambda)

lambda2 <- growth.rate(proj2, vis = TRUE, rem_year = 5)
summary(lambda2$lambda)

lambda3 <- growth.rate(proj3, vis = TRUE, rem_year = 5)
summary(lambda3$lambda)

# Stable stage distributions
ssd0 <- ssd(proj0, vis = TRUE, cols = col_vec)
summary(ssd0$stageMat)


ssd1 <- ssd(proj1, vis = TRUE, cols = col_vec)
summary(ssd1$stageMat)

ssd2 <- ssd(proj2, vis = TRUE, cols = col_vec)
summary(ssd2$stageMat)

ssd3 <- ssd(proj3, vis = TRUE, cols = col_vec)
summary(ssd3$stageMat)

grid.arrange(ssd0$plot, ssd1$plot, ssd2$plot, ssd3$plot) # dont need keys repeated on plot!




# Repeated projections, comparing lambda, pop size and ssd in scenarios! -----
# Basline projection analysis
rep_proj0 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         method = "random",
                         intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                         remyear = NULL, 
                         rem_strat =NULL ,  # if specified removals, "adults, females, yearling females... 
                         bias = NULL ,
                         return.vec= TRUE, 
                         return.remvec = FALSE, 
                         reps = 100) 
# why extinction?

# Yearling_f  Adult_f Yearling_m  Adult_m
# 0          20  0.00000         10 67.00000   
# 1         NaN 17.01827        NaN 58.28236     # should have reproduction between 17 fems and 58 males? 
# 2           0  0.00000          0  0.00000
        

# first removal scenario = 70% random
rep_proj1 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         method = "random",
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear = 5, 
                         rem_strat = "random" ,  # if specified removals, "adults, females, yearling females... 
                         bias = 0.1 ,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) # looks better !


# Second scenario = 70% male biased 
rep_proj2 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat", 
                         method = "random",
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear = 5, 
                         rem_strat = "male" ,  # if specified removals, "adults, females, yearling females... 
                         bias = 0.1,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) 


# Third scenario = 90% female biased 
rep_proj3 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         params, 
                         stagenames = stages,
                         time = 20, 
                         DDapply="Fmat",
                         method = "random",
                         intensity= 80,  # percentage you want REMOVED from pop at time T=ry
                         remyear = 5, 
                         rem_strat = "female" ,  # if specified removals, "adults, females, yearling females... 
                         bias = 0.15 ,
                         return.vec= TRUE, 
                         return.remvec = TRUE, 
                         reps = 100) 

# averages for each scenario
av_proj0 <- pop.av(rep_proj0,  return.Lambda = TRUE, #  lambda per year
                   return.Mats = TRUE)
# calculating av ssd 
stage_vec <- av_proj0$av_prop
stageDist <- colMeans(stage_vec)  # using this in our repeat proj function

# calculating av pop size across baseline

mean.pop <- function(proj){
 mean(proj$pop)
}

rep_meanPops <- sapply(rep_proj0,mean.pop)  # mean across years for each rep

meanPop <- mean(rep_meanPops)  # use mean and sd in normal dist when generating initial vecs
sdPop <- sd(rep_meanPops)

av_proj0 <- pop.av(rep_proj0,  
                   return.Lambda = TRUE, #  lambda per year
                   return.Mats = TRUE)
av_proj1 <- pop.av(rep_proj1, 
                   rep_proj0, 
                   return.Lambda = TRUE, #  lambda per year
                   return.Mats = TRUE )
av_proj2 <- pop.av(rep_proj2, 
                   rep_proj0, 
                   return.Lambda = TRUE, #  lambda per year
                   return.Mats = TRUE)
av_proj3 <- pop.av(rep_proj3, 
                   rep_proj0, 
                   return.Lambda = TRUE, #  lambda per year
                   return.Mats = TRUE)

proj0_lambda <- av_proj0$av_lambda
proj1_lambda <- av_proj1$av_lambda
proj2_lambda <- av_proj2$av_lambda
proj3_lambda <- av_proj3$av_lambda

# setting up dataframe for plots
av_df <- data.frame(proj0_lambda, proj1_lambda, proj2_lambda, proj3_lambda) # spread so col for projection name, lamb value (later )
av_df %>% 
  pivot_longer(av_df, cols = c(proj0_lambda, proj1_lambda), 
               names_to = "Projection",
              values_to = "av_lambda")



# Meta population projection ---------
# generating initial vec 
rand <- runif(8, 0, 20) # likely not more than 20 individuals of a given stage in a group 
init <- list(rand[1:4], rand[5:8])  

Dmat <- matrix(0, ncol= length(stages), nrow = 2) # row = number patches
colnames(Dmat) <- stages

# prob individ in patch1 remains in patch 1 - rnorm
stay1 <- rnorm(4,mean = 0.5, sd =0.2)  # vec of prob for each stage (4)
Dmat1 <- Dmat
Dmat1[1,] <- stay1 # stay for each class in patch 1 
Dmat1[2,] <- 1- stay1   # for multiple patches, will need multiple values that sum to 1, can't use 1-p

stay2 <- rnorm(4,mean = 0.5, sd =0.2)
Dmat2 <- Dmat
Dmat2[1,] <- stay2 # stay for each class in patch 2 
Dmat2[2,] <- 1 - stay2 

D_list <- list(Dmat1, Dmat2)

meta_proj0 <- meta.proj(Umat,   # vector of stage names
                        params, # adjusted beta = 
                        stagenames = stages,
                        initial = init,   # list of initial abundances
                        Dmat = D_list, 
                        time = 20,
                        return.vec = TRUE,
                        return.mats= FALSE)


