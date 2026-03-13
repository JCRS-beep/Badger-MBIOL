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
source(here("Code/Functions/all_functions.R"))

# data extraction
source(here("Code/data.extraction.R"))  # select options 2 (existing data) and 1 (all data)

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

col_vec <- c("#FF6A6A", "#87CEEB")

(proj0_plot <- dd_plot(proj0, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16))
(proj0_Nplot <- dd_plot(proj0, 
                        y_val= "N", 
                        ylab = "Abundance", 
                        xlab = "Time (t)",
                        mytheme = theme_classic(), 
                        cols= col_vec,    # can be vector of cols
                        legend.pos = "top" ,
                        base_size = 16))


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

(proj1_plot <- dd_plot(proj1, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 5,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16))



# scenario 2 - biased female removals
proj2 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  DDapply="Fmat", 
                  intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                  remyear = 10, 
                  rem_strat = "females" ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = 0.15 ,
                  return.vec= TRUE, 
                  return.remvec = TRUE) 

(proj2_plot <- dd_plot(proj2, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       base_size = 16))


# comparisons 
lambda0 <- growth.rate(proj0, vis = TRUE)
summary(lambda0$lambda)

lambda1 <- growth.rate(proj1, vis = TRUE, rem_year = 10)
summary(lambda1$lambda)

lambda2 <- growth.rate(proj2, vis = TRUE, rem_year = 10)
summary(lambda2$lambda)

ssd0 <- ssd(proj0, vis = TRUE, cols = col_vec)
summary(ssd0$stageMat)
ssd1 <- ssd(proj1, vis = TRUE, cols = col_vec)
ssd2 <- ssd(proj2, vis = TRUE, cols = col_vec)


                           




