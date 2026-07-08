# Removal methods
# 31.12.2025
# 
library(tidyverse)
library(ggplot2)

# testing no removal scenario----
proj0 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  DDapply="Fmat", 
                  intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                  remyear = NULL, 
                  rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL , # strength of bias as percentage (range??) how to ignore in function if null?
                  return.vec= TRUE, 
                  return.remvec = FALSE) 

col_vec <- c("#FF6A6A", "#87CEEB")

(proj0_plot <- dd_plot(proj0,
                       y_val = "Vec",
                       ylab = "Abundance",
                       xlab = "Time (t)",
                       rem_year = NULL,
                       mytheme = theme_classic(),
                       cols = col_vec,
                       legend.pos = "top",
                       base_size = 16))
(N0_plot <- dd_plot(proj0, 
                   y_val= "N", 
                   ylab = "Pop size", 
                   xlab = "Time (t)",
                   mytheme = theme_classic(), 
                   cols= col_vec,    # can be vector of cols
                   legend.pos = "topright",
                   base_size = 16))  



# NEXT STEPS = Multiple removals rem_year1, 2....
# goal - "remove X% of pop every 2 years for 50 years, long term pop growth rate.
# Syntax = remove at remyear = seq(10,30, by=2) 
 
# updating function design - multi removals 
