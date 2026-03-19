# Repetitions and analysis
# 15.03.26

# model projection script to get params and Umat, libraries and model scenarios

#  comparison - run models 100 repetitions, plot mean lambda, mean ssd for scenarios


# Basline projection analysis
test <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
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


testN <- dd_plot(test[[2]], 
                  y_val= "N",   # plot type - N or Vec 
                  ylab = "abundance", 
                  xlab = "time (t)",
                  rem_year = NULL,
                  mytheme = theme_classic(), 
                  cols= col_vec,    # can be vector of 2 cols
                  legend.pos = "top",
                  base_size = 16)

testV <- dd_plot(test[[2]], 
                 y_val= "Vec",   # plot type - N or Vec 
                 ylab = "abundance", 
                 xlab = "time (t)",
                 rem_year = NULL,
                 mytheme = theme_classic(), 
                 cols= col_vec,    # can be vector of 2 cols
                 legend.pos = "top",
                 base_size = 16)


grid.arrange(testN, testV)   # diff projections out of 10


# first removal scenario = 70% random
proj1 <- repeat.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
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

proj1N <- dd_plot(proj1[[1]], 
                  y_val= "N",   # plot type - N or Vec 
                  ylab = "abundance", 
                  xlab = "time (t)",
                  rem_year = 5,
                  mytheme = theme_classic(), 
                  cols= col_vec,    # can be vector of 2 cols
                  legend.pos = "top",
                  base_size = 16)

proj1_plot <- dd_plot(proj1[[1]], 
                 y_val= "Vec",   # plot type - N or Vec 
                 ylab = "abundance", 
                 xlab = "time (t)",
                 rem_year = 5,
                 mytheme = theme_classic(), 
                 cols= col_vec,    # can be vector of 2 cols
                 legend.pos = "top",
                 base_size = 16)


metrics1 <- pop.av(proj1)

# trial comparison = lambda box plots, 
