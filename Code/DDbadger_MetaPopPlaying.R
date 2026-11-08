### Code for building up the metapopulation model.
# Written by Chrissy Hernandez
# July 13 2026
# Last updated: Nov 2 2026

# unequal initial vec
n0_3patch <- round(c(rep(initials[1,]/2, 3)))      # equal patch sizes initially

# multiple patch UMAT ----
threepatch_Umat <- matrix(data=0, nrow=12, ncol=12)
threepatch_Umat[1:4,1:4]<- umat
threepatch_Umat[5:8,5:8]<- umat
threepatch_Umat[9:12, 9:12] <- umat

# no dispersal function testing
output<- proj_meta_nodisp_norem(threepatch_Umat,    
                                n0_3patch,
                                params,
                                stages,
                                npatches = 3,
                                time = 30,
                                DDapply="fertility",
                                return.vec= TRUE)   


# Testing fixed dispersal projection with unequal initial outputs
output <- proj_meta_fixeddisp_norem(threepatch_Umat, n0_3patch, params, stages, 
                                    npatches=3, DDapply = 'fertility', time = 50,
                                    dispersal_prob=0.1, 
                                    dispersal_stages = c("Adult_f", "Adult_m"),
                                    return.vec = TRUE)

# population over time
(meta_fixeddisp_norem_pop <- meta_plot(output,   # output obj of dd.proj
                                      y_val= "N",   # plot type - N or Vec 
                                      ylab = "Population size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16))

(meta_fixeddisp_norem_group <- meta_plot(output,   # output obj of dd.proj
                                      y_val= "Group",   # plot type - N or Vec 
                                      ylab = "Group size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16))
meta_fixeddisp_norem_oldparams <- meta_fixeddisp_norem_pop + meta_fixeddisp_norem_group



# setting new beta value for tighter DD params (explained in other script)
plateau_rick <- exp(-params$b *160) # val when N = 160 (plateau) =  0.2927426  (ricker value at plateau)   
#  want plateau when N = 28 
new_b <- log(plateau_rick)/-30    ## solving for new_b exp(-new_b * 28) = 0.2927426

# updated param set up - tighter DD with smaller beta
new_params <- params 
new_params$b = new_b

# new initial vecs to match smaller group sizes
size <- c(10,20, 16)  #  group sizes
patch <- matrix(0, nrow = length(size), ncol = 4)

for (i in 1:length(size)){   # loop to fill rows of matrix with vector
  patch[i,] <- floor(stagedist*size[i])
}

# tuning this into a vector
n0_3patch <- c(t(patch))
n0_3patch_eq <- rep(round(initials[1,]/5), 3)


# Testing fixed dispersal projection with new params
output <- proj_meta_fixeddisp_norem(threepatch_Umat, n0_3patch, new_params, stages, 
                                    npatches=3, DDapply = 'fertility', time = 100,
                                    dispersal_prob = 0.2, 
                                    dispersal_stages = c("Adult_f", "Adult_m"),
                                    return.vec = TRUE)

# visualise
(meta_fixeddisp_norem_pop <- meta_plot(output,   # output obj of dd.proj
                                      y_val= "N",   # plot type - N or Vec 
                                      ylab = "Population size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16))

(meta_fixeddisp_norem_group <- meta_plot(output,   # output obj of dd.proj
                                        y_val= "Group",   # plot type - N or Vec 
                                        ylab = "Group size", 
                                        xlab = "time (t)",
                                        rem_year = NULL,
                                        mytheme = theme_classic(), 
                                        cols= "black",    # can be vector of 2 cols
                                        legend.pos = "top",
                                        base_size = 16))
meta_fixeddisp_norem_newparams <- meta_fixeddisp_norem_pop + meta_fixeddisp_norem_group

## DENSITY DEPENDENCE ---------
# old params and density dependent movement 
output <- proj_meta_DDdisp_norem(threepatch_Umat, n0_3patch, params, stages, 
                                 npatches=3, DDapply = 'fertility', time = 50,
                                 dispersal_stages = c("Adult_f", "Adult_m"),
                                 return.vec=TRUE)

# population over time
(meta_DDdisp_norem_pop <- meta_plot(output,   # output obj of dd.proj
                                    y_val= "N",   # plot type - N or Vec 
                                    ylab = "Population size", 
                                    xlab = "time (t)",
                                    rem_year = NULL,
                                    mytheme = theme_classic(), 
                                    cols= "black",    # can be vector of 2 cols
                                    legend.pos = "top",
                                    base_size = 16))

(meta_DDdisp_norem_group <- meta_plot(output,   # output obj of dd.proj
                                      y_val= "Group",   # plot type - N or Vec 
                                      ylab = "Group size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16))
meta_DDdisp_norem_oldparams <- meta_DDdisp_norem_pop + meta_DDdisp_norem_group


#new params and density dependent movement 
output <- proj_meta_DDdisp_norem(threepatch_Umat, n0_3patch, new_params, stages, 
                       npatches=3, DDapply = 'fertility', time = 50,
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec=TRUE)

# population over time
(meta_DDdisp_norem_pop <- meta_plot(output,   # output obj of dd.proj
                                      y_val= "N",   # plot type - N or Vec 
                                      ylab = "Population size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16))

(meta_DDdisp_norem_group <- meta_plot(output,   # output obj of dd.proj
                                        y_val= "Group",   # plot type - N or Vec 
                                        ylab = "Group size", 
                                        xlab = "time (t)",
                                        rem_year = NULL,
                                        mytheme = theme_classic(), 
                                        cols= "black",    # can be vector of 2 cols
                                        legend.pos = "top",
                                        base_size = 16))
meta_DDdisp_norem_newparams <- meta_DDdisp_norem_pop + meta_DDdisp_norem_group

# next step is to have plateaus at different sizes (group that moves between 10-20 size, 2-8, 16-26 ect..)
# Can be done by changing beta per patch to represent spatial variation



# NEXT STEPS removals ----
# first proj = old params, no dispersal
test_proj_rem1 <- proj_meta_nodisp_rem(threepatch_Umat, n0_3patch, 
                                      stagedist,  # proportion of each stage
                                      params,     # using old parameters at first
                                      stagenames = stages, 
                                      npatches = 3,
                                      time = 50,     
                                      DDapply="fertility", 
                                      remyear = 10,  # removal year = vector of years 
                                      rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                                      bias = NULL , # strength of bias as percentage (range??)
                                      return.vec= TRUE, 
                                      return.remvec = TRUE)
# string of removal years (subsequent only)
test_proj_rem2 <- proj_meta_nodisp_rem(threepatch_Umat, n0_3patch, 
                                      stagedist,  # proportion of each stage
                                      params,     # using old parameters at first
                                      stagenames = stages, 
                                      npatches = 3,
                                      time = 50,     
                                      DDapply="fertility", 
                                      remyear = c(10:14),  # removal year = vector of years 
                                      rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                                      bias = NULL , # strength of bias as percentage (range??)
                                      return.vec= TRUE, 
                                      return.remvec = TRUE)
# supplementary removal years 
test_proj_rem3 <- proj_meta_nodisp_rem(threepatch_Umat, n0_3patch, 
                                      stagedist,  # proportion of each stage
                                      params,     # using old parameters at first
                                      stagenames = stages, 
                                      npatches = 3,
                                      time = 50,     
                                      DDapply="fertility", 
                                      remyear = c(10:19),  # removal year = vector of years 
                                      rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                                      bias = NULL , # strength of bias as percentage (range??)
                                      return.vec= TRUE, 
                                      return.remvec = TRUE)

group1 <- meta_plot(test_proj_rem1,   # output obj of dd.proj
                                      y_val= "Group",   # plot type - N or Vec 
                                      ylab = "Group size", 
                                      xlab = "time (t)",
                                      rem_year = 10,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16)
group2 <- meta_plot(test_proj_rem2,   # output obj of dd.proj
                    y_val= "Group",   # plot type - N or Vec 
                    ylab = "Group size", 
                    xlab = "time (t)",
                    rem_year = 10,
                    mytheme = theme_classic(), 
                    cols= "black",    # can be vector of 2 cols
                    legend.pos = "top",
                    base_size = 16)
group3 <- meta_plot(test_proj_rem3,   # output obj of dd.proj
                    y_val= "Group",   # plot type - N or Vec 
                    ylab = "Group size", 
                    xlab = "time (t)",
                    rem_year = 10,
                    mytheme = theme_classic(), 
                    cols= "black",    # can be vector of 2 cols
                    legend.pos = "top",
                    base_size = 16)
# plot altogether for quick vis 
removal_summary_plot <- (group1+group2/group3)



# testing no dispersal, removals and new params
meta_nodisp_rem_newparams <- proj_meta_nodisp_rem(threepatch_Umat, n0_3patch, 
                                       stagedist,  # proportion of each stage
                                       new_params,     # using old parameters at first
                                       stagenames = stages, 
                                       npatches = 3,
                                       time = 50,     
                                       DDapply="fertility", 
                                       remyear = 10,  # removal year = vector of years 
                                       rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                                       bias = NULL , # strength of bias as percentage (range??)
                                       return.vec= TRUE, 
                                       return.remvec = TRUE)
meta_nodisp_rem_newplot <- meta_plot(meta_nodisp_rem_newparams,   # output obj of dd.proj
                                       y_val= "Group",   # plot type - N or Vec 
                                       ylab = "Group size", 
                                       xlab = "time (t)",
                                       rem_year = 10,
                                       mytheme = theme_classic(), 
                                       cols= "black",    # can be vector of 2 cols
                                       legend.pos = "top",
                                       base_size = 16)



# DDdisp and old params
meta_DDdisp_rem_oldp <- proj_meta_DDdisp_rem(threepatch_Umat, n0_3patch, 
                                        stagedist,  # proportion of each stage
                                        params, 
                                        stagenames = stages, 
                                        npatches = 3,
                                        time = 50,     
                                        DDapply="fertility", 
                                        dispersal_stages = c("Adult_f", "Adult_m"),
                                        remyear = 10,  # removal year = vector of years 
                                        rem_strat = "random",  # if specified 
                                        bias = NULL, 
                                        return.vec=TRUE, 
                                        return.remvec = TRUE)
meta_DDdisp_rem_oldplot <- meta_plot(meta_DDdisp_rem_oldp,   # output obj of dd.proj
                                     y_val= "Group",   # plot type - N or Vec 
                                     ylab = "Group size", 
                                     xlab = "time (t)",
                                     rem_year = 10,
                                     mytheme = theme_classic(), 
                                     cols= "black",    # can be vector of 2 cols
                                     legend.pos = "top",
                                     base_size = 16)



# testing DDdisp and removals with new params
# ISSUES HERE:
# PLATEAUS AFTER REMOVALS - FAILING TO REBOUND ABOVE A FIXED GROUP SIZE
# MAY BE AN ISSUE IN THE LOOP AFTER REMOVALS?
test_DDdisp_rem_newp <- proj_meta_DDdisp_rem(threepatch_Umat, n0_3patch, 
                                        stagedist,  # proportion of each stage
                                        new_params, 
                                        stagenames = stages, 
                                        npatches = 3,
                                        time = 50,     
                                        DDapply="fertility", 
                                        dispersal_stages = c("Adult_f", "Adult_m"),
                                        remyear = 10,  # removal year = vector of years 
                                        rem_strat = "random",  # if specified 
                                        bias = NULL, 
                                        return.vec=TRUE, 
                                        return.remvec = TRUE)

(meta_DDdisp_rem_pop <- meta_plot(test_DDdisp_rem_newp,   # output obj of dd.proj
                                  y_val= "N",   # plot type - N or Vec 
                                  ylab = "Population size", 
                                  xlab = "time (t)",
                                  rem_year = 10,
                                  mytheme = theme_classic(), 
                                  cols= "black",    # can be vector of 2 cols
                                  legend.pos = "top",
                                  base_size = 16))

(meta_DDdisp_rem_group <- meta_plot(test_DDdisp_rem_newp,   # output obj of dd.proj
                                    y_val= "Group",   # plot type - N or Vec 
                                    ylab = "Group size", 
                                    xlab = "time (t)",
                                    rem_year = 10,
                                    mytheme = theme_classic(), 
                                    cols= "black",    # can be vector of 2 cols
                                    legend.pos = "top",
                                    base_size = 16))

meta_DDdisp_rem_newplot <- meta_DDdisp_rem_pop + meta_DDdisp_rem_group



## NEXT STEPS
# Sett biased removals
