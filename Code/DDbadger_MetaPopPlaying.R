### Code for building up the metapopulation model.
# Written by Chrissy Hernandez
# July 13 2026
# Last updated:

# We're gonna build up from the absolute simplest model. Two patches, no dispersal:
twopatch_Umat<- matrix(data=0, nrow=8, ncol=8)
twopatch_Umat[1:4,1:4]<- umat
twopatch_Umat[5:8,5:8]<- umat

# fertility parameters are set in script 03_parameter_rate_setup.R (copied here for reference)
# params<- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
#                     b= beta,       # calculated from mcdonald 2016
#                     rep_K= rogers_k,          # max litter size (K), 
#                     h= 10)   # harem size per male - assume single male sufficient for groups

# I want the two patches to have different starting populations, so I'm just
# pasting together two of the initial pop vectors from the non-spatial model
n0_2patch<- c(initials[1,], initials[2,]) 

output<- proj_meta_nodisp_norem(twopatch_Umat, n0_2patch, params, stages, 
                         npatches=2, DDapply = "fertility", time=20, 
                         return.vec = TRUE)

# Okay, now two patches with dispersal: The problem with adding dispersal is
# that we want them to disperse before they mate, and then mating is
# density-dependent. So this is why we keep running into issues with the
# megamatrix approach, and switching to loops. I think we also wanted to work
# with discrete individuals, which is also a bit clunky in the matrix setting.
# So I've written a sort of hybrid approach, where we do a patch-based model
# where each patch follows an MPM.

output <- proj_meta_fixeddisp_norem(twopatch_Umat, n0_2patch, params, stages, 
                                   npatches=2, DDapply = 'fertility', time = 20,
                                   dispersal_prob=0.1, 
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec = TRUE)


# multiple patch UMAT 
n0_3patch <- c(initials[1,], initials[2,], initials[3,]) 
threepatch_Umat <-  matrix(data=0, nrow=12, ncol=12)
threepatch_Umat[1:4,1:4]<- umat
threepatch_Umat[5:8,5:8]<- umat
threepatch_Umat[9:12, 9:12] <- umat

# Testing fixed dispersal projection 
output <- proj_meta_fixeddisp_norem(threepatch_Umat, n0_3patch, params, stages, 
                                    npatches=3, DDapply = 'fertility', time = 20,
                                    dispersal_prob=0.1, 
                                    dispersal_stages = c("Adult_f", "Adult_m"),
                                    return.vec = TRUE)
# population over time
meta_fixeddisp_norem_pop <- meta_plot(output,   # output obj of dd.proj
                                      y_val= "N",   # plot type - N or Vec 
                                      ylab = "Population size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16)

meta_fixeddisp_norem_group <- meta_plot(output,   # output obj of dd.proj
                                      y_val= "Group",   # plot type - N or Vec 
                                      ylab = "Group size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16)



# Increasing group number
set.seed(123)  # setting repitition number 
groups <- 20        # how many groups?
limits <- sample(10:28, groups, replace = TRUE)    # max group size = 28, what to choose for min (10 temporary number)? 


# multiple patches and density dependent movement -------
output <- proj_meta_DDdisp_norem(threepatch_Umat, n0_3patch, params, stages, 
                       npatches=3, DDapply = 'fertility', time = 20,
                                   max_group = 28,  # based on reading
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec=TRUE)


# interesting - decreasing pop size due to movement limits. Definitely issues in leavers matrix calculation - impossible numbers
# population over time
meta_DDdisp_norem_pop <- meta_plot(output,   # output obj of dd.proj
                                      y_val= "N",   # plot type - N or Vec 
                                      ylab = "Population size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16)

meta_fixedDD_norem_group <- meta_plot(output,   # output obj of dd.proj
                                        y_val= "Group",   # plot type - N or Vec 
                                        ylab = "Group size", 
                                        xlab = "time (t)",
                                        rem_year = NULL,
                                        mytheme = theme_classic(), 
                                        cols= "black",    # can be vector of 2 cols
                                        legend.pos = "top",
                                        base_size = 16)



## Setting tighter density dependence for groups? How to apply density dependence across entire population, limiting group sizes?
m_params <- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
                       b= 0.012,       # should vary by group if we set different max group sizes
                       rep_K= rogers_k,          # max litter size (K), 
                       h= 10)   # harem size per male - assume single male sufficient for groups

# new initial vecs
# group sizes
size <- c(10,20, 12)
n0_3patch <- matrix(0, nrow = length(size), ncol = 4)

for (i in 1:length(size)){   # loop to fill rows of matrix with vector
  n0_3patch[i,] <- floor(stagedist*size[i])
}


# trying new params with fixed disp
output <- proj_meta_fixeddisp_norem(threepatch_Umat,
                                    n0_3patch,  # if we change max group sizes, must also change starting vecs as there are very large
                                    m_params, 
                                    stages, 
                                    npatches=3, DDapply = 'fertility', time = 20,
                                    dispersal_prob=0.1, 
                                    dispersal_stages = c("Adult_f", "Adult_m"),
                                    return.vec = TRUE)
meta_fixeddisp_norem_pop2 <- meta_plot(output,   # output obj of dd.proj
                                      y_val= "N",   # plot type - N or Vec 
                                      ylab = "Population size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16)

meta_fixeddisp_norem_group2 <- meta_plot(output,   # output obj of dd.proj
                                        y_val= "Group",   # plot type - N or Vec 
                                        ylab = "Group size", 
                                        xlab = "time (t)",
                                        rem_year = NULL,
                                        mytheme = theme_classic(), 
                                        cols= "black",    # can be vector of 2 cols
                                        legend.pos = "top",
                                        base_size = 16)
