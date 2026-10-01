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
# NOTE - to test model outputs I'm using equal starting pops for now
n0_2patch <- round(c((initials[1,]/4), (initials[1,]/4)))  # starting with equal group sizes to test if function performs correctly

output<- proj_meta_nodisp_norem(twopatch_Umat,     # must be a big matrix with diaongal Umats
                                n0_2patch,
                                params,
                                stages,
                                npatches = 2,
                                time = 20,
                                DDapply="fertility",
                                return.vec= TRUE)   


# Okay, now two patches with dispersal: The problem with adding dispersal is
# that we want them to disperse before they mate, and then mating is
# density-dependent. So this is why we keep running into issues with the
# megamatrix approach, and switching to loops. I think we also wanted to work
# with discrete individuals, which is also a bit clunky in the matrix setting.
# So I've written a sort of hybrid approach, where we do a patch-based model
# where each patch follows an MPM.

output <- proj_meta_fixeddisp_norem(twopatch_Umat, n0_2patch, params, stages, 
                                   npatches=2, time = 20,
                                   DDapply = 'fertility', 
                                   dispersal_prob=0.1, 
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec = TRUE)
# debugging completed - projections now equal with equal starting vectors. 

# multiple patch UMAT ----
# unequal initial vec
n0_eq_3patch <- round(c(rep(initials[1,]/4,3)))
n0_3patch <- round(c((initials[1,]/4), (initials[2,]/4), (initials[3,]/4))) # equal patch sizes initially

threepatch_Umat <-  matrix(data=0, nrow=12, ncol=12)
threepatch_Umat[1:4,1:4]<- umat
threepatch_Umat[5:8,5:8]<- umat
threepatch_Umat[9:12, 9:12] <- umat

# Testing fixed dispersal projection with equal initial outputs
output <- proj_meta_fixeddisp_norem(threepatch_Umat, n0_eq_3patch, params, stages, 
                                    npatches=3, DDapply = 'fertility', time = 50,
                                    dispersal_prob=0.1, 
                                    dispersal_stages = c("Adult_f", "Adult_m"),
                                    return.vec = TRUE)
# good - equal group sizes across years


# Testing fixed dispersal projection with unequal initial outputs
output <- proj_meta_fixeddisp_norem(threepatch_Umat, n0_3patch, params, stages, 
                                    npatches=3, DDapply = 'fertility', time = 50,
                                    dispersal_prob=0.1, 
                                    dispersal_stages = c("Adult_f", "Adult_m"),
                                    return.vec = TRUE)
# good - equal group sizes across years


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



# multiple patches and density dependent movement -------
output <- proj_meta_DDdisp_norem(threepatch_Umat, n0_3patch, params, stages, 
                       npatches=3, DDapply = 'fertility', time = 50,
                                   max_group = 28,  # based on reading - carrying capacity 
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec=TRUE)
# High number of leavers (43 in year 1, group 1) due to large group sizes compared to limits set in DD movement function
# may improve once new params introduced

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
# increases above 1500?

meta_fixedDD_norem_group <- meta_plot(output,   # output obj of dd.proj
                                        y_val= "Group",   # plot type - N or Vec 
                                        ylab = "Group size", 
                                        xlab = "time (t)",
                                        rem_year = NULL,
                                        mytheme = theme_classic(), 
                                        cols= "black",    # can be vector of 2 cols
                                        legend.pos = "top",
                                        base_size = 16)
# why such different patterns for the groups?


# looking how reproduction vaies with Nm and Nf in density dependent version
# creating an updated function that tells up 



## Setting tighter density dependence for groups
m_params <- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
                       b= 0.046,  # how to claculate instead of choosing random value? Want groups to max out around 28 
                                  # plot 
                                  # should vary by group if we set different max group sizes
                       rep_K= rogers_k,          # max litter size (K), 
                       h= 10)   # harem size per male - assume single male sufficient for groups

# new initial vecs to match smaller group sizes
size <- c(10,20, 16)  #  group sizes
patch <- matrix(0, nrow = length(size), ncol = 4)

for (i in 1:length(size)){   # loop to fill rows of matrix with vector
  patch[i,] <- floor(stagedist*size[i])
}

# tuning this into a vector
n0_3patch <- c(t(patch))


# trying new params with fixed disp
output <- proj_meta_DDdisp_norem(threepatch_Umat,
                                 n0_3patch,  
                                 m_params, 
                                 stages, 
                                 npatches = 3,
                                 DDapply = 'fertility', 
                                 time = 20,
                                 max_group = 28,
                                 dispersal_stages = c("Adult_f", "Adult_m"),
                                 return.vec = TRUE)

# visualising output
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




## Increasing group number
# 10 patch Umat
ten_patch_umat <- matrix(data=0, nrow=40, ncol=40)

for(i in 1:10){
  x <- (1 +4*(i-1)) : (4*i)      # equation to get row and col boundaries per patch   (4, 8, 12, ...)
  ten_patch_umat[x,x] <- umat
}

 
set.seed(123)  # setting repitition number 
groups <- 10        # how many groups?
limits <- sample(10:20, groups, replace = TRUE)    # max group size = 28, what to choose for min (10 temporary number)? 

# generating initial vecs
many_patch <- matrix(0, nrow = 10, ncol = 4)
for (i in 1:length(limits)){   # loop to fill rows of matrix with vector
  many_patch[i,] <- floor(stagedist*limits[i])
}

n0_many <- c(t(many_patch))

output <- proj_meta_DDdisp_norem(ten_patch_umat,
                                 n0_many,  
                                 m_params, 
                                 stages, 
                                 npatches = 10,
                                 DDapply = 'fertility', 
                                 time = 20,
                                 max_group = 28,
                                 dispersal_stages = c("Adult_f", "Adult_m"),
                                 return.vec = TRUE)


 # using group output from function
g_df <- as.data.frame(output$group)
g_df$Year <- c(0:20)

g_df <- g_df |>
   pivot_longer(
  cols = starts_with("v"),
  names_to = "Group",
  names_prefix = "Group",
  values_to = "Abundance"
  )



plot <- ggplot(data= g_df, aes(x= Year, y=Abundance, colour = Group)) +
  geom_line() +
  geom_point()


many_group_pop <- meta_plot(output,   # output obj of dd.proj
                             y_val= "Pop",   # plot type - N or Vec 
                             ylab = "Population size", 
                             xlab = "time (t)",
                             rem_year = NULL,
                             mytheme = theme_classic(), 
                             cols= "black",    # can be vector of 2 cols
                             legend.pos = "top",
                             base_size = 16)

