### Code for building up the metapopulation model.
# Written by Chrissy Hernandez
# July 13 2026
# Last updated:

# multiple patch UMAT ----
# unequal initial vec
n0_eq_3patch <- round(c(rep(initials[1,]/4,3)))      # equal patch sizes initially
n0_3patch <- round(c((initials[1,]/4), (initials[2,]/4), (initials[3,]/4))) # different patch sizes initially

threepatch_Umat <-  matrix(data=0, nrow=12, ncol=12)
threepatch_Umat[1:4,1:4]<- umat
threepatch_Umat[5:8,5:8]<- umat
threepatch_Umat[9:12, 9:12] <- umat


output<- proj_meta_nodisp_norem(threepatch_Umat,     # must be a big matrix with diaongal Umats
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


# Testing fixed dispersal projection withnew params
output <- proj_meta_fixeddisp_norem(threepatch_Umat, n0_3patch, new_params, stages, 
                                    npatches=3, DDapply = 'fertility', time = 100,
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


## DENSITTY DEPENDENCE ---------
# multiple patches and density dependent movement 
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

# next step is to have plateaus at different sizes (group that moves between 10-20 size, 2-8, 16-26 ect..)
# Can be done by changing beta per patch to represent spatial variation


## Increasing group number-----
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
                                 new_params, 
                                 stages, 
                                 npatches = 10,
                                 DDapply = 'fertility', 
                                 time = 100,
                                 dispersal_stages = c("Adult_f", "Adult_m"),
                                 return.vec = TRUE)


 # using group output from function
g_df <- as.data.frame(output$group)
g_df$Year <- c(0:100)

g_df <- g_df |>
   pivot_longer(
  cols = starts_with("v"),
  names_to = "Group",
  names_prefix = "Group",
  values_to = "Abundance"
  )



plot <- ggplot(data= g_df, aes(x= Year, y=Abundance, colour = Group)) +
  geom_line() +
  geom_point()     # too many groups to see patterns


many_group_pop <- meta_plot(output,   # output obj of dd.proj
                             y_val= "Pop",   # plot type - N or Vec 
                             ylab = "Population size", 
                             xlab = "time (t)",
                             rem_year = NULL,
                             mytheme = theme_classic(), 
                             cols= "black",    # can be vector of 2 cols
                             legend.pos = "top",
                             base_size = 16)


# NEXT STEPS removals (in removal script) ----
test_proj_rem <- proj_meta_nodisp_rem(threepatch_Umat, n0_3patch, 
                                      stagedist,  # proportion of each stage
                                      params, 
                                      stagenames = stages, 
                                      npatches = 3,
                                      time = 50,     
                                      DDapply="fertility", 
                                      #      intensity= 70,  # percentage you want REMOVED from pop at time T=ry. Remove this - always 70
                                      remyear = c(10:16),  # removal year = vector of years 
                                      rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                                      bias = NULL , # strength of bias as percentage (range??)
                                      return.vec= TRUE, 
                                      return.remvec = TRUE)

(meta_nodisp_rem_pop <- meta_plot(test_proj_rem,   # output obj of dd.proj
                                    y_val= "N",   # plot type - N or Vec 
                                    ylab = "Population size", 
                                    xlab = "time (t)",
                                    rem_year = NULL,
                                    mytheme = theme_classic(), 
                                    cols= "black",    # can be vector of 2 cols
                                    legend.pos = "top",
                                    base_size = 16))

(meta_nodisp_rem_group <- meta_plot(test_proj_rem,   # output obj of dd.proj
                                      y_val= "Group",   # plot type - N or Vec 
                                      ylab = "Group size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16))

# incorporating dispersal 
