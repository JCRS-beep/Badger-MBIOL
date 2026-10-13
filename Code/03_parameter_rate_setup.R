# Load all required rates and params 
# Script to source before any models are run


# setting up parameters and vital rates for demographic model --------
stages <- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")

umat <- matrix(0, nrow=4, ncol=4)
rownames(umat) <- stages
colnames(umat) <- stages
umat[2,1]<- bright_survival_vec[1]  # yearling f survival
umat[2,2]<- bright_survival_vec[2]  # adult f survival - could use macdonald 2002 paper values
umat[4,3]<-bright_survival_vec[3]   # yearling m survival
umat[4,4]<- bright_survival_vec[4]  # adult m survival

# extract straight from data?
params<- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
                    b= beta,       # calculated from mcdonald 2016
                    rep_K= rogers_k,          # max litter size (K), 
                    h= 10)   # harem size per male - assume single male sufficient for groups


# creating initial vectors here?
# inital vec for single projection 
stagedist <- c(0.12, 0.43, 0.12, 0.34)
n0 <- stagedist * 100 # vec structure = yf, af, ym, am. if pop size = 100, ssd gives this vec


# initial vecs for repeated projection scenarios
# need to generate initial vecs in repeatable way
set.seed(123)  # setting our seed to get repeatable code 
reps <- 100
pop_sizes <- runif(reps, min = 30, max=240) # min pop size = 25, max = 240. 

initials <- matrix(0, nrow = reps, ncol = 4)

for (t in 1:reps){   # loop to fill rows of matrix with vector
  initials[t,] <- floor(stagedist*pop_sizes[t])
}


# PART II - meta-pop model set up ----
# setting new beta value for tighter DD params (explained in other script)
plateau_rick <- exp(-params$b *160) # val when N = 160 (plateau) =  0.2927426  (ricker value at plateau)   
#  want plateau when N = 28 
new_b <- log(plateau_rick)/-28    ## solving for new_b exp(-new_b * 28) = 0.2927426

# updated param set up - tighter DD with smaller beta
meta_params <- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
                    b= new_b,       # should vary by group if we set different max group sizes
                    rep_K= rogers_k,          # max litter size (K), 
                    h= 10)   # harem size per male - assume single male sufficient for groups



# setting up repititions

nGroups <- 10        # how many groups?
# how to split pop_sizes across 10 groups unequally?
# x <- seq(0.05, 0.3, 0.025) # proportions (group size 1 = 0.05* total pop size, until 30% of total pop). will not sum to total pop size. 
                          # to equal total pop_size, need to sample 10 proportions that sum to 1

# set.seed(123)  # setting repitition number 


# limits <- sample(4:28, groups, replace = TRUE)    # max group size = 28, wht to choose for min - single individual of each stage



