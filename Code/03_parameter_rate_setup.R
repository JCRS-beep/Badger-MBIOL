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
temp_stagedist <- c(0.12, 0.43, 0.12, 0.34) # approx stagedist
n0 <- temp_stagedist * 100 # vec structure = yf, af, ym, am. if pop size = 100, ssd gives this vec


# initial vecs for repeated projection scenarios
# need to generate initial vecs in repeatable way
set.seed(123)  # setting our seed to get repeatable code 
reps <- 100
pop_sizes <- runif(reps, min = 30, max=240) # min pop size = 25, max = 240. 

initials <- matrix(0, nrow = reps, ncol = 4)

for (t in 1:reps){   # loop to fill rows of matrix with vector
  initials[t,] <- floor(temp_stagedist*pop_sizes[t])
}


# PART II - meta-pop model set up ----
threepatch_Umat <- matrix(data=0, nrow=12, ncol=12)
threepatch_Umat[1:4,1:4]<- umat
threepatch_Umat[5:8,5:8]<- umat
threepatch_Umat[9:12, 9:12] <- umat



# setting new beta value for tighter DD params (explained in other script)
# calculating max group size from projection
proj0 <- multi.rem(umat,     
                   initial = n0, 
                   stagedist = temp_stagedist,   # approx stage dist 
                   params = params, 
                   stagenames = stages,
                   time = 100, 
                   DDapply="Fmat", 
                   intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                   remyear = NULL, 
                   rem_strat = NULL ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                   bias = NULL ,
                   return.vec= TRUE, 
                   return.remvec = FALSE) 

ssd0 <- ssd(proj0, vis = FALSE)
stagedist <- round(ssd0$stageMat[100,], 3)   # final dist used as ssd, rounded to 3 sf
maxN <- max(proj0$pop)

plateau_rick <- exp(-params$b * maxN) # val when N = 160 (plateau) =  0.2927426  (ricker value at plateau)   
#  want plateau when N = 28 
new_b <- log(plateau_rick)/-28    ## solving for new_b exp(-new_b * 28) = 0.2927426

# updated param set up - tighter DD with smaller beta
# updated param set up - tighter DD with smaller beta
meta_params <- params 
meta_params$b = new_b


# setting up new initial vecs to match smaller group sizes
size <- c(25,20, 16)  #  group sizes
patch <- matrix(0, nrow = length(size), ncol = 4)

for (i in 1:length(size)){   # loop to fill rows of matrix with vector
  patch[i,] <- floor(stagedist*size[i])
}

# tuning this into a vector
n0_3patch <- c(t(patch))