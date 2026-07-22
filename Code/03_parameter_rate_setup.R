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
set.seed(123)  # setting our number 
reps <- 100
pop.sizes <- runif(reps, min = 30, max=240) # min pop size = 25, max = 240. 

initials <- matrix(0, nrow = reps, ncol = 4)

for (t in 1:reps){   # loop to fill rows of matrix with vector
  initials[t,] <- floor(stagedist*pop.sizes[t])
}


