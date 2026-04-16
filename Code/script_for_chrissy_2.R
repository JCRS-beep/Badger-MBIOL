#  Meta mat set up and projection
# Second script for review - feedback on functions and outputs
# 14/04

# loading required packages ---------------
# library(ggplot2)
# library(gridExtra)
library(tidyr)
library(tidyverse)
library(dplyr)
library(readr)
library(here)

# sourcing scripts--------
# # sourcing functions used throughout 
source(here("Code/Functions/01_all_functions.R"))

# data extraction
source(here("Code/02_data_extraction.R"))  # select options 2 (existing data) and 1 (all data)

# parameter definition and df set ups
source(here("Code/03_parameter_rate_setup.R"))   # enter 2 then 1

# function to create meta mat
source(here("Code/Meta_MPM_creation.R"))   # contains creation function (can review in file)

# testing this function ------
set.seed(1)  # repeatable
# random movement prob per patch, equal across stages and sexes
cname <- c("Adults")
Dmat <- Dmat.create(colNames = cname, nPatches = 3)  #only movement probs per patch, ncol  number moving indivs (just adults or all stages?)

# initial vector
initial <- c(1,4,1,4,2,6,2,6,2,8,2,7)  # 3 patches

meta_mat <- meta.mat(Umat,   # matrix of stage-specific survival/growth transitions (will be the same in every patch)
                     params, # parameters for density-dependent reproduction
                     post_mat = matrix(initial, ncol = 4, byrow = TRUE), # matrix with post movement abundances
                     Dmat, 
                     stagenames = colnames(Umat), 
                     nDim = 4, 
                     nPatches = 3)



# projection function -----
meta.proj <- function(Umat,   # vector of stage names
                      params, # adjusted beta = 
                      stagenames,
                      initial,   # list of initial abundances - string length patches * nStages
                      colNames, # length = nCol Dmat
                      time = 20,
                      return.vec = TRUE,
                      return.group = FALSE){
  
  # Time 
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  
  nDim <- length(stagenames)     # how many stages
  patches = length(initial)/ nDim # how many times are stages repeated?
  # if(is.integer(patches) == FALSE) stop("Initial vector provided is incorrect length, must be nPatches*length(patch_vec)")
  
  # Set up the output - how to set up if multiple patches? Lists within lists...
  
  out <- list(vec = list(matrix(0, ncol = length(stagenames), nrow = time + 1)),    # list of matrices per patch. # I dont like the formatting, why $vec$vec[[i]]?
              group_size = matrix(0, ncol =  patches, nrow = time + 1),    #   each row of matrix is patch size that year
              total_pop = vector(),   # vector for total pop size per year
              Nremoved = vector(), # total removed in rem year, NO LIST, vec length = num patches
              remvec = list())    # vector of length(stages) per patch - list
  
  # renaming before repeating
  colnames(out$vec[[1]]) <- stagenames   # naming cols matrix as stages - cant name lists of obj
  rownames(out$vec[[1]]) <- 0:time   # rows correspond to each year of projection. Row 0 = initial or initial
  
  rownames(out$group_size) <- 0:time
  
  Vec <- rep(out$vec, patches) # matrices to fill with stage abundance.  row= time. empty list of matrices 
  Group <- out$group_size
  Pop <- out$total_pop       # vector to fill with total pop size each year
  Remvec <- out$remvec # list of vectors per patch, only 1 year 
  
  # filling initial year for outputs
  abun <- matrix(initial, ncol = 4, byrow = TRUE)   # initial matrix with vec per row, nrow = patches
  
  for (p in 1:patches){  # filling first row of each list vec 
    Vec[[p]][1,] <- abun[p,]          # for each list in Vec, first row is initial 
  }
  
  Group[1,] <- rowSums(abun)   # summing each row of matrix which contains abundance vec. Write in terms of vec?
  Pop[1] <- sum(Group[1,]) # first pop entry is sum of group sizes
  
  # Loop = matrix proj each year  --------
  for (i in 1:time) {   # repeat for as many years as we have
    # meta creation function to get large matrix for this year
    # combining entries for list vector format to long vec
    
    # Calculating random movement this year per patch
    thisDmat <- Dmat.create(colNames, patches) # creates random p(move) per patch, either single prob per patch or varies by stage/sex
    
    # calculating post move vec for use in meta_create
    # current vec as string
    this_mat <- t(sapply(Vec, function(x) x[i,]))    # getting initial vec as a matrix with this years abundance per patch (patch = row)
    movers_mat <- this_mat  
    movers_mat[,1] <- 0    # non moving individuals as 0 - following assumes only adults move
    movers_mat[,3] <- 0 
    
    move_mat <- matrix(0 , nrow = nrow(movers_mat), ncol = ncol(movers_mat))
    new_mat <- move_mat
    for (p in 1:patches){    # all adults equal p moving? Adult f and m?
      move_mat[p,] <- floor(movers_mat[p,] * (thisDmat[p,]/ (patches -1)))    # only works if Dmat 1 row (equal across sexes). 
      # numbers moving from 1 into any other
      ###CH: Why are you dividing by the number of other patches? The number of
      ###individuals leaving patch p doesn't depend on how many other patches
      ###there are.
      
      # moving INTO patch?
      new_vec <- this_mat[p,] - floor(movers_mat[p, ] * (thisDmat[p,])) + colSums(move_mat[-p,])   # those that remain in patch || vec - move out 
      # adding indivs that move in from other patches
      # if any negatives, set to 0
      
      new_mat[p, ] <- new_vec
    }
    
    # convert back to vec form? or just update meta create to accomodate matrix entry?
    
    if (any(new_mat < 0, na.rm = TRUE) || any(is.na(new_mat))) {
      new_mat[new_mat] <- 0
      new_mat[is.na(new_mat)] <- 0
    }
    
    ###CH: If no one dies while moving, then the column sums of new_mat should
    ###be the same as the column sums of this_mat. When I tested this, it seems
    ###like a bunch of badgers are dying off while moving.
    
    # creating meta mat using this years abundances
    this_meta <- meta.mat(Umat,   # vector of stage names
                          params,  # stronger params here to prevent group sizes reaching too high?
                          new_mat, # matrix with post movement abundances
                          thisDmat, 
                          stagenames, 
                          nDim, 
                          patches)
    
    # If the projection matrix has any negative values in it, stop iterating and
    # return the projection up until this point.
    if (any(sapply(this_meta, function(m) any(m < 0, na.rm = TRUE)))) {
      warning(paste("Negative values in projection matrix at time step", i,
                    "- setting negative entries to 0 and continuing."))
      meta_mat <- lapply(this_meta, function(mat) { mat[mat < 0] <- 0; mat}) 
    }
    
    # multiply by vec to project
    this_vec <- c(t(this_mat))  # turning back into long vec for multiplication
    next_vec <- floor(as.numeric(this_meta %*% this_vec))    # vector abundance result following year
    # check if any values < 0 
    if (any(next_vec < 0, na.rm = TRUE) || any(is.na(next_vec))) {
      next_vec[next_vec] <- 0
      next_vec[is.na(next_vec)] <- 0
    }      
    
    # filling next row of vec in each list
    for(p in 1:patches){
      Vec[[p]][i + 1,] <- next_vec[((p*nDim)-3) : (p*nDim)]      # the next row for each list vec
    } 
    # if any stage becomes negative, set to zero and continue
    # issues here, syntax wrong?
    if (any(sapply(Vec, function(mat) any(mat < 0, na.rm = TRUE)))) {
      warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
      Vec <- lapply(Vec, function(mat) { mat[mat < 0] <- 0; mat })
    }
    
    for(p in 1:patches){
      # calculating group size per patch this year
      Group[i + 1,] <- sapply(Vec, function(x) sum(x[i +1, ]))  # list of group size vectors 
    }
    # if any stage becomes negative, set to zero and continue
    if (any(sapply(Group, function(mat) any(mat < 0, na.rm = TRUE)))) {
      warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
      Group <- apply(Group, function(mat) { mat[mat < 0] <- 0; mat })
    }
    
    Pop[i+1] <- sum(Group[i + 1,])
    if (Pop[i] <= 0) {       # if pop size <= 0, stop and return
      warning(paste("Projection stopped at time step", i, "because total pop size reached 0 or below"))
      break
    }
  }
  
  # out objects
  out$total_pop <- Pop 
  
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  }
  if (isTRUE(return.group)) {
    out$group_size <- Group        
  }
  
  return(out)
  
}





# testing function projection with same initial vector ----
proj_test <- meta.proj(Umat,   # vector of stage names
                       params, # adjusted beta = 
                       stages,
                       initial,   # list of initial abundances - string length patches * nStages
                       cname, # length = nCol Dmat
                       time = 20,
                       return.vec = TRUE,
                       return.group = TRUE)
# group sizes and pop sizes reaching too high - max pop = 250 and max pop ~ 30?


# trying greater beta val must be greater than in non-spatial model as max group sizes are around 1/10 as large
params2 <- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
                     b=0.05,       # temp value 0.04 - 0.05 seems about right
                     rep_K= rogers_k,          # max litter size (K), 
                     h= 6)

proj_test <- meta.proj(Umat,   # vector of stage names
                       params2, # adjusted beta = 
                       stages,
                       initial,   # list of initial abundances - string length patches * nStages
                       cname, # length = nCol Dmat
                       time = 20,
                       return.vec = TRUE,
                       return.group = TRUE)


# greater number of patches might make more realistic

# generating vec for 15 patches
stagedist <- c(0.110, 0.445, 0.110, 0.335)
groups <- floor(runif(15, min = 3, max = 30))  # vector of group sizes per patch
initial_mat <- matrix(0, nrow = 15, ncol = 4)

for (i in 1: length(groups)){
initial_mat[i,] <- floor(stagedist * groups[i])
}

large_vec <- c(t(initial_mat))

proj_test <- meta.proj(Umat,   # vector of stage names
                       params2, # adjusted beta = 
                       stages,
                       large_vec,   # list of initial abundances - string length patches * nStages
                       cname, # length = nCol Dmat
                       time = 20,
                       return.vec = TRUE,
                       return.group = TRUE)

# to discuss - how to set limits for 'large and small' group sizes? At the
# moment, all finish around same size, so 10 groups of 20 more realistic, 10
# groups with sizes ranging 5 - 30? Set params per patch?



