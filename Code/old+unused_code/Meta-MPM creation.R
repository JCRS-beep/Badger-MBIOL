# meta-population matrix creation
# 10/04/26

# Creating large matrix, setting up dmat with custom function. 
# ONLY FOR ADULTS MOVING - will need editing is movement differs between sexes



# function to create large matrix - call within projection. Acts as our apply DD function to create Amat each year
meta.mat <- function(Umat,   # vector of stage names
                     params,
                     initial_vec, # big vec length = 4* npatches
                     Dmat){      # random movement probs for males and females or per patch?
  
  stagenames <- c(colnames(Umat))
  nDim <- length(stagenames)  # number of stages
  nPatches <- nrow(Dmat)
  
  # setting up out obj
  # meta_out <- list(meta_mat = matrix())  # what do we want from this?
  
  # setting base matrix dims
  mat <- matrix(0, nrow = nrow(Umat), ncol = ncol(Umat))
  
  # initial vec vs new vec after movements
  initial <- matrix(initial_vec, ncol = 4, byrow = TRUE)
  movers_mat <- initial   # we want non moving individuals as 0 - following assumes only adults move
  movers_mat[,1] <- 0
  movers_mat[,3] <- 0 
  
  move_mat <- matrix(0 , nrow = nrow(movers_mat), ncol = ncol(movers_mat))
  new_mat <- move_mat
   for (p in 1:nPatches){    # all adults equal p moving? Adult f and m?
     move_mat[p,] <- floor(movers_mat[p,] * (Dmat[p,]/ (nPatches -1)))    # only works if Dmat 1 row (equal across sexes). 
    # numbers moving from 1 into any other
   }
  
  for ( p in 1:nPatches){
  # moving INTO patch?
    new_vec <- initial[p,] - floor(movers_mat[p, ] * (Dmat[p,])) + colSums(move_mat[-p,])   # those that remain in patch || initial - move out 
    # adding indivs that move in from other patches
    # if any negatives, set to 0
    
    new_mat[p, ] <- new_vec
  }
  
  # Amat creation for more patches
  a_list <- list()
  
  
  for (p in 1:nPatches){
    thisN <- rowSums(new_mat)[p]   # using new mat for fertility calc
    thisNf <- new_mat[p, 2]
    thisNm <- new_mat[p,4]
  
    this_amat <- apply.DD(params, 
                          Umat, 
                          thisN,   # yearling and adults
                          DDapply="fertility", 
                          stagenames,   # Stages in life cycle graph 
                          thisNf,        # Adult female abundance
                          thisNm         # Adult males
    )
    
    a_list[[p]] <- this_amat    # output list length = nPatches
  }
  
  # large mat set up
  rep_cb <- function(mat, times) do.call(cbind, rep(list(mat), times))   # function to cbind repeated number of times
  meta <- rep_cb(mat, nPatches)
  
  rep_rb <- function(mat, times) do.call(rbind, rep(list(mat), times))  # function to repeatedly stack mat vertically with rbind
  meta_mat <- rep_rb(meta, nPatches)   # adds on 1 dim each loop
  
  # filling diags with Amats from list
  for (p in 1:nPatches){
    meta_mat[((p*nDim)-3):(p*nDim),((p*nDim)-3):(p*nDim)] <- a_list[[p]]
  }
  
  # More difficult to id matrix positions once cols are named
  # colnames(meta_mat) <- rep(stagenames, nPatches)  # rep times = npatches
  # rownames(meta_mat) <- rep(stagenames, nPatches)
  
  # staying probs in relevant matrix entries 
  for (p in 1 : nPatches){   # filling female and male remaining probs 
    meta_mat[,((p*nDim)-2)] <- meta_mat[,((p*nDim)-2)] * (1 - Dmat[p,])  # females remaining in patch
    meta_mat[,(p * nDim)] <- meta_mat[,(p * nDim)] * (1 - Dmat[p,])   # male staying in patch
  }
  
  # base matrix to fill with survival rates of moving indivs if all inidivuals are moving, bmat = umat
  bmat <- mat
  bmat[2,2] <- Umat[2,2]   # adult fem survival 
  bmat[4,4] <- Umat[4,4]   # adult m survival
  
  for (p in 1:nPatches){    # for each patch, work in each col of matrix
    col_lims <- ((p * nDim)-3) :(p * nDim)  # defines col boundaries for each patch 
    
    this_bmat <- bmat
    # fertility vals needed
    thisA <- alist
    this_bmat[1,2] <- a_list[[p]][1,2]   
    this_bmat[3,2] <- this_bmat[1,2]
    # dividing movement prob across remaining patches - assumes moves equally likely
    this_bmat[,2] <- this_bmat[,2] * (Dmat[p,]/ (nPatches-1))      # adult fem row * movement prob in this patch
    this_bmat[,4] <- this_bmat[,4] * (Dmat[p,]/ (nPatches-1))      # male survival * movement prob 
    
    
    # number of bmats above amat = p-1  (when p = 2, 1 bmat above Amat)
    # number of bmats below = nPatches - p   (when p = 2 and nPatches = 3, one below)
    # Define boundaries of Amat position, bmats until and after that within same cols
    
    top <- ((p * nDim) - 3)  # upper row of this Amat
    topN <- p - 1    # how many above? what if 0?
    
    bot <- p * nDim   # bottom row of this Amat
    botN <- nPatches - p   # what if 0?
    if(topN > 0){
      meta_mat[1:(top - 1), col_lims] <- rep_rb(this_bmat, topN)
    }
    if(botN > 0){
      meta_mat[(bot +1) : nrow(meta_mat), col_lims] <- rep_rb(this_bmat, botN)
    }
    
  }
  return(meta_mat)
}   # almost, first and final cols correct but bmat not included in middle. When run individually in lines, this builds correctly. Issue with loop set up?


# changing Dmat set up - only movement probs per patch, ncol  number moving indivs (just adults or all stages?)

set.seed(1)  # replicable
# alternative, random movement per patch
Dmat <- matrix(rnorm(3,mean = 0.5, sd =0.1))

# testing function 
init <- c(1,4,1,4,2,6,2,6,2,8,2,7)  # 3 patches
meta_mat <- meta.mat(Umat,   # vector of stage names
         params,
         init, # big vec length = 4* npatches
         Dmat)
         # SUCCESS!



# Defining movement function for Dmat
# base logistic equation = L / (1 - e^(-k(x-x0)))
# where L = 1 (max value, prob bound between 0 and 1)
# K = steepness of curve (vis using desmos). We want pop size as some density val?
# x0 = 0.5 (midpoint) 

# Idea = group_size calculated as proportion of some carrying capacity (max group size = )
Dmat.create <- function(colnames, nPatches) {
  Dmat <- matrix(0, ncol= length(colnames), nrow = nPatches) # row = number patches
  colnames(Dmat) <- colnames
  
  for (p in 1:nPatches){
    move <- rnorm(length(colnames), mean = 0.5, sd =0.2) # eq to calculate move prob value given vars, for now rnorm
    Dmat[p,] <- move
  }
  return(Dmat)
}

  

# how to set up more patches 
set.seed(123)  # setting our number 
reps <- 10
sizes <- runif(reps, min= 3, max=30) 
# multiply each of these pop sizes with stage dist

init_mat <- matrix(0, ncol = 4, nrow = reps) # fill each row with an inital vec per patch 
for(i in 1:reps){
  init_mat[i, ] <- floor(sizes[i]* stagedist)
}

initial <- c(t(init_mat))



# FUTURE - adding density dependent dispersal
ddDmat.create <- function(colnames, nPatches, 
                        group_size, sex_ratio = NULL, 
                        k) {
  Dmat <- matrix(0, ncol= length(colnames), nrow = nPatches) # row = number patches
  colnames(Dmat) <- colnames
  
  max_group_size <- 28 # based on papers recording 26 and 27 
  
  # logistic equation - total group size
  for (p in 1:nPatches){
    x <- group_size/ max_group_size    # proportion of max size
    move <- 1 / (1 - exp((-k(x-0.5)))) 
    Dmat[p,] <- move   # equal for all indivs or variable between sexes? vary K value for males and females?
  }
  
  if(!is.null(sex_ratio) == FALSE){
    for (p in 1:nPatches){
      x <- sex_ratio    # proportion of max size
      move <- 1 / (1 - exp((-k[2](x-0.5)))) # second value for k needed
      Dmat[p,2] <- move   # male movement prob 
    }
  }
  return(Dmat)
}