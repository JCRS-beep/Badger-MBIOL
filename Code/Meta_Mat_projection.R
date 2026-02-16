# meta projection
# 5.02.26
params<- data.frame(fmax= 0.84360,   # F fecundity max (max cubs per adult female) 
                    Sc_max=0.76,   # max cub survival (equal for sexes)
                    b=0.004,       # temp value- must be calculated from provided datasets
                    rep_K= 2.299,          #litter size (K)
                    h= 6   # harem size per male
)

# create blank base matrix
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
Umat <- matrix(0, nrow=4, ncol=4)
rownames(Umat) <- stages
colnames(Umat) <- stages
# input vital rates
Umat[2,1]<- 0.8509136    # yearling f survival
Umat[2,2]<- 0.8034251    # adult f survival
Umat[4,3]<- 0.8093391    # yearling m survival
Umat[4,4]<- 0.7490891  # adult m survival

# stating with 2 patches, dispersal between at different rates ----
n0 <- list(c(10, 5, 3, 2), c(3, 2, 2, 2))

# creating Fmat for site 1 based on abundance ----
Fmat1_out <- mating.func(params,     # density dependent parameters
                         stages,   # Stages in life cycle graph (single sex)
                         Nf = 5,        # Adult and yearling females  
                         Nm = 2,             # Adult and yearling males
                         Mfunction= "min", 
                         return.mat= TRUE)    

Fmat2_out <- mating.func(params,     
                         stages,   
                         Nf = 2,        
                         Nm = 2,           
                         Mfunction= "min",  
                         return.mat= TRUE)    # both Fmats identical!

Fmat1 <- Fmat1_out$Fmat # same lol
Fmat2 <- Fmat2_out$Fmat
# no mating in dispersing individuals: count -> breed -> move

#  streamline this section with function for Dmat creation based on distance and patch numbers
Dmat <- matrix(0, ncol= length(stages), nrow= 2) # row = number patches
colnames(Dmat) <- stages

# prob individ in patch1 remains in patch 1 - rnorm
stay1 <- rnorm(4,mean = 0.5, sd =0.2) 
Dmat1 <- Dmat
Dmat1[1,] <- stay1 # stay (and disperse) for each class in patch 1 
Dmat1[2,] <- 1- stay1

stay2 <- rnorm(4,mean = 0.5, sd =0.2)
Dmat2 <- Dmat
Dmat2[1,] <- stay2 # stay (and disperse) for each class in patch 1 
Dmat2[2,] <- 1 - stay2 

# generating metamat with 2 patches = SKIP. just for visualisation
metamat <- meta.mat(Umat,   # vector of stage names
                    Fmat1,
                    Dmat1,
                    Fmat2, 
                    Dmat2,
                    return.mats= TRUE)


# projection of 2 sub pops with this - outputs needed are vector per patch, abundance per patch, total abundance
meta.proj <- function(Umat,   # vector of stage names
                      params, # created each year instead?
                      stagenames,
                      initial,   # list of initial abundances
                      Dmat1,  # generated each year instead
                      Dmat2,
                      time = 30,
                      return.mats= TRUE){
  # Time 
  if (time<=1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))   # renaming for simplicity
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  if(is.list(initial)){   # when initial abundances is a list of vectors, we have more than 1 patch 
  patches = length(initial)
  } else if(is.list(initial) == FALSE ){  # if a vector/ not list, only 1 patch
    patches = 1
    stop("only one vector of abundances provided. Provide abundances per patch as a list, or use rem.proj for single patch projections")
  }
 # if (sum(length(n0[1]), length(n0[2])) != length(stages)) stop("initial pop vector must equal length(stages)")
 # depending on patch number, will not sum to length stages

  nStages <- length(stages)/2      # how many stages
  
  # Calculating Unions with mating function
  Nf_v <- n0[[1:patches]][nStages]     # pulls ONLY adult female entry from initial vector - vec of Nfs per patch
  Nm_v <- n0[[1:patches]][2*nStages]      # adult male in vector (final entry) - vec of Nms per patch
  
  
 
  # fmat creation first year - how to repeat for each patch
  mate <- list(matrix(0, ncol = ncol(Umat), nrow= nrow(Umat)))
  mating_out <- rep(mate, patches)   # blank Fmats to fill
  for (i in 1:patches) { # repeat for each patch
    fmat <- mating.func(params,     # density dependent parameters
                        stages,   # Stages in life cycle graph (single sex)
                        Nf_v[i],        # Adult and yearling females  
                        Nm_v[i],             # Adult and yearling males
                        Mfunction= "min", 
                        return.mat= TRUE)
    
    mating_out[i] <- fmat$Fmat   # adding appropriate Fmat to list 
    return(mating_out)  # returns list of Fmats. length= no patches
  } 
 
  # Set up the output - how to set up if multiple patches? Lists within lists...
  out <- list(pop = list(rep(NA, (time + 1))),    #  list of vectors per patch
              vec = list(matrix(0, ncol = length(stages), nrow = time + 1)))    # list of matrices per patch - can't list!

  Vec <- rep(vec, patches)  # matrices to fill with stage abundance.  row= time, col= stage. can't list matrix
  Pop <- rep(pop, patches)       # vector to fill with total pop size each year

  
  colnames(Vec) <- stages   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  
  for (m in 1:patches){  # filling first row of each list vec and pop size 
  Vec[[m]][1,] <- n0[[m]]          # for each list in Vec, first row is initial         
  
  Pop[[m]][1] <- lapply(n0, sum)  # each list is sum of n0
  }
  
  
  # Loop = matrix proj each year
  for (i in 1:time) {   # repeat for as many years as we have

        # Mating Fmat creation 
    thisNf <- as.vector(rep(NA, length(stages)))
    for (n in 1:patches){
    Nf <- sum(Vec[[i]][i,nStages-1], Vec[[i]][i,nStages])    # Nf sums yearling and adult fems in Vec matrix 
    thisNf[[i]] <- Nf  # inputting calculated value to appr place in vec
    return(thisNf)  # returns a vector of no. fems per patch
  }
  
  thisNm <- as.vector(rep(NA, length(stages)))
  for (n in 1:patches) {
    Nm <- sum(Vec[[i]][i,(2*nStages-1)], Vec[[i]][i,2*nStages])    # Nf sums yearling and adult fems in Vec matrix 
  thisNm[[i]] <- Nm  # inputting calculated value to appr place in vec
  return(thisNm)
  }


    # apply mating func to calculate pairs per patch 
  thisMating <- list(matrix(0, ncol = ncol(Umat), nrow= nrow(Umat)))
  thisMating <- rep(thisMating, patches)        # list of repeated matrix 
  for (m in 1:patches) { # repeat for each patch
    fmat <- mating.func(params,     # density dependent parameters
                        stages,   # Stages in life cycle graph (single sex)
                        thisNf[m],        # Adult and yearling females  
                        thisNm[m],             # Adult and yearling males
                        Mfunction= "min", 
                        return.mat= TRUE)
    
    thisMating[m] <- fmat$Fmat   # adding appropriate Fmat to list 
    return(thisMating)  # returns list of Fmats. length= no patches
  } 

  
  # UPDATES NEEDED -----------
    # ricker density dependence each year - pop size 
   N <- list(rep(NA, 11))
   Nlist <- rep(N, 2)
   
   for (k in 1:patches) {
    n <- rowSums(list[[k]])
    Nlist[[k]] <- n 
  
  return(Nlist) 
    }
   thisN <- Nlist[i]
    
   # loop needed to calculate Amat each year - will vary by patch 
   thisAmat <- list(matrix(0, ncol = ncol(Umat), nrow= nrow(Umat)))
   thisAmat <- rep(thisAmat, patches)        # list of repeated matrix 
    for (p in 1:patches){
      amat <- apply.DD(params, thisMating[p], Umat, thisN[p], "fertility")  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
      thisAmat[p] <- amat
      return(thisAmat)
    }
    
    # If the projection matrix has any negative values in it, stop iterating but
    # return the projection up until this point.
    if (sum(thisAmat<0)>0){
      warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
      break
    }
    # if pop size <= 0, stop and return
    if(Pop[i]<= 0){
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0)) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
    }
    
   # ------ needs updating
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
    Pop[i + 1] <- sum(Vec[(i + 1), ])
  }

}  


# easier?
meta.proj <- function(Umat,   # MAX SURVIVAL
                      initial, # vector of patch abundances
                      params, 
                      # stagenames, 
                      time, 
                      memberN=NULL,  # which individuals contribute to pop size? (as vec)
                      DDapply="fertility", 
                      Mfunction= "Min",
                      intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                      remyear = NULL,  # removal year = decrease from following year
                      return.vec= TRUE, 
                      return.remvec = FALSE,
                      Dmat1,  # generated each year instead
                      Dmat2){
  # loop proj for as many patches as we have, just change intial[i] - need to multiply by Dmat at some point!. multiply Umat by  
  
 for(i in length(intial)){
   out_obj 
   
 patch1 <- rem.proj(Umat,   # MAX SURVIVAL 
                    initial[i], 
                    params, 
                    time, 
                    memberN,  # which individuals contribute to pop size? (as vec)
                    DDapply="fertility", 
                    Mfunction= "Min",
                    intensity,  # percentage you want REMOVED from pop at time T=ry
                    remyear,  # removal year = decrease from following year
                    return.vec= TRUE, 
                    return.remvec = FALSE) # for patch 1, no removals
 }
  
 
 patch2 <- rem.proj(Umat,   # MAX SURVIVAL
                    initial[2], 
                    params, 
                    time, 
                    memberN,  # which individuals contribute to pop size? (as vec)
                    DDapply="fertility", 
                    Mfunction= "Min",
                    intensity,  # percentage you want REMOVED from pop at time T=ry
                    remyear,  # removal year = decrease from following year
                    return.vec= TRUE, 
                    return.remvec = FALSE)# for patch 2
  
  movement <-   
  
}
