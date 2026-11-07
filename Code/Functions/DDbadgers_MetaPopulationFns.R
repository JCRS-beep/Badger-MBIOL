## Functions for simple projections of any MPM, and metapopulation functions
# Written by Chrissy Hernandez
# July 2026

### Simple projection function for a metapopulation with no dispersal and no removals:
proj_meta_nodisp_norem<- function(Umat,     # must be a big matrix with diaongal Umats
                                  initial,
                                  params,
                                  stagenames,
                                  npatches,
                                  time,
                                  DDapply="fertility",
                                  return.vec=TRUE){
  # input checks
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
  
  # Check that the n_stages*n_patches = n_col of Umat:
  nStages <- length(stagenames)
  if (ncol(Umat)!=nStages*npatches){
    stop("Umat dimensions do not match the number of stages multiplied by the 
         number of patches.")
  }
  
  # Check that they used "Adult_f" and "Adult_m" in the stage names, since we
  # use a mating function that looks for those.
  I_Af<- which(stagenames=="Adult_f")
  I_Am<- which(stagenames=="Adult_m")
  I_Yf<- which(stagenames=="Yearling_f")
  I_Ym<- which(stagenames=="Yearling_m")
  if (is.null(I_Af) | is.null(I_Am) | is.null(I_Yf) | is.null(I_Ym)){
    stop("To use the mating function and density dependence, the function call 
         expects the stagenames to include the following four stages, in any order: 
         'Yearling_f', 'Adult_f', 'Yearling_m' and 'Adult_m'")
  }
  
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year

  
  colnames(Vec) <- rep(stagenames, npatches)  # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals       
  Pop[1] <- as.numeric(sum(n0))
  
  for(i in 1:npatches){
    x <- (1 +4*(i-1)) : (4*i)   # boundaries containing each patch abundances
    Group[1,i] <- sum(Vec[1, x])     # filling first row of group matrix
  }
  
  # Loop = density dependent matrix application for each year and patch
  for (i in 1:time) {   #  projection until end time
    # each timestep, reset the size of thisAmat to match Umat.
    thisAmat <- Umat
    
    for (p in 1:npatches){
      # extract the Umat for this patch:
      thisUmat <- Umat[(1+4*(p-1)):(4*p), (1+4*(p-1)):(4*p)]
      
      # Nf calculation
      subNf <- sum(Vec[i, (4*(p-1)+I_Af)], Vec[i,(4*(p-1)+I_Yf)])    # Nf sums yearling and adult fems in Vec matrix
      subNm <- sum(Vec[i, (4*(p-1)+I_Am)], Vec[i,(4*(p-1)+I_Ym)])   # Nm 
      
      # ricker density dependence each year
      subN <- sum(Vec[i,(1+4*(p-1)):(4*p)])  # pop sizes sums row i for cols corresponding to group size in patch p
      subAmat <- apply.DD(params, thisUmat, subN, DDapply, stagenames,   
                          subNf,        
                          subNm)    
      
      # If the projection matrix has negative or NA values, return message, replace with 0, and continue
      if (any(subAmat < 0, na.rm = TRUE)) {
        warning(paste("Negative values in projection matrix at time step", i,
                      "- setting negative entries to 0 and continuing."))
        subAmat[subAmat < 0] <- 0
        subAmat[is.na(subAmat)] <- 0
      }
      
      # Insert the subAmat into the correct place in the matrix (on the diagonal block):
      thisAmat[(1+4*(p-1)):(4*p), (1+4*(p-1)):(4*p)]<- subAmat
      
    }
    # multiplying to project
    Vec[(i + 1), ] <- floor(thisAmat %*% Vec[i, ])  # amat values multiplied by vec, round down for integers  
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
      Vec[is.na(Vec)] <- 0
    }
    
    # group sizes next year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[i + 1, (1+4*(p-1)):(4*p)])   # previous indexing issue - vec[i,] keeps calculating this years size
    }
    
    # filling in group size this year
    Group[i + 1,] <- sizes
    
    Pop[i + 1] <- sum(Vec[(i + 1), ])
    
    # if pop size <= 0, stop and return message with year
    if(Pop[i]<= 0 || is.na(Pop[i])) {
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
  }
  # output objects
  if (isTRUE(return.vec)) {
    out <- list(pop = vector(), 
                vec = matrix())
    out$pop <- Pop
    out$group <- Group
    out$vec <- Vec 
    
    
    return(out)
  } else {
    return(Pop)
  }
}

### Simple projection function for a metapopulation with fixed dispersal
### probability and no removals. 
# Note, we are making the following assumptions: (1) a fixed probability of
# dispersing, (2) dispersing occurs before any survival, growth, or
# reproduction, (3) dispersers are equally likely to arrive in any of the other
# patches besides the one they started in.
proj_meta_fixeddisp_norem <- function(Umat,
                                     initial,
                                     params,
                                     stagenames,
                                     npatches,
                                     time,
                                     DDapply="fertility",
                                     dispersal_prob = 0.1,
                                     dispersal_stages = c("Adult_f", "Adult_m"),
                                     return.vec=TRUE){
  # input checks
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
  
  # Check that the n_stages*n_patches = n_col of Umat:
  nStages<- length(stagenames)
  if (ncol(Umat)!=nStages*npatches){
    stop("Umat dimensions do not match the number of stages multiplied by the 
         number of patches.")
  }
  
  # Check that they used "Adult_f" and "Adult_m" in the stage names, since we
  # use a mating function that looks for those.
  I_Af<- which(stagenames=="Adult_f")
  I_Am<- which(stagenames=="Adult_m")
  I_Yf<- which(stagenames=="Yearling_f")
  I_Ym<- which(stagenames=="Yearling_m")
  if (is.null(I_Af) | is.null(I_Am) | is.null(I_Yf) | is.null(I_Ym)){
    stop("To use the mating function and density dependence, the function call 
         expects the stagenames to include the following four stages, in any order: 
         'Yearling_f', 'Adult_f', 'Yearling_m' and 'Adult_m'")
  }
  
  # Check that the stages identified as dispersers are found in stagenames:
  I_leavers<- which(stagenames %in% dispersal_stages)
  if (length(I_leavers)!= length(dispersal_stages)){
    stop("Check names of dispersal stages - the number of stagenames identified 
         as dispersers does not match the number of dispersal stages specified.")
  }
  
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year
  Leavers <- matrix(0, ncol = ncol(Vec), nrow = time)   # matrix to fill with the number of individuals leaving each patch each year (by each stage)
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming cls matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals  
  Pop[1] <- as.numeric(sum(n0))
 
  for(p in 1:npatches){
    x <- (1 +4*(p-1)) : (4*p)   # boundaries containing each patch abundances
    Group[1,p] <- sum(Vec[1, x])     # filling first row of group matrix
  }
  
  # Loop = density dependent matrix application for each year and patch
  # storing dispersal history in array to track origin and destination
  for (i in 1:time) {   #  projection until end time
    # setting up array to store stages moving, origin and destination
    moves <- array(0, dim = c(npatches, npatches, length(stagenames)),
                   dimnames = list(from = 1:npatches, to = 1:npatches,
                                   stage = stagenames))
    
    leavers <- rep(0, ncol(Vec))    # individuals leaving patches - to go in output
    arrivers <- rep(0, ncol(Vec))   # individuals arriving in other patch
    
    # dispersal loop 
    for (p in 1:npatches){
      Ipatch <- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      subVec <- as.vector(Vec[i, Ipatch])    # have to specify row (year) we want, to remove names format as numeric vector
      
      
      # apply the probability of leaving to the stages that leave:
      subleavers<- rep(0,nStages)    # leavers from this patch
      subleavers[I_leavers]<- round(dispersal_prob*subVec[I_leavers])   # round instead of floor - does not underestimate moving individuals
      leavers[Ipatch] <- subleavers
      
      arrival_patches <- setdiff(1:npatches, p)   # setdiff compares differences between dataframes, removes p from all other patches
      
      for (s in I_leavers) {
        n <- subleavers[s]
        if (n > 0) {
          # index into arrival_patches 
          dest <- arrival_patches[sample.int(length(arrival_patches), n, replace = TRUE)]
          moves[p, , s] <- tabulate(dest, nbins = npatches)   # from p to  dest
        }
      }
    }
    
    # arrivers = sum over source patches -> (to x stage), flattened to match Vec layout
    to_by_stage <- apply(moves[, , ], c(2, 3), sum)    # syntaxing errs here
    arrivers <- as.vector(t(to_by_stage))
    
    # check : everyone who left must arrive somewhere
    if(sum(arrivers) != sum(leavers)) stop("Numbers of leavers and arrivers not equal, issues in dispersal calculation")
    
       # DEBUGGING : closing bracket in wrong place - reproduction was calculated before dispersal vec 
        # was completed for all patches
    
    postdispVec<- Vec[i,] - leavers + arrivers    # works in long vec - dont need to loop to calculate for each patch
    Leavers[i,] <- leavers # assigning leavers vec to relevant row in Leavers matrix
    
    # now we do another loop over the patches to calculate reproduction:
    for (p in 1:npatches){
      Ipatch<- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      # extract the Umat for this patch:
      thisUmat<- Umat[Ipatch, Ipatch]
      
      # Nf calculation
      subNf <- sum(postdispVec[(4*(p-1)+I_Af)], postdispVec[(4*(p-1)+I_Yf)])    # Nf sums yearling and adult fems in Vec matrix
      subNm <- sum(postdispVec[(4*(p-1)+I_Am)], postdispVec[(4*(p-1)+I_Ym)])   # Nm 
      
      # ricker density dependence
      subN <- sum(postdispVec[Ipatch])  # pop sizes sums row i for cols corresponding to patch p
      subAmat <- apply.DD(params, thisUmat, subN, DDapply, stagenames,   
                          subNf,        
                          subNm)    
      
      # If the projection matrix has negative or NA values, return message, replace with 0, and continue
      if (any(subAmat < 0, na.rm = TRUE)) {
        warning(paste("Negative values in projection matrix at time step", i,
                      "- setting negative entries to 0 and continuing."))
        subAmat[subAmat < 0] <- 0
        subAmat[is.na(subAmat)] <- 0
      }
      
      # Project just this patch forward:
      Vec[(i+1), Ipatch]<- floor(subAmat %*% postdispVec[Ipatch]) # amat values multiplied by vec, round down for integers 
                                            # DEBUGGING:
                                            # should multiply subAmat by POST DISP vec (not vec) - rep values based on new abundances, multiply rep values by 
                                            # post mvoement vec
    }
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
      Vec[is.na(Vec)] <- 0
    }
    
    # group sizes next year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[i + 1, (1+4*(p-1)):(4*p)])
    }
    
    # filling in group size this year
    Group[i + 1,] <- sizes
    
    # pop sizes by year
    Pop[i + 1] <- sum(Vec[(i + 1), ])
    
    # if pop size <= 0, stop and return message with year
    if(Pop[i]<= 0 || is.na(Pop[i])) {
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
   }
  
  # output objects
  if (isTRUE(return.vec)) {
    out <- list(pop = vector(), 
                vec = matrix(), 
                group = matrix(),
                leavers = matrix()
                )
    
    out$pop <- Pop
    out$vec <- Vec 
    out$group <- Group
    out$leavers <- Leavers
    
    return(out)
  } else {
    return(Pop)
  }
  }

   

# DD dispersal probability matrix ----
# social groups range 1 - 27, not observed larger than this,
# Dispersal prob calculated with logistic between 0 - 0.3, 
# Dispersal equal for males and females in each patch
ddDmat <- function(stagenames,dispersal_stages, 
                   npatches, 
                   group_size) {    # what is the carrying capacity (can be vector (k per patch))
  
  # getting the index of each stage
  I_Af<- which(stagenames=="Adult_f")
  I_Am<- which(stagenames=="Adult_m")
  I_Yf<- which(stagenames=="Yearling_f")
  I_Ym<- which(stagenames=="Yearling_m")
  
  if (is.null(I_Af) | is.null(I_Am) | is.null(I_Yf) | is.null(I_Ym)){
    stop("stagenames must include the following four stages, in any order: 
         'Yearling_f', 'Adult_f', 'Yearling_m' and 'Adult_m'")
  }
  
  if(length(group_size) != npatches){
    stop("You must provided the sizes for each group")
  }
  
  I_leavers <- which(stagenames %in% dispersal_stages)
  
  Dmat <- matrix(0, ncol = length(stagenames), nrow = npatches) # row = number patches
  colnames(Dmat) <- stagenames
  
  max_move <- 0.3  # greatest proportion of moving individuals?
  xhalf <- 15
  steepness = 0.5
  for (p in 1:npatches){
    p_move<- round(max_move/(1+exp(-steepness*(group_size-xhalf))), 3)       
    
    p_move[p_move<0] <- 0   # setting any negatives to 0
    p_move[p_move >1] <- 1
    Dmat[p,I_leavers] <- p_move[p]   # equal for all indivs or variable between sexes? vary K value for males and females?
  }
  
  return(Dmat)
} 


# DD moving projection function -----
# incorporating density dependent movement probs in projection function
proj_meta_DDdisp_norem <- function(Umat,
                                   initial,
                                   params,
                                   stagenames,
                                   npatches,
                                   time,
                                   DDapply="fertility",
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec=TRUE){
  # input checks
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
  
  # Check that the n_stages*n_patches = n_col of Umat:
  nStages<- length(stagenames)
  if (ncol(Umat)!=nStages*npatches){
    stop("Umat dimensions do not match the number of stages multiplied by the 
         number of patches.")
  }
  
  # Check that they used "Adult_f" and "Adult_m" in the stage names, since we
  # use a mating function that looks for those.
  I_Af<- which(stagenames=="Adult_f")
  I_Am<- which(stagenames=="Adult_m")
  I_Yf<- which(stagenames=="Yearling_f")
  I_Ym<- which(stagenames=="Yearling_m")
  if (is.null(I_Af) | is.null(I_Am) | is.null(I_Yf) | is.null(I_Ym)){
    stop("To use the mating function and density dependence, the function call 
         expects the stagenames to include the following four stages, in any order: 
         'Yearling_f', 'Adult_f', 'Yearling_m' and 'Adult_m'")
  }
  
  # Check that the stages identified as dispersers are found in stagenames:
  I_leavers<- which(stagenames %in% dispersal_stages)
  if (length(I_leavers)!= length(dispersal_stages)){
    stop("Check names of dispersal stages - the number of stagenames identified 
         as dispersers does not match the number of dispersal stages specified.")
  }
  
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year
  Leavers <- matrix(0, ncol = ncol(Vec), nrow = time)   # matrix to fill with the number of individuals leaving each patch each year (by each stage)
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming cls matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals  
  Pop[1] <- as.numeric(sum(n0))
  
  for(p in 1:npatches){
    x <- (1 +4*(p-1)) : (4*p)   # boundaries containing each patch abundances
    Group[1,p] <- sum(Vec[1, x])     # filling first row of group matrix
  }
  
  # Loop = density dependent matrix application for each year and patch
  for (i in 1:time) {   #  projection until end time
    # group sizes this year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[i , (1+4*(p-1)):(4*p)])
    }
    
    
    dmat <- ddDmat(stagenames,dispersal_stages, # names of all stages, names of only dispersing stages
                   npatches, # number of patches
                   group_size = sizes)
    
    # record of individuals moving from group p to p'
    moves <- array(0, dim = c(npatches, npatches, length(stagenames)),
                   dimnames = list(from = 1:npatches, to = 1:npatches,
                                   stage = stagenames))
    
    leavers <- rep(0, ncol(Vec))    # individuals leaving patches - to go in output
    arrivers <- rep(0, ncol(Vec))   # individuals arriving in other patch
    
    # dispersal loop 
    for (p in 1:npatches){
      Ipatch <- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      subVec <- as.vector(Vec[i, Ipatch])    # have to specify row (year) we want, to remove names format as numeric vector
      
      
      # apply the probability of leaving to the stages that leave:
      subleavers<- rep(0,nStages)    # leavers from this patch
      
      # DD MATRIX APPLICATION HERE
      subleavers[I_leavers]<- round(dmat[p,I_leavers]*subVec[I_leavers])   # round instead of floor - does not underestimate moving individuals
      
      
      leavers[Ipatch] <- subleavers
      
      arrival_patches <- setdiff(1:npatches, p)   # setdiff compares differences between dataframes, removes p from all other patches
      
      for (s in I_leavers) {
        n <- subleavers[s]
        if (n > 0) {
          # index into arrival_patches 
          dest <- arrival_patches[sample.int(length(arrival_patches), n, replace = TRUE)]
          moves[p, , s] <- tabulate(dest, nbins = npatches)   # from p to  dest
        }
      }
    }  # Brackets now in correct places
    
    # arrivers = sum over source patches -> (to x stage), flattened to match Vec layout
    to_by_stage <- apply(moves[, , ], c(2, 3), sum)    # syntaxing errs here
    arrivers <- as.vector(t(to_by_stage))
    
    # check : everyone who left must arrive somewhere
    if(sum(arrivers) != sum(leavers)) stop("Numbers of leavers and arrivers not equal, issues in dispersal calculation")
    
    postdispVec<- Vec[i,] - leavers + arrivers    # works in long vec - dont need to loop to calculate for each patch
    Leavers[i,] <- leavers # assigning leavers vec to relevant row in Leavers matrix
    
    # now we do another loop over the patches to calculate reproduction:
    for (p in 1:npatches){
      Ipatch<- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      # extract the Umat for this patch:
      thisUmat<- Umat[Ipatch, Ipatch]
      
      # Nf calculation
      subNf <- sum(postdispVec[(4*(p-1)+I_Af)], postdispVec[(4*(p-1)+I_Yf)])    # Nf sums yearling and adult fems in Vec matrix
      subNm <- sum(postdispVec[(4*(p-1)+I_Am)], postdispVec[(4*(p-1)+I_Ym)])   # Nm 
      
      # ricker density dependence
      subN <- sum(postdispVec[Ipatch])  # pop sizes sums row i for cols corresponding to patch p
      subAmat <- apply.DD(params, thisUmat, subN, DDapply, stagenames,   
                          subNf,        
                          subNm)    
      
      # If the projection matrix has negative or NA values, return message, replace with 0, and continue
      if (any(subAmat < 0, na.rm = TRUE)) {
        warning(paste("Negative values in projection matrix at time step", i,
                      "- setting negative entries to 0 and continuing."))
        subAmat[subAmat < 0] <- 0
        subAmat[is.na(subAmat)] <- 0
      }
      
      # Project just this patch forward:
      Vec[(i+1), Ipatch]<- floor(subAmat %*% postdispVec[Ipatch]) # amat values multiplied by vec, round down for integers 
      # DEBUGGING:
      # should multiply subAmat by POST DISP vec (not vec) - rep values based on new abundances, multiply rep values by 
      # post mvoement vec
    }
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
      Vec[is.na(Vec)] <- 0
    }
    
    # group sizes next year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[i + 1, (1+4*(p-1)):(4*p)])
    }
    
    # filling in group size this year
    Group[i + 1,] <- sizes
    
    # pop sizes by year
    Pop[i + 1] <- sum(Vec[(i + 1), ])
    
    # if pop size <= 0, stop and return message with year
    if(Pop[i]<= 0 || is.na(Pop[i])) {
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
  }
  
  # output objects
  if (isTRUE(return.vec)) {
    out <- list(pop = vector(), 
                vec = matrix(), 
                group = matrix(),
                leavers = matrix()
    )
    
    out$pop <- Pop
    out$vec <- Vec 
    out$group <- Group
    out$leavers <- Leavers
    
    return(out)
  } else {
    return(Pop)
  }
}


### REMOVAL FUNCTIONS--
# truncated normal dist function to sample from 
# Sample n integers between min and max with a bell-shaped (truncated normal) distribution
# mean = midpoint of min and max 
# sd defaults to (max - min)/4, so min and max sit ~2 SDs from the mean
trunc.norm <- function(n = 1, min, max, mean = (min + max)/2, sd = (max - min)/3) {
  
  min <- max(0, min)       # can never remove fewer than 0 individuals
  max <- max(min, max)     # guards against max < min
  #  mean = (min + max)/2
  #  sd = (max - min)/4
  if (max == min) return(rep(round(min), n))   # no range to sample from
  
  out <- numeric(n)
  for (k in 1:n) {
    x <- rnorm(1, mean, sd)
    # resample until the draw falls in range (true truncated normal)
    while (x < min || x > max) {
      x <- rnorm(1, mean, sd)
    }
    out[k] <- x
  }
  
  round(out)   # integer number of individuals
}

# first removal function without dispersal
proj_meta_nodisp_rem <- function(Umat,   
                                 initial, 
                                 stagedist,  # proportion of each stage
                                 params, 
                                 stagenames, 
                                 npatches,
                                 time,     
                                 DDapply="fertility", 
                                 remyear = NULL,  # removal year = vector of years 
                                 rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                                 bias = NULL , # strength of bias as percentage (range??)
                                 return.vec= TRUE, 
                                 return.remvec = FALSE){
  # input checks
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
  
  nStages <- length(stagenames)
  
  # Check that the n_stages*n_patches = n_col of Umat:
  if (ncol(Umat)!=nStages*npatches){
    stop("Umat dimensions do not match the number of stages multiplied by the 
         number of patches.")
  }
  
  
  # Check that they used "Adult_f" and "Adult_m" in the stage names, since we
  # use a mating function that looks for those.
  I_Af<- which(stagenames=="Adult_f")
  I_Am<- which(stagenames=="Adult_m")
  I_Yf<- which(stagenames=="Yearling_f")
  I_Ym<- which(stagenames=="Yearling_m")
  if (is.null(I_Af) | is.null(I_Am) | is.null(I_Yf) | is.null(I_Ym)){
    stop("To use the mating function and density dependence, the function call 
         expects the stagenames to include the following four stages, in any order: 
         'Yearling_f', 'Adult_f', 'Yearling_m' and 'Adult_m'")
  }
  
  # removal set up
  if(is.numeric(remyear) == TRUE){
    # ensure remyear are integers and between 1 - time
    remyear <- unique(as.integer(remyear))
    remyear <- remyear[remyear >= 1 & remyear <= time] 
    if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
    
    nYears <- length(remyear)   # over how many years are removals taking place?  remove year_goal for each value of rem_index
    # Map removal year -> index into Remvec / Nremoved:
    rem_index <- setNames(seq_along(remyear), remyear)  # takes length rem year, names assigned as seq of integers 
  }
  else (rem_index <- NULL)
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year
  
  if (is.numeric(remyear) == TRUE) {
    Nremoved <- numeric(length(remyear))      # how many removed total each year
    Remgroup <- matrix(0, nrow = length(remyear), ncol = npatches)    # how many removed per group per year (col per group)
    Remvec <- matrix(0, nrow = length(remyear), ncol =nStages*npatches)  
  }
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals       
  Pop[1] <- as.numeric(sum(n0))
  
  for(p in 1:npatches){
    x <- (1 +4*(p-1)) : (4*p)   # boundaries containing each patch abundances
    Group[1,p] <- sum(Vec[1, x])     # filling first row of group matrix
  }
  
  
  # Loop = density dependent matrix application for each year and patch
  for (i in 1:time) {   #  projection until end time
    # each timestep, reset the size of thisAmat to match Umat.
    thisAmat <- Umat
    
    for (p in 1:npatches){
      # extract the Umat for this patch:
      thisUmat <- Umat[(1+4*(p-1)):(4*p), (1+4*(p-1)):(4*p)]
      
      # Nf calculation
      subNf <- sum(Vec[i, (4*(p-1)+I_Af)], Vec[i,(4*(p-1)+I_Yf)])    # Nf sums yearling and adult fems in Vec matrix
      subNm <- sum(Vec[i, (4*(p-1)+I_Am)], Vec[i,(4*(p-1)+I_Ym)])   # Nm 
      
      # ricker density dependence each year
      subN <- sum(Vec[i,(1+4*(p-1)):(4*p)])  # pop sizes sums row i for cols corresponding to group size in patch p
      subAmat <- apply.DD(params, thisUmat, subN, DDapply, stagenames,   
                          subNf,        
                          subNm)    
      
      # If the projection matrix has negative or NA values, return message, replace with 0, and continue
      if (any(subAmat < 0, na.rm = TRUE)) {
        warning(paste("Negative values in projection matrix at time step", i,
                      "- setting negative entries to 0 and continuing."))
        subAmat[subAmat < 0] <- 0
        subAmat[is.na(subAmat)] <- 0
      }
      
      # Insert the subAmat into the correct place in the matrix (on the diagonal block):
      thisAmat[(1+4*(p-1)):(4*p), (1+4*(p-1)):(4*p)]<- subAmat
    }  
    
    # multiplying to project
    Vec[(i + 1), ] <- round(thisAmat %*% Vec[i, ])  # amat values multiplied by vec, round down for integers  
    # CHANGED FROM FLOOR TO ROUND TO FIX RECOVERY ISSUES IN PROJECTION
    # overwritten when this is a remyear
    
    # checking for removal years --------
    if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
      idx <- rem_index[as.character(i)]     # idx = 1: remyear. What removal year are we 
      
      rem <- vector(length = nStages*npatches)       # blank vector to fill as removal vector
      
      # stage distribution this year to sample probability of removal
      thisDist <- matrix(Vec[i-1,]/Pop[i-1], nrow = npatches, ncol = nStages)     # getting a stagedist for this year - needs averaging across groups
      avStageDist <- round(colMeans(thisDist), 2) # place in table and average across cols
      
      if(rem_strat == "random"){
        dist <- avStageDist
      }
      
      
      ## bias vector for different removal strategies
      # sex biased
 #     if(rem_strat == )
      
      # stage biased
      
      # sett biased
      
      # FIRST YEAR
      if(idx == 1){       # when idx = 1, remove 70% with goal 
        preN <- Pop[remyear[1]]  # population size pre cull # move until i
        Nrem <- round(preN * 0.7)   # goal to remove  - 70% of pop size before first rem  
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        group_remvec <- rep(round(Nrem/ npatches), npatches)
        
        #using stagedist to remove group totals across stages 
        for (p in 1:length(group_remvec)){   
          rem[(1 +4*(p-1)) : (4*p)] <- group_remvec[p] * dist
        }    
        
        thisRemvec <- round(rem)   # actual removed = rounded abundance * dist
        
        # SUBSEQUENT YEARS - very low numbers removed
      }  else {
        if(idx > 1 && idx <= 5){    # if supplementary, change strat to calculate max and min
          estN <- 0.8 * preN - sum(Nremoved)     # ISSUE : negative estN?
          # Use estimate from year before culls minus the culled N since. reduce by 1/5 (from published methods)
          estN[estN < 0] <- 0                                               
          # set min and max removals
          base <- 0.36 * estN                # base using estimated pop size
          min = base - (0.25*estN)          # difference between min and max = 25% of pre-cull total
          max = base + (0.25*estN)
        }
        
        # SUPPLEMENTARY CULLS - ISSUE this part removes too many
        else if(idx > 5){    # if SBC, change strat to calculate max and min. Needs to be much lower intensity?
          base <- 0.1 * Nremoved[1]   # 36% of first year cull total - SBC. CHanging this to 10% because of low population sizes
          # set min and max removals
          min = base - (0.125*preN)          # difference between min and max = 25% of pre-cull total
          max = base                         # set maximum as as base to avoid large removals freqs?
          
        }
        # need to prevent minimum falling bellow 0 - removing negative individuals
        min[min < 0] <- 0
        
        # set seed for repeatability and use a normal dist instead
        Nrem <- trunc.norm(1, min, max) # how many to remove this year using truncated norm dist between min and max values
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
       
        # sampling the patches that we remove from
        rem_mat <- matrix(0, ncol = (nStages*npatches), nrow = Nrem)    # mat that tracks each removal stage and patch for each individual removed
        
        group_rem <- vector(length = npatches)     # vec of how many removals from each group 
        stage_rem <- matrix(0, nrow = npatches, ncol = nStages)    # our removal matrix, tracking which stage is removed and from what group
        
        for (x in 1:nrow(rem_mat)){
          a <- sample(c(1:npatches), 1)        # patch sample
          s <- sample(c(1:4), 1, prob = avStageDist)        # stage sample
          rem_mat[x, (s + 4*(a-1))] <- 1    
        }
        
        thisRemvec <- colSums(rem_mat)
        
      }
      
      # calculating actual number removed from each stage for output
      
      
      group_rem <- vector()       # previous indexing issue - vec[i,] keeps calculating this years size
      for (p in 1:npatches){
        group_rem[p] <- sum(thisRemvec[(1+4*(p-1)):(4*p)])      # incorrect indexing here. We want to sum every 4 entries in the vector
      }# removed per group
      
      
      #  into outputs
      Nremoved[idx] <- sum(thisRemvec)  # total number removed   
      Remvec[idx,] <- thisRemvec
      Remgroup[idx,] <- group_rem
      
      # new population size following removals
      Vec[(i + 1), ] <- floor(thisAmat %*% Vec[i, ]) - thisRemvec  # year after remyear filled with new stage vec. Why is this not filling correct year?
    }
    
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
      Vec[is.na(Vec)] <- 0
    }
    
    # group sizes next year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[i + 1, (1+4*(p-1)):(4*p)])  
    }
    
    # filling in group size this year
    Group[i + 1,] <- sizes
    
    Pop[i + 1] <- sum(Vec[(i + 1), ])
    
    # if pop size <= 0, stop and return message with year
    if(Pop[i]<= 0 || is.na(Pop[i])) {
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
  }  
  
  # output objects
  if (isTRUE(return.vec)) {    # if we want remvec only
    out <- list(pop = vector(), 
                group = matrix(),
                vec = matrix())
    out$pop <- Pop
    out$group <- Group
    out$vec <- Vec 
    
    if(isTRUE(return.vec) && isTRUE(return.remvec)){   # remvec and 
      out$Nremoved <- Nremoved
      out$Group_remvec <- Remgroup
      out$Remvec <- Remvec
    }
    
    return(out)
  } else {
    return(Pop)
  }
}


# incorporating density dependent dispersal into removal function
proj_meta_DDdisp_rem <- function(Umat,
                                 initial,
                                 stagedist,  # proportion of each stage
                                 params,
                                 stagenames,
                                 npatches,
                                 time,
                                 DDapply="fertility",
                                 dispersal_stages = c("Adult_f", "Adult_m"),
                                 remyear = NULL,  # removal year = vector of years 
                                 rem_strat = "random",  # if specified 
                                 bias = NULL, 
                                 return.vec=TRUE, 
                                 return.remvec = FALSE){
  # input checks
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
  
  # Check that the n_stages*n_patches = n_col of Umat:
  nStages<- length(stagenames)
  if (ncol(Umat)!=nStages*npatches){
    stop("Umat dimensions do not match the number of stages multiplied by the 
         number of patches.")
  }
  
  # Check that they used "Adult_f" and "Adult_m" in the stage names, since we
  # use a mating function that looks for those.
  I_Af<- which(stagenames=="Adult_f")
  I_Am<- which(stagenames=="Adult_m")
  I_Yf<- which(stagenames=="Yearling_f")
  I_Ym<- which(stagenames=="Yearling_m")
  if (is.null(I_Af) | is.null(I_Am) | is.null(I_Yf) | is.null(I_Ym)){
    stop("To use the mating function and density dependence, the function call 
         expects the stagenames to include the following four stages, in any order: 
         'Yearling_f', 'Adult_f', 'Yearling_m' and 'Adult_m'")
  }
  
  # Check that the stages identified as dispersers are found in stagenames:
  I_leavers<- which(stagenames %in% dispersal_stages)
  if (length(I_leavers)!= length(dispersal_stages)){
    stop("Check names of dispersal stages - the number of stagenames identified 
         as dispersers does not match the number of dispersal stages specified.")
  }
  
  # removal set up
  if(is.numeric(remyear) == TRUE){
    # ensure remyear are integers and between 1 - time
    remyear <- unique(as.integer(remyear))
    remyear <- remyear[remyear >= 1 & remyear <= time] 
    if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
    
    nYears <- length(remyear)   # over how many years are removals taking place?  remove year_goal for each value of rem_index
    # Map removal year -> index into Remvec / Nremoved:
    rem_index <- setNames(seq_along(remyear), remyear)  # takes length rem year, names assigned as seq of integers 
  }
  else (rem_index <- NULL)
  
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year
  Leavers <- matrix(0, ncol = ncol(Vec), nrow = time)   # matrix to fill with the number of individuals leaving each patch each year (by each stage)
  
  if (is.numeric(remyear) == TRUE) {
    Nremoved <- numeric(length(remyear))      # how many removed total each year
    Remgroup <- matrix(0, nrow = length(remyear), ncol = npatches)    # how many removed per group per year (col per group)
    Remvec <- matrix(0, nrow = length(remyear), ncol =nStages*npatches)  
  }
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming cls matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals  
  Pop[1] <- as.numeric(sum(n0))
  
  for(p in 1:npatches){
    x <- (1 +4*(p-1)) : (4*p)   # boundaries containing each patch abundances
    Group[1,p] <- sum(Vec[1, x])     # filling first row of group matrix
  }
  
  # Loop = density dependent matrix application for each year and patch
  for (i in 1:time) {   #  projection until end time
    # each timestep, reset the size of thisAmat to match Umat.
    thisAmat <- Umat
    
    # group sizes this year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[i , (1+4*(p-1)):(4*p)])
    }
    
    
    dmat <- ddDmat(stagenames,dispersal_stages, # names of all stages, names of only dispersing stages
                   npatches, # number of patches
                   group_size = sizes)
    
    # record of individuals moving from group p to p'
    moves <- array(0, dim = c(npatches, npatches, length(stagenames)),
                   dimnames = list(from = 1:npatches, to = 1:npatches,
                                   stage = stagenames))
    
    leavers <- rep(0, ncol(Vec))    # individuals leaving patches - to go in output
    arrivers <- rep(0, ncol(Vec))   # individuals arriving in other patch
    
    # dispersal loop 
    for (p in 1:npatches){
      Ipatch <- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      subVec <- as.vector(Vec[i, Ipatch])    # have to specify row (year) we want, to remove names format as numeric vector
      
      
      # apply the probability of leaving to the stages that leave:
      subleavers<- rep(0,nStages)    # leavers from this patch
      
      # DD MATRIX APPLICATION HERE
      subleavers[I_leavers]<- round(dmat[p,I_leavers]*subVec[I_leavers])   # round instead of floor - does not underestimate moving individuals
      
      
      leavers[Ipatch] <- subleavers
      
      arrival_patches <- setdiff(1:npatches, p)   # setdiff compares differences between dataframes, removes p from all other patches
      
      for (s in I_leavers) {
        n <- subleavers[s]
        if (n > 0) {
          # index into arrival_patches 
          dest <- arrival_patches[sample(c(1:length(arrival_patches)), size = n, replace = TRUE)]
          moves[p, , s] <- tabulate(dest, nbins = npatches)   # from p to  dest
        }
      }
    }  # Brackets now in correct places
    
    # arrivers = sum over source patches -> (to x stage), flattened to match Vec layout
    to_by_stage <- apply(moves[, , ], c(2, 3), sum)    # syntaxing errs here
    arrivers <- as.vector(t(to_by_stage))
    
    # check : everyone who left must arrive somewhere
    if(sum(arrivers) != sum(leavers)) stop("Numbers of leavers and arrivers not equal, issues in dispersal calculation")
    
    postdispVec<- Vec[i,] - leavers + arrivers    # works in long vec - dont need to loop to calculate for each patch
    Leavers[i,] <- leavers # assigning leavers vec to relevant row in Leavers matrix
    
    # another loop over the patches to calculate reproduction:
    for (p in 1:npatches){
      Ipatch<- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      # extract the Umat for this patch: why do we need this step?
      thisUmat<- Umat[Ipatch, Ipatch]
      
      # Nf calculation
      subNf <- sum(postdispVec[(4*(p-1)+I_Af)], postdispVec[(4*(p-1)+I_Yf)])    # Nf sums yearling and adult fems in Vec matrix
      subNm <- sum(postdispVec[(4*(p-1)+I_Am)], postdispVec[(4*(p-1)+I_Ym)])   # Nm 
      
      # ricker density dependence
      subN <- sum(postdispVec[Ipatch])  # pop sizes sums row i for cols corresponding to patch p
      subAmat <- apply.DD(params, thisUmat, subN, DDapply, stagenames,   
                          subNf,        
                          subNm)    
      
      # If the projection matrix has negative or NA values, return message, replace with 0, and continue
      if (any(subAmat < 0, na.rm = TRUE)) {
        warning(paste("Negative values in projection matrix at time step", i,
                      "- setting negative entries to 0 and continuing."))
        subAmat[subAmat < 0] <- 0
        subAmat[is.na(subAmat)] <- 0
      }
      
      # Insert the subAmat into the correct place in the matrix (on the diagonal block):
      thisAmat[Ipatch, Ipatch]<- subAmat
    }  
    
    
    # Project just this patch forward:
    Vec[(i+1), ] <- floor(thisAmat %*% postdispVec) # amat values multiplied by vec, round down for integers 
    
    # checking for removal years --------
    # checking for removal years --------
    if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
      idx <- rem_index[as.character(i)]     # idx = 1: remyear. What removal year are we 
      
      rem <- vector(length = nStages*npatches)       # blank vector to fill as removal vector
      
      # FIRST YEAR
      if(idx == 1){       # when idx = 1, remove 70% with goal 
        preN <- Pop[remyear[1]]  # population size pre cull # move until i
        Nrem <- round(preN * 0.7)   # goal to remove  - 70% of pop size before first rem  
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        group_remvec <- rep(round(Nrem/ npatches), npatches)
        
        #using stagedist to remove group totals across stages 
        for (p in 1:length(group_remvec)){   
          rem[(1 +4*(p-1)) : (4*p)] <- group_remvec[p] * stagedist
        }    
        
        thisRemvec <- round(rem)   # actual removed = rounded abundance * dist
        
        # SUBSEQUENT YEARS - very low numbers removed
      }  else {
        if(idx > 1 && idx <= 5){    # if supplementary, change strat to calculate max and min
          estN <- 0.8 * preN - sum(Nremoved)     # ISSUE : negative estN?
          # Use estimate from year before culls minus the culled N since. reduce by 1/5 (from published methods)
          estN[estN < 0] <- 0                                               
          # set min and max removals
          base <- 0.36 * estN                # base using estimated pop size
          min = base - (0.25*estN)          # difference between min and max = 25% of pre-cull total
          max = base + (0.25*estN)
        }
        
        # SUPPLEMENTARY CULLS - ISSUE this part removes too many
        else if(idx > 5){    # if SBC, change strat to calculate max and min. Needs to be much lower intensity?
          base <- 0.1 * Nremoved[1]   # 36% of first year cull total - SBC. CHanging this to 10% because of low population sizes
          # set min and max removals
          min = base - (0.125*preN)          # difference between min and max = 25% of pre-cull total
          max = base                         # set maximum as as base to avoid large removals freqs?
          
        }
        # need to prevent minimum falling bellow 0 - removing negative individuals
        min[min < 0] <- 0
        
        # set seed for repeatability and use a normal dist instead
        Nrem <- trunc.norm(1, min, max) # how many to remove this year using truncated norm dist between min and max values
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
        # NEEDS WORK - HOW TO TRACK WHICH PATCH AND STAGE WE REMOVE FROM TO MAKE REMVEC------
        thisDist <- matrix(Vec[i-1,]/Pop[i-1], nrow = npatches, ncol = nStages)     # getting a stagedist for this year - needs averaging across groups
        avDist <- round(colMeans(thisDist), 2) # place in table and average across cols
        
        # sampling the patches that we remove from
        rem_mat <- matrix(0, ncol = (nStages*npatches), nrow = Nrem)    # mat that tracks each removal stage and patch for each individual removed
        
        group_rem <- vector(length = npatches)     # vec of how many removals from each group 
        stage_rem <- matrix(0, nrow = npatches, ncol = nStages)    # our removal matrix, tracking which stage is removed and from what group
        
        for (x in 1:nrow(rem_mat)){
          a <- sample(c(1:npatches), 1)        # patch sample
          s <- sample(c(1:4), 1, prob = avDist)        # stage sample
          rem_mat[x, (s + 4*(a-1))] <- 1    
        }
        
        thisRemvec <- colSums(rem_mat)
        
      }
      
      # calculating actual number removed from each stage for output
      
      
      group_rem <- vector()       # previous indexing issue - vec[i,] keeps calculating this years size
      for (p in 1:npatches){
        group_rem[p] <- sum(thisRemvec[(1+4*(p-1)):(4*p)])      # incorrect indexing here. We want to sum every 4 entries in the vector
      }# removed per group
      
      
      #  into outputs
      Nremoved[idx] <- sum(thisRemvec)  # total number removed   
      Remvec[idx,] <- thisRemvec
      Remgroup[idx,] <- group_rem
      
      # new population size following removals
      Vec[(i + 1), ] <- round(thisAmat %*% Vec[i, ]) - thisRemvec       # year after remyear filled with new stage vec. 
      # CHANGED FORM FLOOR TO ROUND
    }
    
    
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
      Vec[is.na(Vec)] <- 0
    }
    
    # group sizes next year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[i + 1, (1+4*(p-1)):(4*p)])
    }
    
    # filling in group size this year
    Group[i + 1,] <- sizes
    
    # pop sizes by year
    Pop[i + 1] <- sum(Vec[(i + 1), ])
    
    # if pop size <= 0, stop and return message with year
    if(Pop[i]<= 0 || is.na(Pop[i])) {
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
  }
  
  # output objects
  if (isTRUE(return.vec)) {
    out <- list(pop = vector(), 
                vec = matrix(), 
                group = matrix(),
                leavers = matrix()
    )
    
    out$pop <- Pop
    out$vec <- Vec 
    out$group <- Group
    out$leavers <- Leavers
    
    if(isTRUE(return.vec) && isTRUE(return.remvec)){   # remvec and 
      out$Nremoved <- Nremoved
      out$Group_remvec <- Remgroup
      out$Remvec <- Remvec
    }
    
    return(out)
  } else {
    return(Pop)
  }
}






# plot functions  - have not done vec 
meta_plot <- function(out,   # output obj of dd.proj
                               y_val= "N",   # plot type - N or Vec 
                               ylab = "abundance", 
                               xlab = "time (t)",
                               rem_year = NULL,
                               mytheme = theme_classic(), 
                               cols= "black",    # can be vector of 2 cols
                               legend.pos = "top",
                               base_size = 16) {
  # loading required libraries
  require(tidyr)
  require(ggplot2)
  # creating time vector for n years
  t <- nrow(out$vec) - 1                 # t= n years (0-t = t+1 entries)
  time <- as.numeric(c(0:t))       # vector 0:t
  
  # for pop size over time graph - must be as df
  if(y_val %in% c("N", "Pop Size", "pop size", "Pop", "pop")) {
    pop_df <- data.frame(Year = time,
                         Pop  = out$pop )# converting pop to a df
    
    # creating pop projection plot
    base.plot <- ggplot(pop_df, aes(x = Year, y = Pop)) +
      geom_line(size = 1, colour = "grey30") +
      geom_point(size = 1.5, colour = "black") +
      labs(
        x = xlab,
        y = ylab) +
      scale_y_continuous(
        name = "Population size",
        breaks = scales::pretty_breaks(n = 5),
        expand = expansion(mult = c(0.1, 0.2)))+  # 5% below, 30% above
      mytheme +
      theme(
        text = element_text(size = base_size),
        plot.title = element_text(size = base_size + 2, face = "bold"),
        axis.title = element_text(size = base_size),
        axis.text = element_text(size = base_size - 2),
        legend.position = legend.pos
      ) 
    
    if (!is.null(rem_year)){
      plot <- base.plot +
        geom_vline(xintercept = rem_year,                       # Adding a line to show removal year
                   colour = "red3", 
                   linetype = "dashed", linewidth = 1, 
                   alpha = 0.5) 
      
    } else if (is.null(rem_year)){
      plot <- base.plot 
    }
    
    
    # Group size by year plot  - update from vis file
   } else if(y_val %in% c("Group", "group")) {
     df <- as.data.frame(out$group)
     colnames(df) <- paste0("Group", seq_len(ncol(df)))
     df$Year <- time
     group_df <- pivot_longer(df, cols = starts_with("Group"),
                                         names_to = "Group", values_to = "Value") 
     
      # creating group projection plot
     base.plot <- ggplot(group_df, aes(Year, Value, colour = Group)) +
        geom_line(size = 1) +
        geom_point(size = 1.5) +
        labs(
          x = xlab,
          y = ylab) +
        scale_y_continuous(
          name = "Group size",
          breaks = scales::pretty_breaks(n = 5),
          expand = expansion(mult = c(0.1, 0.2)))+  # 5% below, 30% above
        mytheme +
        theme(
          text = element_text(size = base_size),
          plot.title = element_text(size = base_size + 2, face = "bold"),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 2),
          legend.position = legend.pos
        ) 
      if (!is.null(rem_year)){
        plot <- base.plot +
          geom_vline(xintercept = rem_year,                       # Adding a line to show removal year
                     colour = "red3", 
                     linetype = "dashed", linewidth = 1, 
                     alpha = 0.5) 
        
      } else if (is.null(rem_year)){
        plot <- base.plot 
      }
  }
  
  return(plot)
  
}   