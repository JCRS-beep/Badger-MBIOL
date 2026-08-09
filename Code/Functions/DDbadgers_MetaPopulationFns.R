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
  
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals              
  Pop[1] <- as.numeric(sum(n0))
  
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
proj_meta_fixeddisp_norem<- function(Umat,
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
  Leavers <- matrix(0, ncol = ncol(Vec), nrow = time) # matrix to fill with the number of individuals leaving each patch each year (by each stage)
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming cls matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals              
  Pop[1] <- as.numeric(sum(n0))
  
#  colnames(Leavers) <- seq(1, npatches)   # patch names 
 # rownames(Leavers) <- 0:(time) 
  
  # Loop = density dependent matrix application for each year and patch
  for (i in 1:time) {   #  projection until end time
    
    Leavers <- rep(0, ncol(Vec))    # individuals leaving patch p - to go in output
    arrivers <- rep(0, ncol(Vec))   # individuals arriving in patch p+1 (other patch)
    
    # dispersal loop (NOTE that for this simple version with a fixed amount of
    # dispersal, we don't need to do a loop, but I'm pointing you in the
    # direction of the next steps)
    for (p in 1:npatches){
      Ipatch <- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      subVec <- as.vector(Vec[p, Ipatch])    # have to specify row (year) we want, to remove names format as numeric vector
      
      # apply the probability of leaving to the stages that leave:
      subleavers<- rep(0,4)
      subleavers[I_leavers]<- round(dispersal_prob*subVec[I_leavers])   # round instead of floor - does not underestimate moving individuals
      Leavers[Ipatch] <- subleavers

      # distribute the leavers equally among the other patches, as arrivers, in
      # the same stage as they left:
      if (npatches==2){
        subarrivers<- subleavers      
        arrivers[-Ipatch]<- subarrivers
        
      } else {     # divide movers across other patches
        # random sample of patches to decide where individuals go (later replace with some distance equation?)
        nleavers <- sum(subleavers) # number of moving individuals
        x <- seq(npatches)
        arrival_patches <- x[x != p]   # all patches excluding current patch
        
        remaining_subleavers <- subleavers   # how many leavers from this patch remain
        mat_subarrival <- matrix(0, nrow = length(arrival_patches), ncol = nStages)
        rownames(mat_subarrival) <- arrival_patches
        
        # option - loop for each stage in disp stage - incorrectly assigning - overly complicated, simplify ---
        for(j in 1:length(I_leavers)){    # why is j giong to 4?
          stageI <- I_leavers[j]
          stage_I_arrival <- sample(arrival_patches, subleavers[stageI], replace = TRUE)   # samples patches for the stage j in dispersal stages
          
          for (a in arrival_patches){  
            thisIpatch <- c(1+4*(a-1)):(4*a)          
            matrow <- which(rownames(mat_subarrival) == a)
            
          if(matrow == 0) {
            arrivers[thisIpatch] <- rep(0, 4)
          }
         
          else {mat_subarrival[matrow,stageI] <- sum(stage_I_arrival == a)   # arrivers for this patch
        
          arrivers[thisIpatch] <- c(mat_subarrival[matrow,])
          }
        }
        }
        
    }
    
    postdispVec<- Vec[i,] - Leavers + arrivers    # works in long vec - dont need to loop to calculate for each patch
    
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
      Vec[(i+1), Ipatch]<- floor(subAmat %*% Vec[i, Ipatch]) # amat values multiplied by vec, round down for integers 
      
    }
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
      Vec[is.na(Vec)] <- 0
    }
    
    Pop[i + 1] <- sum(Vec[(i + 1), ])
    
    # if pop size <= 0, stop and return message with year
    if(Pop[i]<= 0 || is.na(Pop[i])) {
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
   }
  }
  # output objects
  if (isTRUE(return.vec)) {
    out <- list(pop = vector(), 
                vec = matrix(), 
                leavers = matrix())
    out$pop <- Pop
    out$vec <- Vec 
    out$leavers <- Leavers
    return(out)
  } else {
    return(Pop)
  }
  }





    
   

# Linking density with probability of leaving
# Using absolute number of badgers to calculate density compared to max group size 
# social groups range 1 - 27, not observed larger than this.
# Dispersal prob calculated with quadratic between 0 - 1, equal for males and females (dispersal per patch)
ddDmat <- function(stagenames,dispersal_stages, # names of all stages, names of only dispersing stages
                 npatches, # number of patches
                 group_size, 
                 max_group) {    # what is the carrying capacity (can be vector of K per patch)
  
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
  
  max_group_size <- 28 # based on papers recording 26 and 27
  x <- group_size/ max_group    # proportion each group of max size
  
  for (p in 1:npatches){
    move <- round(x^2, 2)       # quadratic equation for now - rounding to 2 dec places
    move[move<0] <- 0   # setting any negatives to 0
    Dmat[p,I_leavers] <- move[p]   # equal for all indivs or variable between sexes? vary K value for males and females?
  }
  
  return(Dmat)
} 


proj_meta_DDdisp_norem <- function(Umat,
                                   initial,
                                   params,
                                   stagenames,
                                   npatches,
                                   time,
                                   DDapply="fertility",
                                   max_group, 
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
  Leavers <- matrix(0, ncol = ncol(Vec), nrow = time) # matrix to fill with the number of individuals leaving each patch each year (by each stage)
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming cls matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals              
  Pop[1] <- as.numeric(sum(n0))
  
  
  # Loop = density dependent matrix application for each year and patch
  for (i in 1:time) {   #  projection until end time
    
    Leavers <- rep(0, ncol(Vec))    # individuals leaving patch p - to go in output
    arrivers <- rep(0, ncol(Vec))   # individuals arriving in patch p+1 (other patch)
    
    # group sizes this year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
    sizes[p] <- sum(Vec[i, (1+4*(p-1)):(4*p)])
    }
    
    dmat <- ddDmat(stagenames,dispersal_stages, # names of all stages, names of only dispersing stages
                   npatches, # number of patches
                   group_size = sizes, 
                   max_group)
    
    # dispersal loop (NOTE that for this simple version with a fixed amount of
    # dispersal, we don't need to do a loop, but I'm pointing you in the
    # direction of the next steps)
    for (p in 1:npatches){
      Ipatch <- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      subVec <- as.vector(Vec[p, Ipatch])    # have to specify row (year) we want, to remove names format as numeric vector
      
      # apply the probability of leaving to the stages that leave:
      subleavers<- rep(0,4)
      subleavers[I_leavers]<- round(dmat[p,I_leavers]*subVec[I_leavers])   # round instead of floor - does not underestimate moving individuals
      Leavers[Ipatch] <- subleavers
      
      # distribute the leavers equally among the other patches, as arrivers, in
      # the same stage as they left:
      if (npatches==2){
        subarrivers<- subleavers      
        arrivers[-Ipatch]<- subarrivers
        
      } else {     # divide movers across other patches
        # random sample of patches to decide where individuals go (later replace with some distance equation?)
        nleavers <- sum(subleavers) # number of moving individuals
        x <- seq(npatches)
        arrival_patches <- x[x != p]   # all patches excluding current patch
        
        remaining_subleavers <- subleavers   # how many leavers from this patch remain
        mat_subarrival <- matrix(0, nrow = length(arrival_patches), ncol = nStages)
        rownames(mat_subarrival) <- arrival_patches
        
        # option - loop for each stage in disp stage - incorrectly assigning - overly complicated, simplify ---
        for(j in 1:length(I_leavers)){    # why is j giong to 4?
          stageI <- I_leavers[j]
          stage_I_arrival <- sample(arrival_patches, subleavers[stageI], replace = TRUE)   # samples patches for the stage j in dispersal stages
          
          for (a in arrival_patches){  
            thisIpatch <- c(1+4*(a-1)):(4*a)          
            matrow <- which(rownames(mat_subarrival) == a)
            
            if(matrow == 0) {
              arrivers[thisIpatch] <- rep(0, 4)
            }
            
            else {mat_subarrival[matrow,stageI] <- sum(stage_I_arrival == a)   # arrivers for this patch
            
            arrivers[thisIpatch] <- c(mat_subarrival[matrow,])
            }
          }
        }
        
      }
      
      postdispVec<- Vec[i,] - Leavers + arrivers    # works in long vec - dont need to loop to calculate for each patch
      
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
        Vec[(i+1), Ipatch]<- floor(subAmat %*% Vec[i, Ipatch]) # amat values multiplied by vec, round down for integers 
        
      }
      # if any stage becomes negative, set to zero and continue
      if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
        warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
        Vec[Vec < 0] <- 0
        Vec[is.na(Vec)] <- 0
      }
      
      Pop[i + 1] <- sum(Vec[(i + 1), ])
      
      # if pop size <= 0, stop and return message with year
      if(Pop[i]<= 0 || is.na(Pop[i])) {
        warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
        break
      }
    }
  }
  # output objects
  if (isTRUE(return.vec)) {
    out <- list(pop = vector(), 
                vec = matrix(), 
                leavers = matrix())
    out$pop <- Pop
    out$vec <- Vec 
    out$leavers <- Leavers
    return(out)
  } else {
    return(Pop)
  }
}