# Meta removal trial
# 8/10/2026
# 

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


# updating my removal functions to remove per patch
# starting simple - how to remove a determined number at a determined year
# ISSUES : it looks like a very small proportion are removed in subsequent culls, instead of 36% closer to 4%. Needs fixing in both functions

proj_meta_nodisp_rem<- function(Umat,   
                                initial, 
                                stagedist,  # proportion of each stage
                                params, 
                                stagenames, 
                                npatches,
                                time,     
                                DDapply="fertility", 
                                remyear = integer(0),  # removal year = vector of years 
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
  # ensure remyear are integers and between 1 - time
  remyear <- unique(as.integer(remyear))
  remyear <- remyear[remyear >= 1 & remyear <= time] 
  if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
  
  nYears <- length(remyear)   # over how many years are removals taking place?  remove year_goal for each value of rem_index
  # Map removal year -> index into Remvec / Nremoved:
  rem_index <- if (length(remyear) > 0) setNames(seq_along(remyear), remyear)  # takes length rem year, names assigned as seq of integers 
  else integer(0)   # otherwise remyear = 0
  
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year
  
  if (length(remyear) > 0) {
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
    Vec[(i + 1), ] <- floor(thisAmat %*% Vec[i, ])  # amat values multiplied by vec, round down for integers  
                                                    # overwritten when this is a remyear
    
    # checking for removal years --------
    if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
      idx <- rem_index[as.character(i)]     # idx = 1: remyear. What removal year are we 
      
      # FIRST YEAR
      if(idx == 1){       # when idx = 1, remove 70% with goal or proportion to generate rem (number of each sex and stage removed)
        # setting removal goal
        real_int <- rnorm(1, mean = 0.7, sd = 0.05)   # variation in numbers removed # for repeatability, just use intensity?
        
        preN <- Pop[remyear[1]]  # population size pre cull # move until i
        Nrem <- round(preN * real_int)   # goal to remove  - pop size before first rem  
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
        # SUBSEQUENT YEARS - very low numbers removed
      }  else if(idx > 1 && idx <= 4){    # if supplementary, change strat to calculate max and min
        estN <- 0.8 * preN - sum(Nremoved)     # ISSUE : estimated pop size- cannot be tracked directly. This leads to negative estN
        # Use estimate from year before culls minus the culled N since. 2/5 multiplier from published methods
        estN[estN < 0] <- 0                                               
        # set min and max removals
        base <- 0.36 * estN                # base using estimated pop size
        min = base - (0.25*estN)          # difference between min and max = 25% of pre-cull total
        max = base + (0.25*estN)
        # need to prevent minimum falling bellow 0 - removing negative individuals
        min[min < 0] <- 0
        
        # set seed for repeatability and use a normal dist instead
        Nrem <- trunc.norm(1, min, max) # how many to remove this year using truncated norm dist between min and max values
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
        # SUPPLEMENTARY CULLS - ISSUE this part removes too many
      }  else if(idx > 4){    # if SBC, change strat to calculate max and min. Needs to be much lower intensity?
        base <- 0.2 * Nremoved[1]   # 36% of first year cull total - SBC. CHanging this to 10% because of low population sizes
        # set min and max removals
        min = base - (0.125*preN)          # difference between min and max = 25% of pre-cull total
        max = base + (0.125*preN)       # set maximum as as base to avoid large removals freqs?
        min[min < 0] <- 0
        
        # dist between min and max with mean = baseline?
        # set seed for repeatability and use a normal dist instead
        Nrem <- trunc.norm(1, min, max)  # how many to remove this year
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
      }
      
      # distributing this total across the number of groups - sett based removals or random
      group_remvec <- floor(rep(Nrem/npatches, npatches))     # for now, since groups approx equal, divide equally across groups   
                                                             # how to round so that vec sums to Nrem? 
      rem <- rep(NA, nStages*npatches)       # blank vector to fill as removal vector
      
      
      # pop removal -------------------------
      if(rem_strat == "random"){
        if (is.null(bias) == FALSE) paste("ignoring bias value as removal is random across ages and sexes")
        #generating the distribution - varies with rem strat
        dist <- stagedist
        
        for (p in 1:length(group_remvec)){    # getting distribution per patch using loop 
          rem[(1 +4*(p-1)) : (4*p)] <- group_remvec[p] * dist
        }    # where does variation come in?
        
        # stage biased
      } else if(rem_strat != "random" && is.numeric(bias) == FALSE){   # for biased rems that are NOT index specific...
        bias_vec <- rep(bias, 4)    # adults add bias, y remove bias
        
          if (rem_strat %in% c("adults", "Adults", "adult", "Adult")){
            # how to bias removals for classes? 
            bias_vec * c(-1, 1, -1, 1)   # bias is removed from yearlings, added to adults
            
          }
          else if (rem_strat %in% c("yearlings", "Yearlings", "yearling", "Yearlings")){
            bias_vec * c(1, -1, 1, -1)   # bias is added to yearlings, subtracted from adults
          }
          
          # sex biased
          if(rem_strat %in% c("females", "Females", "female", "Female")){
            bias_vec * c(1, 1, -1, -1)   # bias is added to fems, subtracted from males
            
          } else if(rem_strat %in% c("males", "Males", "male", "Male")){
            bias_vec * c(-1, -1, 1, 1)   # bias is subtracted from fems, added to males
          }
          dist <- stagedist + bias_vec    # combining into new removal distribution
          for (p in 1:length(group_remvec)){    # getting distribution per patch
            rem[(1 +4*(p-1)) : (4*p)] <- group_remvec[p] * dist
          }    # number to remove per stage
          
          if(rem_strat %in% c("Sett", "sett", "Group", "group")){
            
          }
        }
       
      # if numeric?
      
      # calculating number removed from each stage
      thisRemvec <- round(rem)   # actual removed = rounded abundance * dist
      
       grem <- vector()       # previous indexing issue - vec[i,] keeps calculating this years size
        for (p in 1:npatches){
          grem[p] <- sum(thisRemvec[(1+4*(p-1)):(4*p)])      # incorrect indexing here. We want to sum every 4 entries in the vector
        }# removed per group
      
      
      #  into outputs
      Nremoved[idx] <- sum(thisRemvec)  # total number removed   
      Remvec[idx,] <- thisRemvec
      Remgroup[idx,] <- grem
      
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


# testing patch based removals
n0_3patch <- round(c(rep(initials[1,]/4,3))) 

test_proj <- proj_meta_nodisp_rem(threepatch_Umat, n0_3patch, 
                                        stagedist,  # proportion of each stage
                                        params, 
                                        stagenames = stages, 
                                        npatches = 3,
                                        time = 50,     
                                        DDapply="fertility", 
                                        remyear = integer(0),  # removal year = vector of years 
                                        rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                                        bias = NULL , # strength of bias as percentage (range??)
                                        return.vec= TRUE, 
                                        return.remvec = FALSE)


test_proj_rem <- proj_meta_nodisp_rem(threepatch_Umat, n0_3patch, 
                                      stagedist,  # proportion of each stage
                                      params, 
                                      stagenames = stages, 
                                      npatches = 3,
                                      time = 50,     
                                      DDapply="fertility", 
                                      remyear = c(10:16),  # removal year = vector of years 
                                      rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                                      bias = NULL , # strength of bias as percentage (range??)
                                      return.vec= TRUE, 
                                      return.remvec = TRUE)
# no removals in subsequent cull years



# incorporating dispersal into function
proj_meta_DDdisp_rem <- function(Umat,
                                 initial,
                                 stagedist,  # proportion of each stage
                                 params,
                                 stagenames,
                                 npatches,
                                 time,
                                 DDapply="fertility",
                                 dispersal_stages = c("Adult_f", "Adult_m"),
                                 remyear = integer(0),  # removal year = vector of years 
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
  # ensure remyear are integers and between 1 - time
  remyear <- unique(as.integer(remyear))
  remyear <- remyear[remyear >= 1 & remyear <= time] 
  if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
  
  nYears <- length(remyear)   # over how many years are removals taking place?  remove year_goal for each value of rem_index
  # Map removal year -> index into Remvec / Nremoved:
  rem_index <- if (length(remyear) > 0) setNames(seq_along(remyear), remyear)  # takes length rem year, names assigned as seq of integers 
  else integer(0)   # otherwise remyear = 0
  
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year
  Leavers <- matrix(0, ncol = ncol(Vec), nrow = time)   # matrix to fill with the number of individuals leaving each patch each year (by each stage)
  
  if (length(remyear) > 0) {
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
    if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
      idx <- rem_index[as.character(i)]     # idx = 1: remyear. What removal year are we 
      
      # FIRST YEAR
      if(idx == 1){       # when idx = 1, remove 70% with goal or proportion to generate rem (number of each sex and stage removed)
        # setting removal goal
        real_int <- rnorm(1, mean = 0.7, sd = 0.05)   # variation in numbers removed # for repeatability, just use intensity?
        
        preN <- Pop[remyear[1]]  # population size pre cull # move until i
        Nrem <- round(preN * real_int)   # goal to remove  - pop size before first rem  
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
        # SUBSEQUENT YEARS
      }  else if(idx > 1 && idx <= 4){    # if supplementary, change strat to calculate max and min
        estN <- 0.8 * preN - sum(Nremoved)     # ISSUE : estimated pop size- cannot be tracked directly. This leads to negative estN
        # Use estimate from year before culls minus the culled N since. 2/5 multiplier from published methods
        estN[estN < 0] <- 0                                               
        # set min and max removals
        base <- 0.36 * estN                # base using estimated pop size
        min = base - (0.25*estN)          # difference between min and max = 25% of pre-cull total
        max = base 
        # need to prevent minimum falling bellow 0 - removing negative individuals
        min[min < 0] <- 0
        
        # set seed for repeatability and use a normal dist instead
        Nrem <- trunc.norm(1, min, max) # how many to remove this year using truncated norm dist between min and max values
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
        # SUPPLEMENTARY CULLS - ISSUE this part removes too many
      }  else if(idx > 4){    # if SBC, change strat to calculate max and min. Needs to be much lower intensity?
        base <- 0.2 * Nremoved[1]   # 36% of first year cull total - SBC. CHanging this to 10% because of low population sizes
        # set min and max removals
        min = base - (0.125*preN)          # difference between min and max = 25% of pre-cull total
        max = base + (0.125*preN)       # set maximum as as base to avoid large removals freqs?
        min[min < 0] <- 0
        
        # dist between min and max with mean = baseline?
        # set seed for repeatability and use a normal dist instead
        Nrem <- trunc.norm(1, min, max)  # how many to remove this year
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
      }
      
      # distributing this total across the number of groups - sett based removals or random
      group_remvec <- floor(rep(Nrem/npatches, npatches))     # for now, since groups approx equal, divide equally across groups   
      # how to round so that vec sums to Nrem? 
      rem <- rep(NA, nStages*npatches)       # blank vector to fill as removal vector
      
      
      # pop removal -------------------------
      if(rem_strat == "random"){
        if (is.null(bias) == FALSE) paste("ignoring bias value as removal is random across ages and sexes")
        #generating the distribution - varies with rem strat
        dist <- stagedist
        
        for (p in 1:length(group_remvec)){    # getting distribution per patch using loop 
          rem[(1 +4*(p-1)) : (4*p)] <- group_remvec[p] * dist
        }    # where does variation come in?
        
        # stage biased
      } else if(rem_strat != "random" && is.numeric(bias) == FALSE){   # for biased rems that are NOT index specific...
        bias_vec <- rep(bias, 4)    # adults add bias, y remove bias
        
        if (rem_strat %in% c("adults", "Adults", "adult", "Adult")){
          # how to bias removals for classes? 
          bias_vec * c(-1, 1, -1, 1)   # bias is removed from yearlings, added to adults
          
        }
        else if (rem_strat %in% c("yearlings", "Yearlings", "yearling", "Yearlings")){
          bias_vec * c(1, -1, 1, -1)   # bias is added to yearlings, subtracted from adults
        }
        
        # sex biased
        if(rem_strat %in% c("females", "Females", "female", "Female")){
          bias_vec * c(1, 1, -1, -1)   # bias is added to fems, subtracted from males
          
        } else if(rem_strat %in% c("males", "Males", "male", "Male")){
          bias_vec * c(-1, -1, 1, 1)   # bias is subtracted from fems, added to males
        }
        dist <- stagedist + bias_vec    # combining into new removal distribution
        for (p in 1:length(group_remvec)){    # getting distribution per patch
          rem[(1 +4*(p-1)) : (4*p)] <- group_remvec[p] * dist
        }    # number to remove per stage
        
        if(rem_strat %in% c("Sett", "sett", "Group", "group")){
          
        }
      }
      
      # if numeric?
      
      # calculating number removed from each stage
      thisRemvec <- round(rem)   # actual removed = rounded abundance * dist
      
      grem <- vector()       # previous indexing issue - vec[i,] keeps calculating this years size
      for (p in 1:npatches){
        grem[p] <- sum(thisRemvec[(1+4*(p-1)):(4*p)])      # incorrect indexing here. We want to sum every 4 entries in the vector
      }# removed per group
      
      
      #  into outputs
      Nremoved[idx] <- sum(thisRemvec)  # total number removed   
      Remvec[idx,] <- thisRemvec
      Remgroup[idx,] <- grem
      
      # new population size following removals
      # we don't use Amat * vec anymore, need to use postdispersal vec and large Amat?
      Vec[(i + 1), ] <- floor(thisAmat %*% postdispVec) - thisRemvec  # year after remyear filled with new stage vec. Why is this not filling correct year?
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


# testing DDdisp and removals
test_DDdisp_rem <- proj_meta_DDdisp_rem(threepatch_Umat, n0_3patch, 
                                        stagedist,  # proportion of each stage
                                        params, 
                                        stagenames = stages, 
                                        npatches = 3,
                                        time = 50,     
                                        DDapply="fertility", 
                                        dispersal_stages = c("Adult_f", "Adult_m"),
                                        remyear = c(10:14),  # removal year = vector of years 
                                        rem_strat = "random",  # if specified 
                                        bias = NULL, 
                                        return.vec=TRUE, 
                                        return.remvec = TRUE)

# Error in proj_meta_DDdisp_rem(threepatch_Umat, n0_3patch, stagedist, params,  : 
# object 'thisAmat' not found


