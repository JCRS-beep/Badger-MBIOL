# Meta removal trial
# 8/10/2026

# updating my removal functions to remove per patch
# starting simple - how to remove a determined number at a determined year
proj_meta_nodisp_rem<- function(Umat,     
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
  
  if (length(remyear) > 0) {
    # ensure remyear are integers and between 1 - time
    
    remyear <- unique(as.integer(remyear))
    remyear <- remyear[remyear >= 1 & remyear <= time]
    if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
  }
  
  # Check that the n_stages*n_patches = n_col of Umat:
  if (ncol(Umat)!=nStages*npatches){
    stop("Umat dimensions do not match the number of stages multiplied by the 
         number of patches.")
  }
  nStages<- length(stagenames)
  
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
  nYears <- length(remyear)   # over how many years are removals taking place?  remove year_goal for each value of rem_index
  # Map removal year -> index into Remvec / Nremoved:
  rem_index <- if (length(remyear) > 0) setNames(seq_along(remyear), remyear)  # takes length rem year, names assigned as seq of integers 
  else integer(0)   # otherwise remyear = 0
  
  
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year
  Nremoved <- numeric(length(remyear))     
  Remgroup <- matrix()    # how many removed per group per year?
  Remvec <- matrix(0, nrow = length(remyear), ncol =4)  
  
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
    
    # checking for removal years --------
    if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
      idx <- rem_index[as.character(i)]     # idx = 1: remyear. What removal year are we 
      
      # FIRST YEAR
      if(idx == 1){       # when idx = 1, remove 70% with goal or proportion to generate rem (number of each sex and stage removed)
        # setting removal goal
        real_int <- rnorm(1, mean = intensity/100, sd = 0.05)   # variation in numbers removed # for repeatability, just use intensity?
        
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
      group_rems <- vector(NA, npatches)     
      
      # pop removal -------------------------
      if(rem_strat == "random"){
        if (is.null(bias) == FALSE) paste("ignoring bias value as removal is random across ages and sexes")
        #generating the distribution - varies with rem strat
        dist <- stagedist
        rem <-  round(Nrem * dist)    # where does variation come in?
        
        # stage biased
      } else if(rem_strat != "random" && is.numeric(bias) == FALSE){   # for biased rems that are NOT index specific...
        bias_vec <- rep(bias, 4)    # adults add bias, y remove bias
        
        if (is.numeric(intensity) && rem_strat %in% c("adult", "Adult", "yearling", "Yearling")){   # want to specify age and sex prob  - adult male, yearling fem.. 
          if (rem_strat %in% c("adults", "Adults", "adult", "Adult")){
            # how to bias removals for classes? 
            bias_vec * c(-1, 1, -1, 1)   # bias is removed from yearlings, added to adults
            
          }
          else if (rem_strat %in% c("yearlings", "Yearlings", "yearling", "Yearlings")){
            bias_vec * c(1, -1, 1, -1)   # bias is added to yearlings, subtracted from adults
          }
          
          # sex biased
        } else if (is.numeric(intensity) && rem_strat %in% c("females", "Females", "female", "Female", "males", "Males", "male", "Male")){  
          if(rem_strat %in% c("females", "Females", "female", "Female")){
            bias_vec * c(1, 1, -1, -1)   # bias is added to fems, subtracted from males
            
          } else if(rem_strat %in% c("males", "Males", "male", "Male")){
            bias_vec * c(-1, -1, 1, 1)   # bias is subtracted from fems, added to males
          }
          dist <- stagedist + bias_vec    # combining into new removal distribution
          rem <-  Nrem * dist     # number to remove per stage
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



