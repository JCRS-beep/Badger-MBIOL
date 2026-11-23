## Badger metapopulation function, rewritten as individual-based
# Written by Chrissy Hernandez 
# 11 November 2026

# DD arrival (Here) probability matrix ----
# individuals 'attracted' to low density patches
# P arrival decreases as group size increases
# Same layout as ddDmat: rows = patches, cols = stages.
# Hmat[p, s] = probability that a dispersing individual of stage s chooses
# patch p as its destination. Low-density patches are more attractive.
# Each dispersing-stage column sums to 1 across patches.
ddHmat <- function(stagenames, dispersal_stages,
                   npatches,
                   group_size,                 # current size of every group
                   max_group,                  # carrying capacity (single value or vector)
                   steepness = 0.3) {          # higher = sharper preference for low density
  
  if (!all(dispersal_stages %in% stagenames)) {
    stop("dispersal_stages must all be present in stagenames.")
  }
  if (length(group_size) != npatches) {
    stop("You must provide the size of each group.")
  }
  
  I_arrivers <- which(stagenames %in% dispersal_stages)
  xhalf = 0.3 * max_group    # group size where attractiveness = 0.5, varies depending on what we set max group
  
  # attractiveness weight per patch: ~1 when empty, 0.5 at xhalf, -> 0 when large
  weight <- 1 / (1 + exp(steepness * (group_size - xhalf)))  # this is not probability - needs to sum to 1
  weight[group_size > max_group] <- 0            # full patches are not attractive
  
  # convert weights to probabilities that sum to 1 across patches
  if (sum(weight) < 1e-12) {
    p_arrive <- rep(1 / npatches, npatches)  # all saturated: equal chance
  } else {
    p_arrive <- round(weight / sum(weight), 3)
  }
  
  Hmat <- matrix(0, nrow = npatches, ncol = length(stagenames))
  colnames(Hmat) <- stagenames
  rownames(Hmat) <- paste0("patch_", seq_len(npatches))
  
  for (s in I_arrivers) {
    Hmat[, s] <- p_arrive    # equal across stages for now; could vary by sex later
  }
  
  return(Hmat)
}

# DD moving projection function as IBM -----
# incorporating density dependent movement probs in projection function
proj_meta_DDdisp_norem_IBM <- function(Umat,
                                       initial,
                                       params,
                                       stagenames,
                                       npatches,
                                       time,
                                       max_group = 28,
                                       DDapply="fertility",
                                       dispersal_stages = c("Adult_f", "Adult_m"),
                                       return.vec=TRUE){
  # input checks
  if (time <= 1) {stop("Time must be a positive integer")
    } else {time <- as.integer(time)}
  
  if (is.null(initial)) {
    stop("You must provide initial population vector with each stage abundance")
  } else {n0 <- as.numeric(initial)}
  
  if(is.null(stagenames) || length(stagenames) == 0) {
    stop("stagenames must be provided for correct matrix dimensions")
  }
  
  if (is.null(params)) {
    stop("Please provide parameters for selected density-dependent function")
  }
  
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
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming columns of matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals  
  Pop[1] <- as.numeric(sum(n0))
  
  for (p in 1:npatches){
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
    
    ## dispersal process:
    prob_leaving <- ddDmat(stagenames, dispersal_stages, # names of all stages, names of only dispersing stages
                           npatches, # number of patches
                           group_size = sizes)
    # reshape prob_leaving into a vector with stages within patches:
    prob_leaving<- as.vector(t(prob_leaving))
    
    # individuals leaving patches
    leavers <- rbinom(n = length(Vec[i,]), 
                      size = Vec[i,],
                      prob = prob_leaving)
    
    #### Arrival process: Needs updating to density dependence! for now, let's
    #say probability of arrival is equal for all patches other than the
    #departure patch
    prob_arrival <- ddHmat(stagenames, dispersal_stages, # names of all stages, names of only dispersing stages
                           npatches, # number of patches
                           group_size = sizes,
                           max_group,                  # carrying capacity (single value or vector)
                           steepness = 0.3)
    
    # individuals arriving to patches
    arrivers <- matrix(0, nrow=npatches, ncol=length(stagenames))   # individuals arriving in other patches
    for (p in 1:npatches){
      Ipatch <- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      subleavers<- leavers[Ipatch]
      if (sum(subleavers)<1){ next }   # only continues if there are individuals leaving current patch
      # possible destinations
      arrival_patches <- setdiff(1:npatches, p)   # setdiff compares differences between dataframes, removes p from all other patches
      # destinations for the number of individuals (currently, we sum all the
      # individuals leaving a given location. Would need to update if
      # probability of arrival was stage-specific)
      subarrivers<- matrix(0, nrow=npatches, ncol=length(stagenames))
      for (s in I_leavers){
        destinations<- sample(x = arrival_patches, 
                              size = subleavers[s],
                              replace = TRUE,
                              prob = prob_arrival[arrival_patches])
        subarrivers[,s]<- tabulate(destinations, nbins = npatches)
      }
      arrivers<- arrivers + subarrivers
    }
    # change arrivers into a row vector with stages within patches:
    arrivers<- as.vector(t(arrivers))
    
    # Check that the number of leavers is equal to the number of arrivers.
    if (sum(leavers) != sum(arrivers)){
      stop('Code bug: the number of leavers does not equal the number of arrivers during dispersal step.')
    }
    
    postdispVec<- Vec[i,] - leavers + arrivers    
    Leavers[i,] <- leavers # assigning leavers vec to relevant row in Leavers matrix
    
    # Next, individuals reproduce in their new patches:
    births<- rep(0, length(postdispVec))
    for (p in 1:npatches){
      Ipatch<- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
     
      # Nf calculation
      subNf <- sum(postdispVec[(4*(p-1)+I_Af)], postdispVec[(4*(p-1)+I_Yf)])    # Nf sums yearling and adult fems in Vec matrix
      subNm <- sum(postdispVec[(4*(p-1)+I_Am)], postdispVec[(4*(p-1)+I_Ym)])   # Nm 
      
      # ricker density dependence
      subN <- sum(postdispVec[Ipatch])  # pop sizes sums row i for cols corresponding to patch p
      repro <- badger_repro_IBM(params, subN, subNf, subNm)
      # Each female has a Poisson-distributed clutch size with mean size equal to repro:
      subbirths<- rpois(postdispVec[Ipatch[I_Af]], repro)
      # Those offspring are female with binomial probability 0.5:
      female_births<- rbinom(n=1, size=sum(subbirths), prob=0.5)
      male_births<- sum(subbirths) - female_births
      
      births[Ipatch[I_Yf]]<- female_births
      births[Ipatch[I_Ym]]<- male_births
      }
    
    # Next, we apply survival and growth to the postdispVec:
    prob_surv<- colSums(Umat) # survival probabilities for all stages/patches
    ### Note that the input params could be updated to be a survival probability vector instead of a matrix
    postsurvVec<- rbinom(n = length(postdispVec),
                         size = postdispVec,
                         prob = prob_surv)
    # when yearlings survive, they become adults. Adults that survive stay in their class.
    yearlingfems<- seq(from=I_Yf, by=length(stagenames), to=length(Vec[i,]))
    adultfems<- seq(from=I_Af, by=length(stagenames), to=length(Vec[i,]))
    postsurvVec[adultfems]<- postsurvVec[adultfems]+postsurvVec[yearlingfems]
    postsurvVec[yearlingfems]<- 0
    
    yearlingmales<- seq(from=I_Ym, by=length(stagenames), to=length(Vec[i,]))
    adultmales<- seq(from=I_Am, by=length(stagenames), to=length(Vec[i,]))
    postsurvVec[adultmales]<- postsurvVec[adultmales]+postsurvVec[yearlingmales]
    postsurvVec[yearlingmales]<- 0
    
    # Finally, we add the postsurvival vector and the births to get next year's counts:
    Vec[(i+1),]<- postsurvVec + births
    
    # group sizes next year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[(i + 1), (1+4*(p-1)):(4*p)])
    }
    
    # filling in group size this year
    Group[(i + 1),] <- sizes
    
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


# projection functions converted to IBM
# Converted by Claude
# updated claude conversion, preventing creation of smaller functions
# fns coonverted - nodisp no removal, fixeddisp no rem, dddisp rem
proj_meta_nodisp_norem_IBM <- function(Umat,     # must be a big matrix with diagonal Umats
                                       initial,
                                       params,
                                       stagenames,
                                       npatches,
                                       time,
                                       DDapply="fertility",
                                       return.vec=TRUE){
  # input checks
  if (time <= 1) {stop("Time must be a positive integer")
  } else {time <- as.integer(time)}
  
  if (is.null(initial)) {
    stop("You must provide initial population vector with each stage abundance")
  } else {n0 <- as.numeric(initial)}
  
  if(is.null(stagenames) || length(stagenames) == 0) {
    stop("stagenames must be provided for correct matrix dimensions")
  }
  
  if (is.null(params)) {
    stop("Please provide parameters for selected density-dependent function")
  }
  
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
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming columns of matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals  
  Pop[1] <- as.numeric(sum(n0))
  
  for (p in 1:npatches){
    x <- (1 +4*(p-1)) : (4*p)   # boundaries containing each patch abundances
    Group[1,p] <- sum(Vec[1, x])     # filling first row of group matrix
  }
  
  # Loop = individual-based projection for each year and patch
  for (i in 1:time) {   #  projection until end time
    
    # no dispersal: individuals stay in their starting patch
    postdispVec<- Vec[i,]
    
    # Individuals reproduce in their patches:
    births<- rep(0, length(postdispVec))
    for (p in 1:npatches){
      Ipatch<- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      
      # Nf calculation
      subNf <- sum(postdispVec[(4*(p-1)+I_Af)], postdispVec[(4*(p-1)+I_Yf)])    # Nf sums yearling and adult fems in Vec matrix
      subNm <- sum(postdispVec[(4*(p-1)+I_Am)], postdispVec[(4*(p-1)+I_Ym)])   # Nm 
      
      # ricker density dependence
      subN <- sum(postdispVec[Ipatch])  # pop sizes sums row i for cols corresponding to patch p
      repro <- badger_repro_IBM(params, subN, subNf, subNm)
      # Each female has a Poisson-distributed clutch size with mean size equal to repro:
      subbirths<- rpois(postdispVec[Ipatch[I_Af]], repro)
      # Those offspring are female with binomial probability 0.5:
      female_births<- rbinom(n=1, size=sum(subbirths), prob=0.5)
      male_births<- sum(subbirths) - female_births
      
      births[Ipatch[I_Yf]]<- female_births
      births[Ipatch[I_Ym]]<- male_births
    }
    
    # Next, we apply survival and growth to the postdispVec:
    prob_surv<- colSums(Umat) # survival probabilities for all stages/patches
    postsurvVec<- rbinom(n = length(postdispVec),
                         size = postdispVec,
                         prob = prob_surv)
    # when yearlings survive, they become adults. Adults that survive stay in their class.
    yearlingfems<- seq(from=I_Yf, by=length(stagenames), to=length(Vec[i,]))
    adultfems<- seq(from=I_Af, by=length(stagenames), to=length(Vec[i,]))
    postsurvVec[adultfems]<- postsurvVec[adultfems]+postsurvVec[yearlingfems]
    postsurvVec[yearlingfems]<- 0
    
    yearlingmales<- seq(from=I_Ym, by=length(stagenames), to=length(Vec[i,]))
    adultmales<- seq(from=I_Am, by=length(stagenames), to=length(Vec[i,]))
    postsurvVec[adultmales]<- postsurvVec[adultmales]+postsurvVec[yearlingmales]
    postsurvVec[yearlingmales]<- 0
    
    # Finally, we add the postsurvival vector and the births to get next year's counts:
    Vec[(i+1),]<- postsurvVec + births
    
    # group sizes next year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[(i + 1), (1+4*(p-1)):(4*p)])
    }
    
    # filling in group size this year
    Group[(i + 1),] <- sizes
    
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
                group = matrix()
    )
    
    out$pop <- Pop
    out$vec <- Vec 
    out$group <- Group
    
    return(out)
  } else {
    return(Pop)
  }
}



proj_meta_fixeddisp_norem_IBM <- function(Umat,
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
  if (time <= 1) {stop("Time must be a positive integer")
  } else {time <- as.integer(time)}
  
  if (is.null(initial)) {
    stop("You must provide initial population vector with each stage abundance")
  } else {n0 <- as.numeric(initial)}
  
  if(is.null(stagenames) || length(stagenames) == 0) {
    stop("stagenames must be provided for correct matrix dimensions")
  }
  
  if (is.null(params)) {
    stop("Please provide parameters for selected density-dependent function")
  }
  
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
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming columns of matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals  
  Pop[1] <- as.numeric(sum(n0))
  
  for (p in 1:npatches){
    x <- (1 +4*(p-1)) : (4*p)   # boundaries containing each patch abundances
    Group[1,p] <- sum(Vec[1, x])     # filling first row of group matrix
  }
  
  # fixed probability of leaving, only for the dispersing stages
  prob_leaving<- rep(0, nStages)
  prob_leaving[I_leavers]<- dispersal_prob
  # reshape prob_leaving into a vector with stages within patches:
  prob_leaving<- rep(prob_leaving, npatches)
  
  # Loop = individual-based projection for each year and patch
  for (i in 1:time) {   #  projection until end time
    
    ## dispersal process:
    # individuals leaving patches
    leavers <- rbinom(n = length(Vec[i,]), 
                      size = Vec[i,],
                      prob = prob_leaving)
    
    #### Arrival process: probability of arrival is equal for all patches other
    # than the departure patch
    # individuals arriving to patches
    arrivers <- matrix(0, nrow=npatches, ncol=length(stagenames))   # individuals arriving in other patches
    for (p in 1:npatches){
      Ipatch <- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      subleavers<- leavers[Ipatch]
      if (sum(subleavers)<1){ next }
      # possible destinations
      arrival_patches <- setdiff(1:npatches, p)   # setdiff compares differences between dataframes, removes p from all other patches
      subarrivers<- matrix(0, nrow=npatches, ncol=length(stagenames))
      for (s in I_leavers){
        # sample.int indexes into arrival_patches: sample() on a single patch number
        # (npatches = 2) would sample from 1:x and send leavers to the wrong patch
        destinations<- arrival_patches[sample.int(n = length(arrival_patches), 
                                                  size = subleavers[s],
                                                  replace = TRUE)]
        subarrivers[,s]<- tabulate(destinations, nbins=npatches)
      }
      arrivers<- arrivers + subarrivers
    }
    # change arrivers into a row vector with stages within patches:
    arrivers<- as.vector(t(arrivers))
    
    # Check that the number of leavers is equal to the number of arrivers.
    if (sum(leavers) != sum(arrivers)){
      stop('Code bug: the number of leavers does not equal the number of arrivers during dispersal step.')
    }
    
    postdispVec<- Vec[i,] - leavers + arrivers    
    Leavers[i,] <- leavers # assigning leavers vec to relevant row in Leavers matrix
    
    # Next, individuals reproduce in their new patches:
    births<- rep(0, length(postdispVec))
    for (p in 1:npatches){
      Ipatch<- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      
      # Nf calculation
      subNf <- sum(postdispVec[(4*(p-1)+I_Af)], postdispVec[(4*(p-1)+I_Yf)])    # Nf sums yearling and adult fems in Vec matrix
      subNm <- sum(postdispVec[(4*(p-1)+I_Am)], postdispVec[(4*(p-1)+I_Ym)])   # Nm 
      
      # ricker density dependence
      subN <- sum(postdispVec[Ipatch])  # pop sizes sums row i for cols corresponding to patch p
      repro <- badger_repro_IBM(params, subN, subNf, subNm)
      # Each female has a Poisson-distributed clutch size with mean size equal to repro:
      subbirths<- rpois(postdispVec[Ipatch[I_Af]], repro)
      # Those offspring are female with binomial probability 0.5:
      female_births<- rbinom(n=1, size=sum(subbirths), prob=0.5)
      male_births<- sum(subbirths) - female_births
      
      births[Ipatch[I_Yf]]<- female_births
      births[Ipatch[I_Ym]]<- male_births
    }
    
    # Next, we apply survival and growth to the postdispVec:
    prob_surv<- colSums(Umat) # survival probabilities for all stages/patches
    postsurvVec<- rbinom(n = length(postdispVec),
                         size = postdispVec,
                         prob = prob_surv)
    # when yearlings survive, they become adults. Adults that survive stay in their class.
    yearlingfems<- seq(from=I_Yf, by=length(stagenames), to=length(Vec[i,]))
    adultfems<- seq(from=I_Af, by=length(stagenames), to=length(Vec[i,]))
    postsurvVec[adultfems]<- postsurvVec[adultfems]+postsurvVec[yearlingfems]
    postsurvVec[yearlingfems]<- 0
    
    yearlingmales<- seq(from=I_Ym, by=length(stagenames), to=length(Vec[i,]))
    adultmales<- seq(from=I_Am, by=length(stagenames), to=length(Vec[i,]))
    postsurvVec[adultmales]<- postsurvVec[adultmales]+postsurvVec[yearlingmales]
    postsurvVec[yearlingmales]<- 0
    
    # Finally, we add the postsurvival vector and the births to get next year's counts:
    Vec[(i+1),]<- postsurvVec + births
    
    # group sizes next year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[(i + 1), (1+4*(p-1)):(4*p)])
    }
    
    # filling in group size this year
    Group[(i + 1),] <- sizes
    
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



proj_meta_DDdisp_rem_IBM <- function(Umat,
                                     initial,
                                     stagedist,  # proportion of each stage (no longer used - individuals are sampled)
                                     params,
                                     stagenames,
                                     npatches,
                                     time,
                                     DDapply="fertility",
                                     dispersal_stages = c("Adult_f", "Adult_m"),
                                     remyear = NULL,  # removal year = vector of years 
                                     rem_strat = "random",  # not used
                                     bias = NULL,           # not used
                                     return.vec=TRUE, 
                                     return.remvec = FALSE){
  # input checks
  if (time <= 1) {stop("Time must be a positive integer")
  } else {time <- as.integer(time)}
  
  if (is.null(initial)) {
    stop("You must provide initial population vector with each stage abundance")
  } else {n0 <- as.numeric(initial)}
  
  if(is.null(stagenames) || length(stagenames) == 0) {
    stop("stagenames must be provided for correct matrix dimensions")
  }
  
  if (is.null(params)) {
    stop("Please provide parameters for selected density-dependent function")
  }
  
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
    remyear <- sort(unique(as.integer(remyear)))
    remyear <- remyear[remyear >= 1 & remyear <= time] 
    if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
    
    nYears <- length(remyear)   # over how many years are removals taking place?
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
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming columns of matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals  
  Pop[1] <- as.numeric(sum(n0))
  
  for (p in 1:npatches){
    x <- (1 +4*(p-1)) : (4*p)   # boundaries containing each patch abundances
    Group[1,p] <- sum(Vec[1, x])     # filling first row of group matrix
  }
  
  # Loop = individual-based projection for each year and patch
  for (i in 1:time) {   #  projection until end time
    # group sizes this year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[i , (1+4*(p-1)):(4*p)])
    }
    
    ## dispersal process:
    prob_leaving <- ddDmat(stagenames, dispersal_stages, # names of all stages, names of only dispersing stages
                           npatches, # number of patches
                           group_size = sizes)
    # reshape prob_leaving into a vector with stages within patches:
    prob_leaving<- as.vector(t(prob_leaving))
    
    # individuals leaving patches
    leavers <- rbinom(n = length(Vec[i,]), 
                      size = Vec[i,],
                      prob = prob_leaving)
    
    #### Arrival process: Needs updating to density dependence! for now, let's
    #say probability of arrival is equal for all patches other than the
    #departure patch
    # individuals arriving to patches
    arrivers <- matrix(0, nrow=npatches, ncol=length(stagenames))   # individuals arriving in other patches
    for (p in 1:npatches){
      Ipatch <- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      subleavers<- leavers[Ipatch]
      if (sum(subleavers)<1){ next }
      # possible destinations
      arrival_patches <- setdiff(1:npatches, p)   # setdiff compares differences between dataframes, removes p from all other patches
      subarrivers<- matrix(0, nrow=npatches, ncol=length(stagenames))
      for (s in I_leavers){
        # sample.int indexes into arrival_patches: sample() on a single patch number
        # (npatches = 2) would sample from 1:x and send leavers to the wrong patch
        destinations<- arrival_patches[sample.int(n = length(arrival_patches), 
                                                  size = subleavers[s],
                                                  replace = TRUE)]
        subarrivers[,s]<- tabulate(destinations, nbins=npatches)
      }
      arrivers<- arrivers + subarrivers
    }
    # change arrivers into a row vector with stages within patches:
    arrivers<- as.vector(t(arrivers))
    
    # Check that the number of leavers is equal to the number of arrivers.
    if (sum(leavers) != sum(arrivers)){
      stop('Code bug: the number of leavers does not equal the number of arrivers during dispersal step.')
    }
    
    postdispVec<- Vec[i,] - leavers + arrivers    
    Leavers[i,] <- leavers # assigning leavers vec to relevant row in Leavers matrix
    
    # Next, individuals reproduce in their new patches:
    births<- rep(0, length(postdispVec))
    for (p in 1:npatches){
      Ipatch<- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      
      # Nf calculation
      subNf <- sum(postdispVec[(4*(p-1)+I_Af)], postdispVec[(4*(p-1)+I_Yf)])    # Nf sums yearling and adult fems in Vec matrix
      subNm <- sum(postdispVec[(4*(p-1)+I_Am)], postdispVec[(4*(p-1)+I_Ym)])   # Nm 
      
      # ricker density dependence
      subN <- sum(postdispVec[Ipatch])  # pop sizes sums row i for cols corresponding to patch p
      repro <- badger_repro_IBM(params, subN, subNf, subNm)
      # Each female has a Poisson-distributed clutch size with mean size equal to repro:
      subbirths<- rpois(postdispVec[Ipatch[I_Af]], repro)
      # Those offspring are female with binomial probability 0.5:
      female_births<- rbinom(n=1, size=sum(subbirths), prob=0.5)
      male_births<- sum(subbirths) - female_births
      
      births[Ipatch[I_Yf]]<- female_births
      births[Ipatch[I_Ym]]<- male_births
    }
    
    # Next, we apply survival and growth to the postdispVec:
    prob_surv<- colSums(Umat) # survival probabilities for all stages/patches
    postsurvVec<- rbinom(n = length(postdispVec),
                         size = postdispVec,
                         prob = prob_surv)
    # when yearlings survive, they become adults. Adults that survive stay in their class.
    yearlingfems<- seq(from=I_Yf, by=length(stagenames), to=length(Vec[i,]))
    adultfems<- seq(from=I_Af, by=length(stagenames), to=length(Vec[i,]))
    postsurvVec[adultfems]<- postsurvVec[adultfems]+postsurvVec[yearlingfems]
    postsurvVec[yearlingfems]<- 0
    
    yearlingmales<- seq(from=I_Ym, by=length(stagenames), to=length(Vec[i,]))
    adultmales<- seq(from=I_Am, by=length(stagenames), to=length(Vec[i,]))
    postsurvVec[adultmales]<- postsurvVec[adultmales]+postsurvVec[yearlingmales]
    postsurvVec[yearlingmales]<- 0
    
    # Finally, we add the postsurvival vector and the births to get next year's counts:
    Vec[(i+1),]<- postsurvVec + births
    
    # checking for removal years --------
    # removals are taken from the projected population (row i+1), the same point
    # as in proj_meta_DDdisp_rem where the removal vector was subtracted from year i+1
    if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
      idx <- rem_index[as.character(i)]     # idx = 1: remyear. What removal year are we 
      
      projN <- sum(Vec[(i + 1), ])   # number of individuals available to remove
      
      # FIRST YEAR
      if(idx == 1){       # when idx = 1, remove 70% with goal 
        preN <- Pop[remyear[1]]  # population size pre cull
        Nrem <- round(preN * 0.7)   # goal to remove  - 70% of pop size before first rem  
        Nrem <- min(Nrem, projN)    # cannot remove more individuals than are present
        
        # equal number removed from each group
        group_remvec <- rep(Nrem %/% npatches, npatches)
        extra <- Nrem %% npatches     # remainder goes to randomly chosen groups
        if (extra > 0){
          extra_patches <- sample.int(npatches, extra)
          group_remvec[extra_patches] <- group_remvec[extra_patches] + 1
        }
        
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
        
        # SUPPLEMENTARY CULLS
        else if(idx > 5){    # if SBC, change strat to calculate max and min. Needs to be much lower intensity?
          base <- 0.1 * Nremoved[1]   # 36% of first year cull total - SBC. CHanging this to 10% because of low population sizes
          # set min and max removals
          min = base - (0.125*preN)          # difference between min and max = 25% of pre-cull total
          max = base                         # set maximum as as base to avoid large removals freqs?
          
        }
        # need to prevent minimum falling bellow 0 - removing negative individuals
        min[min < 0] <- 0
        
        # set seed for repeatability and use a normal dist instead
        Nrem <- round(trunc.norm(1, min, max)) # how many to remove this year using truncated norm dist between min and max values
        Nrem <- min(Nrem, projN)    # cannot remove more individuals than are present
        
        # each removed individual is assigned to a random group
        group_remvec <- tabulate(sample.int(npatches, size = Nrem, replace = TRUE), nbins = npatches)
      }
      
      # removing individuals: drawn without replacement from each group, so the
      # stage removed is proportional to stage abundance and counts can never go negative
      thisRemvec <- rep(0, nStages*npatches)
      for (p in 1:npatches){
        Ipatch <- (1+4*(p-1)):(4*p)
        subVec <- Vec[(i + 1), Ipatch]
        n_p <- min(group_remvec[p], sum(subVec))    # cannot remove more than are in the group
        if (n_p > 0){
          inds <- rep(1:nStages, times = subVec)    # one entry per individual in the group (its stage)
          drawn <- inds[sample.int(length(inds), size = n_p, replace = FALSE)]
          thisRemvec[Ipatch] <- tabulate(drawn, nbins = nStages)
        }
      }
      
      # if a group could not supply its share, remove the shortfall from remaining individuals
      shortfall <- sum(group_remvec) - sum(thisRemvec)
      if (shortfall > 0){
        remaining <- Vec[(i + 1), ] - thisRemvec
        shortfall <- min(shortfall, sum(remaining))
        if (shortfall > 0){
          inds <- rep(1:length(remaining), times = remaining)
          drawn <- inds[sample.int(length(inds), size = shortfall, replace = FALSE)]
          thisRemvec <- thisRemvec + tabulate(drawn, nbins = length(remaining))
        }
      }
      
      # calculating actual number removed from each group for output
      group_rem <- vector()
      for (p in 1:npatches){
        group_rem[p] <- sum(thisRemvec[(1+4*(p-1)):(4*p)])      # sum every 4 entries in the vector
      }# removed per group
      
      #  into outputs
      Nremoved[idx] <- sum(thisRemvec)  # total number removed   
      Remvec[idx,] <- thisRemvec
      Remgroup[idx,] <- group_rem
      
      # new population size following removals
      Vec[(i + 1), ] <- Vec[(i + 1), ] - thisRemvec   # removals taken from the projected population in year i+1
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
      sizes[p] <- sum(Vec[(i + 1), (1+4*(p-1)):(4*p)])
    }
    
    # filling in group size this year
    Group[(i + 1),] <- sizes
    
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