## Badger metapopulation function, rewritten as individual-based
# Written by Chrissy Hernandez 
# 11 November 2026


# DD moving projection function as IBM -----
# incorporating density dependent movement probs in projection function
proj_meta_DDdisp_norem_IBM <- function(Umat,
                                       initial,
                                       params,
                                       stagenames,
                                       npatches,
                                       time,
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
    prob_arrival <- rep(1/(npatches-1), npatches)
    
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