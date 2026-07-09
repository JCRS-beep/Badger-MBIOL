# meta projection
# 5.02.26
# projecting the large matrix forward given an initial matrix

# Dmat creation function - can be customised to how many moving individuals we have 
Dmat.create <- function(colNames, nPatches) { 
  Dmat <- matrix(0, ncol = length(colNames), nrow = nPatches) # row = number patches
  colnames(Dmat) <- colNames
  
  for (p in 1:nPatches){
    move <- rnorm(length(colNames), mean = 0.2, sd = 0.05) # rnorm to calculate pmove from each patch # IMPROVEMENT = eq to calculate move prob value given vars
    move[move < 0] <- 0    # if any values lower than 0, set to 0
    Dmat[p,] <- move
  }
  return(Dmat)
} 

dtest <- Dmat.create(colNames = c("Yearling_female", "Adult_female", "Yearling_male",  "Adult_male"), 
                     nPatches = 3)

# load Umat, params, stages and meta script
# meta.projection - ORIGINAL using large meta-matrix
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
  patches = length(initial) / nDim # how many times are stages repeated?
  # if(is.integer(patches) == FALSE) stop("Initial vector provided is incorrect length, must be nPatches *length(patch_vec)")
  
  # Set up the output - how to set up if multiple patches? Lists within lists...
  
  out <- list(vec = list(matrix(0, ncol = length(stagenames), nrow = time + 1)),    # list of matrices per patch. # I dont like the formatting, why $vec$vec[[i]]?
              group_size = matrix(0, ncol =  patches, nrow = time + 1),    #   each row of matrix is patch size that year
              total_pop = vector(),   # vector for total pop size per year
              Nremoved = vector(), # total removed in rem year, vec length = num patches
              patch_rem = list(), # number removed per patch
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
    # meta creation function to get large matrix for this year - don't need large mat if working in loops
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
  
      # moving INTO patch?
      new_vec <- this_mat[p,] - floor(movers_mat[p, ] * (thisDmat[p,])) + colSums(move_mat[-p,])   # those that remain in patch || vec - move out 
      # adding indivs that move in from other patches
      # if any negatives, set to 0
      
      new_mat[p, ] <- new_vec
    }
    
    # convert back to vec form? or just update meta create to accomodate matrix entry?
    
     if (any(new_mat < 0, na.rm = TRUE) || any(is.na(new_mat))) {
       new_mat[new_mat < 0] <- 0
       new_mat[is.na(new_mat)] <- 0
    }
    
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



# No meta matrix vers - FORGET BIG MAT, WORK IN LOOPS 
meta.proj <- function(Umat,   # vector of stage names
                      params, # adjusted beta = tighter density constraints per group
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
              patch_rem = list(), # number removed per patch
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
   
    # Calculating random movement this year per patch
    thisDmat <- Dmat.create(colNames, patches) # creates random p(move) per patch, either single prob per patch or varies by stage/sex
    
    
    # calculating post move vec for use in meta_create
    # NEED TO FIX ROUNDING ISSUE
    # current vec as string
    this_mat <- t(sapply(Vec, function(x) x[i,]))    # getting initial vec as a matrix with this years abundance per patch (patch = row)
    movers_mat <- this_mat  
    movers_mat[,1] <- 0    # non moving individuals as 0 -  assumes only adults move
    movers_mat[,3] <- 0 
    
    move_mat <- matrix(0 , nrow = nrow(movers_mat), ncol = ncol(movers_mat))
    post_mat <- move_mat
    
   # MOVEMENT -----
      # dividing moving individuals across patches to fix rounding issues
      patch <- seq(1:patches)       # vector representing patches. how to exclude movement into own patch
      arrival_mat <- matrix(0, nrow = patches, ncol = patches)   # want this for each sex
      dimnames(arrival_mat) <- list(origin = patch, destination = patch)  # doesn't distinguish sexes
      colnames(arrival_mat) <- patch
      rownames(arrival_mat) <- patch
      
      arrival_matF <- arrival_mat
      arrival_matM <- arrival_mat
      
      for (p in 1:patches){   # how many males and females leaving given patch, where they arrive
      move_mat[p,] <- floor(movers_mat[p,] * thisDmat[p,])   # leaving each patch
        
      move_nf <- move_mat[p,2]  # number fems leaving group
      move_nm <- move_mat[p,4]  # males leaving group
      
      if(sum(move_mat[p,] > 0, na.rm = TRUE)) {   # only run if we have moving individuals
      dest_f <- sample(patch[-p], size = sum(move_nf), replace = TRUE)  # which patches receive these new fems?
      dest_m <- sample(patch[-p], size = sum(move_nm), replace = TRUE)  # which patches receive these new fems?
      
      # for col in matrix, mark how many times it occurs in vec
      f_counts <- tabulate(match(dest_f, patch), nbins = patches)   # matches individuals move to destination patch
      m_counts <- tabulate(match(dest_m, patch), nbins = patches)
      
      arrival_matF[p,] <- f_counts    # filling this row with where individuals arrive
      arrival_matM[p,] <- m_counts
      }      
      
      # moving INTO patch?
      post_vec <- this_mat[p,] - move_mat[p,]       # remaining = initial - move out 
      # adding females that move in from other patches
      post_vec[2] <- post_vec[2] + sum(arrival_matF[,p]) # adding any females that are destined for this patch
      post_vec[4] <- post_vec[4] + sum(arrival_matM[,p]) # adding any females that are destined for this patch
      

    if (any(post_vec < 0, na.rm = TRUE) || any(is.na(post_vec))) {
      post_vec[post_vec < 0] <- 0
      post_vec[is.na(post_vec)] <- 0
    }
      }
    
    # SURVIVAL AND REPRODUCTION AT PATCH LEVEL
    # Amat creation for patches
    a_list <- list()

    for (p in 1:patches){
      thisN <- sum(post_vec)   # post move abundances for fertiltiy calc
      thisNf <- post_vec[2]
      thisNm <- post_vec[4]
      
      this_amat <- apply.DD(params, 
                            Umat, 
                            thisN,   # yearling and adults POST MOVEMENT 
                            DDapply ="fertility", 
                            stagenames,   # Stages in life cycle graph 
                            thisNf,        # Adult female abundance
                            thisNm         # Adult males
      )
    # If the projection matrix has any negative values in it, stop iterating and
    # return the projection up until this point.
    if (any(sapply(this_amat, function(m) any(m < 0, na.rm = TRUE)))) {
      warning(paste("Negative values in projection matrix at time step", i,
                    "- setting negative entries to 0 and continuing."))
      meta_mat <- lapply(this_amat, function(mat) { mat[mat < 0] <- 0; mat}) 
      
      a_list[[p]] <- this_amat    # output list length = nPatches
      
    }
    
    # multiply by vec to project
    next_vec <- floor(as.numeric(this_amat %*% post_vec))    # vector abundance result following year
    # check if any values < 0 
    if (any(next_vec < 0, na.rm = TRUE) || any(is.na(next_vec))) {
      next_vec[next_vec] <- 0
      next_vec[is.na(next_vec)] <- 0
    }      
    }
    
    # filling next row of vec in each list
    for(p in 1:patches){
      Vec[[p]][i + 1,] <- next_vec[((p*nDim)-3) : (p*nDim)]      # the next row for each list vec
    } 
    # if any stage becomes negative, set to zero and continue
    if (any(sapply(Vec, function(mat) any(mat < 0, na.rm = TRUE || is.na(mat))))) {
      warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
      Vec <- lapply(Vec, function(mat) { mat[mat < 0 || is.na(mat)] <- 0; mat })
    }

    for(p in 1:patches){
      # calculating group size per patch this year
      Group[i + 1,] <- sapply(Vec, function(x) sum(x[i +1, ]))  # list of group size vectors 
    }
    # if any stage becomes negative, set to zero and continue
    if (any(sapply(Group, function(mat) any(mat < 0, na.rm = TRUE|| is.na(mat))))) {
      warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
      Group <- apply(Group, function(mat) { mat[mat < 0] <- 0; mat })
    }
    
    Pop[i+1] <- sum(Group[i + 1,])
    if (Pop[i] <= 0 || is.na(Pop[i])) {       # if pop size <= 0, stop and return
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

# initial vec generation - 3 patches
sizes <- sample(2:20, 3, replace = TRUE)   # group sizes to generate
vec <- matrix(0, nrow = length(sizes), ncol = 4)
 for( i in 1:length(sizes)){
  vec[i,] <- round(sizes[i] *stagedist)
}

patch_n0 <- as.vector(t(vec))

meta_test <- meta.proj(Umat,   # vector of stage names
                       params, # adjusted beta = tighter density constraints per group
                       stagenames = stages,
                       initial = patch_n0,   # list of initial abundances - string length patches * nStages
                       colNames = c("females", "males"), # length = nCol Dmat
                       time = 20,
                       return.vec = TRUE,
                       return.group = TRUE)
# ISSUES - still getting NA values in vec, group and total pop. 


# meta removals
# setting removal goal across population, then dividing across number of patches. How to account for patches of variable size?
# patch proportions - proportion each patch contributes to pop size. if N = 100 and patch = 20, prop = 0.2. Minus 0.2 of goal?
meta.rem <- function(Umat,   # vector of stage names
                      params, # adjusted beta = 
                      stagenames,
                      initial,   # list of initial abundances - string length patches * nStages
                      stagedist,  # proportion of each stage
                      colNames, # length = nCol Dmat
                      time = 20,
                      intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                      remyear = integer(0),  # removal year = vector of years 
                      rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                      bias = NULL , # strength of bias as percentage (range??)
                      return.vec= TRUE, 
                      return.group = FALSE,
                      return.group_remvec = TRUE){
  
  # Time 
  if (time <= 1){ stop("Time must be a positive integer")
  }else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  if (length(remyear) > 0) {
    # ensure remyear are integers and between 1 - time
    
    remyear <- unique(as.integer(remyear))
    remyear <- remyear[remyear >= 1 & remyear <= time]
    if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
  }
  
  nDim <- length(stagenames)     # how many stages
  patches <-  length(initial)/ nDim # how many times are stages repeated?
  nYears <- length(remyear)   # over how many years are removals taking place?  
  
  # if(is.integer(patches) == FALSE) stop("Initial vector provided is incorrect length, must be nPatches*length(patch_vec)")
  
  # Set up the output - how to set up if multiple patches? Lists within lists...
  
  out <- list(vec = list(matrix(0, ncol = length(stagenames), nrow = time + 1)),    # list of matrices per patch. # I dont like the formatting, why $vec$vec[[i]]?
              group_size = matrix(0, ncol =  patches, nrow = time + 1),    #   each row of matrix is patch size that year
              total_pop = vector(),   # vector for total pop size per year
              Nremoved = vector(), # total removed in rem year, NO LIST, vec length = num patches
              patch_rem = matrix(0, nrow = nYears, ncol = patches), # number removed per patch
              remvec = list(matrix()))    # vector of length(stages) per patch - list
  
  # renaming before repeating
  colnames(out$vec[[1]]) <- stagenames   # naming cols matrix as stages - cant name lists of obj
  rownames(out$vec[[1]]) <- 0:time   # rows correspond to each year of projection. Row 0 = initial or initial
  
  rownames(out$group_size) <- 0:time
  
  Vec <- rep(out$vec, patches) # matrices to fill with stage abundance.  row= time. empty list of matrices 
  Group <- out$group_size
  Pop <- out$total_pop       # vector to fill with total pop size each year
  Nrem <- out$Nremoved
  Prem <- out$patch_rem # removed per patch?
  Remvec <- list(matrix(0, ncol = 4, nrow = nYears)) # list of vectors per patch, only 1 year 
  
  # filling initial year for outputs
  abun <- matrix(initial, ncol = 4, byrow = TRUE)   # initial matrix with vec per row, nrow = patches
  
  for (p in 1:patches){  # filling first row of each list vec 
    Vec[[p]][1,] <- abun[p,]          # for each list in Vec, first row is initial 
  }
  
  Group[1,] <- rowSums(abun)   # summing each row of matrix which contains abundance vec. Write in terms of vec?
  Pop[1] <- sum(Group[1,]) # first pop entry is sum of group sizes
  
  # Map removal year -> index into Remvec / Nremoved:
  rem_index <- if(length(remyear) > 0) {setNames(seq_along(remyear), remyear)  # takes length rem year, names assigned as seq of integers 
  }else integer(0)   # otherwise remyear = 0 
  # unexpected else?
  
  # Loop = matrix proj each year  --------
  for (i in 1:time) {   # repeat for as many years as we have
    # Calculating random movement this year per patch
    thisDmat <- Dmat.create(colNames, patches) # creates random p(move) per patch, either single prob per patch or varies by stage/sex
    
    # calculating post move vec for use in meta_create
    # current vec as string
    this_mat <- t(sapply(Vec, function(x) x[i,]))    # getting initial vec as a matrix with this years abundance per patch (patch = row)
    movers_mat <- this_mat  
    movers_mat[,1] <- 0    # non moving individuals as 0 - following assumes only adults move
    movers_mat[,3] <- 0 
    
    move_mat <- matrix(0 , nrow = nrow(movers_mat), ncol = ncol(movers_mat))
    post_mat <- move_mat
    for (p in 1:patches){    # all adults equal p moving? Adult f and m?
      move_mat[p,] <- floor(movers_mat[p,] * (thisDmat[p,]/ (patches -1)))    # only works if Dmat 1 row (equal across sexes). 
      # numbers moving from 1 into any other
      
      # moving INTO patch?
      post_vec <- this_mat[p,] - floor(movers_mat[p, ] * (thisDmat[p,])) + colSums(move_mat[-p,])   # those that remain in patch || vec - move out 
      # adding indivs that move in from other patches
      # if any negatives, set to 0
      
      post_mat[p, ] <- post_vec
    }
    
    # convert back to vec form? or just update meta create to accomodate matrix entry?
    
    if (any(post_mat < 0, na.rm = TRUE) || any(is.na(post_mat))) {
      post_mat[post_mat<0] <- 0
      post_mat[is.na(post_mat)] <- 0
    }
    
    # creating meta mat using this years abundances
    this_meta <- meta.mat(Umat,   # vector of stage names
                          params,  # stronger params here to prevent group sizes reaching too high?
                          post_mat, # matrix with post movement abundances
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
    if (any(sapply(Vec, function(mat) any(mat < 0, na.rm = TRUE)))) {
      warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
      Vec <- lapply(Vec, function(mat) { mat[mat < 0] <- 0; mat })
    }
    
    for(p in 1:patches){
      # calculating group size per patch this year
      Group[i + 1,] <- sapply(Vec, function(x) sum(x[i +1, ]))  # list of group size vectors 
    }
    # if any group size becomes negative, set to zero and continue
    if (any(sapply(Group, function(mat) any(mat < 0, na.rm = TRUE)))) {
      warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
      Group <- apply(Group, function(mat) { mat[mat < 0] <- 0; mat })
    }
    
    Pop[i+1] <- sum(Group[i + 1,])  # filling pop size this year
    if (Pop[i] <= 0) {       # if pop size <= 0, warn and return year
      warning(paste("Projection stopped at time step", i, "because total pop size reached 0 or below"))
      break
    }
    
    # checking if a rem year
    if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
      idx <- rem_index[as.character(i)]  # value of rem_vec entry
      
      # setting removal goals
      goal <- round(Pop[remyear[1]] * (intensity/100))   # goal to remove  - pop size before first rem  
      year_goal <- round(goal/nYears)    # how many removed per year
      
      # how many removed per patch - must change if group removals
      patch_dist <-  rep(1/patches, patches)                # vector of rems across patches  - how to remove the year goal across patches? If we round, might not add to 1
      patch_goal <- year_goal * patch_dist     # vec of removal goal per patch, define dist or divide equally across patches?
     
      
      # pop removal -------------------------
      if(rem_strat == "random"){
        if (is.null(bias) == FALSE) paste("ignoring bias value as removal is random across ages and sexes")
        #generating the distribution - varies with rem strat
        dist <- stagedist
        rem <- matrix(0, ncol = 4, nrow = patches)
        for (p in 1:patches){
        rem[p,] <-  patch_goal[p] * dist    # patch goal is a vector - must create a list or matrix for rem for each patch
        } 
      #  patch_rem <- rem * dist    # if removal across patches equal, the remvec per patch is equal
      
        # biases -----
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
          rem <-  year_goal * dist     # number to remove per stage
        }
        
        
        # specified
      } else if (is.numeric(intensity) && is.numeric(rem_strat)){
        
        nbi <- length(rem_strat)  # how many elements provided?
        
        # WARNING - only works if rem_strat(length = 1)
        dist <- stagedist 
        dist[rem_strat] <- dist[rem_strat] + bias   # increasing specified element
        dist[-rem_strat] <- dist[-rem_strat] - bias/3    # bias removed from others must divide by 3
        
        rem <-  year_goal * dist    # stages removed per year is total removed per year * stage props
      }
      
      # ------------------
      
      # calculating number removed from each stage
      thisRemvec <- round(rem)   # actual removed = rounded abundance * dist
      thisPrem <- rowSums(thisRemvec)  # vec of how many per patch
      thisRem <- sum(thisRemvec)  # total number removed 
  
      
      #  into outputs
      Nrem[idx] <- thisRem  # nth vector entry 
      Prem[idx,] <- thisPrem
      Remvec[[idx]] <- thisRemvec   # this remvec is a matrix - vector for each patch. Current method means all remvecs identical
      
      # new population size following removals
      for(p in 1:patches){
        Vec[[p]][i + 1,] <- Vec[[p]][i ,]  - thisRemvec[p,]  # year after remyear = 2 rows later filled with new stage vec. need to repeat for each list
      } 
      # if any stage becomes negative, set to zero and continue
      if (any(sapply(Vec, function(mat) any(mat < 0, na.rm = TRUE)))) {
        warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
        Vec <- lapply(Vec, function(mat) { mat[mat < 0] <- 0; mat })
      }
      
      for(p in 1:patches){
        # calculating group size per patch this year
        Group[i + 1,] <- sapply(Vec, function(x) sum(x[i + 1, ]))  # list of group size vectors 
      }
      # if any group size becomes negative, set to zero and continue
      if (any(sapply(Group, function(mat) any(mat < 0, na.rm = TRUE)))) {
        warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
        Group <- apply(Group, function(mat) { mat[mat < 0] <- 0; mat })
      }
      
      Pop[i+1] <- sum(Group[i + 1,])  # filling pop size this year
      if (Pop[i] <= 0) {       # if pop size <= 0, warn and return year
        warning(paste("Projection stopped at time step", i, "because total pop size reached 0 or below"))
        break
        }
      }
    }

  
  # out objects
  out$total_pop <- Pop 
  out$Nremoved <- Nrem
  out$remvec <- Remvec
  
  if (isTRUE(return.group_remvec)) {
    out$patch_rem <- Prem  # number removed per patch
  } 

  if (isTRUE(return.group)) {
    out$group_size <- Group        
  }
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 
  
  return(out)
}


# testing this function
rem_test <- meta.rem(Umat,   # vector of stage names
                     params, # adjusted beta = 
                     stagenames = stages,
                     initial,   # list of initial abundances - string length patches * nStages
                     stagedist,  # proportion of each stage
                     colNames = c("Adults"), # length = nCol Dmat
                     time = 20,
                     intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                     remyear = 5,  # removal year = vector of years 
                     rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                     bias = NULL , # strength of bias as percentage (range??)
                     return.vec= TRUE, 
                     return.group = TRUE,
                     return.group_remvec = TRUE)

# NOTES when b = 0.007 used, recovery more likely. When stricter dd, extinction common. How to define this patch based beta? 


# also needed = prop vers so we can do continuous rems