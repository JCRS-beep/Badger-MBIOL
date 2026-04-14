# meta projection
# 5.02.26
# projecting the large matrix forward given an initial matrix

# load (set up script) Umat, params, stages and meta script

# updating function to include dispersal - only works for 2 patches
# meta.projection
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
  Pop[1] <- sum(Group[1]) # first pop entry is sum of group sizes
  
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
      if (any(sapply(group_size, function(mat) any(mat < 0, na.rm = TRUE)))) {
        warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
       # group_size <- group_size, function(mat) { mat[mat < 0] <- 0; mat })
      }
      
    # if pop size <= 0, stop and return
    Pop[i+1] <- sum(group_size[i + 1,])
    if (Pop[i] <= 0) {
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


# test
params2<- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
                     b=0.007,       # temp value- must be greater than in non-spatial model as max group sizes are around 1/10 as large
                     rep_K= rogers_k,          # max litter size (K), 
                     h= 6)
                    
# 3 patch initial vec
init <- c(4,5,2,6,2,9,2,6,2,8,2,7)  # 3 patches

# first projection attempt
proj_test <- meta.proj(Umat,   # vector of stage names
                       params2, # adjusted beta = 
                       stages,
                       init,   # list of initial abundances - string length patches * nStages
                       colNames, # length = nCol Dmat
                       time = 20,
                       return.vec = TRUE,
                       return.group = TRUE)