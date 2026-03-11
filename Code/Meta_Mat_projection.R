# meta projection
# 5.02.26
params<- data.frame(fmax= 0.84360,   # F fecundity max (max cubs per adult female) 
                    Sc_max=0.76,   # max cub survival (equal for sexes)
                    b=0.01,       # temp value- must be calculated from provided datasets
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





# projection of sub pops - outputs vector per patch, abundance per patch, total abundance --------
meta.proj <- function(Umat,   # vector of stage names
                      params2, # adjusted beta = 
                      stagenames,
                      initial,   # list of initial abundances
                      Dmat = list(), 
                      time = 20,
                      return.vec = TRUE,
                      return.mats= FALSE){
  # Time 
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  if(is.list(initial)){   # when initial abundances is a list of vectors, we have more than 1 patch 
    patches = length(initial)
  } else if(is.list(initial) == FALSE ){  # if a vector/ not list, only 1 patch
    patches = 1
    stop("only one vector of abundances provided. Provide abundances per patch as a list, or use rem.proj for single patch projections")
  }
  
  nStages <- length(stagenames)/2      # how many stages
  
  Nf_v <- rep(NA, patches) # vector length patches to fill with adult fem abundance by year
  Nm_v <- rep(NA, patches)
  
  
  for (i in 1:patches) {           # filling each list within vector with Nf and Nm at that patch initially
    Nf_v[i] <- initial[[i]][nStages]   
    Nm_v[i] <- initial[[i]][2*nStages]      # adult male in vector (final entry) - vec of Nms per patch
  }
  
  
  # f calculation first year - how to repeat for each patch
  f_v <- rep(NA, patches)   # blank vec to fill  before adding fertility values straight to amat - f values per patch f[p] 
  
  for (j in 1:patches) { # repeat for each patch
    mat <- mating.func(params,     # density dependent parameters
                       stagenames,   # Stages in life cycle graph (single sex)
                       Nf_v[j],        # Adult and yearling females  
                       Nm_v[j],             # Adult and yearling males
                       Mfunction= "min", 
                       return.mat= FALSE) # don't return mat, just f value
    
    f_v[j] <- mat$f   # adding appropriate fertility rates to list, where entry = patch  
  } 
  
  # Set up the output - how to set up if multiple patches? Lists within lists...
  # differs from AI suggested - change to simplify?
  #{
  out <- list(pop = list(rep(NA, (time + 1))),    #  list of vectors per patch
              vec = list(matrix(0, ncol = length(stagenames), nrow = time + 1)),    # list of matrices per patch
              Nremoved = (rep(NA, 2*nStages)), # total removed in rem year, NO LIST, vec length = num patches
              remvec = list(rep(NA, 2*nStages)))    # vector of length(stages) per patch - list
  
  # renaming before repeating
  colnames(out$vec[[1]]) <- stagenames   # naming cols matrix as stages - cant name lists of obj
  rownames(out$vec[[1]]) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or initial
  
  
  Vec <- rep(out$vec, patches)  # matrices to fill with stage abundance.  row= time, col= stage. out$vec or just vec?
  Pop <- rep(out$pop, patches)       # vector to fill with total pop size each year
  remVec <- rep(out$remvec, patches) # list of vectors per patch, only 1 year 
  
  # filling first row of pop size and abundance vecs
  pop_sizes <- sapply(initial, sum)  # gives total pop size per patch in initial year
  
  for (m in 1:patches){  # filling first row of each list vec and pop size 
    Vec[[m]][1,] <- initial[[m]]          # for each list in Vec, first row is initial         
    
    Pop[[m]][1] <- pop_sizes[m]   # each list row 1 is sum of n0 - issues in plugging into pop
  }
  # } all simplified into fewer lines
  
  # Loop = matrix proj each year  --------
  for (i in 1:time) {   # repeat for as many years as we have
    
    # Mating Fmat creation - f values per patch
    thisNf <- as.vector(rep(NA, patches))  # vec of current Nf for each patch
    thisNm <- as.vector(rep(NA, patches))
    
    for (n in 1:patches){  # loop for each patch this year
      thisNf[n] <- Vec[[n]][i,nStages]  # inputting calculated value to appr place in vec. Issues with using both n and i in this loop 
      thisNm[n] <-  Vec[[n]][i,2*nStages]                                  #repeat for each patch, but as we loop through years will have to update row  of vec selected
    }
    
    # apply mating func to calculate pairs per patch 
    thisMating <- rep(NA, patches)        # list of repeated matrix 
    for (p in 1:patches) { # repeat for each patch
      mat_out <- mating.func(params,     # density dependent parameters
                             stagenames,   # Stages in life cycle graph (single sex)
                             thisNf[p],        # Adult and yearling females  
                             thisNm[p],             # Adult and yearling males
                             Mfunction= "min", 
                             return.mat= FALSE)
      
      thisMating[p] <- mat_out$f  # adding appropriate Fmat to list 
    } 
    
    # ricker density dependence each year - total pop size by patch
    
    thisN <- rep(NA, patches)  # vector length patches to fill with current pop size each patch 
    for (p in 1:patches){
      thisN[p] <- sum(Vec[[p]][i, ])   # length patches, entries = pop size that year. this N entry [1] is pop size(row = year, col = p)
    }
    # Amat each year will vary by patch 
    thisAmat <- list(matrix(0, ncol = ncol(Umat), nrow= nrow(Umat)))
    thisAmat <- rep(thisAmat, patches)        # list of repeated matrix 
    
    for (p in 1:patches){
      amat <- apply.DD(params, 
                       Umat, 
                       thisN[p],   # yearling and adults
                       DDapply="fertility", 
                       stagenames,   # Stages in life cycle graph 
                       thisNf[p],        # Adult female abundance
                       thisNm[p],             # Adult males
                       Mfunction= "min",       # mating function applied
                       return.mat= FALSE)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
      
      thisAmat[[p]] <- amat
    } 
    
    
    # If the projection matrix has any negative values in it, stop iterating and
    # return the projection up until this point.
    
    for (p in 1:patches) {
      if (any(sapply(thisAmat, function(m) any(m < 0, na.rm = TRUE)))) {
        warning(paste("Projection stopped at time step", i, "because the density-dependent matrix has negative values."))
        year_end <- year
        break
      }
    }
    
    # projection for each patch
    for (p in 1:patches){
      Vec[[p]][(i + 1), ] <- thisAmat[[p]] %*% Vec[[p]][i, ] # following year stage vector is this Amat*this abundance * p(stay) (row 1 of Dmat.p) in this patch - loop for patch 1 then 2?
      Pop[[p]][i + 1] <- sum(Vec[[p]][(i + 1), ])   # lapply or sapply? issues running this section - no Pop obj returned
    }
    
    # if pop size <= 0, stop and return
    total_pop<- sum(sapply(Pop, function(x) x[i + 1]))
    if (total_pop <= 0) {
      warning(paste("Projection stopped at time step", i, "because total pop size reached 0 or below"))
      break
    }
    # if any stage becomes negative, set to zero and continue
    if (any(sapply(Vec, function(mat) any(mat < 0, na.rm = TRUE)))) {
      warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
      Vec <- lapply(Vec, function(mat) { mat[mat < 0] <- 0; mat })
    }
  }
  # out objects
  out$pop <- Pop # total group size per patch, do we also want total pop size of pop?
  
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  }
  
  return(out)
}

# SUCCESS!!!!!!

# testing function
params2<- data.frame(fmax= 0.84360,   # F fecundity max (max cubs per adult female) 
                    Sc_max=0.76,   # max cub survival (equal for sexes)
                    b=0.7,       # temp value- must be calculated from provided datasets
                    rep_K= 2.299,          #litter size (K)
                    h= 6   # harem size per male
)
n3 <- list(c(1:4), c(3:6), c(7:10))

test <- meta.proj(Umat,   # vector of stage names
                  params2, # should beta vary by patch?
                  stagenames = stages,
                  initial = n3,   # list of initial abundances
                  time = 30,
                  return.vec = TRUE,
                  return.mats= FALSE)     # group sizes should cap at around 50 individuals max
 




# NEXT STEPS WITH FUNCTION = adding dispersal  ----
# GOAL = Dmats for each patch with movement function. multiply amats by dmat entries then project. calculate number of individuals added to new patch
# for now, p(move) randomly generated

# function to generate dispersal probs given number of patches, sex ratio and group size 

# stating with 2 patches, dispersal between at different rates ----
initial <- list(c(10, 5, 3, 2), c(3, 2, 2, 2))  
# no mating in dispersing individuals: count -> breed -> move. So multiply by Umat only, exclude Fmat

Dmat <- matrix(0, ncol= length(stages), nrow = 2) # row = number patches
colnames(Dmat) <- stages

# prob individ in patch1 remains in patch 1 - rnorm
stay1 <- rnorm(4,mean = 0.5, sd =0.2)  # vec of prob for each stage (4)
Dmat1 <- Dmat
Dmat1[1,] <- stay1 # stay for each class in patch 1 
Dmat1[2,] <- 1- stay1   # for multiple patches, will need multiple values that sum to 1, can't use 1-p

stay2 <- rnorm(4,mean = 0.5, sd =0.2)
Dmat2 <- Dmat
Dmat2[1,] <- stay2 # stay for each class in patch 2 
Dmat2[2,] <- 1 - stay2 


# updating function to include dispersal - only works for 2 patches
meta.proj2 <- function(Umat,   # vector of stage names
                      params, # adjusted beta = 
                      stagenames,
                      initial,   # list of initial abundances
                      Dmat = list(), 
                      time = 20,
                      return.vec = TRUE,
                      return.mats= FALSE){
 
    # Time 
    if (time <= 1) stop("Time must be a positive integer")
    else(time <- as.integer(time))
    
    if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
    
    if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
    
    if(is.list(initial)){   # when initial abundances is a list of vectors, we have more than 1 patch 
      patches = length(initial)
    } else if(is.list(initial) == FALSE ){  # if a vector/ not list, only 1 patch
      patches = 1
      stop("only one vector of abundances provided. Provide abundances per patch as a list, or use rem.proj for single patch projections")
    }
    
    nStages <- length(stagenames)/2      # how many stages
    
    Nf_v <- rep(NA, patches) # vector length patches to fill with adult fem abundance by year
    Nm_v <- rep(NA, patches)
    
    
    for (i in 1:patches) {           # filling each list within vector with Nf and Nm at that patch initially
      Nf_v[i] <- initial[[i]][nStages]   
      Nm_v[i] <- initial[[i]][2*nStages]      # adult male in vector (final entry) - vec of Nms per patch
    }
    
    
    # f calculation first year - how to repeat for each patch
    f_v <- rep(NA, patches)   # blank vec to fill  before adding fertility values straight to amat - f values per patch f[p] 
    
    for (j in 1:patches) { # repeat for each patch
      mat <- mating.func(params,     # density dependent parameters
                         stagenames,   # Stages in life cycle graph (single sex)
                         Nf_v[j],        # Adult and yearling females  
                         Nm_v[j],             # Adult and yearling males
                         Mfunction= "min", 
                         return.mat= FALSE) # don't return mat, just f value
      
      f_v[j] <- mat$f   # adding appropriate fertility rates to list, where entry = patch  
    } 
    
    # Set up the output - how to set up if multiple patches? Lists within lists...
    # differs from AI suggested - change to simplify?
    #{
    out <- list(pop = list(rep(NA, (time + 1))),    #  list of vectors per patch
                vec = list(matrix(0, ncol = length(stagenames), nrow = time + 1)),    # list of matrices per patch
                Nremoved = (rep(NA, 2*nStages)), # total removed in rem year, NO LIST, vec length = num patches
                remvec = list(rep(NA, 2*nStages)))    # vector of length(stages) per patch - list
    
    # renaming before repeating
    colnames(out$vec[[1]]) <- stagenames   # naming cols matrix as stages - cant name lists of obj
    rownames(out$vec[[1]]) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or initial
    
    
    Vec <- rep(out$vec, patches)  # matrices to fill with stage abundance.  row= time, col= stage. out$vec or just vec?
    Pop <- rep(out$pop, patches)       # vector to fill with total pop size each year
    remVec <- rep(out$remvec, patches) # list of vectors per patch, only 1 year 
    
    # filling first row of pop size and abundance vecs
    pop_sizes <- sapply(initial, sum)  # gives total pop size per patch in initial year
    
    for (m in 1:patches){  # filling first row of each list vec and pop size 
      Vec[[m]][1,] <- initial[[m]]          # for each list in Vec, first row is initial         
      
      Pop[[m]][1] <- pop_sizes[m]   # each list row 1 is sum of n0 - issues in plugging into pop
    }
    # } all simplified into fewer lines
    
    # Loop = matrix proj each year  --------
    for (i in 1:time) {   # repeat for as many years as we have
      
      # Mating Fmat creation - f values per patch
      thisNf <- as.vector(rep(NA, patches))  # vec of current Nf for each patch
      thisNm <- as.vector(rep(NA, patches))
      
      for (n in 1:patches){  # loop for each patch this year
        thisNf[n] <- Vec[[n]][i,nStages]  # inputting calculated value to appr place in vec. Issues with using both n and i in this loop 
        thisNm[n] <-  Vec[[n]][i,2*nStages]                                  #repeat for each patch, but as we loop through years will have to update row  of vec selected
      }
      
      # apply mating func to calculate pairs per patch 
      thisMating <- rep(NA, patches)        # list of repeated matrix 
      for (p in 1:patches) { # repeat for each patch
        mat_out <- mating.func(params,     # density dependent parameters
                               stagenames,   # Stages in life cycle graph (single sex)
                               thisNf[p],        # Adult and yearling females  
                               thisNm[p],             # Adult and yearling males
                               Mfunction= "min", 
                               return.mat= FALSE)
        
        thisMating[p] <- mat_out$f  # adding appropriate Fmat to list 
      } 
      
      # ricker density dependence each year - total pop size by patch
      
      thisN <- rep(NA, patches)  # vector length patches to fill with current pop size each patch 
      for (p in 1:patches){
        thisN[p] <- sum(Vec[[p]][i, ])   # length patches, entries = pop size that year. this N entry [1] is pop size(row = year, col = p)
      }
      # Amat each year will vary by patch 
      thisAmat <- list(matrix(0, ncol = ncol(Umat), nrow= nrow(Umat)))
      thisAmat <- rep(thisAmat, patches)        # list of repeated matrix 
      
      for (p in 1:patches){
        amat <- apply.DD(params, 
                         Umat, 
                         thisN[p],   # yearling and adults
                         DDapply="fertility", 
                         stagenames,   # Stages in life cycle graph 
                         thisNf[p],        # Adult female abundance
                         thisNm[p],             # Adult males
                         Mfunction= "min",       # mating function applied
                         return.mat= FALSE)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
        
        thisAmat[[p]] <- amat
      } 
      
      
      # If the projection matrix has any negative values in it, stop iterating and
      # return the projection up until this point.
      
      for (p in 1:patches) {
        if (any(sapply(thisAmat, function(m) any(m < 0, na.rm = TRUE)))) {
          warning(paste("Projection stopped at time step", i, "because the density-dependent matrix has negative values."))
          year_end <- year
          break
        }
      }
      
      remains <- vector("list", length(stages))   # list of remaining individual abundances 
      moves <- vector("list", length(stages))     # list of moves individual abundances
      
      # projection for each patch before dispersal
      for (p in 1:patches) {
        current_vec <- as.numeric(Vec[[p]][i, ])         # current stage-abundance vector
        
        pStay <- as.numeric(Dmat[[p]][1, ])        # stage prob of staying in patch p
        pMove <- as.numeric(Dmat[[p]][2, ])        # stage p(leave) patch p
        
        # abundance BEFORE  survival/rep
        remains[[p]] <- current_vec * pStay        #  number who remain to breed
        moves[[p]]  <- current_vec * pMove       # number who leave (no breeding)
      }
      
      # Project stayers using Amat (they reproduce in original patch)
      remains_proj <- vector("list", patches)
      for (p in 1:patches) {
        remains_proj[[p]] <- as.numeric(thisAmat[[p]] %*% remains[[p]])
      }
      
      # moves:  leave BEFORE mating, only survive (Umat)
      moves_surv <- vector("list", patches)
      for (p in 1:patches) {
        moves_surv[[p]] <- as.numeric(Umat %*% moves[[p]])
      }
      
      # Assemble next-year vectors: each patch gets its own stayers_proj + incoming movers from the other patch
      for (p in 1:patches) {
        other <- ifelse(p == 1, 2, 1)
        
        next_vec <- remains_proj[[p]] + moves_surv[[other]]  # combine remaining abundances with those who leave other patch
        
        # guard against numeric noise/negatives
        next_vec[next_vec < 0] <- 0
        
        Vec[[p]][i + 1, ] <- next_vec
        Pop[[p]][i + 1] <- sum(Vec[[p]][i + 1, ])
      }
      
      # if pop size <= 0, stop and return
      total_pop<- sum(sapply(Pop, function(x) x[i + 1]))
      if (total_pop <= 0) {
        warning(paste("Projection stopped at time step", i, "because total pop size reached 0 or below"))
        break
      }
      # if any stage becomes negative, set to zero and continue
      if (any(sapply(Vec, function(mat) any(mat < 0, na.rm = TRUE)))) {
        warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
        Vec <- lapply(Vec, function(mat) { mat[mat < 0] <- 0; mat })
      }
    }
    # out objects
    out$pop <- Pop # total group size per patch, do we also want total pop size of pop?
    
    if (isTRUE(return.vec)) {
      out$vec <- Vec        
    }
    
    return(out)
}

n0 <- list(c(10,14,6,8),c(6,10,5,2))
dis.trial <- meta.proj(Umat,   # vector of stage names
                        params, # adjusted beta = 
                        stagenames = stages,
                        initial =n0,   # list of initial abundances
                        Dmat, 
                        time = 20,
                        return.vec = TRUE,
                        return.mats= FALSE)

# now we can project with movement between 2 patches, what about greater patch numbers?
Dmat.create <- function(patches, group_size, sex_ratio) { # inputs = number of patches, current sex ratio, current group size
  
}