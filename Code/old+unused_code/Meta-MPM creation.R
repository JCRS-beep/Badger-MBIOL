# meta-population matrix creation
# 10/04/26

# changing Dmat set up - only movement probs per patch
# still 2 patch specific 
# OBLY FOR ADULTS MOVING
 

# next steps = adapting this function for multiple patch inputs = only once i create Dmat automatically from distance and Nf and Nm (logit link)
# changing Dmat set up - only movement probs per patch
Dmat <- matrix(0, ncol= 2, nrow= 3) # row = number patches
colnames(Dmat) <- c("Females", "Males")

# prob individ in patch1 remains in patch 1 = rnorm
set.seed(1)  # replicable
move1 <- rnorm(2,mean = 0.5, sd =0.05)  # vec of prob for each stage (4)
# make sure does not fall bind between 0 and 1
move1[move1<0] <- 0 
move1[move1>1] <- 1
Dmat[1,] <- move1 # stay for each class in patch 1 

set.seed(2)
move2 <- rnorm(2,mean = 0.3, sd =0.1)
Dmat[2,] <- move2 # stay for each class in patch 2 
     
move3 <- rnorm(2,mean = 0.6, sd =0.1)
Dmat[3,] <- move3 # stay for each class in patch 2 

initial_vec <- c(1,4,1,4,2,5,2,6,2,8,3,7)  # 3 patch vec

# function to create large matrix - call within projection. Acts as our apply DD function to create Amat each year
# OBLY FOR ADULTS MOVING
meta.mat <- function(Umat,   # vector of stage names
                     initial_vec, # big vec length = 4* npatches
                     Dmat,
                     return.mats = FALSE){    
  
  stagenames <- c(colnames(Umat))
  nDim <- length(stagenames)  # number of stages
  nPatches <- nrow(Dmat)
  
  # setting up out obj
  # meta_out <- list(meta_mat = matrix())  # what do we want from this?
  
  
  # Amat creation for more patches
  A_list <- list()
   for (p in 1:nPatches){
     thisN <- sum(initial_vec[((p*nDim)-3): (p * nDim) ])
     thisNf <- (p*nDim)-2
     thisNm <- p * nDim
     
  this_amat <- apply.DD(params, 
                        Umat, 
                        thisN,   # yearling and adults
                        DDapply="fertility", 
                        stagenames,   # Stages in life cycle graph 
                        thisNf,        # Adult female abundance
                        thisNm         # Adult males
                        )
  
   A_list[[p]] <- this_amat    # output list length = nPatches
  }
  
  # large mat set up
  mat <- matrix(0, nrow = nrow(Umat), ncol = ncol(Umat))   # setting base matrix dims
  meta <- cbind(mat, mat)  # creates 2 patch mat
  
  if (nPatches > 2){   
  for (i in 1:(nPatches-2)){   # need as many reps as there are patches, above row adds first 2 why does it rep too many times?
  meta <- cbind(meta, mat)   # adds on 1 dim each loop
  }
  }
  
  meta_mat <- rbind(meta, meta)
  for (i in 1:(nPatches-2)){   # need as many reps as there are patches, above row adds first 2 why does it rep too many times?
    meta_mat <- rbind(meta_mat, meta)   # adds on 1 dim each loop
  }

  
  for (p in 1:nPatches){
    meta_mat[((p*nDim)-3):(p*nDim),((p*nDim)-3):(p*nDim)] <- A_list[[p]]
  }
  
  colnames(meta_mat) <- rep(stagenames, nPatches)  # rep times = npatches
  rownames(meta_mat) <- rep(stagenames, nPatches)
  
  # staying probs in relevant matrix entries - cols 2 or 4
  for (p in 1 : nPatches){   # filling female and male remaining probs 
  meta_mat[,((p*nDim)-2)] <- meta_mat[,((p*nDim)-2)] * (1 - Dmat[p,1])  # females remaining in patch
  meta_mat[,(p * nDim)] <- meta_mat[,(p * nDim)] * (1 - Dmat[p,2])   # male staying in patch
  }
  
  # creating another base matrix?
  # if all inidivuals are moving, bmat = umat
  bmat <- mat
  bmat[2,2] <- Umat[2,2]
  bmat[4,4] <- Umat[4,4]
  
  for (p in 1:nPatches){
   this_bmat <- bmat
   this_bmat <- this_bmat[,2] * (Dmat[p,1]/ (nPatches-1))    # dividing prob across remaining patches - assumes all equal 
   this_bmat <- this_bmat[,4] * (Dmat[p,2]/ (nPatches-1))
     
  if(nPatches == 2){
   # for 2 patches
   # when p = 1
   meta_mat[5:8, 1:4] <- this_bmat
   
   # when p = 2
   meta_mat[1:4, 5:8] <- this_bmat
   
   } else if(nPatches > 2){
    # when p =1  (patch 1), all rows below = bmat. col lims = 1:4 = p*nDim - 3 : p * nDim
    meta_mat[(nDim + 1):nrow(meta_mat), 1:4] <- bmat   # all rows below (5;end)
     
   }
   

  }
  
  # survival * move probs for indivs moving patch 1 to 2
  meta_mat[6,2] <- Saf *  Dmat[1,1]  # below patch 1 [((p * nStages) + 2), nStages -2] = Af survival * Dmat[p, 1 = F]
  meta_mat[8,4] <- Umat[4,4] * Dmat[1,2]  # Am survival * move from 1 to 2
  

    meta_mat[,((p*nDim)-2)] <- meta_mat[,((p*nDim)-2)] * (1 - Dmat[p,1])  # females moving 
    meta_mat[,(p * nDim)] <- meta_mat[,(p * nDim)] * (1 - Dmat[p,2])   # male moving
  
  
  # for movement * survival in off diags - 
  
  meta_mat[8,4] <- Umat[4,4] * Dmat[1,2]  # Am survival * move from 1 to 2
  meta_mat[4,8] <- Umat[4,4] * Dmat[2,2]  # survival * 
  
  # for fems, multiiply whole cols before adding fem survival
  meta_mat[,2] <- meta_mat[,2] * (1 - Dmat[1,1])
  meta_mat[6,2] <- Umat[2,2] * (1 - Dmat[1,1])  # Af survival * 
  
  meta_mat[,6] <- meta_mat[,6] * Dmat[2,1]
  meta_mat[2,6] <- Umat[2,2] * (1 - Dmat[2,1])
  
  return(meta_mat)
}
   

meta.mat <- function(Umat,   # vector of stage names
                     params,
                     initial_vec, # big vec length = 4* npatches
                     Dmat){    
  
  stagenames <- c(colnames(Umat))
  nDim <- length(stagenames)  # number of stages
  nPatches <- nrow(Dmat)
  
  # setting up out obj
  # meta_out <- list(meta_mat = matrix())  # what do we want from this?
  
  # setting base matrix dims
  mat <- matrix(0, nrow = nrow(Umat), ncol = ncol(Umat))
  
  # Amat creation for more patches
  A_list <- list()
  for (p in 1:nPatches){
    thisN <- sum(initial_vec[((p*nDim)-3): (p * nDim) ])
    thisNf <- (p*nDim)-2
    thisNm <- p * nDim
    
    this_amat <- apply.DD(params, 
                          Umat, 
                          thisN,   # yearling and adults
                          DDapply="fertility", 
                          stagenames,   # Stages in life cycle graph 
                          thisNf,        # Adult female abundance
                          thisNm         # Adult males
    )
    
    A_list[[p]] <- this_amat    # output list length = nPatches
  }
  
  # large mat set up
  rep_cb <- function(mat, times) do.call(cbind, rep(list(mat), times))   # function to cbind repeated number of times
  meta <- rep_cb(mat, nPatches)
  
  rep_rb <- function(mat, times) do.call(rbind, rep(list(mat), times))  # function to repeatedly stack mat vertically with rbind
  meta_mat <- rep_rb(meta, nPatches)   # adds on 1 dim each loop
  
  # filling diags with Amats from list
  for (p in 1:nPatches){
    meta_mat[((p*nDim)-3):(p*nDim),((p*nDim)-3):(p*nDim)] <- A_list[[p]]
  }
  
 # colnames(meta_mat) <- rep(stagenames, nPatches)  # rep times = npatches
 # rownames(meta_mat) <- rep(stagenames, nPatches)
  
  # staying probs in relevant matrix entries 
  for (p in 1 : nPatches){   # filling female and male remaining probs 
    meta_mat[,((p*nDim)-2)] <- meta_mat[,((p*nDim)-2)] * (1 - Dmat[p,1])  # females remaining in patch
    meta_mat[,(p * nDim)] <- meta_mat[,(p * nDim)] * (1 - Dmat[p,2])   # male staying in patch
  }
  
  # base matrix to fill if all inidivuals are moving, bmat = umat
  bmat <- mat
  bmat[2,2] <- Umat[2,2]
  bmat[4,4] <- Umat[4,4]
  
  for (p in 1:nPatches){
    this_bmat <- bmat
    this_bmat[,2] <- this_bmat[,2] * (Dmat[p,1]/ (nPatches-1))    # dividing prob across remaining patches - assumes all equal 
    this_bmat[,4] <- this_bmat[,4] * (Dmat[p,2]/ (nPatches-1))
    
    col_lims <- ((p * nDim)-3) :(p * nDim)  # defines cols for each patch 
    
    # for each patch, bmat identical (survival constant, movement divided equally). SO each col either bmst or Amat
    # Amat position in meta mat = patch number (when p = 2, second quadrant from top)
    # number of bmats above amat = p-1  (when p = 2, 1 bmat above Amat)
    # number of bmats below = nPatches - p   (when p = 2 and nPatches = 3, one below)
    # Define boundaries of Amat position, bmats until and after that within same cols
    
    if(p == 1){   # for first entry, A at top, everything below bmat
      # won't work, need to rbind with custom func
      times <- nPatches - p    # defining times 
      meta_mat[5:nrow(meta_mat), col_lims] <- rep_rb(this_bmat, times) # rep for all quadrants below - NEED TO MAKE SURE THEY ARE STACKED VERTICALLY
      
    } else if(1 < p && p > nPatches){   # middle rows - formula to work out bmats above and below
      top <- ((p * nDim) - 3)  # upper row of this Amat
      bot <- p * nDim   # bottom row of this Amat

      meta_mat[1:(top - 1), col_lims] <- this_bmat
      meta_mat[(bot +1) :nrow(meta_mat), col_lims] <- this_bmat
      
    } else if(p == nPatches){  # for final loop, Amat at bottom, everything above bmat
      meta_mat[1:((p-1)*nDim), col_lims] <- this_bmat
      
    }
  }
  return(meta_mat)
}

# testing function 
init <- c(1,4,1,4,2,6,2,6,2,8,2,7)  # 3 patches
meta.mat(Umat,   # vector of stage names
         params,
         init, # big vec length = 4* npatches
         Dmat)
         # off diags not filled correct


# testing matrix sytax
big_mat <- matrix(0, nrow = 12, ncol = 12)
a <- seq(1:16)
mini_mat <- matrix(a, nrow = 4, ncol = 4, byrow= TRUE)
big_mat[1:8, 1:4] <- mini_mat

rep_block <- function(mat, times) do.call(rbind, rep(list(mat), times))   # AI PMO it solves things so quick

big_mat[1:8, 1:4] <- rep_block(mini_mat, 2)
  
  

