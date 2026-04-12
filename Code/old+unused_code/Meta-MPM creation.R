# meta-population matrix creation
# 10/04/26

# Creating large matrix, setting up dmat with custom function. 
# ONLY FOR ADULTS MOVING
 

# changing Dmat set up - only movement probs per patch, nrow  nPatches, ncol  number moving indivs (just adults or all stages?)
Dmat <- matrix(0, ncol= 2, nrow= 3) # row = number patches
colnames(Dmat) <- c("Females", "Males")

# prob individ in patch1 remains in patch 1 = rnorm
set.seed(1)  # replicable
move1 <- rnorm(2,mean = 0.5, sd =0.05)  # vec of prob for each stage (4)
# make sure does not fall bind between 0 and 1
move1[move1<0] <- 0 
move1[move1>1] <- 1
Dmat[1,] <- move1 # stay for each class in patch 1 

set.seed(1)
move2 <- rnorm(2,mean = 0.3, sd =0.1)
Dmat[2,] <- move2 # stay for each class in patch 2 
     
move3 <- rnorm(2,mean = 0.6, sd =0.1)
Dmat[3,] <- move3 # stay for each class in patch 2 

initial_vec <- c(1,4,1,4,2,5,2,6,2,8,3,7)  # 3 patch vec

# function to create large matrix - call within projection. Acts as our apply DD function to create Amat each year

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
  
  # More difficult to id matrix positions once cols are named
  # colnames(meta_mat) <- rep(stagenames, nPatches)  # rep times = npatches
  # rownames(meta_mat) <- rep(stagenames, nPatches)
  
  # staying probs in relevant matrix entries 
  for (p in 1 : nPatches){   # filling female and male remaining probs 
    meta_mat[,((p*nDim)-2)] <- meta_mat[,((p*nDim)-2)] * (1 - Dmat[p,1])  # females remaining in patch
    meta_mat[,(p * nDim)] <- meta_mat[,(p * nDim)] * (1 - Dmat[p,2])   # male staying in patch
  }
  
  # base matrix to fill with survival rates of moving indivs if all inidivuals are moving, bmat = umat
  bmat <- mat
  bmat[2,2] <- Umat[2,2]   # adult fem survival 
  bmat[4,4] <- Umat[4,4]   # adult m survival
  
  for (p in 1:nPatches){    # for each patch, work in each col of matrix
    col_lims <- ((p * nDim)-3) :(p * nDim)  # defines col boundaries for each patch 
    
    this_bmat <- bmat
    # dividing movement prob across remaining patches - assumes moves equally likely
    this_bmat[,2] <- this_bmat[,2] * (Dmat[p,1]/ (nPatches-1))      # female survival * movement prob in this patch
    this_bmat[,4] <- this_bmat[,4] * (Dmat[p,2]/ (nPatches-1))      # male survival * movement prob 
    
    
    # number of bmats above amat = p-1  (when p = 2, 1 bmat above Amat)
    # number of bmats below = nPatches - p   (when p = 2 and nPatches = 3, one below)
    # Define boundaries of Amat position, bmats until and after that within same cols
    
    top <- ((p * nDim) - 3)  # upper row of this Amat
    topN <- p - 1    # how many above? what if 0?
    
    bot <- p * nDim   # bottom row of this Amat
    botN <- nPatches - p   # what if 0?
    if(topN > 0){
      meta_mat[1:(top - 1), col_lims] <- rep_rb(this_bmat, topN)
    }
    if(botN > 0){
      meta_mat[(bot +1) : nrow(meta_mat), col_lims] <- rep_rb(this_bmat, botN)
    }
    
  }
  return(meta_mat)
}   # almost, first and final cols correct but bmat not included in middle. When run individually in lines, this builds correctly. Issue with loop set up?


# testing function 
init <- c(1,4,1,4,2,6,2,6,2,8,2,7)  # 3 patches
meta_mat <- meta.mat(Umat,   # vector of stage names
         params,
         init, # big vec length = 4* npatches
         Dmat)
         # SUCCESS!


  
  
