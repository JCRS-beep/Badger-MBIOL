# meta-population matrix creation
# 10/04/26

# Creating large matrix, setting up dmat with custom function. 
# ONLY FOR ADULTS MOVING - will need editing is movement differs between sexes



# function to create large matrix - call within projection. Acts as our apply DD function to create Amat each year
# Input post movement vecs to calculate dd fertility after individuals have moved into this years group to mate.
meta.mat <- function(Umat,   # vector of stage names
                     params,
                     post_mat, # matrix with post movement abundances
                     Dmat, 
                     stagenames, 
                     nDim, 
                     nPatches){      # random movement probs for males and females or per patch?
    # setting base matrix dims
  mat <- matrix(0, nrow = nDim, ncol = nDim)
  
  # Amat creation for more patches
  a_list <- list()
  
  for (p in 1:nPatches){
    thisN <- rowSums(post_mat)[p]   # post move abundances for fertiltiy
    thisNf <- post_mat[p, 2]
    thisNm <- post_mat[p,4]
  
    this_amat <- apply.DD(params, 
                          Umat, 
                          thisN,   # yearling and adults
                          DDapply="fertility", 
                          stagenames,   # Stages in life cycle graph 
                          thisNf,        # Adult female abundance
                          thisNm         # Adult males
    )
    
    a_list[[p]] <- this_amat    # output list length = nPatches
  }
  
  # large mat set up
  rep_cb <- function(mat, times) do.call(cbind, rep(list(mat), times))   # function to cbind repeated number of times
  meta <- rep_cb(mat, nPatches)
  
  rep_rb <- function(mat, times) do.call(rbind, rep(list(mat), times))  # function to repeatedly stack mat vertically with rbind
  meta_mat <- rep_rb(meta, nPatches)   # adds on 1 dim each loop
  
  # filling diags with Amats from list
  for (p in 1:nPatches){
    meta_mat[((p*nDim)-3):(p*nDim),((p*nDim)-3):(p*nDim)] <- a_list[[p]]
  }
  
  # More difficult to id matrix positions once cols are named
   colnames(meta_mat) <- rep(stagenames, nPatches)  # rep times = npatches
   rownames(meta_mat) <- rep(stagenames, nPatches)
  
  # staying probs in relevant matrix entries 
  for (p in 1 : nPatches){   # filling female and male remaining probs 
    meta_mat[,((p*nDim)-2)] <- meta_mat[,((p*nDim)-2)] * (1 - Dmat[p,])  # females remaining in patch
    meta_mat[,(p * nDim)] <- meta_mat[,(p * nDim)] * (1 - Dmat[p,])   # male staying in patch
  }
  
  # base matrix to fill with survival rates of moving indivs if all inidivuals are moving, bmat = umat
  bmat <- mat
  bmat[2,2] <- Umat[2,2]   # adult fem survival 
  bmat[4,4] <- Umat[4,4]   # adult m survival
  
  for (p in 1:nPatches){    # for each patch, work in each col of matrix
    col_lims <- ((p * nDim)-3) :(p * nDim)  # defines col boundaries for each patch 
    
    this_bmat <- bmat
    
    # fertility vals needed
    thisA <- alist
    this_bmat[1,2] <- a_list[[p]][1,2]   
    this_bmat[3,2] <- this_bmat[1,2]
    
    # dividing movement prob across remaining patches - assumes moves equally likely
    this_bmat[,2] <- this_bmat[,2] * (Dmat[p,]/ (nPatches-1))      # adult fem row * movement prob in this patch
    this_bmat[,4] <- this_bmat[,4] * (Dmat[p,]/ (nPatches-1))      # male survival * movement prob 
    
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
}  

# Dmat creation function - can be customised to how many moving individuals we have 
Dmat.create <- function(colNames, nPatches) {
  Dmat <- matrix(0, ncol= length(colNames), nrow = nPatches) # row = number patches
      colnames(Dmat) <- colNames
  
  for (p in 1:nPatches){
    move <- rnorm(length(colNames), mean = 0.5, sd =0.2) # eq to calculate move prob value given vars, for now rnorm
    Dmat[p,] <- move
  }
  return(Dmat)
} 
###CH: I think the movement probability might be too high. I would suggest
###setting the mean to 0.2 or something, and then put a line of code that
###catches any values less than 0 (either set them to 0 or take the absolute
###value)






