# Density Dependent dispersal update
# Changing my dispersal equation to logistic

# logistic eq 
# p_move<- max_move/(1+exp(-steepness*(x-xhalf)))
# max_move = 0.3
# steepness?
# xhalf = group_size when p_move = 0.15 (half max_move)

x<- seq(0, 50, 1)    # group size 
steepness<- 0.5
max_move<- 0.3
xhalf<- 15
p_move<- max_move/(1+exp(-steepness*(x-xhalf)))
plot(x, p_move, type='l')


# the maximum value is not enforced - want to max out 0.3 but may rise over
# occurs when groups shoots over max size, (x > 1)


# DD moving projection function -----
# incorporating density dependent movement probs in projection function
MOD_meta_DDdisp_norem <- function(Umat,
                                   initial,
                                   params,
                                   stagenames,
                                   npatches,
                                   time,
                                   DDapply="fertility",
                                   max_group,        # vector of single value for carying capacity of each group
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec=TRUE){
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
  
  # Set up  output
  Vec <- matrix(0, ncol = length(stagenames)*npatches, nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Group <- matrix(0, ncol = npatches, nrow = time+1)   # group sizes each year
  Leavers <- matrix(0, ncol = ncol(Vec), nrow = time)   # matrix to fill with the number of individuals leaving each patch each year (by each stage)
  
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
    # group sizes this year
    sizes <- vector()
    # sum every 4 entries
    for (p in 1:npatches){
      sizes[p] <- sum(Vec[i , (1+4*(p-1)):(4*p)])
    }
    
    
    dmat <- ddDmat(stagenames,dispersal_stages, # names of all stages, names of only dispersing stages
                   npatches, # number of patches
                   group_size = sizes, 
                   max_group)
    
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
    
    # now we do another loop over the patches to calculate reproduction:
    for (p in 1:npatches){
      Ipatch<- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      # extract the Umat for this patch:
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
      
      # Project just this patch forward:
      Vec[(i+1), Ipatch]<- floor(subAmat %*% postdispVec[Ipatch]) # amat values multiplied by vec, round down for integers 
      # DEBUGGING:
      # should multiply subAmat by POST DISP vec (not vec) - rep values based on new abundances, multiply rep values by 
      # post mvoement vec
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
    
    return(out)
  } else {
    return(Pop)
  }
}


# test 
test_proj <- MOD_meta_DDdisp_norem(threepatch_Umat, n0_small_3patch, new_params, stages, 
                                              npatches=3, DDapply = 'fertility', time = 50,
                                              max_group = 28,  # based on reading - carrying capacity 
                                              dispersal_stages = c("Adult_f", "Adult_m"),
                                              return.vec=TRUE)

meta_fixedDD_norem_group <- meta_plot(test_proj,   # output obj of dd.proj
                                      y_val= "Group",   # plot type - N or Vec 
                                      ylab = "Group size", 
                                      xlab = "time (t)",
                                      rem_year = NULL,
                                      mytheme = theme_classic(), 
                                      cols= "black",    # can be vector of 2 cols
                                      legend.pos = "top",
                                      base_size = 16)
# looking better!