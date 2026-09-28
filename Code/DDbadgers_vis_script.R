# modified functions to visualise outputs

MOD_meta_fixeddisp_norem<- function(Umat,
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
  Fert <- matrix(NA, nrow = time, ncol = npatches)   # fertility entry used to project year i -> i+1
  DDval <- matrix(NA, nrow = time, ncol = npatches)  # density dependent multiplier year by year (ricker output)
  
  colnames(Vec) <- rep(stagenames, npatches)  # naming cls matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals  
  Pop[1] <- as.numeric(sum(n0))
  
  for(i in 1:npatches){
    x <- (1 +4*(i-1)) : (4*i)   # boundaries containing each patch abundances
    Group[1,i] <- sum(Vec[1, x])     # filling first row of group matrix
  }
  
  #  colnames(Leavers) <- seq(1, npatches)   # patch names 
  # rownames(Leavers) <- 0:(time) 
  
  # Loop = density dependent matrix application for each year and patch
  for (i in 1:time) {   #  projection until end time
    
    leavers <- rep(0, ncol(Vec))    # individuals leaving patches - to go in output#
    arrivers <- rep(0, ncol(Vec))   # individuals arriving in patch p+1 (other patch)
    
    # dispersal loop 
    for (p in 1:npatches){
      Ipatch <- (1+4*(p-1)):(4*p) # indices of the patch rows/cols in Vec and Umat
      subVec <- as.vector(Vec[i, Ipatch])    # have to specify row (year) we want, to remove names format as numeric vector
      
      
      # apply the probability of leaving to the stages that leave:
      subleavers<- rep(0,4)    # leavers from this patch
      subleavers[I_leavers]<- round(dispersal_prob*subVec[I_leavers])   # round instead of floor - does not underestimate moving individuals
      leavers[Ipatch] <- subleavers
      
      if(sum(subleavers) == 0) { # if none leave, can skip all further steps and set arrivers as the same value
        arrivers <- arrivers
      }
      else {
      # distribute the leavers equally among the other patches, as arrivers, in
      # the same stage as they left:
      if (npatches==2){
        subarrivers<- subleavers      
        arrivers[-Ipatch]<- subarrivers     # arrivers into other patch 
        
      } else {     # divide movers across other patches
        # random sample of patches to decide where individuals go (later replace with some distance equation?)
        nleavers <- sum(subleavers) # number of moving individuals
        x <- seq(npatches)
        arrival_patches <- x[x != p]   # all patches excluding current patch
        
        remaining_subleavers <- subleavers   # how many leavers from this patch remain
        mat_subarrival <- matrix(0, nrow = length(arrival_patches), ncol = nStages)
        rownames(mat_subarrival) <- arrival_patches
        
        # loop for each stage in disp stage 
        for(j in 1:length(I_leavers)){    # why is j giong to 4?
          stageI <- I_leavers[j]
          stage_I_arrival <- sample(arrival_patches, subleavers[stageI], replace = TRUE)   # samples patches for the stage j in dispersal stages
          
          for (a in arrival_patches){  # for each possible arrival patch
            thisIpatch <- c(1+4*(a-1)):(4*a)    # first arrival patch index          
            matrow <- which(rownames(mat_subarrival) == a)
            
            if(matrow == 0) {
              arrivers[thisIpatch] <- rep(0, 4)
            }
            
            else {mat_subarrival[matrow,stageI] <- sum(stage_I_arrival == a)   # arrivers for this patch
            
            arrivers[thisIpatch] <- c(mat_subarrival[matrow,])
            }
          }
          }
        }
      }
    }   # DEBUGGING : closing bracket in wrong place - reproduction was calculated before dispersal vec 
    # was completed for all patches
    
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
      
      Fert[i, p] <- subAmat[1, 2]   # storing our fertility value for this group this year (multiply by 2 for 2 sexes)
      DDval[i,p] <- ricker(params, subN)  # storing the ricker value for this group this year
      
      # Project just this patch forward:
      Vec[(i+1), Ipatch]<- floor(subAmat %*% postdispVec[Ipatch]) # amat values multiplied by POST DISP vec at current patch.
                                                            # round down for integers 
      
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
                leavers = matrix(),
                fert = matrix(), 
                DDval = matrix())
    
    out$pop <- Pop
    out$vec <- Vec 
    out$group <- Group
    out$leavers <- Leavers
    out$fert <- Fert
    out$DDval <- DDval
    
    
    return(out)
  } else {
    return(Pop)
  }
}


# running projection to visualise
output <- MOD_meta_fixeddisp_norem(threepatch_Umat, n0_3patch, params, stages, 
                                   npatches=3, DDapply = 'fertility', time = 30,
                                   dispersal_prob=0.1, 
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec = TRUE)



# function linking reproductivve vals with group size (in future Nm and Nf)
plot_patch_matrix <- function(mat,                  # e.g. output$fert (Amat) or output$DDval
                              ylab = "Value",       
                              years = NULL,
                              size_mat = NULL,      # e.g. output$group; if given, x = group size
                              facet = FALSE,        # one panel per group
                              path = TRUE) {        # join points in time order
  require(tidyr); require(ggplot2)
  
  n <- nrow(mat)   # either fert or DDval
  df <- as.data.frame(mat)
  colnames(df) <- paste0("Group", seq_len(ncol(mat)))   # naming group 1 to p
  df$Year <- if (is.null(years)) seq_len(n) - 1 else years    # assigning years
  df_long <- pivot_longer(df, cols = starts_with("Group"),
                          names_to = "Group", values_to = "Value")    
  
  if (is.null(size_mat)) {
    # original behaviour: value against year
    p <- ggplot(df_long, aes(Year, Value, colour = Group)) +
      geom_line(linewidth = 1) + geom_point() +
      labs(x = "Time (t)", y = ylab)
  } else {
    # value against group size
    size_mat <- size_mat[seq_len(n), , drop = FALSE]   # group has time+1 rows, fert has time
    sz <- as.data.frame(size_mat)
    colnames(sz) <- paste0("Group", seq_len(ncol(size_mat)))
    sz$Year <- df$Year
    sz_long <- pivot_longer(sz, cols = starts_with("Group"),
                            names_to = "Group", values_to = "Size")
    
    df_long <- merge(df_long, sz_long, by = c("Year", "Group"))
    df_long <- df_long[order(df_long$Group, df_long$Year), ]   # keep time order for geom_path
    
    p <- ggplot(df_long, aes(Size, Value, colour = Group)) +
      geom_point(size = 2) +
      labs(x = "Group size", y = ylab)
    if (path) p <- p + geom_path(alpha = 0.5)
  }
  
  if (facet) p <- p + facet_wrap(~ Group)
  p + theme_classic(base_size = 16)
}

# plots
plot_patch_matrix(output$fert, "Fertility (Adult_f -> Yearling_f)", years = 1:30)  # fertility vs year
plot_patch_matrix(output$fert, "Fertility", size_mat = output$group)   # fertility vs group size

plot_patch_matrix(output$DDval, "ricker value", size_mat = output$group)   # ricker vs group size


# one panel per group (better with many groups)
plot_patch_matrix(output$fert, "Fertility", size_mat = output$group, facet = TRUE)


# what next? we want the exact same curve shape, but instead of max size 150, max size 30?
