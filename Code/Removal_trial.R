# Removal methods
# 31.12.2025
# 
library(tidyverse)
library(ggplot2)

# vec structure = yf, af, ym, am (25, 10, 25, 10)
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
n0 <- c(25, 10, 25, 10)

# first time using params from extraction 
# Creating my control - population projected over 100 years
Umat <- matrix(0, nrow=4, ncol=4)
rownames(Umat) <- stages
Umat[2,1]<- 0.851 # yearling f survival
Umat[2,2]<- 0.803  # adult f survival
Umat[4,3]<- 0.809   # yearling m survival
Umat[4,4]<- 0.749   # adult m survival

params<- data.frame(fmax= 0.8436,   # F fecundity max (max cubs per adult female) 
                    Sc_max=0.76,   # max cub survival (equal for sexes)
                    b=0.004,       # temp value- must be calculated from provided datasets
                    rep_K= 2.299,          #litter size (K)
                    h= 6   # harem size per male
                    )


# reaches plateau around 20 years - removal at year = 25?

#removal funcs similar to DDproj - WHAT ABOUT BIASED REMOVALS?
rem.proj <- function(Umat,   # MAX SURVIVAL
                     initial, 
                     params, 
                     stagenames, 
                     time, 
                     memberN=NULL,  # which individuals contribute to pop size? (as vec)
                     DDapply="matrix", 
                     Mfunction= "Min",
                     intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                     remyear = NULL,  # removal year = decrease from following year
                     rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                     bias = NULL , # strength of bias as percentage (range??)
                     return.vec= TRUE, 
                     return.remvec = FALSE) 
{
  # Time 
  if (time<=1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  if (rem_strat != "random" && 1+bias > 1/(intensity/100)) stop("Bias value provided leads to p(removal) >1. Re-enter below", 1/intensity, "and repeat")
  
  nStages <- length(stagenames)/2      # how many stages
  
  
  # population size
  if (is.null(memberN)){
    memberN <- 1:length(stagenames)   # NULL = all members contribute to pop size
    
  } else if(length(memberN) >= 1){
    memberN <- c(memberN) 
  }
  
  # Calculating Unions with mating function
  Nf <- n0[nStages]     # pulls ONLY adult female entry from initial vector
  Nm <- n0[2*nStages]      # adult male in vector (final entry)
  
  # mating func gives initial Fmat for first year
  mating.out <- mating.func(params, stagenames, Nf, Nm, Mfunction,  return.mat= TRUE)  # Npairs and fmat given 
  U <- mating.out$U 
  Fmat <- mating.out$Fmat
  
  # Set up the output
  out <- list(pop = vector(), 
              vec = matrix(), 
              mat = matrix(), 
              Nremoved = numeric(), # how to leave blank if no removals?
              remvec = vector())  # including number removed and removals from each stage
  
  Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- rep(NA, (time + 1))       # vector to fill with total pop size each year
  Nremoved <- 0
  remvec <- rep(NA, length(stagenames))
  
  colnames(Vec) <- stagenames   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- n0                   
  Pop[1] <- sum(n0)    
  
  # defining rem year = if no removals occur, set to time
  if(is.numeric(intensity)){   # removal scenarios only
    ry <- as.numeric(remyear)           # shortening name for future use
  } else if(is.null(intensity)){
    ry <- time                   # no removals = run until time
  }
  # Loop = density dependent matrix application for each year
  for (i in 1:ry) {   # normal projection until remyear (if removal at t=5, project until t=5, final entry inserted to row 6) remember vec[5,] holds entries for year=4 (0:time)
    
    # Mating Fmat creation 
    thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
    thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
    
    # apply mating func to calculate pairs
    thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
    
    thisFmat <- thisMating$Fmat
    thisU <- thisMating$U
    
    # ricker density dependence each year
    thisN <- sum(Vec[i,memberN])  # pop sizes sums row i for cols included in N
    thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
    
    # If the projection matrix has any negative values in it, stop iterating but
    # return the projection up until this point.
    if (sum(thisAmat<0)>0){
      warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
      break
    }
    # if pop size <= 0, stop and return
    if(Pop[i]<= 0){
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0)) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
    }
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
    Pop[i + 1] <- sum(Vec[(i + 1), ])
  }
  
  # pop removal
  # ------------------------------------------------
  if (is.numeric(intensity)) {
    if(rem_strat == "random"){
  if (is.null(bias) == FALSE) paste("ignoring bias value since removal is random across ages and sexes")
    #generating the distribution - varies with rem strat
    prop <- rnorm(length(stagenames), mean = intensity/100, sd= intensity/1000)  # 4 samples from dist mean 0.5, sd 0.05
 
  } else if (is.numeric(intensity) && rem_strat %in% c("adults", "Adults", "adult", "Adult", "yearlings", "Yearlings", "yearling", "Yearlings")){   # want to specify age and sex prob  - adult male, yearling fem.. 
    if (is.null(bias) == TRUE) stop("please provide strength of bias as value 0-")
    if (rem_strat %in% c("adults", "Adults", "adult", "Adult")){
    y_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd = intensity/1000) 
    a_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd = intensity/1000)    
    }
    else if (rem_strat %in% c("yearlings", "Yearlings", "yearling", "Yearlings")){
      y_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd = intensity/1000) 
      a_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd = intensity/1000) 
    }
    # col binding so each row represents stage
    bind <- cbind(y_rem, a_rem)
    prop <- c(bind[1,], bind[2,])   # 4 proportions to remove in correct order (yf, af, ym, am)
 
  } else if (is.numeric(intensity) && rem_strat %in% c("females", "Females", "female", "Female", "males", "Males", "male", "Male")){  
    if (is.null(bias) == TRUE) stop("please provide strength of bias as value 0-")
    if(rem_strat %in% c("females", "Females", "female", "Female")){
     f_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd = intensity/1000) # bias applied to females
     m_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd = intensity/1000) 
  
       } else if(rem_strat %in% c("males", "Males", "male", "Male")){
         f_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd = intensity/1000) # bias applied to females
         m_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd = intensity/1000) 
       }
    bind <- cbind(f_rem, m_rem)
    prop <- c(bind[,1], bind[,2])   # 4 proportions to remove in correct order (yf, af, ym, am)
 
  } else if (is.numeric(intensity) && is.numeric(rem_strat)){
    if (is.null(bias) == TRUE) stop("please provide strength of bias as value 0-")
    # how to specify is specific element biased?  numeric vector or integer 1:4, then apply bias to element
    # ex rem strat = 3 (ym)
    bi <- length(rem_strat)  # how many elements provided

    rem <- rnorm((1- bi), mean = ((1+bias)*intensity)/100, sd = intensity/1000)   # unbiased, 1- number of biased samples needed
    rembi <- rnorm(bi, mean = ((1+bias)*intensity)/100, sd = intensity/1000)
    
    # how to order? biased p pos matched to bias stage? 
    prop <- rep(NA, length(stagenames))
    prop[rem_strat] <- rembi
    prop[-rem_strat] <- rem    # all non biased stages, add unbiased probs
      
    }

    
  rem_vec <- Vec[ry+1,] * prop    # rem vec is how many individuals we remove from the population 
    
    Vec[ry + 2,] <- Vec[ry + 1,]  - rem_vec  # year after remyear = 2 rows later filled with new stage vec
    Pop[ry + 2] <- sum(Vec[ry + 2]) # filling in total pop size
    
    # repeat loop after removal year
    for (j in (ry + 2):time) {    # starts filling from 27th row (26th year)
      
      # Mating Fmat creation 
      thisNf <- sum(Vec[j,nStages-1], Vec[j,nStages-1])    # Nf is mid col in Vec matrix
      thisNm <- sum(Vec[j,2*nStages-1],Vec[j,2*nStages])   # Nm 
      
      # apply mating func to calculate pairs
      thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
      
      thisFmat <-thisMating$Fmat
      
      # ricker density dependence each year
      thisN <- sum(Vec[j,memberN])  # pop sizes sums row i for cols included in N
      thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
      
      # If the projection matrix has any negative values in it, stop iterating but
      # return the projection up until this point.
      if (sum(thisAmat<0)>0){
        warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
        break
      }
      # if pop size <= 0, stop and return
      if(Pop[i]<= 0){
        warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
        break
      }
      # if any stage becomes negative, set to zero and continue
      if (any(Vec < 0)) {
        warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
        Vec[Vec < 0] <- 0
      }
      
      Vec[(j + 1), ] <- thisAmat %*% Vec[j, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
      Pop[j + 1] <- sum(Vec[(j + 1), ])
    }
    }
  
  # out objects
  out$pop <- Pop
  if(is.numeric(intensity)){   
    out$Nremoved <- sum(rem_vec)    
    if (isTRUE(return.remvec)) {
      out$remvec <- rem_vec        # returns blank?
    } }
  
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 
  out$mat <- thisAmat
  
  return(out)
}



# testing no removal scenario----
proj0 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="Fmat", 
                  Mfunction= "Min",
                  intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                  remyear = NULL, 
                  rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL , # strength of bias as percentage (range??) how to ignore in function if null?
                  return.vec= TRUE, 
                  return.remvec = FALSE) 

col_vec <- c("#FF6A6A", "#87CEEB")

(proj0_plot <- dd_plot(proj0, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       cex.legend = 0.8))
N0_plot <- dd_plot(proj0, 
                   y_val= "N", 
                   ylab = "Pop size", 
                   xlab = "Time (t)",
                   theme = theme_classic(), 
                   cols= col_vec,    # can be vector of cols
                   legend.pos = "topright",
                   cex.legend = 0.8)  # issues plotting

proj2 <- rem.proj(Umat,   # MAX SURVIVAL
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="Fmat", 
                  Mfunction= "Min",
                  intensity= 70,  # percentage you want REMOVED from pop at time t=ry
                  remyear = 10,  # removal year = decrease from following year
                  rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL , # strength of bias as percentage (range??) how to ignore in function if null?
                  return.vec= TRUE, 
                  return.remvec = TRUE) 

(proj2_plot <- dd_plot(proj2, 
                     y_val= "Vec", 
                     ylab = "Abundance", 
                     xlab = "Time (t)",
                     rem_year= 10,
                     theme = theme_classic(), 
                     cols= col_vec,    # can be vector of cols
                     legend.pos = "topright",
                     cex.legend = 0.8))
# NOTES - population very quickly returns to plateau - some environmental stochasticity needed?
#         Update plot function to draw red dotted line at remyear


proj3 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages, 
                  time = 30, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="Fmat", 
                  Mfunction= "Min",
                  intensity= 90,  # percentage you want REMOVED from pop at time T=ry
                  remyear = 10, 
                  rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL , # strength of bias as percentage (range??) how to ignore in function if null?
                  return.vec= TRUE) 

(proj3_plot <- dd_plot(proj3, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       cex.legend = 0.8))
N3_plot <- dd_plot(proj3, 
                  y_val= "N", 
                  ylab = "Pop size", 
                  xlab = "Time (t)",
                  rem_year = 10,
                  theme = theme_classic(), 
                  cols= col_vec,    # can be vector of cols
                  legend.pos = "topright",
                  cex.legend = 0.8)
# almost works fine - ' Error in if (sum(thisAmat < 0) > 0) { : 
# missing value where TRUE/FALSE needed ' when values below 0?


# lambda and SSD 
lambda0 <- growth.rate(proj0, vis = TRUE)
summary(lambda0$lambda)
lambda2 <- growth.rate(proj2, vis = TRUE, rem_year = 10)
summary(lambda2$lambda)


lambda3 <- growth.rate(proj3, vis = TRUE) # lambda value of 24 following rem year?

ssd0 <- ssd(proj0, vis = TRUE, cols = col_vec)
ssd2 <- ssd(proj2, vis = TRUE, cols = col_vec)
# how to compare? sapply doesn't work on matrix
ssd0_df <- as.data.frame(ssd0$stageMat)
sapply(ssd0_df, mean, 1)

ssd2_df <- as.data.frame(ssd2$stageMat)
sapply(ssd2_df, mean, 1)


# testing biased removals
proj_bi <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages, 
                  time = 30, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="Fmat", 
                  Mfunction= "Min",
                  intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                  remyear = 10, 
                  rem_strat = "female",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = 0.2 , # strength of bias as percentage (range??) how to ignore in function if null?
                  return.vec= TRUE, 
                  return.remvec = TRUE) # not returned?

(bias_plot <- dd_plot(proj_bi, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       cex.legend = 0.8))



# Multiple removals? rem_year 1, rem_year 2 
multi.rem.proj <- function(Umat,   # MAX SURVIVAL
                     initial, 
                     params, 
                     stagenames, 
                     time, 
                     memberN=NULL,  # which individuals contribute to pop size? (as vec)
                     DDapply="matrix", 
                     Mfunction= "Min",
                     intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                     remyears = c(NULL),  # removal year = vector of years to create removals
                     return.vec= TRUE, 
                     return.remvec = FALSE) 
{
  # Time 
  if (time<=1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  
  nStages <- length(stagenames)/2      # how many stages
  
  
  # population size
  if (is.null(memberN)){
    memberN <- 1:length(stagenames)   # NULL = all members contribute to pop size
    
  } else if(length(memberN) >= 1){
    memberN <- c(memberN) 
  }
  
  # Calculating Unions with mating function
  Nf <- n0[nStages]     # pulls ONLY adult female entry from initial vector
  Nm <- n0[2*nStages]      # adult male in vector (final entry)
  
  # mating func gives initial Fmat for first year
  mating.out <- mating.func(params, stagenames, Nf, Nm, Mfunction,  return.mat= TRUE)  # Npairs and fmat given 
  U <- mating.out$U 
  Fmat <- mating.out$Fmat
  
  # Set up the output
  out <- list(pop = vector(), 
              vec = matrix(), 
              mat = matrix(), 
              Nremoved = numeric(), # how to leave blank if no removals?
              remvec = vector())  # including number removed and removals from each stage
  
  Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- rep(NA, (time + 1))       # vector to fill with total pop size each year
  Nremoved <- 0
  remvec <- rep(NA, length(stagenames))
  
  colnames(Vec) <- stagenames   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- n0                   
  Pop[1] <- sum(n0)    
  
  # is multiple removals, create multiple ry
 if(length(rem_year) > 1){   # removal scenarios only
#    for (b in 1:length(rem_year))  # repeat for as many removal years listed
    
      ry <- as.numeric(remyear)       # create individual objects (ry1, ry2..)?
 
     } else if(length(rem_year) = 1){
    ry <- as.numeric(remyear)
  }
  
  # Loop = density dependent matrix application for each year
  for (i in 1:ry[1]) {   # normal projection until first entry of remyear (if removal at t=5, project until t=5, final entry inserted to row 6) remember vec[5,] holds entries for year=4 (0:time)
    
    # Mating Fmat creation 
    thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
    thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
    
    # apply mating func to calculate pairs
    thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
    
    thisFmat <- thisMating$Fmat
    thisU <- thisMating$U
    
    # ricker density dependence each year
    thisN <- sum(Vec[i,memberN])  # pop sizes sums row i for cols included in N
    thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
    
    # If the projection matrix has any negative values in it, stop iterating but
    # return the projection up until this point.
    if (sum(thisAmat<0)>0){
      warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
      break
    }
    # if pop size <= 0, stop and return
    if(Pop[i]<= 0){
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0)) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
    }
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
    Pop[i + 1] <- sum(Vec[(i + 1), ])
  }
  
  # pop removal
  # -------------------------
  
    #generating the distribution
    prop <- rnorm(length(stagenames), mean = intensity/100, sd= intensity/1000)  # 4 samples from dist mean 0.5, sd 0.05
    rem_vec <- Vec[ry[1]+1,] * prop    # remaining <- Vec[ry+1] - rem_vec
    
    Vec[ry[1] + 2,] <- Vec[ry[1] + 1,]  - rem_vec  # year after remyear = 2 rows later filled with new stage vec
    Pop[ry[1] + 2] <- sum(Vec[ry[1] + 2]) # filling in total pop size
    
  # -----------------------  
    # repeat loop after removal year
    # if only one removal, project until time, if other wise, project until ry[2]
    if(length(ry) = 1){
      for (j in (ry + 2):time)
        # Mating Fmat creation 
        thisNf <- sum(Vec[j,nStages-1], Vec[j,nStages-1])    # Nf is mid col in Vec matrix
      thisNm <- sum(Vec[j,2*nStages-1],Vec[j,2*nStages])   # Nm 
      
      # apply mating func to calculate pairs
      thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
      
      thisFmat <-thisMating$Fmat
      
      # ricker density dependence each year
      thisN <- sum(Vec[j,memberN])  # pop sizes sums row i for cols included in N
      thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
      
      # If the projection matrix has any negative values in it, stop iterating but
      # return the projection up until this point.
      if (sum(thisAmat<0)>0){
        warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
        break
      }
      # if pop size <= 0, stop and return
      if(Pop[i]<= 0){
        warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
        break
      }
      # if any stage becomes negative, set to zero and continue
      if (any(Vec < 0)) {
        warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
        Vec[Vec < 0] <- 0
      }
      
      Vec[(j + 1), ] <- thisAmat %*% Vec[j, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
      Pop[j + 1] <- sum(Vec[(j + 1), ])
        
    } else if(length(ry) > 1)
      for (l in 1:length(ry))   # rep for 1:length(ry)  = amat creation ry times
    for (j in (ry[1] + 2):ry[2]) {    # starts filling from 27th row (26th year)
      
      # Mating Fmat creation 
      thisNf <- sum(Vec[j,nStages-1], Vec[j,nStages-1])    # Nf is mid col in Vec matrix
      thisNm <- sum(Vec[j,2*nStages-1],Vec[j,2*nStages])   # Nm 
      
      # apply mating func to calculate pairs
      thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
      
      thisFmat <-thisMating$Fmat
      
      # ricker density dependence each year
      thisN <- sum(Vec[j,memberN])  # pop sizes sums row i for cols included in N
      thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
      
      # If the projection matrix has any negative values in it, stop iterating but
      # return the projection up until this point.
      if (sum(thisAmat<0)>0){
        warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
        break
      }
      # if pop size <= 0, stop and return
      if(Pop[i]<= 0){
        warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
        break
      }
      # if any stage becomes negative, set to zero and continue
      if (any(Vec < 0)) {
        warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
        Vec[Vec < 0] <- 0
      }
      
      Vec[(j + 1), ] <- thisAmat %*% Vec[j, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
      Pop[j + 1] <- sum(Vec[(j + 1), ])
    }
 
  # out objects
  out$pop <- Pop
  if(is.numeric(intensity)){   
    out$Nremoved <- sum(rem_vec)    
    if (isTRUE(return.remvec)) {
      out$remvec <- rem_vec        # returns blank?
    } }
  
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 
  out$mat <- thisAmat
  
  return(out)
}