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
  

# testing no removal scenario----
proj0 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  DDapply="Fmat", 
                  intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                  remyear = NULL, 
                  rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL , # strength of bias as percentage (range??) how to ignore in function if null?
                  return.vec= TRUE, 
                  return.remvec = FALSE) 

col_vec <- c("#FF6A6A", "#87CEEB")

(proj0_plot <- dd_plot(proj0,
                       y_val = "Vec",
                       ylab = "Abundance",
                       xlab = "Time (t)",
                       rem_year = NULL,
                       mytheme = theme_classic(),
                       cols = col_vec,
                       legend.pos = "top",
                       base_size = 16))
(N0_plot <- dd_plot(proj0, 
                   y_val= "N", 
                   ylab = "Pop size", 
                   xlab = "Time (t)",
                   mytheme = theme_classic(), 
                   cols= col_vec,    # can be vector of cols
                   legend.pos = "topright",
                   base_size = 16))  

proj2 <- rem.proj(Umat,   # MAX SURVIVAL
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  DDapply="Fmat", 
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
                     mytheme = theme_classic(), 
                     cols= col_vec,    # can be vector of cols
                     legend.pos = "topright",
                     base_size = 16))


proj3 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages, 
                  time = 30, 
                  DDapply="Fmat", 
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
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16))
N3_plot <- dd_plot(proj3, 
                  y_val= "N", 
                  ylab = "Pop size", 
                  xlab = "Time (t)",
                  rem_year = 10,
                  mytheme = theme_classic(), 
                  cols= col_vec,    # can be vector of cols
                  legend.pos = "top",
                  base_size = 16)



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
                  DDapply="Fmat", 
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
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       base_size = 16))











# NEXT STEPS = Multiple removals rem_year1, 2....
# goal - "remove X% of pop every 2 years for 50 years, long term pop growth rate.
# Syntax = remove at remyear = seq(10,30, by=2) 

# function design - project pop until ry, remove until reach timestep. use this as rem.proj instead of seperating?

multi.rem <- function(Umat,   # MAX SURVIVAL
                      initial, 
                      params, 
                      stagenames, 
                      time, 
                      DDapply="matrix", 
                      intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                      remyears = integer(0),  # removal year = vector of years 
                      return.vec= TRUE, 
                      return.remvec = TRUE) 
{
  # input checks
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
  
  
  nStages <- length(stagenames)/2      # how many stages
  ry <- remyears                 
  
  # ensure remyears are integers and between 1 - time
  if (length(remyears) > 0) {
    remyears <- unique(as.integer(remyears))
    remyears <- remyears[remyears >= 1 & remyears <= time]
    if (length(remyears) == 0) warning("No valid remyears in 1:time after filtering")
  }
  
  
  # Set up  output
  out <- list(pop = vector(), 
              vec = matrix(), 
              Nremoved = numeric(length(remyears)), # must be a vector - length = number of removal years
              remvec = vector("list", length(remyears)))  # including number removed and removals from each stage - must be a list - length number of removal years
  
  Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- rep(NA, (time + 1))       # vector to fill with total pop size each year
  Nremoved <- numeric(length(ry))
  Remvec <- vector("list", length(ry))
  
  
  colnames(Vec) <- stagenames   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- n0                   
  Pop[1] <- sum(n0)    
  
  # Map removal year -> index into Remvec / Nremoved:
  rem_index <- if (length(remyears) > 0) setNames(seq_along(remyears), remyears) # index contains number of rem years, links each remyear to a value in vector 
  else integer(0)
  
  
  # Loop = density dependent matrix application for each year UNTIL first rem year
    for (i in 1:time) {   # normal projection 
      # Nf calculation
      thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
      thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
      
      # ricker density dependence each year
      thisN <- sum(Vec[i,])  # pop sizes sums row i for cols included in N
      thisAmat <- apply.DD(params, Umat, thisN, DDapply, stagenames,   
                           thisNf,        
                           thisNm,            
                           return.mat= FALSE)    
      
      
      # If the projection matrix has any negative values in it, stop iterating and
      # return projection up until this point.
      if (sum(thisAmat<0, na.rm = TRUE)>0){
        warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
        break
      }
      # if pop size <= 0, stop and return
      if(Pop[i]<= 0, na.rm = TRUE){
        warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
        break
      }
      # if any stage becomes negative, set to zero and continue
      if (any(Vec < 0, na.rm = TRUE)) {
        warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
        Vec[Vec < 0] <- 0
      }
      
      # multiplying to project
      Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]  
      # set any negatives to 0
      Vec[i + 1, ][Vec[i + 1, ] < 0] <- 0
      
      Pop[i + 1] <- sum(Vec[(i + 1), ])
    
  if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index, remove prob and add to following year
    idx <- rem_index[as.character(i)]  
    
    # pop removal -------------------------
    # haven't added bias yet - next once running

    thisProp <- rnorm(length(stagenames), mean = intensity/100, sd= 0.05)  # 4 samples from dist mean 0.5, sd 0.05
    thisProp <- pmax(0, pmin(1, thisProp))  # clip to [0,1]
    
    thisRemvec <- Vec[i,] * thisProp   
    thisRem <- sum(thisRemvec)  # total number removed 
    
    #  into outputs
    Nremoved[idx] <- thisRem  # nth vector entry 
    Remvec[[idx]] <- thisRemvec
    
    # new population size following removals
    Vec[i + 1 ,] <- Vec[i ,]  - thisRemvec  # year after remyear = 2 rows later filled with new stage vec
    Vec[i + 1, ][Vec[i + 1, ] < 0] <- 0
    Pop[i + 1] <- sum(Vec[i + 1,]) # filling in total pop size
      }
    }
  # -----------------------  
  
  # out objects
  out$pop <- Pop
  out$Nremoved <- Nremoved    
    if (isTRUE(return.remvec)) {
      out$remvec <- Remvec        # returns blank?
    } 
  
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 
    return(out)
}

multi.trial <- multi.rem(Umat,   # MAX SURVIVAL
                         initial = n0, 
                         params, 
                         stagenames = stages, 
                         time = 20, 
                         DDapply="matrix", 
                         intensity= 50,  # percentage you want REMOVED from pop at time T=ry
                         remyears = c(5, 9, 11),  # removal year = vector of years 
                         return.vec= TRUE, 
                         return.remvec = TRUE) 

(trial.plot <- dd_plot(multi.trial, 
                      y_val= "Vec", 
                      ylab = "Abundance", 
                      xlab = "Time (t)",
                      rem_year = c(5,9,11),  # adding line to year - must adapt to incorporate vec inputs
                      mytheme = theme_classic(), 
                      cols= col_vec,    # can be vector of cols
                      legend.pos = "topright",
                      base_size = 16))

trialN <- dd_plot(multi.trial, 
                  y_val= "N", 
                  ylab = "Abundance", 
                  xlab = "Time (t)",
                  rem_year = c(5,9,11),  
                  mytheme = theme_classic(), 
                  cols= col_vec,    # can be vector of cols
                  legend.pos = "topright",
                  base_size = 16)

