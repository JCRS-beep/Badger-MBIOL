# Updated removal function
# 8.7.2026
# 
library(tidyverse)
library(ggplot2)

# NEXT STEPS = what to do for low level continuous killings?

# updating function design - multi removals 
multi.rem <- function(Umat,   # MAX SURVIVAL
                      initial, # initial vec
                      stagedist,  # proportion of each stage
                      params, 
                      stagenames, 
                      time,     
                      DDapply="fertility", 
                      intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                      remyear = integer(0),  # removal year = vector of years 
                      rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                      bias = NULL , # strength of bias as percentage (range??)
                      return.vec= TRUE, 
                      return.remvec = FALSE) 
{
  # input checks
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  
  if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
  # if(rem_strat != c("random", NULL) && is.null(bias) == TRUE){   # is non random but no bias val provided... this leads to err when both null!
  #  stop("please provide strength of bias for no random removal strategies")
  #  }  
  
  if (length(remyear) > 0) {
    # ensure remyear are integers and between 1 - time
    
    remyear <- unique(as.integer(remyear))
    remyear <- remyear[remyear >= 1 & remyear <= time]
    if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
  }
  
  nStages <- length(stagenames)/2      # how many stages
  nYears <- length(remyear)   # over how many years are removals taking place?  remove year_goal for each value of rem_index
  
  # Set up  output
  out <- list(pop = vector(), 
              vec = matrix(), 
              Nremoved = numeric(length(remyear)), # must be a vector - length = number of removal years
              remvec = matrix())  # including number removed and removals from each stage -matrix nrow = number of removal years
  
  Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
  Nremoved <- numeric(length(remyear))
  Remvec <- matrix(0, nrow = length(remyear), ncol =4)  
  
  
  colnames(Vec) <- stagenames   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals              
  Pop[1] <- as.numeric(sum(n0))
  
  # Map removal year -> index into Remvec / Nremoved:
  rem_index <- if (length(remyear) > 0) setNames(seq_along(remyear), remyear)  # takes length rem year, names assigned as seq of integers 
  else integer(0)   # otherwise remyear = 0
  
  
  # Loop = density dependent matrix application for each year UNTIL first rem year
  for (i in 1:time) {   #  projection until end time
    # if removal at t=5, project until t=5, final entry inserted to row 6 (vec[5,] holds entries for year=4 (rows 0:time))
    
    # Nf calculation
    thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
    thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
    
    # ricker density dependence each year
    thisN <- sum(Vec[i,])  # pop sizes sums row i for cols included in N
    thisAmat <- apply.DD(params, Umat, thisN, DDapply, stagenames,   
                         thisNf,        
                         thisNm)    
    
    
    # If the projection matrix has negative or NA values, return message, replace with 0, and continue
    
    if (any(thisAmat < 0, na.rm = TRUE)) {
      warning(paste("Negative values in projection matrix at time step", i,
                    "- setting negative entries to 0 and continuing."))
      thisAmat[thisAmat < 0] <- 0
      thisAmat[is.na(thisAmat)] <- 0
    }
    
    
    # multiplying to project
    Vec[(i + 1), ] <- floor(thisAmat %*% Vec[i, ])  # amat values multiplied by vec, round down for integers  
    
    # checking for removal years --------
    if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
      idx <- rem_index[as.character(i)]  
      # target intensity to remove has variation - can be above or below what is really achieved
      real_int <- rnorm(1, mean = intensity/100, sd = 0.05)
      preN <- Pop[remyear[1]]  # population size pre cull
      
      if(idx ==1){       # when idx = 1, remove 70% with goal or proportion to generate rem (number of each sex and stage removed)
      # setting removal goal
      Nrem <- round(preN * real_int)   # goal to remove  - pop size before first rem  # setting removal goal
      Nrem <- round(preN * real_int)   # goal to remove  - pop size before first rem  
     
      }  else if(idx >=1){    # if a subsequent year, supplementary culls
        base <- 0.36*Nremoved[1]   # 36% of first year cull total
        # set min and max removals
        min = base - (0.125*preN)          # difference between min and max = 25% of pre-cull total
        max = base - (0.125*preN)
        
        # dist between min and max with mean = baseline?
        x <- floor(min:max)  # sample using a sequence
        NRem <- sample(x, 1) # how many to remove this year
      }
      
      # pop removal -------------------------
      if(rem_strat == "random"){
        if (is.null(bias) == FALSE) paste("ignoring bias value as removal is random across ages and sexes")
        #generating the distribution - varies with rem strat
        dist <- stagedist
        rem <-  Nrem * dist    # where does variation come in?
        
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
          rem <-  Nrem * dist     # number to remove per stage
        }
        
        
        # specified
      } else if (is.numeric(intensity) && is.numeric(rem_strat)){
        
        nbi <- length(rem_strat)  # how many elements provided?
        
        # WARNING - only works if rem_strat(length = 1)
        dist <- stagedist 
        dist[rem_strat] <- dist[rem_strat] + bias   # increasing specified element
        dist[-rem_strat] <- dist[-rem_strat] - bias/3    # bias removed from others must divide by 3
        
        rem <-  Nrem * dist    # stages removed per year is total removed per year * stage props
      }
      
      # calculating number removed from each stage
      thisRemvec <- round(rem)   # actual removed = rounded abundance * dist
      thisRem <- sum(thisRemvec)  # total number removed 
      
      #  into outputs
      Nremoved[idx] <- thisRem  # nth vector entry 
      Remvec[idx,] <- thisRemvec
      
      # new population size following removals
      Vec[i + 1 ,] <- Vec[i ,]  - thisRemvec  # year after remyear filled with new stage vec
      Vec[i + 1, ][Vec[i + 1, ] < 0] <- 0
      Pop[i + 1] <- sum(Vec[i + 1,]) # filling in total pop size
    }
    # section to put only before output ---------------
    # set any negatives to 0
    Vec[i + 1, ][Vec[i + 1, ] < 0] <- 0
    Vec[i + 1, ][is.na(Vec[i + 1, ])] <- 0  # setting na values to 0
    Pop[i + 1] <- sum(Vec[(i + 1), ])
    
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
      Vec[is.na(Vec)] <- 0
    }
    
    # if pop size <= 0,  and return message with year
    if(Pop[i]<= 0 || is.na(Pop[i])) {
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
  }
  
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





# testing no removal scenario----
proj0 <- multi.rem(Umat,   # MAX SURVIVAL
                   initial = n0, # initial vec
                   stagedist,  # proportion of each stage
                   params, 
                   stagenames = stages, 
                   time = 25,     
                   DDapply="fertility", 
                   intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                   remyear = integer(0),  # removal year = vector of years 
                   rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                   bias = NULL , # strength of bias as percentage (range??)
                   return.vec= TRUE, 
                   return.remvec = FALSE) 

col_vec <- c("#FF6A6A", "#87CEEB")

(proj0_plot <- dd.plot(proj0,
                       y_val = "Vec",
                       ylab = "Abundance",
                       xlab = "Time (t)",
                       rem_year = NULL,
                       mytheme = theme_classic(),
                       cols = col_vec,
                       legend.pos = "top",
                       base_size = 16))
(N0_plot <- dd.plot(proj0, 
                   y_val= "N", 
                   ylab = "Pop size", 
                   xlab = "Time (t)",
                   mytheme = theme_classic(), 
                   cols= col_vec,    # can be vector of cols
                   legend.pos = "topright",
                   base_size = 16))  


# random removal
proj1 <- multi.rem(Umat,   # MAX SURVIVAL
                   initial = n0, # initial vec
                   stagedist,  # proportion of each stage
                   params, 
                   stagenames = stages, 
                   time = 25,     
                   DDapply="fertility", 
                   intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                   remyear = 10,  # removal year = vector of years 
                   rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                   bias = NULL , # strength of bias as percentage (range??)
                   return.vec= TRUE, 
                   return.remvec = TRUE) 
# year after removal doesn't updte Nf and Nm?