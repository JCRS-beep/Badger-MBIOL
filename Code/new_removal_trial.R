# Updated removal function
# 8.7.2026
# 

# Intensive culls - 70% reduction first year, subsequent years..., MINIMUM 4 years (remyears set 5:8)
# Supplementary  culling - only if initial cull successful.  2 year licence. Can be culled up to 5 years
# setting min and max used for supplementary badger culling

# truncated normal dist function to sample from 
# Sample n integers between min and max with a bell-shaped (truncated normal) distribution
# mean = midpoint of min and max 
# sd defaults to (max - min)/4, so min and max sit ~2 SDs from the mean
trunc.norm <- function(n = 1, min, max, mean = (min + max)/2, sd = (max - min)/3) {
  
  min <- max(0, min)       # can never remove fewer than 0 individuals
  max <- max(min, max)     # guards against max < min
  #  mean = (min + max)/2
  #  sd = (max - min)/4
  if (max == min) return(rep(round(min), n))   # no range to sample from
  
  out <- numeric(n)
  for (k in 1:n) {
    x <- rnorm(1, mean, sd)
    # resample until the draw falls in range (true truncated normal)
    while (x < min || x > max) {
      x <- rnorm(1, mean, sd)
    }
    out[k] <- x
  }
  
  round(out)   # integer number of individuals
}


# updating function design - multi removals 
multi.rem <- function(umat,   # MAX SURVIVAL
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
    thisAmat <- apply.DD(params, umat, thisN, DDapply, stagenames,   
                         thisNf,        
                         thisNm)    
    
    # If the projection matrix has negative or NA values, return message, replace with 0, and continue
    if (any(thisAmat < 0, na.rm = TRUE)) {
      warning(paste("Negative values in projection matrix at time step", i,
                    "- setting negative entries to 0 and continuing."))
      thisAmat[thisAmat < 0] <- 0
      thisAmat[is.na(thisAmat)] <- 0
    }
    Vec[(i + 1), ] <- floor(thisAmat %*% Vec[i, ])  # multiply to project
    
    # checking for removal years --------
    if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
      idx <- rem_index[as.character(i)]     # idx = 1: remyear. What removal year are we 
      
      # FIRST YEAR
      if(idx == 1){       # when idx = 1, remove 70% with goal or proportion to generate rem (number of each sex and stage removed)
      # setting removal goal
        real_int <- rnorm(1, mean = intensity/100, sd = 0.05)   # variation in numbers removed # for repeatability, just use intensity?
        
        preN <- Pop[remyear[1]]  # population size pre cull # move until i
        Nrem <- round(preN * real_int)   # goal to remove  - pop size before first rem  
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)

        # SUBSEQUENT YEARS
      }  else if(idx > 1 && idx <= 4){    # if supplementary, change strat to calculate max and min
        estN <- 0.8 * preN - sum(Nremoved)     # ISSUE : estimated pop size- cannot be tracked directly. This leads to negative estN
                                                       # Use estimate from year before culls minus the culled N since. 2/5 multiplier from published methods
        estN[estN < 0] <- 0                                               
        # set min and max removals
        base <- 0.36 * estN                # base using estimated pop size
        min = base - (0.25*estN)          # difference between min and max = 25% of pre-cull total
        max = base 
        # need to prevent minimum falling bellow 0 - removing negative individuals
        min[min < 0] <- 0

        # set seed for repeatability and use a normal dist instead
        Nrem <- trunc.norm(1, min, max) # how many to remove this year using truncated norm dist between min and max values
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
        # SUPPLEMENTARY CULLS - ISSUE this part removes too many
      }  else if(idx > 4){    # if SBC, change strat to calculate max and min. Needs to be much lower intensity?
        base <- 0.2 * Nremoved[1]   # 36% of first year cull total - SBC. CHanging this to 10% because of low population sizes
        # set min and max removals
        min = base - (0.125*preN)          # difference between min and max = 25% of pre-cull total
        max = base + (0.125*preN)       # set maximum as as base to avoid large removals freqs?
        min[min < 0] <- 0
        
        # dist between min and max with mean = baseline?
        # set seed for repeatability and use a normal dist instead
        Nrem <- trunc.norm(1, min, max)  # how many to remove this year
        Nrem <- min(Nrem, Pop[i])    # if Nrem greater than the current pop size, then remove pop size (shows up correct in Nremoved vec)
        
      }
      
      # pop removal -------------------------
      if(rem_strat == "random"){
        if (is.null(bias) == FALSE) paste("ignoring bias value as removal is random across ages and sexes")
        #generating the distribution - varies with rem strat
        dist <- stagedist
        rem <-  round(Nrem * dist)    # where does variation come in?
        
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
      Vec[(i + 1), ] <- floor(thisAmat %*% Vec[i, ]) - thisRemvec  # year after remyear filled with new stage vec
    }
    # section to put only before output ---------------
    # set any negatives to 0
    Vec[i + 1, ][is.na(Vec[i + 1, ])] <- 0  # setting na values to 0

        # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
      Vec[is.na(Vec)] <- 0
    }
    
    Pop[i + 1] <- sum(Vec[(i + 1), ])

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
proj0 <- multi.rem(umat,   # MAX SURVIVAL
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


# random removal at year 10 - 13 
proj1 <- multi.rem(umat,   # MAX SURVIVAL
                   initial = n0, # initial vec
                   stagedist,  # proportion of each stage
                   params, 
                   stagenames = stages, 
                   time = 50,     
                   DDapply="fertility", 
                   intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                   remyear = c(10:13),  # removal year = vector of years 
                   rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                   bias = NULL , # strength of bias as percentage (range??)
                   return.vec= TRUE, 
                   return.remvec = TRUE) 
# year after removal doesn't updte Nf and Nm?
# Doesn not remove smaller intensity after first rem year
(proj1_vec <- dd.plot(proj1,
                       y_val = "Vec",
                       ylab = "Abundance",
                       xlab = "Time (t)",
                       rem_year = 10,
                       mytheme = theme_classic(),
                       cols = col_vec,
                       legend.pos = "top",
                       base_size = 16))
(proj1_pop <- dd.plot(proj1, 
                    y_val= "N", 
                    ylab = "Pop size", 
                    xlab = "Time (t)",
                    rem_year = 10,
                    mytheme = theme_classic(), 
                    cols= col_vec,    # can be vector of cols
                    legend.pos = "topright",
                    base_size = 16))  


# random removal at year 10 - 19 - SBC leads to too many removed - need to reduce removal target 
proj2 <- multi.rem(umat,   # MAX SURVIVAL
                   initial = n0, # initial vec
                   stagedist,  # proportion of each stage
                   params, 
                   stagenames = stages, 
                   time = 50,     
                   DDapply="fertility", 
                   intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                   remyear = c(10:19),  # removal year = vector of years 
                   rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                   bias = NULL , # strength of bias as percentage (range??)
                   return.vec= TRUE, 
                   return.remvec = TRUE) 
# year after removal doesn't updte Nf and Nm?
# Doesn not remove smaller intensity after first rem year
(proj2_vec <- dd.plot(proj2,
                      y_val = "Vec",
                      ylab = "Abundance",
                      xlab = "Time (t)",
                      rem_year = 10,
                      mytheme = theme_classic(),
                      cols = col_vec,
                      legend.pos = "top",
                      base_size = 16))
(proj2_pop <- dd.plot(proj2, 
                      y_val= "N", 
                      ylab = "Pop size", 
                      xlab = "Time (t)",
                      rem_year = 10,
                      mytheme = theme_classic(), 
                      cols= col_vec,    # can be vector of cols
                      legend.pos = "topright",
                      base_size = 16))  