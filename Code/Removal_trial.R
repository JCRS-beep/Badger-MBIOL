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

# removal function (see all function folder)
rem.proj <- function(Umat,   # MAX SURVIVAL
                     initial, 
                     params,
                     stagenames, # needed for mating func
                     time, 
                     DDapply="fertility", 
                     intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                     remyear = NULL,  # removal year = decrease from following year
                     rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                     bias = NULL , # strength of bias as percentage (range??)
                     return.vec= TRUE, 
                     return.remvec = FALSE) {
  # Time 
  if (time<=1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  if (!is.null(bias) && 1+bias > 1/(intensity/100)) stop("Bias value provided leads to p(removal) >1. Re-enter below", 1/intensity, "and repeat")
  
  
  nStages <- length(stagenames)/2      # how many stages

  # Set up the output
  out <- list(pop = vector(), 
              vec = matrix(), 
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
  
  
  if(is.numeric(intensity)){   # removal scenarios only
    ry <- as.numeric(remyear)           # shortening name for future use
  } else if(is.null(intensity)){
    ry <- time                   # no removals = run until time
  }
  # Loop = density dependent matrix application for each year
  for (i in 1:ry) {   # normal projection until remyear (if removal at t=5, project until t=5, final entry inserted to row 6) remember vec[5,] holds entries for year=4 (0:time)
    
    # f value per year creation
    thisNf <- Vec[i,nStages]    # Nf = adult fems in Vec matrix
    thisNm <- Vec[i,2*nStages]  # Nm 
    
    thisN <- sum(Vec[i,])  # pop sizes sums row i for cols included in N
    
    # ricker density dependence each year
    thisAmat <- apply.DD(params, Umat, thisN, DDapply, stagenames,   
                         thisNf,        
                         thisNm,            
                         Mfunction= "min",       
                         return.mat= FALSE)  
    
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
  
  
  if (is.numeric(intensity)) {
    # pop removal
    # ------------------------------------------------
    if(rem_strat == "random"){
      if (is.null(bias) == FALSE) paste("ignoring bias value since removal is random across ages and sexes")
      #generating the distribution - varies with rem strat
      prop <- rnorm(length(stagenames), mean = intensity/100, sd= 0.05)  # 4 samples from dist mean 0.5, sd 0.05
      
    } else if (is.numeric(intensity) && rem_strat %in% c("adults", "Adults", "adult", "Adult", "yearlings", "Yearlings", "yearling", "Yearlings")){   # want to specify age and sex prob  - adult male, yearling fem.. 
      if (is.null(bias) == TRUE) stop("please provide strength of bias as value 0-")
      if (rem_strat %in% c("adults", "Adults", "adult", "Adult")){
        y_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd = 0.05) 
        a_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd = 0.05)    
      }
      else if (rem_strat %in% c("yearlings", "Yearlings", "yearling", "Yearlings")){
        y_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd = 0.05) 
        a_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd = 0.05) 
      }
      # col binding so each row represents stage
      bind <- cbind(y_rem, a_rem)
      prop <- c(bind[1,], bind[2,])   # 4 proportions to remove in correct order (yf, af, ym, am)
      
    } else if (is.numeric(intensity) && rem_strat %in% c("females", "Females", "female", "Female", "males", "Males", "male", "Male")){  
      if (is.null(bias) == TRUE) stop("please provide strength of bias as value 0-")
      if(rem_strat %in% c("females", "Females", "female", "Female")){
        f_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd = 0.05) # bias applied to females
        m_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd = 0.05) 
        
      } else if(rem_strat %in% c("males", "Males", "male", "Male")){
        f_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd = 0.05) # bias applied to females
        m_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd = 0.05) 
      }
      bind <- cbind(f_rem, m_rem)
      prop <- c(bind[,1], bind[,2])   # 4 proportions to remove in correct order (yf, af, ym, am)
      
    } else if (is.numeric(intensity) && is.numeric(rem_strat)){
      if (is.null(bias) == TRUE) stop("please provide strength of bias as value 0-")
      # how to specify is specific element biased?  numeric vector or integer 1:4, then apply bias to element
      # ex rem_strat = 3 (ym)
      bi <- length(rem_strat)  # how many elements provided
      
      rem <- rnorm((1- bi), mean = ((1+bias)*intensity)/100, sd = 0.05)   # unbiased, 1- number of biased samples needed
      rembi <- rnorm(bi, mean = ((1+bias)*intensity)/100, sd = 0.05)
      
      # how to order? biased p pos matched to bias stage? 
      prop <- rep(NA, length(stagenames))
      prop[rem_strat] <- rembi
      prop[-rem_strat] <- rem    # all non biased stages, add unbiased probs
    }
    # ------------------
    # removals
    rem_vec <- Vec[ry+1,] * prop    # rem vec is how many individuals we remove from the population 
    Rem <- sum(rem_vec)
    
    Vec[ry + 2,] <- Vec[ry + 1,]  - rem_vec  # year after remyear = 2 rows later filled with new stage vec
    Pop[ry + 2] <- sum(Vec[ry + 2]) # filling in total pop size
    
    
    # repeat loop after removal year
    for (j in (ry + 2):time) {    # starts filling from 27th row (26th year)
      
      # Mating Fmat creation 
      thisNf <- sum(Vec[j,nStages-1], Vec[j,nStages-1])    # Nf is mid col in Vec matrix
      thisNm <- sum(Vec[j,2*nStages-1],Vec[j,2*nStages])   # Nm 
      
      # ricker density dependence each year
      thisN <- sum(Vec[j,])  # pop sizes sums row i for cols included in N
      thisAmat <- apply.DD(params, Umat, thisN, DDapply, stagenames,   
                           thisNf,        
                           thisNm,            
                           Mfunction= "min",       
                           return.mat= FALSE)
      
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
    out$Nremoved <- Rem    # where do we define total N removed
    if (isTRUE(return.remvec)) {
      out$remvec <- rem_vec        # returns blank?
    } }
  
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 

  return(out)
}


# ChatGPT updated ddplot code
dd_plot <- function(out,
                    y_val = "N",
                    ylab = "Abundance",
                    xlab = "Time (t)",
                    rem_year = NULL,
                    theme = theme_classic(),
                    cols = "black",
                    legend.pos = "top",
                    base_size = 16) {
  
  require(tidyr)
  require(ggplot2)
  
  # Create time vector
  t <- nrow(out$vec) - 1
  time <- 0:t
  
  # ---------- POPULATION SIZE ----------
  if (y_val %in% c("N", "Pop Size", "pop size", "Pop", "pop")) {
    
    pop_df <- data.frame(
      Year = time,
      Pop  = out$pop
    )
    
    plot <- ggplot(pop_df, aes(x = Year, y = Pop)) +
      geom_line(size = 1.2, colour = "grey30") +
      geom_point(size = 2, colour = "black") +
      labs(title = "Population Size Over Time",
           x = xlab,
           y = ylab) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme +
      theme(
        text = element_text(size = base_size),
        plot.title = element_text(size = base_size + 2, face = "bold"),
        axis.title = element_text(size = base_size),
        axis.text = element_text(size = base_size - 2),
        legend.position = legend.pos
      )
    
    if (!is.null(rem_year)) {
      plot <- plot +
        geom_vline(xintercept = rem_year,
                   colour = "red3",
                   linetype = "dashed",
                   linewidth = 1,
                   alpha = 0.6)
    }
    
    # ---------- STAGE STRUCTURE ----------
  } else if (y_val %in% c("Vec", "Pop Structure", "Stages", "vec")) {
    
    x_val <- ncol(out$vec)
    
    df <- as.data.frame(out$vec)
    df$Year <- time
    
    df_long <- gather(df, key = "Stage", value = "Abundance", 1:x_val)
    df_long <- separate(df_long,
                        col = "Stage",
                        into = c("Stage", "Sex"),
                        sep = "_")
    
    plot <- ggplot(df_long,
                   aes(x = Year,
                       y = Abundance,
                       colour = Sex,
                       linetype = Stage)) +
      geom_line(size = 1.1) +
      geom_point(position= "jitter", alpha=0.8, size = 2.5) + 
      scale_colour_manual(values = cols,
                          labels = c("Female", "Male")) +
      labs(title = "Stage Abundance Over Time",
           x = xlab,
           y = ylab,
           colour = "Sex",
           linetype = "Stage", 
           shape = "Stage") +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme +
      theme(
        text = element_text(size = base_size),
        plot.title = element_text(size = base_size + 2, face = "bold"),
        axis.title = element_text(size = base_size),
        axis.text = element_text(size = base_size - 2),
        legend.position = legend.pos,
        # journal-style compact legend:
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key = element_rect(fill = NA, colour = NA),
        legend.key.size = unit(0.8, "lines"),
        legend.title = element_text(face = "bold", size = base_size),
        legend.text = element_text(size = base_size - 2),
        legend.spacing.x = unit(0.2, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
      ) +
      # make legend compact and show titles above keys (journal style)
      guides(
        colour = guide_legend(title.position = "top",
                              title.hjust = 0.5,
                              nrow = 1,
                              byrow = TRUE,
                              override.aes = list(size = 3, linetype = 1, shape = 16)),
        linetype = guide_legend(title.position = "top",
                                title.hjust = 0.5,
                                nrow = 1,
                                byrow = TRUE,
                                override.aes = list(size = 1.2)),
        shape = guide_legend(title.position = "top",
                             title.hjust = 0.5,
                             nrow = 1,
                             byrow = TRUE,
                             override.aes = list(size = 3))
      )
      
    if (!is.null(rem_year)) {
      plot <- plot +
        geom_vline(xintercept = rem_year,
                   colour = "red3",
                   linetype = "dashed",
                   linewidth = 1,
                   alpha = 0.6)
    }
  }
  
  return(plot)
}
  

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
                       theme = theme_classic(),
                       cols = col_vec,
                       legend.pos = "top",
                       base_size = 16))
(N0_plot <- dd_plot(proj0, 
                   y_val= "N", 
                   ylab = "Pop size", 
                   xlab = "Time (t)",
                   theme = theme_classic(), 
                   cols= col_vec,    # can be vector of cols
                   legend.pos = "topright",
                   base_size = 16))  s

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
                     theme = theme_classic(), 
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
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       base_size = 16))
N3_plot <- dd_plot(proj3, 
                  y_val= "N", 
                  ylab = "Pop size", 
                  xlab = "Time (t)",
                  rem_year = 10,
                  theme = theme_classic(), 
                  cols= col_vec,    # can be vector of cols
                  legend.pos = "topright",
                  base_size = 16)
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
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       base_size = 16))











# NEXT STEPS = Multiple removals rem_year1, 2....
# goal - "remove X% of pop every 2 years for 50 years, long term pop growth rate.
# Syntax = remove at remyear = seq(10,30, by=2) 

# function design - project pop until ry, remove, then ry2, remove, ry 3, remove until reach timestep

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
    for (i in 1:time) {   # normal projection until first entry of remyear - not needed?
      # (if removal at t=5, project until t=5, final entry inserted to row 6 (remember vec[5,] holds entries for year=4 (rows 0:time))
      
      # Nf calculation
      thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
      thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
      
      # ricker density dependence each year
      thisN <- sum(Vec[i,])  # pop sizes sums row i for cols included in N
      thisAmat <- apply.DD(params, Umat, thisN, DDapply, stagenames,   
                           thisNf,        
                           thisNm,            
                           Mfunction= "min",       
                           return.mat= FALSE)    
      
      
      # If the projection matrix has any negative values in it, stop iterating and
      # return projection up until this point.
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

trial.plot <- dd_plot(multi.trial, 
                      y_val= "Vec", 
                      ylab = "Abundance", 
                      xlab = "Time (t)",
                      rem_year = c(5,9,11),  # adding line to year - must adapt to incorporate vec inputs
                      theme = theme_classic(), 
                      cols= col_vec,    # can be vector of cols
                      legend.pos = "topright",
                      base_size = 16)

trialN <- dd_plot(multi.trial, 
                  y_val= "N", 
                  ylab = "Abundance", 
                  xlab = "Time (t)",
                  rem_year = c(5,9,11),  
                  theme = theme_classic(), 
                  cols= col_vec,    # can be vector of cols
                  legend.pos = "topright",
                  base_size = 16)


