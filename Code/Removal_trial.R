# Removal methods
# 31.12.2025
# 
library(tidyverse)
library(ggplot2)

# vec structure = yf, af, ym, am (25, 10, 25, 10)
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
n0 <- c(5, 20, 5, 20)

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
                    b=0.007,       # temp value- must be calculated from provided datasets
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




# NEXT STEPS = Multiple removals rem_year1, 2....
# goal - "remove X% of pop every 2 years for 50 years, long term pop growth rate.
# Syntax = remove at remyear = seq(10,30, by=2) 

# function design - project pop until ry, remove until reach timestep. use this as rem.proj instead of seperating? --------

 # model set ups
 # Single projections = useful in visualisations only. Could combine into rep proj as first year rep?
 # Scenario 1 = 70% removal trial at year 10
 proj1 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                   initial = n0, 
                   params, 
                   stagenames = stages,
                   time = 20, 
                   DDapply="Fmat", 
                   intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                   remyear = 5, 
                   rem_strat = "random" ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                   bias = NULL ,
                   return.vec= TRUE, 
                   return.remvec = TRUE) 
 
 
 
 # scenario 2 - biased male removals
 proj2 <- rem.proj(Umat,      
                   initial = n0, 
                   params, 
                   stagenames = stages,
                   time = 20, 
                   DDapply="Fmat", 
                   intensity= 70,  # percentage you want REMOVED from pop at time=ry
                   remyear = 5, 
                   rem_strat = "males" ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                   bias = 0.15 ,
                   return.vec= TRUE, 
                   return.remvec = TRUE) 
 
 
 # Scenario 3 - biased female adult removals
 proj3 <- rem.proj(Umat,      
                   initial = n0, 
                   params, 
                   stagenames = stages,
                   time = 20, 
                   DDapply="Fmat", 
                   intensity= 70,  # percentage you want REMOVED from pop at time=ry
                   remyear = 5, 
                   rem_strat = 2 ,  #2nd in list = adult fems
                   bias = 0.15,  # bias too strong?
                   return.vec= TRUE, 
                   return.remvec = TRUE) 
 
 
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
                      return.remvec = TRUE) 
 {
   # input checks
   if (time <= 1) stop("Time must be a positive integer")
   else(time <- as.integer(time))
   
   if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
   else(n0 <- as.numeric(initial))
   
   
   if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
   
   if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
   if(rem_strat != "random" && is.null(bias) == TRUE){   # is non random but no bias val provided...
      stop("please provide strength of bias for no random removal strategies")
   }  
   
   if (length(remyear) > 0) {
     # ensure remyear are integers and between 1 - time
      
     remyear <- unique(as.integer(remyear))
     remyear <- remyear[remyear >= 1 & remyear <= time]
     if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
   }
   
   nStages <- length(stagenames)/2      # how many stages
   
   # Set up  output
   out <- list(pop = vector(), 
               vec = matrix(), 
               Nremoved = numeric(length(remyear)), # must be a vector - length = number of removal years
               remvec = vector("list", length(remyear)))  # including number removed and removals from each stage - must be a list - length number of removal years
   
   Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
   Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
   Nremoved <- numeric(length(remyear))
   Remvec <- vector("list", length(remyear))
   
   
   colnames(Vec) <- stagenames   # naming cols matrix as stages 
   rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
   Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals              
   Pop[1] <- as.numeric(sum(n0))
   
   # Map removal year -> index into Remvec / Nremoved:
   rem_index <- if (length(remyear) > 0) setNames(seq_along(remyear), remyear)  # takes length rem year, names assigned as seq of integers 
   else integer(0)   # otherwise remyear = 0
      
   
   # Loop = density dependent matrix application for each year UNTIL first rem year
   for (i in 1:time) {   #  projection until end time
     # if removal at t=5, project until t=5, final entry inserted to row 6 (remember vec[5,] holds entries for year=4 (rows 0:time))
     
     # Nf calculation
     thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
     thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
     
     # ricker density dependence each year
     thisN <- sum(Vec[i,])  # pop sizes sums row i for cols included in N
     thisAmat <- apply.DD(params, Umat, thisN, DDapply, stagenames,   
                          thisNf,        
                          thisNm,            
                          return.mat= FALSE)    
     
     
     # If the projection matrix has negative or NA values, return message, replace with 0, and continue
     
     if (any(thisAmat < 0, na.rm = TRUE)) {
       warning(paste("Negative values in projection matrix at time step", i,
                     "- setting negative entries to 0 and continuing."))
       thisAmat[thisAmat < 0] <- 0
       thisAmat[is.na(thisAmat)] <- 0
     }
     
     
     # multiplying to project
     Vec[(i + 1), ] <- floor(thisAmat %*% Vec[i, ])  # amat values multiplied by vec, round down for integers  
     # set any negatives to 0
     Vec[i + 1, ][Vec[i + 1, ] < 0] <- 0
     Vec[i + 1, ][is.na(Vec[i + 1, ])] <- 0  # setting na values to 0
     
     Pop[i + 1] <- sum(Vec[(i + 1), ])
     
     # if pop size <= 0,  and return message with year
     if(Pop[i]<= 0 || is.na(Pop[i])) {
       warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
       break
     }
     # if any stage becomes negative, set to zero and continue
     if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
       warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
       Vec[Vec < 0] <- 0
       Vec[is.na(Vec)] <- 0
     }
     
     # checking for removal years --------
     if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index
        # setting removal goals
        idx <- rem_index[as.character(i)]  
        nYears <- length(remyear)   # over how many years are removals taking place?  remove year_goal for each value of rem_index
        # setting removal goal
        goal <- round(Pop[remyear[1]] * (intensity/100))   # goal to remove  - pop size before first rem
        year_goal <- round(goal/nYears)    # how many removed per year
        
       # pop removal -------------------------
       if(rem_strat == "random"){
         if (is.null(bias) == FALSE) paste("ignoring bias value as removal is random across ages and sexes")
         #generating the distribution - varies with rem strat
         dist <- stagedist
         rem <-  year_goal * dist    # where does variation come in?
         
         # stage biased
       } else if(rem_strat != "random" && is.numeric(bias) == FALSE){   # for biased rems that are NOT index specific...
       bias_vec <- rep(bias, 4)    # adults add bias, y remove bias
       
       if (is.numeric(intensity) && rem_strat %in% c("adult", "Adult", "yearling", "Yearling")){   # want to specify age and sex prob  - adult male, yearling fem.. 
         if (rem_strat %in% c("adults", "Adults", "adult", "Adult")){
            # how to bias removals for classes? 
           bias_vec *c(-1, 1, -1, 1)   # bias is removed from yearlings, added to adults
           
         }
         else if (rem_strat %in% c("yearlings", "Yearlings", "yearling", "Yearlings")){
            bias_vec *c(1, -1, 1, -1)   # bias is added to yearlings, subtracted from adults
         }
       
         # sex biased
       } else if (is.numeric(intensity) && rem_strat %in% c("females", "Females", "female", "Female", "males", "Males", "male", "Male")){  
         if(rem_strat %in% c("females", "Females", "female", "Female")){
            bias_vec *c(1, 1, -1, -1)   # bias is added to fems, subtracted from males
           
         } else if(rem_strat %in% c("males", "Males", "male", "Male")){
            bias_vec *c(-1, -1, 1, 1)   # bias is subtracted from fems, added to males
         }
          dist <- stagedist + bias_vec    # combining into new removal distribution
          rem <-  year_goal * dist     # number to remove per stage
       }
        
         
         # specified
       } else if (is.numeric(intensity) && is.numeric(rem_strat)){
         
         nbi <- length(rem_strat)  # how many elements provided?
         
         # WARNING - only works if rem_strat(length = 1)
         dist <- stagedist 
         dist[rem_strat] <- dist[rem_strat] + bias   # increasing specified element
         dist[-rem_strat] <- dist[-rem_strat] - bias/3    # bias removed from others must divide by 3
         
         rem <-  year_goal * dist    # stages removed per year is total removed per year * stage props
       }

       # ------------------
       
       # calculating number removed from each stage
       thisRemvec <- round(rem)   # actual removed = rounded abundance * dist
       thisRem <- sum(thisRemvec)  # total number removed 
       
       #  into outputs
       Nremoved[idx] <- thisRem  # nth vector entry 
       Remvec[[idx]] <- thisRemvec
       
       # new population size following removals
       Vec[i + 1 ,] <- Vec[i ,]  - thisRemvec  # year after remyear = 2 rows later filled with new stage vec
       Vec[i + 1, ][Vec[i + 1, ] < 0] <- 0
       Pop[i + 1] <- sum(Vec[i + 1,]) # filling in total pop size
       
       # if pop size <= 0, stop and return
       if(Pop[i]<= 0 || is.na(Pop[i])) {
         stop(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
         break
       }
       # if any stage becomes negative, set to zero and continue
       if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
         warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
         Vec[Vec < 0] <- 0
         Vec[is.na(Vec)] <- 0
       }
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

# FUNCTION OUTLINE - 
# project until rem year, get pop size, multiply by intensity to get removal goal. Year goal by dividing by nYears
# Rem vec is year goal * a certain distribution (random or biased) 
# Remove rem vec, put into next year
# Uncertainty - specific biases? 

single_test <- multi.rem(Umat,   # MAX SURVIVAL
                          initial = n0, # initial vec
                          stagedist = c(0.12, 0.4, 0.12, 0.36),  # proportion of each stage
                          params, 
                          stagenames = stages, 
                          time = 20,     
                          DDapply="fertility", 
                          intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                          remyear = 5,  # removal year = vector of years 
                          rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                          bias = NULL , # strength of bias as percentage (range??)
                          return.vec= TRUE, 
                          return.remvec = TRUE) 

double_test <- multi.rem(Umat,   # MAX SURVIVAL
                         initial = n0, # initial vec
                         stagedist = c(0.12, 0.4, 0.12, 0.36),  # proportion of each stage
                         params, 
                         stagenames = stages, 
                         time = 20,     
                         DDapply="fertility", 
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear = c(5,6),  # removal year = vector of years 
                         rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                         bias = NULL , # strength of bias as percentage (range??)
                         return.vec= TRUE, 
                         return.remvec = TRUE) 

multi_test <- multi.rem(Umat,   # MAX SURVIVAL
                         initial = n0, # initial vec
                         stagedist = c(0.12, 0.4, 0.12, 0.36),  # proportion of each stage
                         params, 
                         stagenames = stages, 
                         time = 20,     
                         DDapply="fertility", 
                         intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                         remyear = c(5,6, 7, 8, 9),  # removal year = vector of years 
                         rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                         bias = NULL , # strength of bias as percentage (range??)
                         return.vec= TRUE, 
                         return.remvec = TRUE) 

bias.test <-  multi.rem(Umat,   # MAX SURVIVAL
                        initial = n0, # initial vec
                        stagedist = c(0.12, 0.4, 0.12, 0.36),  # proportion of each stage
                        params, 
                        stagenames = stages, 
                        time = 20,     
                        DDapply="fertility", 
                        intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                        remyear = c(5,6,7,8,9),  # removal year = vector of years 
                        rem_strat = 2,  # if specified removals, "adult, females, yearling, males, yearling females, 
                        bias = 0.03 , # strength of bias as percentage (range??)
                        return.vec= TRUE, 
                        return.remvec = TRUE) 


# old rem.proj function
rem.proj <- function(Umat,   # MAX SURVIVAL
                     initial, 
                     params, 
                     stagenames, 
                     time, 
                     DDapply="fertility", 
                     intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                     remyear = integer(0),  # removal year = vector of years 
                     rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                     bias = NULL , # strength of bias as percentage (range??)
                     return.vec= TRUE, 
                     return.remvec = TRUE) 
{
   # input checks
   if (time <= 1) stop("Time must be a positive integer")
   else(time <- as.integer(time))
   
   if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
   else(n0 <- as.numeric(initial))
   
   
   if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
   
   if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
   
   # if(rem_strat != "random" && is.null(bias)){         # not defined properly - runs err
   #  stop("please provide strength of bias for no random removal strategies")
   
   #}
   
   nStages <- length(stagenames)/2      # how many stages
   ry <- remyear                 
   
   # ensure remyear are integers and between 1 - time
   if (length(remyear) > 0) {
      remyear <- unique(as.integer(remyear))
      remyear <- remyear[remyear >= 1 & remyear <= time]
      if (length(remyear) == 0) warning("No valid removal year in 1:time after filtering")
   }
   
   # Set up  output
   out <- list(pop = vector(), 
               vec = matrix(), 
               Nremoved = numeric(length(remyear)), # must be a vector - length = number of removal years
               remvec = vector("list", length(remyear)))  # including number removed and removals from each stage - must be a list - length number of removal years
   
   Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
   Pop <- numeric(length= (time + 1))       # vector to fill with total pop size each year
   Nremoved <- numeric(length(ry))
   Remvec <- vector("list", length(ry))
   
   
   colnames(Vec) <- stagenames   # naming cols matrix as stages 
   rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
   Vec[1, ] <- floor(n0)     # makes sure this is as an integer, no decimals              
   Pop[1] <- as.numeric(sum(n0))
   
   # Map removal year -> index into Remvec / Nremoved:
   rem_index <- if (length(remyear) > 0) setNames(seq_along(remyear), remyear)  # takes length rem year, names assigned as seq of integers 
   else integer(0)   # otherwise remyear = 0
   
   
   # Loop = density dependent matrix application for each year UNTIL first rem year
   for (i in 1:time) {   #  projection until end time
      # if removal at t=5, project until t=5, final entry inserted to row 6 (remember vec[5,] holds entries for year=4 (rows 0:time))
      
      # Nf calculation
      thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
      thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
      
      # ricker density dependence each year
      thisN <- sum(Vec[i,])  # pop sizes sums row i for cols included in N
      thisAmat <- apply.DD(params, Umat, thisN, DDapply, stagenames,   
                           thisNf,        
                           thisNm,            
                           return.mat= FALSE)    
      
      
      # If the projection matrix has negative or NA values, return message, replace with 0, and continue
      
      if (any(thisAmat < 0, na.rm = TRUE)) {
         warning(paste("Negative values in projection matrix at time step", i,
                       "- setting negative entries to 0 and continuing."))
         thisAmat[thisAmat < 0] <- 0
         thisAmat[is.na(thisAmat)] <- 0
      }
      
      
      # multiplying to project
      Vec[(i + 1), ] <- floor(thisAmat %*% Vec[i, ])  # amat values multiplied by vec, round down for integers  
      # set any negatives to 0
      Vec[i + 1, ][Vec[i + 1, ] < 0] <- 0
      Vec[i + 1, ][is.na(Vec[i + 1, ])] <- 0  # setting na values to 0
      
      Pop[i + 1] <- sum(Vec[(i + 1), ])
      
      # if pop size <= 0,  and return message with year
      if(Pop[i]<= 0 || is.na(Pop[i])) {
         warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
         break
      }
      # if any stage becomes negative, set to zero and continue
      if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
         warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
         Vec[Vec < 0] <- 0
         Vec[is.na(Vec)] <- 0
      }
      
      if (as.character(i) %in% names(rem_index)) {  # if year i is present in remyear index, remove prob and add to following year
         idx <- rem_index[as.character(i)]  
         
         # pop removal -------------------------
         
         if(rem_strat == "random"){
            if (is.null(bias) == FALSE) paste("ignoring bias value since removal is random across ages and sexes")
            #generating the distribution - varies with rem strat
            thisProp <- rnorm(length(stagenames), mean = intensity/100, sd= 0.05)  # 4 samples from dist mean 0.5, sd 0.2
            thisProp <- pmax(0, pmin(1, thisProp))  # clip to [0,1]       
            
            # stage biased
         } else if (is.numeric(intensity) && rem_strat %in% c("adults", "Adults", "adult", "Adult", "yearlings", "Yearlings", "yearling", "Yearlings")){   # want to specify age and sex prob  - adult male, yearling fem.. 
            if (rem_strat %in% c("adults", "Adults", "adult", "Adult")){
               y_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd= 0.05) 
               a_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd= 0.05)    
            }
            else if (rem_strat %in% c("yearlings", "Yearlings", "yearling", "Yearlings")){
               y_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd= 0.05) 
               a_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd= 0.05) 
            }
            # col binding so each row represents stage
            bind <- cbind(y_rem, a_rem)
            thisProp <- c(bind[1,], bind[2,])   # 4 proportions to remove in correct order (yf, af, ym, am)
            
            
            # sex biased
         } else if (is.numeric(intensity) && rem_strat %in% c("females", "Females", "female", "Female", "males", "Males", "male", "Male")){  
            if(rem_strat %in% c("females", "Females", "female", "Female")){
               f_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd= 0.05) # bias applied to females
               m_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd= 0.05) 
               
            } else if(rem_strat %in% c("males", "Males", "male", "Male")){
               f_rem <- rnorm(nStages, mean = ((1-bias)*intensity)/100, sd= 0.05) # bias applied to females
               m_rem <- rnorm(nStages, mean = ((1+bias)*intensity)/100, sd= 0.05) 
            }
            bind <- cbind(f_rem, m_rem)
            thisProp <- c(bind[,1], bind[,2])   # 4 proportions to remove in correct order (yf, af, ym, am)
            
            
            # specified   - THIS MATH INCORRECT - leads to lower intensity than other strats
         } else if (is.numeric(intensity) && is.numeric(rem_strat)){
            # how to specify is specific element biased?  numeric vector or integer 1:4, then apply bias to element
            bi <- length(rem_strat)  # how many elements provided
            
            rem <- rnorm((length(stagenames) - bi), mean = ((1 - bias)*intensity)/100, sd= 0.05)   # unbiased, 1- number of biased samples needed
            rembi <- rnorm(bi, mean = ((1 + bias)*intensity)/100, sd= 0.05)
            
            # how to order? biased p pos matched to bias stage? 
            thisProp <- numeric(length = length(stagenames))
            thisProp[rem_strat] <- rembi   # biased value entry gets biased proportion - only works if length = 1?
            thisProp[-rem_strat] <- rem    # all non biased stages, add unbiased probs
         }
         # ------------------
         
         # calculating number removed from each stage
         thisRemvec <- floor(Vec[i,] * thisProp)   # round removals down or up? If we remove 13.6 badgers, 13 or 14?
         thisRem <- sum(thisRemvec)  # total number removed 
         
         #  into outputs
         Nremoved[idx] <- thisRem  # nth vector entry 
         Remvec[[idx]] <- thisRemvec
         
         # new population size following removals
         Vec[i + 1 ,] <- Vec[i ,]  - thisRemvec  # year after remyear = 2 rows later filled with new stage vec
         Vec[i + 1, ][Vec[i + 1, ] < 0] <- 0
         Pop[i + 1] <- sum(Vec[i + 1,]) # filling in total pop size
         
         # if pop size <= 0, stop and return
         if(Pop[i]<= 0 || is.na(Pop[i])) {
            stop(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
            break
         }
         # if any stage becomes negative, set to zero and continue
         if (any(Vec < 0, na.rm = TRUE) || any(is.na(Vec))) {
            warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
            Vec[Vec < 0] <- 0
            Vec[is.na(Vec)] <- 0
         }
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

