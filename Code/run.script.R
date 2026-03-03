# RUN script
# 23/02/25
# Can be run enturely in order

# loading required packages
library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(readr)

# data extraction
mcdonald_demo <- as.data.frame(read.csv("Data/mcdonald_2016_supinfo.csv"))  # posterior estimates from IPM

dens_posterior <- mcdonald_demo$pop_size   # posterior estimate from model
dens_posterior_mean <- mean(dens_posterior)
dens_posterior_sd <- sd(dens_posterior)

beta_reported <- 0.239  # from paper 

beta <- beta_reported/dens_posterior_sd  # explain maths in appendix!


# required functions ------
# Mating systems function ----  when Nm or Nf = 0 U = 0
mating.func <- function(params,     # density dependent parameters
                        stagenames,   # Stages in life cycle graph 
                        Nf,        # Adult and yearling females
                        Nm,             # Adult and yearling males
                        Mfunction= "min",       # mating function applied
                        return.mat= FALSE) {       # Fmat output?
  
  if(is.null(stagenames) || length(stagenames)== 0 & return.mat==TRUE) stop("stagenames must be provided for correct matrix dimension calculations")
  
  # naming objects 
  K <- params$rep_K      # MAX litter size
  h <- params$h          # harem size
  S<- params$Sc_max    # assumes different survival for male and female cubs
  
  
  # creating the output obj 
  out<- list(f=numeric(), Fmat=matrix())
  f <- NA    # blank 
  Fmat <- matrix(0, ncol= length(stagenames), nrow= length(stagenames)) # blank matrix
  rownames(Fmat) <- stagenames
  colnames(Fmat) <- stagenames
  
  # Harmonic mean mating function
  if(Mfunction %in% c( "HM", "HMMF", "Harmonic Mean Mating Function", "harmonic mean mating function", "Harmonic Mean", "harmonic mean")){
    U <- (2* Nf * Nm*h)/(Nf + Nm* h) # females in unions given pop structure
    f <- K*Nm/ (Nm + Nf*h^-1)  # fertility coefficient f= cubs produced by adult female
    
    # Minimum mating function
  } else if(Mfunction %in% c("Minimum Mating Function", "minimum mating function", "Minimum", "minimum", "min", "Min")){
    # defining number of pairs formed (U)
    if(Nf < Nm*h){
      U <- Nf
    } else if(Nm*h < Nf){
      U <- Nm*h
    } 
    f<- (K*U)/Nf   # fertility coefficient f= cubs produced by adult female
    
    # Mod Harmonic mean mating function
  } else if(Mfunction %in% c("Modified Hamonic Mean Mating Function", "modified harmonic mean mating function", "ModHarmonic", "modharmonic"))   {
    if(Nf < (2* Nf*Nm*h) / (Nf +Nm*h)) {
      U <- Nf
    } else if((2* Nf*Nm*h) / (Nf +Nm*h) < Nf){
      U <- (2* Nf*Nm*h) / (Nf+ Nm)
    }
    f <- (K*U)/Nf    # fertility coefficient f= cubs produced by adult female
  }
  
  
  if(return.mat==TRUE) {
    nStages <- length(stagenames)/2 # number stages for each sex
    Fmat[1,nStages] <- 0.5* f* S    # female cub production by females
    Fmat[nStages+1, nStages] <- 0.5* f* S  # male cub production by females
    out$Fmat <- Fmat
  }
  
  out$f <- f
  
  return(out) 
}


#  Ricker function and matrix application - must return functional Amat for proj
apply.DD <- function(params, 
                     Umat, 
                     N,   # yearling and adults
                     DDapply="fertility", 
                     stagenames,   # Stages in life cycle graph 
                     Nf,        # Adult female abundance
                     Nm,             # Adult males
                     Mfunction= "min",       # mating function applied
                     return.mat= FALSE) {   # apply ricker to whole matrix, survival or fertility
  
  ricker <- function(params, N){   # Using params (b) and population size input
    dd_fun <- exp(-params$b*N)
    return(dd_fun)      # returns value of multiplier
  }
  rick <- ricker(params, N)    # ricker function embedded - only this part needed?
  # Function ricker
  # params, incl b = strength of density dependence
  # N= population size (can be total, NAdults, other, but explain in comments) 
  # Use: returns dd_fun multiplier.
  
  # mating func embedded
  mating <- mating.func(params,     # density dependent parameters
              stagenames,   # Stages in life cycle graph 
              Nf,        # Adult and yearling females
              Nm,             # Adult and yearling males
              Mfunction= "min",       # mating function applied
              return.mat= FALSE) 
  
  f <- mating$f 
  
  # constructing Amat
  Amat <- Umat
  S <- params$Sc_max
  Amat[1,2] <- 0.5* f* S 
  Amat[3,2] <- 0.5* f* S 
  
  # applying based on method
  if(DDapply %in% c("matrix", "Matrix"))  {
    Amat_N <- Amat*rick  
    
  } else if(DDapply %in% c("Survival", "survival", "Umat")){
    Umat_N <- Umat*rick
    Amat_N <- Amat + Umat_N # how to combine?
    
  } else if(DDapply  %in% c("Fertility", "fertility", "Fmat")) {
    # embedd mating func?
    f_N <- f*rick
    Amat_N <- Amat
    Amat_N[1,2] <- 0.5 * f_N * S 
    Amat_N[3,2] <- 0.5 * f_N * S 
  }
  
  return(Amat_N) 
}
# Function name= applyDD
# Inputs:  
#  params: incl b for ricker
#  N: effective pop size (can be total, NAdults, other, but explain in comments) 
#  Fmat: matrix with reproductive params (can be created by mating.func) 
#  Umat: survival matrix
#  DDapply= across which elements ricker is applied - entire matrix (Amat), survival (Umat), Fertility (Fmat) or recruitment (applies twice to fmat - cub survival and females reproducing
# Use:  Creates a density dependent Amat depending on DDapplication to existing Umat, calculates f values to add in amat


# Removal function - allows bias for age and sex
rem.proj <- function(Umat,   # MAX SURVIVAL
                     initial, 
                     params,
                     stagenames, # needed for mating func
                     time, 
                     memberN=NULL,  # which individuals contribute to pop size? (as vec)
                     DDapply="fertility", 
                     Mfunction= "Min",
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
  mating.out <- mating.func(params, stagenames, Nf, Nm, Mfunction,  return.mat= FALSE)  # return f value only  
  f <- mating.out$f
  
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
  
  
  if(is.numeric(intensity)){   # removal scenarios only
    ry <- as.numeric(remyear)           # shortening name for future use
  } else if(is.null(intensity)){
    ry <- time                   # no removals = run until time
  }
  # Loop = density dependent matrix application for each year
  for (i in 1:ry) {   # normal projection until remyear (if removal at t=5, project until t=5, final entry inserted to row 6) remember vec[5,] holds entries for year=4 (0:time)
    
    # f value per year creation
    thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
    thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
    
    thisN <- sum(Vec[i,memberN])  # pop sizes sums row i for cols included in N
    
    # ricker density dependence each year
    thisAmat <- apply.DD(params, Umat, thisN, DDapply= "fertility", stagenames,   
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
      # ex rem_strat = 3 (ym)
      bi <- length(rem_strat)  # how many elements provided
      
      rem <- rnorm((1- bi), mean = ((1+bias)*intensity)/100, sd = intensity/1000)   # unbiased, 1- number of biased samples needed
      rembi <- rnorm(bi, mean = ((1+bias)*intensity)/100, sd = intensity/1000)
      
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
      
      # apply mating func to calculate pairs
      thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
      
      thisFmat <-thisMating$Fmat
      
      # ricker density dependence each year
      thisN <- sum(Vec[j,memberN])  # pop sizes sums row i for cols included in N
      thisAmat <- apply.DD(params, Umat, thisN, DDapply= "fertility", stagenames,   
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
  out$mat <- thisAmat
  
  return(out)
}

dd_plot <- function(out,   # output obj of dd.proj
                    y_val= "N",   # plot type - N or Vec 
                    ylab = "abundance", 
                    xlab = "time (t)",
                    rem_year = NULL,
                    theme = theme_classic(), 
                    cols= "black",    # can be vector of 2 cols
                    legend.pos = "topright",
                    cex.legend = 0.8){
  # loading required libraries
  require(tidyr)
  require(ggplot2)
  # creating time vector for n years
  t <- nrow(out$vec) -1                 # t= n years (0-t = t+1 entries)
  time <- as.numeric(c(0:t))       # vector 0:t
  
  # for pop size over time graph - must be as df
  if(y_val %in% c("N", "Pop Size", "pop size", "Pop", "pop")) {
    pop <- out$pop
    pop_df <- as.data.frame(pop) # converting pop to a df
    
    if (is.null(rem_year)){
      plot <- ggplot(data = pop_df, aes(x= time, y= pop)) +  # start form year = 0
        geom_point() +
        labs(title = "Pop size over time", x = xlab, y = ylab)  
      geom_smooth(alpha= 0.7)+   # doesn't fit so well for removals!
        theme
      
    } else if (is.numeric(rem_year)){
      plot <- ggplot(data = pop_df, aes(x= time, y= pop)) +  # start form year = 0
        geom_point() +
        geom_line(data = pop_df, alpha = 0.7) +
        labs(title = "Pop size oever time", x = xlab, y = ylab) +
        theme +
        geom_vline(aes(xintercept = rem_year),                       # Adding a line to show removal year
                   colour = "red3", linetype = "dashed", size= 1, alpha = 0.5) 
    }
    
  } else if(y_val %in% c("Vec", "Pop Structure", "Stages", "vec")){
    x_val <- ncol(out$vec)   # number of classes and sexes (if nStages = 2 and sex =2, x =4)
    # turning into dataframe
    df <- as.data.frame(out$vec)
    df$Year <- time    # year column from 0 to t years
    # tidy data - converting to long format so each row is a single observation 
    df_long <- gather(df, key= "Stage", value = "Abundance", 1:x_val)   # creating a stage col in df with abundance
    df_long <- separate(df_long, col= "Stage", into= c("Stage", "Sex"), sep='_')   # splliting by sex, seperated by _
    
    if(is.null(rem_year)) {
      # plotting graph with ggplot2 
      plot <- ggplot(data= df_long, 
                     aes(x=Year, y=Abundance, colour= Sex, linetype= Stage, shape= Stage)) +  # sexes diff cols, shapes and lines diff for stages
        geom_point(position= "jitter", alpha=0.8) +  # jitter to avoid overlap of yearlings
        geom_line(data= df_long, alpha=0.7) +
        scale_colour_manual(values= cols,
                            labels=c("Female", "Male")) +
        labs(title = "Stage Abundance over Time", 
             x = xlab, y = ylab) + 
        theme
      
    } else if(is.numeric(rem_year)){
      plot <- ggplot(data= df_long, 
                     aes(x=Year, y=Abundance, colour= Sex, linetype= Stage, shape= Stage)) +  # sexes diff cols, shapes and lines diff for stages
        geom_point(position= "jitter", alpha=0.8) +  # jitter to avoid overlap of yearlings
        geom_line(data= df_long, alpha=0.7) +
        scale_colour_manual(values=cols,
                            labels=c("Female", "Male")) +
        labs(title = "Stage Abundance over Time", 
             x = xlab, y = ylab) + 
        theme +
        geom_vline(aes(xintercept = rem_year),                       # Adding a line to show removal year
                   colour = "red3", linetype = "dashed", size=0.8, alpha = 0.5) 
      
    } 
    
  } 
  return(plot)
  
}

# growth rate calculate with output format in DDproj
growth.rate <- function(out, vis = FALSE, rem_year = NULL){
  N <- out$pop    # isolating pop size vector
  lambda <- N[2:length(N)]/N[1:(length(N)-1)]
  # output proj
  lambda_out <- list(lambda = vector(), lamb = NA)
  lambda_out$lambda <- lambda
  
  if(is.numeric(rem_year)){
    ry <- rem_year
  }
  
  if(!isFALSE(vis)){
    lamb <- plot(x= 1:length(lambda), y= lambda, 
                 xlab = "Year", 
                 ylab = "lambda") 
    lines(lambda, col= "#94D673")     # keeping cols consistent - red = removals, green = growth rates
    
    title("Pop Growth rate per year", )
    
    lambda_out$lamb <- lamb   # issue here - not loading plot as object within list
  }
  
  return(lambda_out)
}
# Function name : growth.rate
# Arguments : 
# dd_proj or rem_proj output
# vis = TRUE or FALSE for whether you want graph plotted
# Use : generates growth rate per year (gr) by dividing pop size by previous year 



# Stage distribution calculation
# Updating function to match output of ddproj obj
ssd <- function(out, vis = FALSE, cols) {     # input projected matrix
  # separating a matrix from out obj
  mat <- out$vec   
  
  stageMat<- matrix(0, ncol=ncol(mat), nrow=nrow(mat))    # empty matrix to fill
  rownames(stageMat) <- rownames(mat)
  colnames(stageMat) <- colnames(mat)
  
  # out obj is a list with matrix and plot
  ssd_out <- list(stageMat = matrix(), plot = NA)
  
  for(i in 1:nrow(mat)) {   # loop for each column 
    stageMat[i,]<- mat[i,]/sum(mat[i,])          # column i of matrix filled with row i divided by col sum
  } 
  ssd_out$stageMat <- stageMat
  
  if(!isFALSE(vis)){
    # stageMat as df
    nStage <- ncol(stageMat)   # number of classes and sexes (if nStages = 2 and sex =2, x =4)
    # turning into dataframe
    df <- as.data.frame(stageMat)
    df$Year <- as.numeric(rownames(stageMat))    # year column from 0 to t years
    
    # tidy data - converting to long format so each row is a single observation 
    df_long <- gather(df, key= "Stage", value = "Proportion", 1:nStage)   # creating a stage col in df with abundance
    df_long <- separate(df_long, col= "Stage", into= c("Stage", "Sex"), sep='_')   # splliting by sex, seperated by _
    
    plot <- ggplot(data= df_long, 
                   aes(x= Year, y = Proportion, colour = Sex,  # qhy has year ordered so weird?
                       linetype = Stage, shape = Stage)) +  # sexes diff cols, shapes and lines diff for stages
      geom_line(data= df_long, position= "jitter") + # why is this not joining as other 
      scale_colour_manual(values=col_vec,
                          labels=c("Female", "Male")) +
      labs(title = "Stage Proportions over Time", 
           x = "xlab", y = "ylab") 
    
    ssd_out$plot <- plot   # issue here - not loading plot as object within list
  }
  
  return(ssd_out)   # returns matrix of each stage as proportion of total pop
  
}
# Function name= prop.stage
# Arguments: 
# out = stage abundance matrix over a time interval, produced by proj function
# vis = whether to print graph
# cols = vector of 2 colours for each sex
# Purpose:  Calculates the proportion of each stage class out of the total pop size in a given year



# projections (model 1) -----
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
n0 <- c(25, 10, 25, 10)  # vec structure = yf, af, ym, am 

# first time using params from extraction 
Umat <- matrix(0, nrow=4, ncol=4)
rownames(Umat) <- stages
Umat[2,1]<- 0.851 # yearling f survival
Umat[2,2]<- 0.803  # adult f survival
Umat[4,3]<- 0.809   # yearling m survival
Umat[4,4]<- 0.749   # adult m survival

params<- data.frame(fmax= 0.8436,   # F fecundity max (max cubs per adult female) - needed?
                    Sc_max=0.76,   # max cub survival (equal for sexes), rogers 1997
                    b= beta,       # calculated from mcdonald 2016
                    rep_K= 2.7,          # max litter size (K), 
                    h= 6   # harem size per male
)

# control plot----
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
                  rem_strat =NULL ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL ,
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


# Scenario 1 = 70% removal trial at year 10
proj1 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="Fmat", 
                  Mfunction= "Min",
                  intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                  remyear = 10, 
                  rem_strat = "random" ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = NULL ,
                  return.vec= TRUE, 
                  return.remvec = TRUE) 

(proj1_plot <- dd_plot(proj1, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       cex.legend = 0.8))

# scenario 2 - biased female removals
proj2 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages,
                  time = 30, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="Fmat", 
                  Mfunction= "Min",
                  intensity= 70,  # percentage you want REMOVED from pop at time T=ry
                  remyear = 10, 
                  rem_strat = "females" ,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                  bias = 0.15 ,
                  return.vec= TRUE, 
                  return.remvec = TRUE) 

(proj2_plot <- dd_plot(proj2, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       cex.legend = 0.8))


# comparisons 
lambda0 <- growth.rate(proj0, vis = TRUE)
summary(lambda0$lambda)

lambda1 <- growth.rate(proj1, vis = TRUE, rem_year = 10)
summary(lambda1$lambda)

lambda2 <- growth.rate(proj2, vis = TRUE, rem_year = 10)
summary(lambda2$lambda)

ssd0 <- ssd(proj0, vis = TRUE, cols = col_vec)
ssd1 <- ssd(proj1, vis = TRUE, cols = col_vec)
ssd2 <- ssd(proj2, vis = TRUE, cols = col_vec)

# next steps - META MAT PROJECTION
# 