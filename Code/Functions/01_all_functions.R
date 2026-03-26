## All required functions
# Run this script before any additional code. (Ctrl + A, Ctrl Ent)
# Includes
# - Mating system function
# - ricker/ dd application function
# - dd projection and removal function
# - dd projection plot function
# - growth rate calculation and plotting
# - stable stage distribution calculation and plotting
# 

##  Mating systems function ----  when Nm or Nf = 0 U = 0
mating.func <- function(params,     # density dependent parameters
                        stagenames,   # Stages in life cycle graph 
                        Nf,        # Adult females
                        Nm,             # Adult males
                        return.mat= FALSE) {       # Fmat output?
  
  if(is.null(stagenames) || length(stagenames)== 0 & return.mat==TRUE) stop("stagenames must be provided for correct matrix dimension calculations")
  
  # naming objects 
  K <- params$rep_K      # MAX litter size
  h <- params$h          # harem size
  S<- params$Sc_max    # cub survival rate (assumes equal for male and female cubs)
  
  
  # creating outputs  
  out<- list(f=numeric(), Fmat=matrix())
  f <- numeric()    # fertility value (cubs produced by adult female) 
  
  
  # handling NA values 
  if (is.na(Nf) || is.na(Nm)) {
    warning("Nf or Nm is NA: setting f = 0")
    f <- 0
  }
  # handling negatives values
  else if (Nf <= 0 || Nm <= 0) {
    # no mating possible if either sex absent or zero
    out$f <- 0
  }
    
  # Minimum mating function to define number of pairs formed (U)
  # U = min(nf, Nm*h)
  U <- min(Nf, Nm * h)
  
  f<- (K*U)/Nf   # fertility coefficient f= cubs produced by adult female
  
  if(return.mat==TRUE) {  # visualising into fertility matrix
    Fmat <- matrix(0, ncol= length(stagenames), nrow= length(stagenames)) # blank matrix
    rownames(Fmat) <- stagenames
    colnames(Fmat) <- stagenames
    
    nStages <- length(stagenames)/2 # number stages for each sex
    Fmat[1,nStages] <- 0.5* f* S    # female cub production by females
    Fmat[nStages+1, nStages] <- 0.5* f* S  # male cub production by females
    out$Fmat <- Fmat
  }
  
  out$f <- f
  
  return(out) 
}
# Function name= Mating.func
# Inputs:
#  params: density dependent parameters, including k (max litter size),  h (harem size), max cub survival 
#  stagenames Stages in life cycle , where length()= rows of matrix
#  Nf: Number of Adult females in population (can incl yearlings if desired)
#  Nm: Number of Adult males
#  Return.mat: whether Fmat is given in results
# 
# Use: Returns max possible number of pairs formed based on male and female abundance. 
# Outputs:
# f = cub production per female

##  Ricker function and matrix application ------
ricker <- function(params, N){   # Using params (b) and population size input
  dd_fun <- exp(-params$b*N)
  return(dd_fun)      # returns value of multiplier
}
# Function ricker
# params, incl b = strength of density dependence
# N= population size (can be total, NAdults, other, but explain in comments) 

apply.DD <- function(params, 
                     Umat, 
                     N,   # yearling and adults
                     DDapply="fertility", 
                     stagenames,   # Stages in life cycle graph 
                     Nf,        # Adult female abundance
                     Nm,             # Adult males
                     return.mat= FALSE) {   # apply ricker to whole matrix, survival or fertility
  
  # Call to ricker() function
  rick <- ricker(params, N)    # rick = multiplier to be applied to fertility later
  
  # Call to mating.func()
  mating <- mating.func(params,     # density dependent parameters
                        stagenames,   # Stages in life cycle graph 
                        Nf,        # Adult  females
                        Nm,             # Adult  males
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
# Use:  Creates a density dependent Amat depending on DDapplication to existing matrices, Uma and Fmat



## IMPROVEMENTS = integer individuals only. proportion must be bound between 0 and 1. 
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
  
  
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames) == 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("Please provide parameters for selected density-dependent function")
  
  
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
       } else if(rem_strat == "random" && is.null(bias) == TRUE){
         stop("please provide strength of bias for no random removal strategies")
      
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
          
          
          # specified
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
# NOTES  Removal function  ------
# Function name - rem.proj
# Inputs - 
#  Umat: max predicted survival rates
#  initial: population vector for abundance of each stage class
#  params: defined in other funcs, importantly b, 
#  stagenames: vec of stage classes and sex (2 stage classes should have stagenames length 4)
#  time: projection interval
#  DDapply= how is density dep applied (matrix, fertility, survival, recruitment)
#  return.vec= abundance stage class matrix returned? defaults FALSE
# Use - project an initial population vector over t years using density dependence at each time step to adjust vital rates in matrix. 
# Outputs - out
#  $pop: vector of total population size each year
#  $vec: matrix of abundance of each stage each year
#  $mat: returns final Amat produced
#  $Nremoved = total of pop removed
#  $rem_vec = optional vector of removals per class






# Projection plotting function ---- increasing text size and ensuring axis labelled 
dd.plot <- function(out,   # output obj of dd.proj
                    y_val= "N",   # plot type - N or Vec 
                    ylab = "abundance", 
                    xlab = "time (t)",
                    rem_year = NULL,
                    mytheme = theme_classic(), 
                    cols= "black",    # can be vector of 2 cols
                    legend.pos = "top",
                    base_size = 16){
  # loading required libraries
  require(tidyr)
  require(ggplot2)
  # creating time vector for n years
  t <- nrow(out$vec) -1                 # t= n years (0-t = t+1 entries)
  time <- as.numeric(c(0:t))       # vector 0:t
  
  # for pop size over time graph - must be as df
  if(y_val %in% c("N", "Pop Size", "pop size", "Pop", "pop")) {
    pop_df <- data.frame(Year = time,
                         Pop  = out$pop )# converting pop to a df
    
    # creating pop projection plot
      base.plot <- ggplot(pop_df, aes(x = Year, y = Pop)) +
        geom_line(size = 1.2, colour = "grey30") +
        geom_point(size = 2, colour = "black") +
        labs(title = "Population Size Over Time",
             x = xlab,
             y = ylab) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
        mytheme +
        theme(
          text = element_text(size = base_size),
          plot.title = element_text(size = base_size + 2, face = "bold"),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 2),
          legend.position = legend.pos
        ) 
      
     if (!is.null(rem_year)){
      plot <- base.plot +
        geom_vline(xintercept = rem_year,                       # Adding a line to show removal year
                   colour = "red3", 
                   linetype = "dashed", linewidth = 1, 
                   alpha = 0.5) 
      
     } else if (is.null(rem_year)){
      plot <- base.plot 
    }
    
      # alternative plot = vec abundance
  } else if(y_val %in% c("Vec", "Pop Structure", "Stages", "vec")){
    x_val <- ncol(out$vec)   # number of classes and sexes (if nStages = 2 and sex =2, x =4)
    # turning into dataframe
    df <- as.data.frame(out$vec)
    df$Year <- time    # year column from 0 to t years
    # long format so each row is a single observation 
    df_long <- gather(df, key= "Stage", value = "Abundance", 1:x_val)   # creating a stage col in df with abundance
    df_long <- separate(df_long, 
                        col= "Stage", 
                        into= c("Stage", "Sex"), 
                        sep='_')   # splitting by sex, separated by _
    
      # plotting graph 
      base.plot <- ggplot(data= df_long, 
                     aes(x=Year, 
                         y=Abundance, 
                         colour= Sex, 
                         linetype= Stage, 
                         shape= Stage)) +  # sexes diff cols, shapes and lines diff for stages
        geom_point(position= "jitter", size = 2, alpha=0.8) +  # jitter to avoid overlap of yearlings
        geom_line(data= df_long, size = 1.2, alpha=0.7) +
        scale_colour_manual(values= cols,
                            labels=c("Female", "Male")) +
        labs(title = "Stage Abundance Over Time",
             x = xlab,
             y = ylab,
             colour = "Sex",
             linetype = "Stage",
             shape = "Stage") +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
        mytheme +
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
      
     if(!is.null(rem_year)){   # if rem years a vector, can we loop geom_line addition?
      plot <- base.plot +
        geom_vline(xintercept = rem_year,                       # Adding a line to show removal year
                   colour = "red3", linetype = "dashed", size=0.8, alpha = 0.3) 
      
     } else if (is.null(rem_year)){
       plot <- base.plot
     }
  }
    
  return(plot)
  
}   

# NOTES ----
# Function name = dd_plot
# Inputs 
#  out: output obj from dd_proj function, including structure vec or popsize (N)
#  y_val= plot type. N = total pop size by year, Vec = stage abundance by year 
#  ylab : defaults "abundance", 
#  xlab : defaults "time (t)"
#  col: "black",    colours used for sexes in graph, should be length 2
#  legend.pos: "topright",   legend position on graph
#  cex.legend: legend size
# Use - Take output from projection and plot with ggplot, [including custom theme if desired]
# Required packages = ggplot2, tidyr




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
           x = "Year", y = "Proportion of population") 
    
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




# creating repeatable projections using proj functions, varying initial pop vector
# Need a way to stop current loop if pop size falls to zero but continue with other loops
repeat.proj <- function(Umat,   # MAX SURVIVAL
                        params,
                        stagenames, # needed for mating func
                        time = 20, 
                        DDapply="fertility", 
                        method = "norm",  # how to generate initial pop sizes?
                        intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                        remyear = NULL,  # removal year = decrease from following year
                        rem_strat = "random",  # if specified removals, "adults, females, yearlings, males, yearling females, 
                        bias = NULL , # strength of bias as percentage (range??)
                        return.vec= TRUE, 
                        return.remvec = TRUE,
                        reps = 10) {
  # set up the ouput = a list, length = number of repeats, each containing the out obj of  projection function
  out <- vector("list", reps) 
  
  # looping projection for as many repetitions desired
  for (t in 1:reps){
    # generating inital vec = want to maintain similar stage dist but vary population size
    thisInitial <- vector(length = length(stagenames))   # list to fill with initial vector for each repetition
    
    stageDist <- c(0.1204889, 0.4170396, 0.1212302, 0.3412414) # calculated from repeated baseline proj average proportions per rep, then using colmeans
    
   # rnorm = aiming for pops around certain size, with wide variation
   # runif vs random starting pop size (how to define min and max?)
   
    # generating pop size
    if(method == "random"){
    size <- floor(runif(1, min = 40, max = 250))  # defining min and max
    
    } else if (method == "norm"){
      size <- floor(rnorm(1, mean = 137, sd  = 50))   # defining mean from repeated baseline means. large sd for variation
    }
    # seperate samples needed for adults and yearlings = 
    
    thisInitial <- size * stageDist  #multiply generated pop size with fixed stage dist
    
    out[[t]] <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                         thisInitial, 
                         params, 
                         stagenames,
                         time, 
                         DDapply,
                         intensity,  # percentage you want REMOVED from pop at time T=ry
                         remyear,  # removal year = decrease from following year
                         rem_strat,  # if specified removals, "adults, females, yearlings, males, yearling females, 
                         bias,
                         return.vec, 
                         return.remvec) 
    # warning of extinction and what iteration this happened in 
    if(any(out[[t]]$pop<= 0)) {
      warning(paste("Projection reached pop size 0 or below in iteration", t))
    }
  }
  return(out)
  
}
# output syntax = out[[rep]]$pop[] or $vec[,]




# calculating average lambda, av ssd in single function
# Splitting calculations into their own functions? Av lambda, relative pop sizes?

# calculating pop size per rep function
N.extract <- function(proj_list, reps) {
N_list <- vector("list", reps)   # list of N vecs for each repeat

for (t in 1:reps){
  N_list[[t]] <- proj_list[[t]]$pop    # isolating pop size vector
}

return(N_list)
}


lamb.av <- function(proj_list, return.Lambda = FALSE){
  reps <- length(proj_list)
  time <- length(proj_list[[1]]$pop - 1) # minus 1 for the initial population at t = 0
  
  N_list <- N.extract(proj_list, reps)
  
  # calculating lambda for each rep lambda calculation as a function applied to pop size vec
  lambda <- function(N){
    N[2:length(N)]/N[1:(length(N)-1)]
  } 
  
  lambda_list <- lapply(N_list , lambda)
  av_lambda <- sapply(lambda_list, mean)
  
  # setting up out obj
  lamb_out <- list(lambda_list = vector(), 
                   av_lambda = vector() )
  
  if(return.Lambda == TRUE){
    lamb_out$lambda <- lambda_list
  }
  lamb_out$av_lambda <- av_lambda
}
# lambda : lambda per year per rep (list of 100 vecs length 20)
# av_lambda : average lambda per rep (vec length 100)


# ssd function  - 
ssd.av <- function(proj_list, return.Mats = FALSE){
  
  reps <- length(proj_list)
  time <- length(proj_list[[1]]$pop - 1)
  
  # separating a matrix from out obj
  abun_mat <- vector("list", reps)   # list of abundances for each repeat
  
  for (t in 1:reps){
    abun_mat[[t]] <- proj_list[[t]]$vec    # isolating pop size vector
  }
  
  mat <- matrix(0, ncol = ncol(abun_mat[[1]]), nrow = nrow(abun_mat[[1]]))    # list of empty matrix to fill with stage props
  
  
  stageMat <- vector("list", reps)   # fill each stage mats list with mat?
  
  # out obj is a list with matrix and av prop 1 row matrix
  
  for (t in 1:reps){
    for(i in 1:nrow(abun_mat[[1]])) {   # loop for each column 
      mat[i,] <- abun_mat[[t]][i,]/sum(abun_mat[[t]][i,])          # column i of matrix filled with row i divided by col sum
    } 
    stageMat[[t]] <- mat # ssd for each year per rep
  }
  
  list_prop <- lapply(stageMat, colMeans) # list per rep, each single vec of av abundances
  
  # combining lists into a single matrix
  av_prop <- matrix(unlist(list_prop), byrow = TRUE, nrow = reps, ncol = 4) # filling each row with list vec
  colnames(av_prop) <- colnames(abun_mat[[1]])
  
  # calculating sex ratio (proportion female in pop)
  fems <- function(x) {
    sum(x[c(1:2)])   # sum entries 1+2 = proportion fems in pop
    }
  sr <- sapply(list_prop, fems)  # list of sex ratios as proportion female
  
  
  # set up out obj
  ssd_out <- list(ssdMat = matrix(), av_prop = matrix(), sex_ratio = vector())
  
  if(return.Mats == TRUE){
    ssd_out$ssdMat <- stageMat # stage proportions by year
  }
  
  ssd_out$av_prop <- av_prop  # single stage proportion for each rep 
  
  ssd_out$sex_ratio <- sr  
  
  return(ssd_out)
}
# Outputs 
# ssdMat : stage dist mat
# av_prop : average proportions of each stage in matrix   




relative.pop <- function(proj_list,   
                         baseline_list = NULL) # baseline for comparison. if empty, no relative population sizes calculated
  {
  # basic checks
  if(length(baseline_list)!= length(proj_list)){
    stop(paste("length of baseline and projection differ - comparison not possible. Ensure repetitions are equal before continuing."))
  }
  
  reps <- length(proj_list)
  time <- length(proj_list[[1]]$pop - 1) # minus 1 for the initial population at t = 0
  
  N_list <- N.extract(proj_list, reps)   # list of N vecs for each repeat
  base_list <- N.extract(baseline_list, reps)
  
  
  # average pop size proportions
  av_N <- sapply(N_list, mean)   # vector of mean pop sizes per rep
  av_base_N <- mean(sapply(base_list, mean))   # SINGLE value for mean pop size across all reps of baseline
  
  if(!is.null(baseline_list)){   # only if there is comparison proj provided
  # calculating final pop size in each scenario  (av?)
  fin_props <- vector(length = reps)   # vector of proportions
  
  for (t in 1:reps){
    fin_props[t] <- N_list[[t]][time]/base_list[[t]][time]     #  list calculate proportions
  }

  av_base_N <- sapply(base_list, mean)  
  relative_meanN <- av_N / av_base_N   # issues = if baseline proj goes to 0 for some reason, causes
  }
  
  # output proj
  pop_out <- list(fin_props = vector(), av_fin = vector(), 
                  relative_meanN= vector()) # number of matrices = rep, rows = time
                                                 
  
  
  pop_out$fin_props <- fin_props
  pop_out$av_fin <- mean(fin_props, na.rm = TRUE)
  pop_out$relative_meanN <- relative_meanN
  
  
  return(pop_out)
}
# NOTES
# Name =  av_pop
# proj_list,   
# baseline_list : the baseline projection to compare projections
# return.Lambda : return lambda per year for each rep? T/F
# return.Mats : return stage distribution per year for each rep? T/F
# 
# Outputs
# fin_props : final year proportions of projection compared to baseline (vec length 100)
# av_fin : average final pop size of proj compared to baseline (numeric val).
# relative_meanN : mean pop size across proj compared to mean pop size of baseline (vec length 100)


### ISSUES WITH COMPARISONS - if there is a year when baseline goes to extinction, this skews that years comparison. 
  # better to have an average pop size for baseline and compare all repeats to this value?


# recombining all functions into single function to apply in df create?








# meta.projection
meta.proj <- function(Umat,   # vector of stage names
                      params, # adjusted beta = 
                      stagenames,
                      initial,   # list of initial abundances
                      Dmat = list(), 
                      time = 20,
                      return.vec = TRUE,
                      return.mats= FALSE){
  
  # Time 
  if (time <= 1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  if(is.list(initial)){   # when initial abundances is a list of vectors, we have more than 1 patch 
    patches = length(initial)
  } else if(is.list(initial) == FALSE ){  # if a vector/ not list, only 1 patch
    patches = 1
    stop("only one vector of abundances provided. Provide abundances per patch as a list, or use rem.proj for single patch projections")
  }
  
  nStages <- length(stagenames)/2      # how many stages
  
  Nf_v <- numeric(length= patches) # vector length patches to fill with adult fem abundance by year
  Nm_v <- numeric(length= patches)
  
  
  for (i in 1:patches) {           # filling each list within vector with Nf and Nm at that patch initially
    Nf_v[i] <- initial[[i]][nStages]   
    Nm_v[i] <- initial[[i]][2*nStages]      # adult male in vector (final entry) - vec of Nms per patch
  }
  
  
  # f calculation first year - how to repeat for each patch
  f_v <- numeric(length= patches)   # blank vec to fill  before adding fertility values straight to amat - f values per patch f[p] 
  
  for (j in 1:patches) { # repeat for each patch
    mat <- mating.func(params,     # density dependent parameters
                       stagenames,   # Stages in life cycle graph (single sex)
                       Nf_v[j],        # Adult females  
                       Nm_v[j],             # Adult males
                       return.mat= FALSE) # don't return mat, just f value
    
    f_v[j] <- mat$f   # adding appropriate fertility rates to list, where entry = patch  
  } 
  
  # Set up the output - how to set up if multiple patches? Lists within lists...
  
  out <- list(pop = list(numeric(length= (time + 1))),    #  list of vectors per patch
              vec = list(matrix(0, ncol = length(stagenames), nrow = time + 1)),    # list of matrices per patch
              Nremoved = (rep(NA, 2*nStages)), # total removed in rem year, NO LIST, vec length = num patches
              remvec = list(rep(NA, 2*nStages)))    # vector of length(stages) per patch - list
  
  # renaming before repeating
  colnames(out$vec[[1]]) <- stagenames   # naming cols matrix as stages - cant name lists of obj
  rownames(out$vec[[1]]) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or initial
  
  
  Vec <- rep(out$vec, patches)  # matrices to fill with stage abundance.  row= time, col= stage. out$vec or just vec?
  Pop <- rep(out$pop, patches)       # vector to fill with total pop size each year
  remVec <- rep(out$remvec, patches) # list of vectors per patch, only 1 year 
  
  # filling first row of pop size and abundance vecs
  pop_sizes <- sapply(initial, sum)  # gives total pop size per patch in initial year
  
  for (p in 1:patches){  # filling first row of each list vec and pop size 
    Vec[[p]][1,] <- initial[[p]]          # for each list in Vec, first row is initial         
    
    Pop[[p]][1] <- pop_sizes[p]   # each list row 1 is sum of n0 - issues in plugging into pop
  }
  
  
  # Loop = matrix proj each year  --------
  for (i in 1:time) {   # repeat for as many years as we have
    
    # Mating Fmat creation - f values per patch
    thisNf <- as.vector(numeric(length= patches))  # vec of current Nf for each patch
    thisNm <- as.vector(numeric(length= patches))
    
    for (p in 1:patches){  # loop for each patch this year
      thisNf[p] <- Vec[[p]][i,nStages]  # inputting calculated value to appr place in vec. Issues with using both n and i in this loop 
      thisNm[p] <-  Vec[[p]][i,2*nStages]                                  #repeat for each patch, but as we loop through years will have to update row  of vec selected
    }
    
    # apply mating func to calculate pairs per patch 
    thisMating <- rep(NA, patches)        # list of repeated matrix 
    for (p in 1:patches) { # repeat for each patch
      mat_out <- mating.func(params,     # density dependent parameters
                             stagenames,   # Stages in life cycle graph (single sex)
                             thisNf[p],        # Adult females  
                             thisNm[p],             # Adult males
                             return.mat= FALSE)
      
      thisMating[p] <- mat_out$f  # adding appropriate Fmat to list 
    } 
    
    # ricker density dependence each year - total pop size by patch
    
    thisN <- numeric(length= patches)  # vector length patches to fill with current pop size each patch 
    for (p in 1:patches){
      thisN[p] <- sum(Vec[[p]][i, ])   # length patches, entries = pop size that year. this N entry [1] is pop size(row = year, col = p)
    }
    # Amat each year will vary by patch 
    thisAmat <- list(matrix(0, ncol = ncol(Umat), nrow= nrow(Umat)))
    thisAmat <- rep(thisAmat, patches)        # list of repeated matrix 
    
    for (p in 1:patches){
      amat <- apply.DD(params, 
                       Umat, 
                       thisN[p],   # yearling and adults
                       DDapply="fertility", 
                       stagenames,   # Stages in life cycle graph 
                       thisNf[p],        # Adult female abundance
                       thisNm[p],             # Adult males
                       return.mat= FALSE)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
      
      thisAmat[[p]] <- amat
      
      
    }
    # If the projection matrix has any negative values in it, stop iterating and
    # return the projection up until this point.
    if (any(sapply(thisAmat, function(m) any(m < 0, na.rm = TRUE)))) {
      warning(paste("Negative values in projection matrix at time step", i,
                    "- setting negative entries to 0 and continuing."))
      thisAmat <- lapply(thisAmat, function(amat) { amat[amat < 0] <- 0; amat}) 
    }
    
    remains <- vector("list", length(stagenames))   # list of remaining individual abundances 
    moves <- vector("list", length(stagenames))     # list of moves individual abundances
    
    # projection for each patch before dispersal
    for (p in 1:patches) {
      current_vec <- as.numeric(Vec[[p]][i, ])         # current stage-abundance vector
      
      pStay <- as.numeric(Dmat[[p]][1, ])        # stage prob of staying in patch p
      pMove <- as.numeric(Dmat[[p]][2, ])        # stage p(leave) patch p
      
      # abundance BEFORE  survival/rep
      remains[[p]] <- current_vec * pStay        #  number who remain to breed
      moves[[p]]  <- current_vec * pMove       # number who leave (no breeding)
    }
    
    # Project stayers using Amat (they reproduce in original patch)
    remains_proj <- vector("list", patches)
    for (p in 1:patches) {
      remains_proj[[p]] <- as.numeric(thisAmat[[p]] %*% remains[[p]])
    }
    
    # moves:  leave BEFORE mating, only survive (Umat)
    moves_surv <- vector("list", patches)
    for (p in 1:patches) {
      moves_surv[[p]] <- as.numeric(Umat %*% moves[[p]])
    }
    
    # Assemble next-year vectors: each patch gets its own stayers_proj + incoming movers from the other patch
    for (p in 1:patches) {
      other <- ifelse(p == 1, 2, 1)
      
      next_vec <- remains_proj[[p]] + moves_surv[[other]]  # combine remaining abundances with those who leave other patch
      
      # guard against numeric noise/negatives
      next_vec[next_vec < 0] <- 0
      
      Vec[[p]][i + 1, ] <- next_vec
      Pop[[p]][i + 1] <- sum(Vec[[p]][i + 1, ])
    }
    
    # if pop size <= 0, stop and return
    total_pop <- sum(sapply(Pop, function(x) x[i + 1]))
    if (total_pop <= 0) {
      warning(paste("Projection stopped at time step", i, "because total pop size reached 0 or below"))
      break
    }
    # if any stage becomes negative, set to zero and continue
    if (any(sapply(Vec, function(mat) any(mat < 0, na.rm = TRUE)))) {
      warning(paste("Negative abundances produced at time step", i, "— setting negatives to 0 and continuing."))
      Vec <- lapply(Vec, function(mat) { mat[mat < 0] <- 0; mat })
    }
  }
  
  # out objects
  out$pop <- Pop # total group size per patch, do we also want total pop size of pop?
  
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  }
  
  return(out)
  
}