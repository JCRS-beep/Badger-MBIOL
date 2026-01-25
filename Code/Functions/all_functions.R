## All required functions
# Run this script before any additional code. (Ctrl + A, Ctrl Ent)
# Includes
# - ricker/ dd application function
# - Mating system function
# - dd projection function
# - dd projection plot function

##  Ricker function and matrix application ------
apply.DD <- function(params, Fmat, Umat, N, DDapply="matrix") {   # apply ricker to whole matrix, survival or fertility
 
  # making sure mats match 
  if (sum(dim(Fmat)==dim(Umat))!=2){
    stop("Reproductive (Fmat) and survival (Umat) matrices 
         must have the same dimensions.")
  }
  
  # Function ricker
  # params, incl b = strength of density dependence
  # N= population size (can be total, NAdults, other, but explain in comments) 
  # Use: returns dd_fun multiplier.
  ricker <- function(params, N){   # Using params (b) and population size input
    dd_fun <- exp(-params$b*N)
    return(dd_fun)      # returns value of multiplier
  }
  rick <- ricker(params, N)    # ricker function embedded
  
  # constructing Amat
  Amat <- Fmat + Umat
  
  # applying based on method
  if(DDapply %in% c("matrix", "Matrix"))  {
    Amat_N <- Amat*rick  
    
  } else if(DDapply %in% c("Survival", "survival", "Umat")){
    Umat_N <- Umat*rick
    Amat_N <- Fmat + Umat_N
    
  } else if(DDapply  %in% c("Fertility", "fertility", "Fmat")) {
    Fmat_N <- Fmat*rick
    Amat_N <- Fmat_N + Umat
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


##  Mating systems function ----
mating.func <- function(params,     # density dependent parameters
                        stagenames,   # Stages in life cycle graph (single sex)
                        Nf,        # Adult and yearling females
                        Nm,             # Adult and yearling males
                        Mfunction= "min",       # mating function applied
                        return.mat= TRUE){       # Fmat output?

  if(is.null(stagenames) || length(stagenames)== 0 & return.mat==TRUE) stop("stagenames must be provided for correct matrix dimension calculations")
  
  # naming objects 
  K <- params$rep_K      # MAX litter size
  h <- params$h          # harem size
  Sf<- params$Sc_f_max    # assumes different survival for male and female cubs
  Sm<- params$Sc_m_max 
  
  
  # creating the output obj 
  out<- list(U=numeric(), Fmat=matrix())
  U <- NA    # blank 
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
    Fmat[1,nStages] <- 0.5* f* Sf    # female cub production by females
    Fmat[nStages+1, nStages] <- 0.5* f* Sm  # male cub production by females
    out$Fmat <- Fmat
  }
  
  out$U <- U
  
  return(out) 
}
# Function name= Mating.func
# Inputs:
#  params: density dependent parameters, including k (max litter size),  h (harem size), max cub survival 
#  stagenames Stages in life cycle , where length()= rows of matrix
#  Nf: Number of Adult females in population (can incl yearlings if desired)
#  Nm: Number of Adult males
#  Mfunc: Mating function, can be Harmonic mean, minimum or mod harmonic mean
#  Return.mat: whether Fmat is given in results
# 
# Use: Returns max possible number of pairs formed based on male and female abundance. 
# Caveats- h may increase with increasing density, then plateau

##  DD projection function

# WARNING: MUST RUN MATING SYSTEM FUNCTION BEFORE THIS FUNCTION

dd.proj <- function(Umat,   # MAX SURVIVAL
                    initial, # pop structure vec
                    params,  # dd params
                    stagenames, 
                    time,     # how far to project
                    memberN=NULL,  # which individuals contribute to pop size? (as vec)
                    DDapply="matrix",  
                    Mfunction= "Min",
                    return.vec= TRUE) {   # return Fmat from mating func
  # Time 
  if(time<=1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if(is.null(initial)) stop("You must provide initial population vector for each stage abundance")
  else(n0 <- as.numeric(initial))
  
  
  if(length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimensions")
  if(is.null(params)) stop("You must provide the parameters that match your selected function for density-dependent reproduction.")
  
  
  nStages<-length(stagenames)/2      # how many stages in lifecycle
  
  # population size
  if (is.null(memberN)){
    memberN <- 1:length(stagenames)   # all members contribute to pop size
    
  } else if(length(memberN) >= 1){ 
    memberN <- c(memberN)       # vector of stages in lifecycle incl (in case, 2,4 = adult f and m)
  }
  
  # Calculating Unions with mating function
  Nf <- n0[nStages]       # pulls adult female entry from initial vector
  Nm <- n0[2*nStages]      # adult male in vector (final entry)
  
  # mating func gives initial Fmat for first year
  mating.out <- mating.func(params, stagenames, Nf, Nm, Mfunction,  return.mat= TRUE)  # Npairs and Fmat given 
  U <- mating.out$U 
  Fmat <- mating.out$Fmat 
  
  # Initialize the output for population vector:
  out<- list(pop=vector(), vec=matrix(), mat=matrix())
  
  Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- rep(NA, (time + 1))       # vector to fill with total pop size each year
  colnames(Vec) <- stagenames   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- n0                   
  Pop[1] <- sum(n0)    
  
  # Loop = density dependent matrix application for each year
  for (i in 1:time) {  
    
    # Mating Fmat creation 
    thisNf <- Vec[i,nStages]    # Nf is mid col in Vec matrix - should sum yearlings and adults!
    thisNm <- Vec[i,2*nStages]  # Nm 
    
    # apply mating func to calculate pairs
    thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
    
    thisFmat <-thisMating$Fmat
    thisU <- thisMating$U
    
    # ricker density dependence each year
    thisN <- sum(Vec[i,memberN])  # pop sizes sums row i for cols included in N
    thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
    Pop[i + 1] <- sum(Vec[(i + 1), ])
  }
  
  # out objects
  out$pop <- Pop
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 
  
  out$mat <- thisAmat
  
  return(out)
}

# Notes----
# Function name - DD.proj
# Inputs - 
#  Umat: max predicted survival rates
#  initial: population vector for abundance of each stage class
#  params: defined in other funcs, importantly b, 
#  stagenames: vec of stage classes and sex (2 stage classes should have stagenames length 4)
#  time: projection interval
#  memberN= which stages included in total pop size count? Defautls NULL, all individuals included. Opts- c(a,b,c) will select those entries of initial vec
#  DDapply= how is density dep applied (matrix, fertility, survival, recruitment)
#  return.vec= abundance stage class matrix returned? defaults FALSE
# Use - project an initial population vector over t years using density dependence at each time step to adjust vital rates in matrix. 
# Outputs - out
#  $pop: vector of total population size each year
#  $vec: matrix of abundance of each stage each year
#  $mat: returns final Amat produced






# Projection plotting function ---- 
dd_plot <- function(out,   # output obj of dd.proj
                    y_val= "N",   # plot type - N or Vec 
                    ylab = "abundance", 
                    xlab = "time (t)",
                    theme = theme_classic(), 
                    cols= "black",    # can be vector of 2 cols
                    legend.pos = "topright",
                    cex.legend = 0.8){
  # loading required libraries
  library(tidyr)
  library(ggplot2)
  
  # creating time vector for n years
  t <- nrow(out$vec) -1                 # t= n years (0-t = t+1 entries)
  time <- as.numeric(c(0:t))       # vector 0:t
  
  # for pop size over time graph
  if(y_val %in% c("N", "Pop Size", "pop size", "Pop", "pop")) {
    plot <- ggplot(data = out, aes(x= time, y= pop)) +  # start form year = 0
                      geom_point() +
                      xlab(xlab) +
                      ylab(ylab) +
                      geom_smooth()+
                      theme_bw()
    
  } else if(y_val %in% c("Vec", "Pop Structure", "Stages", "vec")){
    x_val <- ncol(out$vec)   # number of classes and sexes (if nStages = 2 and sex =2, x =4)
    # turning into dataframe
    df <- as.data.frame(out$vec)
    df$Year <- time    # year column from 0 to t years
    # tidy data - converting to long format so each row is a single observation 
    df_long <- gather(df, key= "Stage", value = "Abundance", 1:x_val)   # creating a stage col in df with abundance
    df_long <- separate(df_long, col= "Stage", into= c("Stage", "Sex"), sep='_')   # splliting by sex, seperated by _
    
    # plotting graph with ggplot2 
    plot <- ggplot(data= df_long, aes(x=Year, y=Abundance, colour= Sex, linetype= Stage, shape= Stage))+  # sexes diff cols, shapes and lines diff for stages
      geom_point(position= "jitter", alpha=0.8)+  # jitter to avoid overlap of yearlings
      geom_line(data= df_long, alpha=0.7) +
      scale_colour_manual(values=cols,
                          labels=c("Female", "Male")) +
      labs(title = "Stage Abundance over Time", 
           x = xlab, y = ylab) + 
      theme
    
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
