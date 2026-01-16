# Removal  function
# 14.01.2025

rem.proj <- function(Umat,   # MAX SURVIVAL
                     initial, 
                     params, 
                     stagenames, 
                     time, 
                     memberN=NULL,  # which individuals contribute to pop size? (as vec)
                     DDapply="matrix", 
                     Mfunction= "Min",
                     intensity= 50,  # percentage you want REMOVED from pop at time T=ry
                     remyear = 25,  # removal year = decrease from following year
                     return.vec= TRUE) 
{
  # Time 
  if (time<=1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  
  nStages<-length(stagenames)/2      # how many stages
  ry <- as.numeric(remyear)           # shortening name for future use
  
  # population size
  if (is.null(memberN)){
    memberN <- 1:length(stagenames)   # NULL = all members contribute to pop size
    
  } else if(length(memberN) >= 1){
    memberN <- c(memberN) 
  }
  
  # Calculating Unions with mating function
  Nf <- n0[nStages]     # pulls yearling and adult female entry from initial vector
  Nm <- n0[2*nStages]      # adult male in vector (final entry)
  
  # mating func gives initial Fmat for first year
  mating.out <- mating.func(params, stagenames, Nf, Nm, Mfunction,  return.mat= TRUE)  # Npairs and fmat given 
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
  for (i in 1:ry) {   # normal projection until remyear (if removal at t=5, project until t=5, final entry inserted to row 6) remember vec[5,] holds entries for year=4 (0:time)
    
    # Mating Fmat creation 
    thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
    thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
    
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
  # ------------------------------------------------
  # removal year 
  Rem <- Pop[ry + 1] * intensity/100   # ry +1 IS ROW REPRESENTING YEAR remyear! Nremoved by calculating latest N, multiplying by percentage 
  
  # rem_vec <- rep(NA, length(stagenames))  # generate number vec length= stages # sum(rem_vec) == Rem
  
  # sample_vec <- c(1:Rem)
  # for (r in 1:Rem){    # sampling loop until numbers generated
  #  rem_vec <- sample(sample_vec, size = length(stagenames), replace = TRUE)
  # }   
  
  # for now, simply divide rem/4 and remove from each class
  rem_vec <- rep(Rem/4, length(stagenames))
  
  Vec[ry + 2,] <- Vec[ry + 1,]  - rem_vec  # year after remyear = 2 rows later filled with new stage vec
  Pop[ry + 2] <- Pop[ry+1] - Rem  # remyear starting pop = prev N - Nremoved. *edit after*
  # -------------------------------------------------
  
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
    
    Vec[(j + 1), ] <- thisAmat %*% Vec[j, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
    Pop[j + 1] <- sum(Vec[(j + 1), ])
  }
  
  # out objects
  out$pop <- Pop
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 
  
  out$mat <- thisAmat
  
  return(out)
}
# Function name - DD.proj
# Inputs - 
#  Umat: max predicted survival rates
#  initial: population vector for abundance of each stage class
#  params: defined in other funcs, importantly b, 
#  stagenames: vec of stage classes and sex (2 stage classes should have stagenames length 4)
#  time: projection interval
#  memberN= which stages included in total pop size count? Defautls NULL, all individuals included. Opts- c(a,b,c) will select those entries of initial vec
#  DDapply= how is density dep applied (matrix, fertility, survival, recruitment)
#  intensity = what percentage are you removing
#  remyear = after how many years is removal occuring?
#  return.vec= abundance stage class matrix returned? defaults FALSE
# Use - project an initial population vector over t years using density dependence at each time step to adjust vital rates in matrix. 
# Outputs - out
#  $pop: vector of total population size each year
#  $vec: matrix of abundance of each stage each year
#  $mat: returns final Amat produced
# SUCESS !!