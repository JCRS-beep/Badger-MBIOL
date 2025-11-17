# DD projection function
# 17/11/2025

DD.proj<-function(Fmat, 
                  Umat, 
                  initial, 
                  params, 
                  stagenames, 
                  time, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  AdultNx= TRUE,  # are adults/ oldest class the only individuals included in Nf and Nm calcs? If elderly incl and yearlings not, tweaks needed
                  DDapply="matrix", 
                  return.vec= FALSE) 
{
  # Time 
  if (time<=1) {stop("Time must be a positive integer")
  } else(time <- as.integer(time))
  
  if (is.null(initial)) {stop("You must provide initial population vector for each stage abundance")
  } else(n0 <- as.numeric(initial))
  
  nStages<-length(stagenames)/2      # how many stages
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)){
    stop("You must provide the parameters that match your selected function for density-dependent reproduction.")
  }
  
  # population size
  if (is.null(memberN)){
    memberN <- 1:length(stagenames)   # all members contribute to pop size
    
  } else if(length(memberN >= 1)){
    memeberN <- c(memberN) 
  }
  
  # Nf and Nm from n0. adult only vs yearlings included
  if(AdultNx==TRUE){ # if adults only, final entries for each sex = Nx
    Nf <- n0[nStages]      
    Nm <- n0[2*nStages]
    
  } else if(AdultNx==FALSE){   # if multiple stages included, sum all stages except cubs (n0[1] and n0[nStages+1])
    Nf <- sum(n0[2:nStages])      # omitting 1 (fem cub), yearling to final. 
    Nm <- sum(n0[nStages+2:2*nStages]) 
    
  }
  
  # creating Amat
  Amat <- Fmat+Umat
  
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
    if(AdultNx == TRUE) {   # if only adults
      thisNf <- Vec[i,nStages]    # Nf is mid col in Vec matrix
      thisNm <- Vec[i,2*nStages]       
    } else if(AdultNx == FALSE){
      thisNf <- sum(Vec[i,2:nStages])   # all fems excluding cubs
      thisNm <- sum(Vec[i,(nStages+2):(2*nStages)]) # all males excl cubs- ISSUES IN THIS SYNTAX, SUBSCRIPT OUT OF BOUNDS
    }
    
    thisFmat <- Mating.Fmat(params, stagenames,thisNf, thisNm)   # applies Fmat creation to calculated Nx 
    
    # ricker density dependence each year
    thisN <- sum(Vec[i,memberN])  # pop sizes sums row i for cols included in N
    thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]  # following year stage vector is this Amat* this year vector
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
# Function name- DD.proj
# Inputs
# Fmat: blank matrix that is filled using mating system function for vital rates based on pop structure
# Umat: max predicted survival rates
# initial: population vector for abundance of each stage class
# params: defined in other funcs, importantly b, 
# stagenames: vec of stage classes and sex (2 stage classes should have stagenames length 6)
# time: projection interval
# memberN= which stages included in total pop size count? Defautls NULL, all individuals included. Opts- c(a,b,c) will select those entries of initial vec
# AdultNx= are adults/ oldest class the only individuals included in Nf and Nm? If elderly incl and yearlings not, tweaks needed. TRUE default includes only oldest age class, FALSE incl all stages except first 
# DDapply= how is density dep applied (matrix, fertility, survival, recruitment)
# return.vec= abundance stage class matrix returned? defualts FALSE