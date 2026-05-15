# Mating systems function
# 16/11/2025


mating.func <- function(params,     # density dependent parameters
                        stagenames,   # Stages in life cycle graph (single sex)
                        Nf,        # Adult and yearling females
                        Nm,             # Adult and yearling males
                        Mfunction= "min", 
                        return.mat= TRUE)        # mating function applied
{
 if(is.null(stagenames) || length(stagenames)== 0 & return.mat==TRUE) stop("stagenames must be provided for correct matrix dimension calculations")

  # naming objects used
  K <- params$rep_K
  h <- params$h
  Sf<- params$Sc_f_max    # assumes different survival for male and female cubs
  Sm<- params$Sc_m_max 

  
  # creating the output obj 
  out<- list(U=numeric(), Fmat=matrix())
  U <- NA    # blank
  Fmat <- matrix(0, ncol= length(stagenames), nrow= length(stagenames)) # blank matrix
  rownames(Fmat) <- stagenames
  colnames(Fmat) <- stagenames
  
  #Harmonic mean mating function
  if(Mfunction %in% c( "HM", "HMMF", "Harmonic Mean Mating Function", "harmonic mean mating function", "Harmonic Mean", "harmonic mean")){
  U <- (2* Nf * Nm*h)/(Nf + Nm* h) # females in unions given pop structure
  f <- K*Nm/ (Nm + Nf*h^-1)
  
  # Minimum mating function
  } else if(Mfunction %in% c("Minimum Mating Function", "minimum mating function", "Minimum", "minimum", "min", "Min")){
   # defining number of pairs formed (U)
    if(Nf < Nm*h){
    U <- Nf
   } else if(Nm*h < Nf){
     U <- Nm*h
   } #fertility coefficient
    f<- (K*U)/Nf  # f= cubs produced by adult female
     
     # Mod Harmonic mean mating function?
  } else if(Mfunction %in% c("Modified Hamonic Mean Mating Function", "modified harmonic mean mating function", "ModHarmonic", "modharmonic"))   {
    if(Nf < (2* Nf*Nm*h) / (Nf +Nm*h)) {
      U <- Nf
    } else if((2* Nf*Nm*h) / (Nf +Nm*h) < Nf){
      U <- (2* Nf*Nm*h) / (Nf+ Nm)
    }
    f <- (K*U)/Nf   }
  
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
# params: density dependent parameters, including k (max litter size) & h (harem size), max cub survival 
# stagenames Stages in life cycle , where length()= rows of matrix
# Nf: Number of Adult females in population
# Nm:Number of Adult males
# Mfunc : Mating function, can be Harmonic mean, minimum or mod harmonic mean
# Return.mat= whether Fmat is given in results
# 
# Use: number of pairs formed based on chosen mating function
# Returns maximum possible number of pairs formed in pop based on male and female abundance. Caveats- h may increase with increasing density, then plateau
# For use in matrix projection- use U as Nf instead of vector input. 
# Can create Fmat for use in projection
