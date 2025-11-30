# Mating systems function
# 16/11/2025


mating.func <- function(params,     # density dependent parameters
                    #   stagenames,   # Stages in life cycle graph (single sex)
                       Nf,        # Adult and yearling females
                       Nm,             # Adult and yearling males
                     Mfunction= "min")        # mating function applied
{
#  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimension calculations")
#  nStages <- length(stagenames) 
  #
#  Fmat <- matrix(0, nrow=(nStages), ncol= (nStages))   # blank matrix with dim of stages
#  rownames(Fmat) <- stagenames
#  colnames(Fmat) <- stagenames
  
  K <- params$rep_K
  h <- params$h
  Sf<- params$Sc_f_max
  Sm<- params$Sc_m_max 
  
  #Harmonic mean mating function
  if(Mfunction %in% c( "HM", "HMMF", "Harmonic Mean Mating Function", "harmonic mean mating function", "Harmonic Mean", "harmonic mean")){
  U <- (2* Nf * Nm*h)/ (Nf + Nm* h) # females in unions given pop structure
  
 #   f <- (K*Nm)/(Nm+Nf*h^-1)    # fertility function= cubs produced per adult female - Necessary?
  
  # Minimum mating function
  } else if(Mfunction %in% c("Minimum Mating Function", "minimum mating function", "Minimum", "minimum", "min", "Min")){
   # defining number of pairs formed (U)
    if(Nf < Nm*h){
    U <- Nf
   } else if(Nm*h < Nf){
     U <- Nm*h
   } #fertility coefficient
#    f<- (K*U)/ Nf  # f= cubs produced by adult female
     
     # Mod Harmonic mean mating function?
  } else if(Mfunction %in% c("Modified Hamonic Mean Mating Function", "modified harmonic mean mating function", "ModHarmonic", "modharmonic"))   {
    if(Nf < (2* Nf*Nm*h) / (Nf +Nm*h)) {
      U <- Nf
    } else if((2* Nf*Nm*h) / (Nf +Nm*h) < Nf){
      U <- (2* Nf*Nm*h) / (Nf +Nm*h)
    }
  }
#  Fmat[1,nStages/2] <- 0.5*f * Sf    # female cub production by females
#  Fmat[(nStages/2)+1, nStages/2] <- 0.5*f* Sm  # male cub production by females

   return(U) 
}
# Function name= Mating.func
# Inputs:
# params: density dependent parameters, including k (litter size) & h (harem size)
# stagenames Stages in life cycle , where length()= rows of matrix
# Nf: Number of Adult females in population
# Nm:Number of Adult males
# Mfunc : Mating function, can be Harmonic mean, minimum or mod harmonic mean
# 
# Use: number of pairs formed based on chosen mating function
# Returns maximum possible number of pairs formed in pop based on male and female abundance. Cavveats- h may increase with increasing density, then plateau
# For use in matrix projection- use U as Nf instead of vector input


# Creating Fmat - must run second
mating.Fmat <- function(params,     # density dependent parameters
                        stagenames,   # Stages in life cycle graph (single sex)
                        Nf,        # Adult and yearling females
                        Nm,             # Adult and yearling males
                        Mfunction= "min") {
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimension calculations")

 nStages <- length(stagenames) 

 Fmat <- matrix(0, nrow=(nStages), ncol= (nStages))   # blank matrix with dim of stages
 rownames(Fmat) <- stagenames
 colnames(Fmat) <- stagenames
  
  K <- params$rep_K
  h <- params$h
  Sf<- params$Sc_max
  Sm<- params$Sc_max 
  
  #Harmonic mean mating function
  if(Mfunction %in% c( "HM", "HMMF", "Harmonic Mean Mating Function", "harmonic mean mating function", "Harmonic Mean", "harmonic mean")) {
    U <- (2* Nf * Nm*h)/ (Nf + Nm* h) # females in unions given pop structure
    
    # fertility function= cubs produced per adult female 
      f <- (K*Nm)/(Nm+Nf*h^-1)    
    
    # Minimum mating function
  } else if(Mfunction %in% c("Minimum Mating Function", "minimum mating function", "Minimum", "minimum", "min", "Min")) {
    # defining number of pairs formed (U)
    if(Nf < Nm*h) {
      U <- Nf
    } else if(Nm*h < Nf){
      U <- Nm*h
      
    } #fertility coefficient
      f<- (K*U)/ (Nm+Nf*h^-1)  # f= cubs produced by adult female
    
    # Mod Harmonic mean mating function?
   } else if(Mfunction %in% c("Modified Hamonic Mean Mating Function", "modified harmonic mean mating function", "ModHarmonic", "modharmonic"))   {
    if(Nf < (2* Nf*Nm*h) / (Nf +Nm*h)) {
      U <- Nf
    } else if((2* Nf*Nm*h) / (Nf +Nm*h) < Nf){
      U <- (2* Nf*Nm*h) / (Nf +Nm*h)
      
    } #fertility coefficient
     f <- (K*U)/(Nm+Nf*h^-1)
  }
   Fmat[1,nStages/2] <- 0.5* f * Sf    # female cub production by females
   Fmat[(nStages/2)+1, nStages/2] <- 0.5*f* Sm  # male cub production by females
  
  return(Fmat)
  
}
