# Mating systems
# 10/11/2025

# Incorporating mating systems 

# Seperating Fmat and Umat
stages<- c( "Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
mat.blank<- matrix(0, ncol=4, nrow=4)
rownames(mat.blank)<- stages
colnames(mat.blank)<- stages

ex_Fmat<- mat.blank   # copying to visualise Fertility matrix (not female mat)
# Assumes females only reproduction
ex_Fmat[1,3]<- "0.5f"      # where f= litter size per female for N males and females (Nm, Nf), multiplied cub survival
ex_Fmat[4,3]<-  "0.5f"


# Inspo from ex functions
params<- data.frame(      # dataframe with density dependent parameters
  fmax= 3*0.65,     # fecundity max (max cubs per adult female=3, 0.65= max cub survival)
  Sc_f_max=0.65,   # 
  Sc_m_max=0.65, 
  b=0.02,
  rep_K= 3,          #litter size (K)
  h= 6   # harem size per male
  )

# Creating mating function.  Uses number of pairs for max number of possible unions (if no opposite sex, no mating) -----
mating.func <- function(params,     # density dependent parameters
                        stagenames,   # Stages in life cycle graph (single sex)
                        Nf,        # Adult and yearling females
                        Nm,             # Adult and yearling males
                        Mfunction= "min", 
                        return.mat= FALSE)        # mating function applied
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
    U <- (2* Nf * Nm*h)/ (Nf + Nm* h) # females in unions given pop structure
    f <- K*Nm/ (Nm + Nf*h^-1)
    
    # Minimum mating function
  } else if(Mfunction %in% c("Minimum Mating Function", "minimum mating function", "Minimum", "minimum", "min", "Min")){
    # defining number of pairs formed (U)
    if(Nf < Nm*h){
      U <- Nf
    } else if(Nm*h < Nf){
      U <- Nm*h
    } #fertility coefficient
    f<- (K*U)/ (Nf+ Nm)  # f= cubs produced by adult female
    
    # Mod Harmonic mean mating function?
  } else if(Mfunction %in% c("Modified Hamonic Mean Mating Function", "modified harmonic mean mating function", "ModHarmonic", "modharmonic"))   {
    if(Nf < (2* Nf*Nm*h) / (Nf +Nm*h)) {
      U <- Nf
    } else if((2* Nf*Nm*h) / (Nf +Nm*h) < Nf){
      U <- (2* Nf*Nm*h) / (Nf+ Nm)
    }
    f <- (K*U)/ (Nf+ Nm)}
  
  if(return.mat==TRUE) {
    nStages <- length(stagenames)/2 # number stages for each sex
    Fmat[1,nStages] <- 0.5* f* Sf    # female cub production by females
    Fmat[nStages+1, nStages] <- 0.5* f* Sm  # male cub production by females
    out$Fmat <- Fmat
  }
  
  out$U <- U
  
  return(out) 
}


# Testing function ----
pairs <- mating.func(params,     # density dependent parameters
                     stagenames= stages,   # Stages in life cycle graph (single sex)
                     Nf=25,        # Adult and yearling females
                     Nm=12,             # Adult and yearling males
                     Mfunction= "min", 
                     return.mat=TRUE)        # mating function applied)  # starting pop = 10 fems, 8 males
pairs   # SUCCESS