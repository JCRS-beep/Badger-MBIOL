# Density dependent recruitment
# 12/11/25

# creating blank base matrix
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
mat.blank<- matrix(0, ncol=4, nrow=4)
rownames(mat.blank)<- stages
colnames(mat.blank)<- stages

# creating Umat, martix with only survival and transitions
ex_Umat<- mat.blank
ex_Umat[2,1]<- "G2f"
ex_Umat[2,2]<- "Sf"
ex_Umat[4,3]<- "G2m"
ex_Umat[4,4]<- "Sm"

Umat<-mat.blank  # creating my working matrix. Values from 
Umat[2,1]<- 0.67  
Umat[2,2]<- 0.78
Umat[3,3]<- 0.65
Umat[4,4]<- 0.72

# loading Parameters for density and Fmat creation 
params<- data.frame(      # dataframe 
  fmax= 3.2*0.65,     # F fecundity max (max cubs per adult female) * cub survival
  Sc_max=0.65,   # max cub survival (equal for sexes)
  b=0.002,       # temp value- must be calculated from provided datasets
  rep_K= 3,          #litter size (K)
  h= 6   # harem size per male
)


# Fmat example for density dependence
ex_Fmat<- mat.blank
ex_Umat[1,2]<- "0.5*e^-2bN*f"    # female reproduction AND cub survival density dependent. Must apply dependence twice
ex_Umat[3,2]<- "0.5*e^-2bN*f"  


# creating Fmat basic, with max reproductive output
Fmat<- mat.blank # Fmat will be blank, each year of projection reproductive rates are calculated
f<- params$fmax
Sf<- params$Sc_max
Sm<- params$Sc_max
Fmat[1,2]<- 0.5*f*Sf   # this Fmat has max values of reproduction (freq and density independent)
Fmat[4,2]<- 0.5*f * Sm

# Mating function -----
# can be loaded from script in Function folder
mating.func <- function(params,     # density dependent parameters
                        #   stagenames,   # Stages in life cycle graph (single sex)
                        Nf,        # Adult and yearling females
                        Nm,             # Adult and yearling males
                        Mfunction= "HM")        # mating function applied
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
  
  # incorporating in projection- using U instead of Nf in population size
  
  return(U) 
}
 
# Ricker function----
# or load from script in Function folder ----
ricker<- function(params, N){   # Using params (b) and population size input
  dd_fun<-exp(-params$b*N)
  return(dd_fun)      # returns value of multiplier
}

# Matrix Application loaded from ricker script----
apply.DD <- function(params, Fmat, Umat, N, DDapply="matrix")   # apply ricker to whole matrix, survival or fertility
{ 
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
    
  } else if(DDapply  %in% c("Ferility", "fertility", "Fmat")) {
    Fmat_N <- Fmat*rick
    Amat_N <- Fmat_N + Umat
    
  } else if(DDapply %in% c("Recruitment", "recruitment")) {
    # Applying to Fmat TWICE (survival and fertility)
    Fmat_N <- Fmat*(rick^2)
    Amat_N <- Fmat_N + Umat
  } 
  
  return(Amat_N) 
}



# Applying to matrix with function applyDD for a static population - 40 adults, 20 females, 20 males -----
Amat<-apply.DD(params, Fmat, Umat, N=40, DDapply="recruitment") # N includes only adults and yearlings, so that N= Nf + Nm
Amat  

# DD projection function ----


# testing projection function----
initial <- c(10, 10, 10, 10)
test_proj <-dd.proj( Fmat, 
                     Umat, 
                     initial, 
                     params, 
                     stagenames = stages, 
                     time = 20, 
                     memberN=NULL,  # which individuals contribute to pop size? (as vec)
                     DDapply= recruitment, 
                     Mfunction= "Min",
                     return.vec= FALSE) 
# SUCCESS!!