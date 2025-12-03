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
Umat[2,1]<- 0.67   # yearling f survival
Umat[2,2]<- 0.78   # adult f survival
Umat[4,3]<- 0.65   # yearling m survival
Umat[4,4]<- 0.72

# loading Parameters for density and Fmat creation 
params<- data.frame(      # dataframe 
  fmax= 3.2,     # F fecundity max (max cubs per adult female) 
  Sc_f_max=0.65,   # max cub survival (equal for sexes)
  Sc_m_max=0.65,
  b=0.002,       # temp value- must be calculated from provided datasets
  rep_K= 4,          #litter size (K)
  h= 6   # harem size per male
)


# Fmat example for density dependence
ex_Fmat<- mat.blank
ex_Umat[1,2]<- "0.5*e^-2bN*f"    # female reproduction AND cub survival density dependent. Must apply dependence twice
ex_Umat[3,2]<- "0.5*e^-2bN*f"  


# creating Fmat basic, with max reproductive output
Fmat<- mat.blank # Fmat will be blank, each year of projection reproductive rates are calculated
f<- params$fmax
Sf<- params$Sc_f_max
Sm<- params$Sc_m_max
Fmat[1,2]<- 0.5*f*Sf   # this Fmat has max values of reproduction (freq and density independent)
Fmat[3,2]<- 0.5*f * Sm

# Mating function -----
# can be loaded from script in Function folder
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

# Creating Fmat for max rates WITH MATING FUNCTIONS
Fmat2 <- mating.func(params, stages, Nf= 20, Nm= 20, Mfunction= "modharmonic", return.mat=TRUE)



# Applying to matrix with function applyDD for a static population - 40 adults, 20 females, 20 males -----
Amat<-apply.DD(params, Fmat2$Fmat, Umat, N=40, DDapply="recruitment") # N includes only adults and yearlings, so that N= Nf + Nm
Amat  

# DD projection function ----


# testing projection function----
initial <- c(10, 10, 10, 10)
test_proj <-dd.proj(Umat, 
                    initial, 
                    params, 
                    stagenames = stages, 
                    time = 20, 
                    memberN=NULL,  # which individuals contribute to pop size? (as vec)
                    DDapply= "recruitment", 
                    Mfunction= "min",
                    return.vec= TRUE) 
# SUCCESS? 

 
# Plotting pop strucutre over time
library(ggplot2)
time_vec <- c(0:20)
(Nplot <- ggplot(NULL, aes(x=time_vec, y=test_proj$pop)) + 
                           xlab("time(years)") +
                           ylab("Population size") +
  geom_point())      # nice plateau, but 300 individuals appears high!

(test_plot_N <-dd_plot(test_proj, 
                      y_val= "N", 
                      ylab = "Population size", 
                      xlab = "time (t)",
                      col= "black",
                      legend.pos = "topleft",
                      cex.legend = 0.8))


(test_plot_vec <- dd_plot(test_proj, 
                    y_val= "Vec", 
                    ylab = "abundance", 
                    xlab = "time (t)",
                    col= c("red", "blue"),
                    legend.pos = "topleft",
                    cex.legend = 1))
# why do yearling m and f have identical projections? Why do adults appear much more abudance than yearlings? 
