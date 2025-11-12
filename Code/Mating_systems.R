# Mating systems
# 10/11/2025

# Incorporating mating systems 

# Seperating Fmat and Umat
stages<- c("Cub_f", "Yearling_f", "Adult_f", "Cub_m", "Yearling_m", "Adult_m")
mat.blank<- matrix(0, ncol=6, nrow=6)
rownames(mat.blank)<- stages
colnames(mat.blank)<- stages

ex_Fmat<- mat.blank   # copying to visualise Fertility matrix (not female mat)
ex_Fmat[1,3]<- "sr*Ff"
ex_Fmat[1,6]<- "sr*Fm"
ex_Fmat[4,3]<-  "(1-sr)*Ff"
ex_Fmat[4,6]<- "(1-sr)*Fm"

# Inspo from ex functions
params<- data.frame(      # dataframe with density dependent parameters
  Ffmax= 1.6,     # F fecundity max (max cubs per adult female)
  Fmmax= 2.1,     # Male fecundity  (max cubs per adult female)
  Sc_f_max=0.65,   # 
  Sc_m_max=0.65, 
  b=0.02,
  rep_K= 3,          #litter size (K)
  h= 6   # harem size per male
  )


# creating Fmat with mating systems
Mating.Fmat<- function(params,     # density dependent parameters
                       nStages,   # Stages in life cycle graph
                       Nf,        # Adult females
                       Nm)        # Adult males
  {
  Fmat<-matrix(0, nrow=(2*nStages), ncol= (2*nStages))   # blank matrix with dim of stages
  p<- Nf/(Nm+Nf)   # sex ratio as prop females
  K<- params$rep_K
  h<- params$h
  Ff<- (K*Nm)/(Nm+Nf*h^-1)    # female fertility function
  Fm<- (K*Nf)/(Nf+Nm*h^-1)    # male fertility function
  Fmat[1,nStages]<-p*Ff      # female cub production by females
  Fmat[1,2*nStages]<- p*Fm   # female cub production by males
  Fmat[(nStages+1),nStages]<- (1-p)*Ff   # male cub production by females
  Fmat[(nStages+1),(2*nStages)]<- (1-p)*Fm  # male cub production by males
  return(Fmat)
  }

# Function name= Mating.Fmat
# Calculates resulting matrix (Fmat) with 

Fmat<- Mating.Fmat(params=params, nStages=3, Nf=20, Nm=18)  # starting pop = 20 fems, 18 males
Fmat   # SUCCESS!

# Do I include Harmonic fertility function and ricker equation to manage matF?




#Projecting this over time steps ---- 

Fmat.proj<-function(initial, params, stagenames, time, memberN, return.vec="FALSE") {
  out<- list(pop=vector(), vec=matrix(), mat=matrix())  # objects returned
  
  # Initial vector:
   n0 <- initial
   
   nStages<-length(stagenames)
   Nf<- n0[3]        # extract Nf and Nm from n0
   Nm<- n0[6]
   
   # Initialize the output for population vector:
   Vec <- matrix(0, ncol = nStages, nrow = time + 1)
   Pop <- rep(NA, (time + 1))
   dimnames(Vec)[[2]] <- stagenames
   dimnames(Vec)[[1]] <- 1:(time+1)
   Vec[1, ] <- n0
   Pop[1] <- sum(n0)
   for (i in 1:time) {  
     thisNf<- Vec[i,(nStages/2)]    
     thisNm<- Vec[i,nStages]        # assumes adult males final entry
     thisFmat<-Mating.Fmat(params, nStages,thisNf, thisNm)   # results in matrix from
     Vec[(i + 1), ] <- thisFmat %*% Vec[i, ]    # new stage abundance vector
     Pop[i + 1] <- sum(Vec[(i + 1), ])        # new total pop size
   }
   
  # out objects
   dim(Pop) <- time + 1
   out$pop <- Pop
   if (return.vec== TRUE) {
     out$vec <- Vec
  } out$mat <- thisFmat
   
   return(out)
     }

# Testing function----
n0=c(10,10,10,10,10,10)  # equal all sexes and ages
