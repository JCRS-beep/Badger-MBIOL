# Density dependence and mating systems
# 10/11/2025

# TRIAL DENSITY DEPENDENCE----
# Incorporating mating systems and DD reproduction

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


# creating Fmat
create.Fmat<- function(params,     # density dependent parameters
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
  
Fmat2<- create.Fmat(params=params, nStages=3, Nf=20, Nm=18)  # starting pop = 20 fems, 18 males
Fmat2   # SUCCESS!

Fmat.d<- create.Fmat(params= params, nStages=3, Nf=110, Nm=60)   # values from L. M. Rogers et al (1997)
Fmat.d  # female skew

Fmat.l<- create.Fmat(params= params, nStages=3, Nf=34, Nm=23)  
Fmat.l

# Ricker function----
ricker<- function(params, N){
  dd_fun<-exp(-params$b*N)
  return(dd_fun)      # returns value of multiplier
}

trial<- data.frame(b= 0.02)
ricker_trial<- ricker(params=trial, N=100)
ricker_trial

# DD SURVIVAL ----
# creating Umat, martix with only survival and transitions
ex_Umat<- matrix(0, nrow=6, ncol= 6)   # visualising matrix
ex_Umat[2,1]<- "e^-bN*G1f"    # only cub survival density dependent
ex_Umat[3,2]<- "G2f"
ex_Umat[3,3]<- "Sf"
ex_Umat[5,4]<- "e^-bN*G1m"   
ex_Umat[6,5]<- "G2m"
ex_Umat[6,6]<- "Sm"

Umat<-mat.blank  # copying blank matrix
Umat[2,1]<- 0.32   # only cub survival density dependent
Umat[3,2]<- 0.67
Umat[3,3]<- 0.86
Umat[5,4]<- 0.27 
Umat[6,5]<- 0.65
Umat[6,6]<- 0.77

gN<- ricker(params=params, N=50)
# multiplying cub survival values
Umat[2,1]<- 0.32*gN 
Umat[5,4]<- 0.27*gN  
Umat  

# generalising into a function
create.Umat<- function(params, 
                       nStages,
                       N, 
                       g2f, 
                       g2m,
                       S) {   # Adult survival vector
Umat<- matrix(0, nrow=(2*nStages), ncol= (2*nStages))
b<- params$b
g1f<- params$Sc_f_max
g1m<-params$Sc_m_max
rick<-exp(-b*N)    # ricker function embedded
Umat[2,1]<- g1f*rick   # only cub survival density dependent
Umat[3,2]<- g2f
Umat[3,3]<- S[1]
Umat[5,4]<- g1m*rick
Umat[6,5]<- g2m
Umat[6,6]<- S[2]
return(Umat)
}


# Combining these into a final matrix (Amat)----
# Theoretical pop of 55f, 45m. near carrying capacity, strong density dependence
params_th<- data.frame(
  Ffmax= 1.6,     # F fecundity max (max cubs per adult female)
  Fmmax= 2.1,     # Male fecundity  (max cubs per adult female)
  Sc_f_max= 0.67,   
  Sc_m_max= 0.65, 
  b=0.025,
  rep_K= 3,          #litter size (K)
  h= 6               # harem size per male
)

(Fmat_th<-create.Fmat(params=params_th, nStages=3, Nf=55, Nm=45))
(Umat_th<-create.Umat(params= params_th, nStages=3, N=100, 
                     g2f=0.7, 
                     g2m=0.67, 
                     S=c(0.86, 0.84)))

Amat_th<- Umat_th+ Fmat_th 

#Projecting this over time steps ---- 
