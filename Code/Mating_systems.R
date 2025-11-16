# Mating systems
# 10/11/2025

# Incorporating mating systems 

# Seperating Fmat and Umat
stages<- c("Cub_f", "Yearling_f", "Adult_f", "Cub_m", "Yearling_m", "Adult_m")
mat.blank<- matrix(0, ncol=6, nrow=6)
rownames(mat.blank)<- stages
colnames(mat.blank)<- stages

ex_Fmat<- mat.blank   # copying to visualise Fertility matrix (not female mat)
# Assumes females only reproduction
ex_Fmat[1,3]<- "0.5f"      # where f= litter size per *adult* female for given number males and females (Nm, Nf)
ex_Fmat[4,3]<-  "0.5f"


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
                       stagenames,   # Stages in life cycle graph (single sex)
                       Nf,        # Adult and yearling females
                       Nm)        # Adult and yearling males
  {
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimension calculations")
   nStages<- length(stagenames) 
     
  Fmat<-matrix(0, nrow=(nStages), ncol= (nStages))   # blank matrix with dim of stages
  rownames(Fmat)<- stagenames
  colnames(Fmat)<- stagenames
  
  K<- params$rep_K
  h<- params$h
  f<- (K*Nm)/(Nm+Nf*h^-1)    # fertility function
  Fmat[1,nStages/2]<-0.5*f      # female cub production by females
  Fmat[(nStages/2+1),nStages/2]<- 0.5*f   # male cub production by females
  return(Fmat)
  }

# Function name= Mating.Fmat
# Calculates resulting matrix (Fmat) using harmonic fertility function  

Fmat<- Mating.Fmat(params=params, stagenames=stages, Nf=20, Nm=18)  # starting pop = 20 fems, 18 males
Fmat   # SUCCESS!



#Can't project without Amat= Fmat + Umat, otherwise extinction ---- 
Fmat.proj<-function(initial, params, stagenames, time, return.vec= FALSE) {
 
   out<- list(pop=vector(), vec=matrix(), mat=matrix())  # objects returned
 
   # Time as integer
  time <- as.integer(time)
  # Initial vector:
  dims<-length(stagenames)      # dims is nrow and ncol of Fmat 
  if (length(initial) != dims) stop("initial pop vector must equal length(stagenames)")
  n0 <- as.numeric(initial)     # no intial pop vector
   
  
#   if (is.null(memberN)){
#     memberN<- 1:length(stagenames)   # all members contribute to pop size
#   }
   Nf<- n0[(dims/2)-1] + n0[(dims/2)]      # extract Nf and Nm from n0. adult and yearlings included, only cubs excluded
   Nm<- n0[dims]
   
   # Initialize the output for population vector:
   Vec <- matrix(0, ncol = dims, nrow = time + 1)  # matrix to fill with stage abundance for each time step (each row)
   Pop <- rep(NA, (time + 1))       # vector to fill with total pop size each year
   colnames(Vec) <- stagenames   # naming cols matrix as stages 
   rownames(Vec) <- 0:(time)   # rows correspond to each year of projection
    Vec[1, ] <- n0                   
    Pop[1] <- sum(n0)    
 
     for (i in 1:time) {  
     thisNf<- Vec[i,dims/2]    # Nf is mid col in Vec matrix
     thisNm<- Vec[i,dims]        # assumes adult males final entry in matrix
    
     thisFmat<-Mating.Fmat(params, dims,thisNf, thisNm)   # applies Fmat creation to 
     Vec[(i + 1), ] <- thisFmat %*% Vec[i, ]    # mat multiplication of Fmat rates and prev year stage abundance gives new stage abundance
     Pop[i + 1] <- sum(Vec[(i + 1), ])        # new total pop size sums stage abundance
   }
   
  # out objects
   dim(Pop) <- time + 1
   out$pop <- Pop
   if (isTRUE(return.vec)) {
     out$vec <- Vec        # issues in this line?
  } 
   out$mat <- thisFmat
   
   return(out)
     }

# Testing function----
initial=c(10,10,10,10,10,10)  # equal all sexes and ages
proj.test<- Fmat.proj(initial, params, stagenames= stages, time=10, return.vec= TRUE)
# 'Error in thisFmat %*% Vec[i, ] : non-conformable arguments'   

