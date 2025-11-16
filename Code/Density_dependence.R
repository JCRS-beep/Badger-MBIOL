# Density dependent recruitment
# 12/11/25

# creating blank base matrix
stages<- c("Cub_f", "Yearling_f", "Adult_f", "Cub_m", "Yearling_m", "Adult_m")
mat.blank<- matrix(0, ncol=6, nrow=6)
rownames(mat.blank)<- stages
colnames(mat.blank)<- stages

# creating Umat, martix with only survival and transitions
ex_Umat<- mat.blank
ex_Umat[2,1]<- "e^-bN*G1f"    # only cub survival density dependent
ex_Umat[3,2]<- "G2f"
ex_Umat[3,3]<- "Sf"
ex_Umat[5,4]<- "e^-bN*G1m"   
ex_Umat[6,5]<- "G2m"
ex_Umat[6,6]<- "Sm"



Umat<-mat.blank  # creating my working matrix. Values from 
Umat[2,1]<- 0.32   # only cub survival density dependent
Umat[3,2]<- 0.67  
Umat[3,3]<- 0.86
Umat[5,4]<- 0.27 
Umat[6,5]<- 0.65
Umat[6,6]<- 0.77

# loading Parameters for density and Fmat creation 
params<- data.frame(      # dataframe 
  fmax= 3.2,     # F fecundity max (max cubs per adult female)
  Sc_f_max=0.65,   # max female  survival
  Sc_m_max=0.65, 
  b=0.002,       # temp value- must be calculated from provided datasets
  rep_K= 3,          #litter size (K)
  h= 6   # harem size per male
)

# creating Fmat basic, with max reproductive output
Fmat<- mat.blank # Fmat will be blank, each year of projection reproductive rates are calculated
f<- params$fmax
Fmat[1,3]<- 0.5*f
Fmat[4,3]<- 0.5*f

Fmat2<- Mating.Fmat(params, stages, Nf=20, Nm=20)  # expected based on N=40
Fmat2

# Ricker function----
# or load from script in Function folder
ricker<- function(params, N){   # Using params (b) and population size input
  dd_fun<-exp(-params$b*N)
  return(dd_fun)      # returns value of multiplier
}

# Matrix Application - can also be loaded from ricker script
apply.DD<- function(params, Fmat, Umat, N, DDapply="matrix")   # apply ricker to whole matrix, survival or fertility
{ 
  # making sure mats match 
  if (sum(dim(Fmat)==dim(Umat))!=2){
    stop("Reproductive (Fmat) and survival (Umat) matrices 
         must have the same dimensions.")
  }
  rick<-ricker(params, N)    # ricker function embedded
  Amat<- Fmat + Umat
  
  # applying based on method
  if(DDapply %in% c("matrix", "Matrix"))  {
    Amat_N<- Amat*rick  
    
  } else if(DDapply %in% c("Survival", "survival", "Umat")){
    Umat_N<- Umat*rick
    Amat_N<- Fmat + Umat_N
    
  } else if(DDapply  %in% c("Ferility", "fertility", "Fmat")) {
    Fmat_N<- Fmat*rick
    Amat_N<- Fmat_N + Umat
    
  } else if(DDapply %in% c("Recruitment", "recruitment")){
    # working out positions of cub in matrix
    nStages<- nrow(Umat)/2  
    # Mcub<- c(nStages+2, nStages+1)         # assumes female cubs are [1,], male cubs [nStages+1, ] 
    Umat_N<- Umat
    Umat_N[2,1]<- Umat[2,1]* rick # female cub survival to yearlings
    Umat_N[(nStages+2), (nStages+1)]<- Umat_N[(nStages+2), (nStages+1)]* rick # male cub survival to yearlings  INCORRECT NUMBER OF DIMS?
    Amat_N<- Fmat + Umat_N
    
    # Applying to Fmat
    Fmat_N<- Fmat*rick
    Amat_N<- Fmat_N + Umat_N
  } 
  
  return(Amat_N) 
}


# Mating function -----
# can be loaded from script
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
  Fmat[(nStages/2)+1, nStages/2]<- 0.5*f   # male cub production by females
  
  return(Fmat)
}

# Applying to matrix with function applyDD for a static population - 40 adults, 20 females, 20 males -----
Amat<-apply.DD(params, Fmat2, Umat, N=40, DDapply="recruitment") # N includes only adults and yearlings, so that N= Nf + Nm
Amat  

# Trialing projection ----
DD.proj<-function(Fmat, 
                  Umat, 
                  initial, 
                  params, 
                  stagenames, 
                  time, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  AdultNx= TRUE,  # are adults the only individuals included in Nf and Nm calcs?
                  return.vec= FALSE) 
  {
  # Time 
  time <- as.integer(time)
  n0 <- as.numeric(initial) 
  
  nStages<-length(stagenames)/2      # dims is nrow and ncol of Fmat 
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimension calculations")
  
  # population size
  if (is.null(memberN)){
      memberN<- 1:length(stagenames)   # all members contribute to pop size
      
  } else if(length(memberN >= 1)){
    memeberN<- c(memberN) 
  }
  
  # Nf and Nm from n0. adult only vs yearlings included
  if(AdultNx==TRUE){ # if adults only, final entries for each sex = Nx
  Nf<- n0[nStages]      
  Nm<- n0[2*nStages]
  
  } else if(AdultNx==FALSE){   # if multiple stages included, sum all stages except cubs (n0[1] and n0[nStages+1])
    Nf<- sum(n0[c(2:nStages)])      
    Nm<- sum(n0[c(nStages+2:2*nStages)]) 
      
  }
  
  # creating Amat
  Amat<- Fmat+Umat
  
  # Initialize the output for population vector:
  Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance. time= row
  Pop <- rep(NA, (time + 1))       # vector to fill with total pop size each year
  colnames(Vec) <- stagenames   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection
  Vec[1, ] <- n0                   
  Pop[1] <- sum(n0)    
  
  # Looping density dependent matrix application for each year
  for (i in 1:time) {  
    # Mating Fmat creation per year
    if(AdultNx == TRUE) {   # if only adults
    thisNf<- Vec[i,nStages]    # Nf is mid col in Vec matrix
    thisNm<- Vec[i,2*nStages]       
    } else if(AdultNx == FALSE){
      thisNf<- sum(Vec[i,c(2:nStages)])   # all fems excluding cubs
      thisNm<- sum(Vec[i,c(nStages+2:2*nStages)]) # all males excl cubs- ISSUES IN THIS SYNTAX, SUBSCRIPT OUT OF BOUNDS
    }
    
    thisFmat<-Mating.Fmat(params, dims,thisNf, thisNm)   # applies Fmat creation to calculated Nx 
    
    # ricker density dependence each year
    thisN<- sum(Vec[i,memberN])  # pop sizes sums row i for cols included in N
    thisAmat<- apply_Dfactor(thisFmat, Umat, thisN, params, apply.DD)
  }
  
  # out objects
  out$pop <- Pop
  if (isTRUE(return.vec)) {
    out$vec <- Vec        # issues in this line?
  } 
  out$mat <- thisFmat
  
  return(out)
}

# testing projection function----
initial<- c(10, 10, 10, 10, 10, 10)
test_proj<-DD.proj(Fmat, 
                   Umat, 
                   initial, 
                   params, 
                   stagenames= stages, 
                   time=5, 
                   memberN=c(2,3,5,6),  # which individuals contribute to pop size? (as vec)
                   AdultNx= FALSE,  # including yearlings and adults
                   return.vec= TRUE) 
  


# test syntax -----
a<- c(1,2,3,1,2,3)
m<- matrix(0, ncol=length(a), nrow=length(a))
m[1,]<- a
b<- sum(m[1,c(1,2,6)])
c<- sum(a[c(1,3)])
c