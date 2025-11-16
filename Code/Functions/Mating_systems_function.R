# Mating systems Fmat function
# 16/11/2025


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
# Function name= Mating.Fmat
# Inputs:
# params: density dependent parameters, including k (litter size) & h (harem size)
# stagenames Stages in life cycle , where length()= rows of matrix
# Nf: Number of Adult and yearling females in population
# Nm:Number of Adult and yearling males
# 
# Use: Creates a fertility matrix of expected cub production per female for male and female cubs. Ignores male reproduction aside from mating contribution
# f= female fertility function (number of cubs per female given Nf females and Nm males in pop)
# 0.5 assumes sex ratio is equal and density independent, so half the cubs produced are female


