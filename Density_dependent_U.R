# Density dependent survival
# 12/11/25

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

# Ricker function----
ricker<- function(params, N){
  dd_fun<-exp(-params$b*N)
  return(dd_fun)      # returns value of multiplier
}



# Applying to matrix with function
apply.DD<- function(params, Fmat, Umat, N, DDapply="matrix")   # apply ricker to matrix, survival or fertility
{ 
  Amat<- Fmat + Umat
  rick<-ricker(params, N)    # ricker function embedded
  if(DDapply %in% c("matrix", "Matrix"))  {
    DDmat<- Amat*rick  
    
  } else if(DDapply %in% c("Survival, survival, Umat")){
    DDmat<- Umat*rick
    
  } else if(DDapply  %in% c("Ferility, fertility, Fmat")) {
    DDmat<- Fmat*rick
  }
  
  return(DDmat) 
}

# Combining these into a final matrix (Amat)----
# Theoretical pop of 55f, 45m. near carrying capacity, strong density dependence


fun.test<-apply.DD(params, Fmat, Umat, 60, DDapply=matrix)
