# Ricker function and matrix application
# 16/11/2025


ricker <- function(params, N){   # Using params (b) and population size input
  dd_fun <- exp(-params$b*N)
  return(dd_fun)      # returns value of multiplier
}

# Function name= ricker
# Inputs: 
# params, incl b = strength of density dependence
# N= population size (can be total, NAdults, other, but explain in comments) 
# Use: Using given population size, returns dd_fun multiplier.



# Applying ricker function density dependence to matrix
apply.DD <- function(params, Fmat, Umat, N, DDapply="matrix")   # apply ricker to whole matrix, survival or fertility
{ 
  # making sure mats match 
  if (sum(dim(Fmat)==dim(Umat))!=2){
    stop("Reproductive (Fmat) and survival (Umat) matrices 
         must have the same dimensions.")
  }
  rick <- ricker(params, N)    # ricker function embedded
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
  
 } else if(DDapply %in% c("Recruitment", "recruitment")){
   # working out positions of cub in matrix
   nStages <- nrow(Umat)/2  
  # Mcub<- c(nStages+2, nStages+1)         # assumes female cubs are [1,], male cubs [nStages+1, ] 
   Umat_N <- Umat
   Umat_N[2,1] <- Umat[2,1]* rick # female cub survival to yearlings
   Umat_N[(nStages+2), (nStages+1)] <- Umat_N[(nStages+2), (nStages+1)]* rick # male cub survival to yearlings 
   Amat_N <- Fmat + Umat_N
  
   # Applying to Fmat
   Fmat_N <- Fmat*rick
   Amat_N <- Fmat_N + Umat_N
 } 

  return(Amat_N) 
}
# Function name= applyDD
# Inputs: 
# params, incl b for ricker
# N for ricker
# Fmat: matrix wihth reproductive params 
# Umat: 
# DDapply= across which elements ricker is applied. Can include entire matrix, survival (Umat only), Fertility (fmat only) or recruitment (cub survival and Fmat)
# N= population size (can be total, NAdults, other, but explain in comments) 
# Use:  Creates a density dependent Amat depending on DDapplication to existing matrices, Uma and Fmat

# issues in application- recruitment subsetting for male survival incorrect?