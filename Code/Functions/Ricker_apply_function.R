# Ricker function and matrix application
# 16/11/2025


# Applying ricker function density dependence to matrix
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