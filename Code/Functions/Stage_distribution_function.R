# Function- Stage distribution
#24/10/25

# function to calculate stage distribution per year
prop.stage<- function(projected)     # input projected matrix
{ 
  stageMat<- matrix(0, ncol=ncol(projected), nrow=nrow(projected))    # empty matrix to fill
  for(i in 1:ncol(projected)) {   # loop for each column 
    stageMat[,i]<- projected[,i]/sum(projected[,i])          # column i of matrix filled with row i divided by col sum
  } 
  return(stageMat)
}

prop.stage(project1)   # SUCCESS!

# Function name= prop.stage
# Arguments: 
# projected= a projected matrix over a time interval, produced by mat.proj
# Purpose:  Calculates the proportion of each stage class out of the total pop size in a given year


# Updating function to match output of ddproj obj
prop.stage<- function(out) {     # input projected matrix
  # separating a matrix from out obj
  mat <- out$vec  
  stageMat<- matrix(0, ncol=ncol(mat), nrow=nrow(mat))    # empty matrix to fill
  for(i in 1:nrow(mat)) {   # loop for each column 
    stageMat[i,]<- mat[i,]/sum(mat[i,])          # column i of matrix filled with row i divided by col sum
  } 
  return(stageMat)   # returns matrix of each stage as proportion of total pop
}