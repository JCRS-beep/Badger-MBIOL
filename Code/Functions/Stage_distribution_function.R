# Function- Stage distribution
#24/10/25

# function to calculate stage distribution per year
prop.stage<- function(projected)     # input projected matrix
{ 
  popStage<- matrix(0, ncol=ncol(projected), nrow=nrow(projected))    # empty matrix to fill
  for(i in 1:ncol(projected)) {   # loop for each column 
    popStage[,i]<- projected[,i]/sum(projected[,i])          # column i of matrix filled with row i divided by col sum
  } 
  return(popStage)
}

prop.stage(project1)   # SUCCESS!

# Function name= prop.stage
# Arguments: 
# projected= a projected matrix over a time interval, produced by mat.proj
# Purpose:  Calculates the proportion of each stage class out of the total pop size in a given year

