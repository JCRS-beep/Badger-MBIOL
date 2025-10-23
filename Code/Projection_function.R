# Function- projection
# 23/10/25

# writing my own func to compute this ----
mat.proj<- function(N0, mat, time_int){
  out.mat<-matrix(0, nrow=time_int+1, ncol=length(N0)) # creates empty matrix to fill
  rownames(out.mat)<- rownames(mat)   
  colnames(out.mat)<- c(0:length(time_int))  
  out.mat[,1]<- N0                            # first col= N0 
  for(t in 2:(time_int+1))
  {
    out.mat[, t]<- mat%*% out.mat[,t-1]
  }
  return(out.mat)
}  
# Function = mat.proj  
# Arguments
# N0       : initial population structure for each age class  
# mat      : vital rate matrix 
# time_int : (time interval) number of years projected
# Purpose= out.mat is the resulting matrix filled with abundance for age classes by column, with number of columns reflecting number of years projected.
# Can be plotted with mat.plot