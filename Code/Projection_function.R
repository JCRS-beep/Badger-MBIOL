# Function- projection and plotting this
# 23/10/25

# Function to compute projection over given time interval ----
mat.proj<- function(N0, mat, time_int){
  out.mat<-matrix(0, nrow=time_int+1, ncol=length(N0)) # creates empty matrix to fill
  rownames(out.mat)<- rownames(mat)   
  colnames(out.mat)<- c(0:length(time_int))   # columns named for years of projection (initial labelled 0)
  out.mat[,1]<- N0                            # first col= N0 
  for(t in 2:(time_int+1))
  {
    out.mat[, t]<- mat%*% out.mat[,t-1]
  }
  return(out.mat)
}  
# Function name = mat.proj  
# Arguments
# N0       : initial population structure for each age class  
# mat      : vital rate matrix 
# time_int : (time interval) number of years projected
# Purpose= out.mat is the resulting matrix filled with abundance for age classes by column, with number of columns reflecting number of years projected.
# Can be plotted with plot.plot function below
 


# Function to plot output of mat.proj----
plot.proj<- function(pmat, vec_col){
matplot(pmat, type="l",  # type l for line
        col=c(vec_col),  # 2 cols for sexes- NEEDS EDITING
        lty= 1:(nrow(pmat))/2,           # line types repeated for sexes to match age classes
        ylab= "abundance", xlab="time (t)")
legend("topright", legend= rownames(pmat),  
       cex=0.65,
       lty= 1:(nrow(pmat))/2,  # line type matches graph
       col=vec_col)            # colours match graph
     }
# Function name = plot.proj  
# Arguments
# pmat   : matrix output of projection (out.mat)
# v_cols : Colour vector length 2, col1 for sex1
# Purpose= plots the projected matrix with a key corresponding to stage and sex. MUST HAVE EQUAL STAGES FOR EACH SEX