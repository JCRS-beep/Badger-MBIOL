# Function- projection and plotting this- DOESN'T RUN
# 23/10/25

# Function to compute projection over given time interval - issues with nrow?
mat.proj<- function(N0, mat, time_int){
  out.mat<-matrix(0, nrow=nrow(mat), ncol= (time_int+1)) # creates empty matrix to fill
  rownames(out.mat)<- rownames(mat)   
  colnames(out.mat)<- c(0:time_int)    # col named for years of projection (initial labelled 0)
  out.mat[,1]<- N0                            # first col= N0 
  for(t in 2:(time_int+1))                # starting from second col, until final 
  {
    out.mat[, t]<- mat %*% out.mat[,t-1]
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

# TESTING FUNCTION
n0<- c(10, 10, 10)
names<- c("cub", "yearling", "adult")
# creating matrix 
Fmat<- matrix(0, nrow=3, ncol=3)# blank female (F) matrix of zeros to fill
Fmat[1,3]<- 0.27   # female reproduction
Fmat[2,1]<- 0.65    # cub to yearling 
Fmat[3,]<- c(0, 0.65, 0.86)  # yearling to adult, adult survival as(n-1/n)mat.proj()

rownames(Fmat)<- names
colnames(Fmat)<- names

proj<- mat.proj(n0, Fmat, 20)
proj  # SUCCESS! 




# Function to plot output of mat.proj
plot.proj<- function(pmat,
                     vec_col= "black", 
                     nSex=1,
                     ylab = "abundance", 
                     xlab = "time (t)",
                     legend.pos = "topright", 
                     cex.legend = 0.8) {
  nStages=(nrow(pmat))
  if(nSex==1) 
  {
    lty_vec<- seq_len(nStages)   # single sex, lty = age classes
  }
  else {
    lty_vec<- rep(seq_len(nStages), times=2) # line types repeated for sexes 
  }        
  
  matplot(t(pmat), # transposing the matrix to plot col on x and abundance on y
        type="l",  # type l for line
        lty= lty_vec,
        col=c(vec_col), 
        ylab = ylab, 
        xlab = xlab)
                 
legend("topright", 
       legend= rownames(pmat),  
       cex=0.65,
       lty= lty_vec,  # line type matches graph
       col=vec_col)            # colours match graph
}

# Function name = plot.proj  
# Arguments
# pmat   : matrix output of projection (out.mat)
# v_cols : Colour vector length 2, col1 for sex1
# Purpose= plots the projected matrix with a key corresponding to stage and sex. MUST HAVE EQUAL STAGES FOR EACH SEX

# TEST
plot<- plot.proj(proj)
