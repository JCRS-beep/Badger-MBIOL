# Matrix creation function- DOESN'T RUN
# 29/10/2025


create.matrix<- function(names,   # vector of stage names    
                         nSex=1,    # 1 or 2 sex matrix (default 1)             
                         r1,             # reproduction sex 1
                         r2= NULL,              # reproduction sex 2
                         g1,   # vector of all growth rates/ transition probabilities
                         g2= NULL, # second sex values for growth
                         S)
                        {         
  nStages<- length(names)                            # number of stages
  mat<- matrix(0, nrow=(nSex*nStages), ncol=(nSex*nStages))     # empty matrix with equal rows and columns as lifestages
  colnames(mat)<- names                                         # naming columns after stages
  rownames(mat)<- names                                         # naming rows after stages

  if(length(r1)==1){    # if only 1 reproductive stage
    mat[1,nStages]<- r1[1]
     
  }
  else{
    mat[1,(nStages-1)]<-  r1[1]   # otherwise assume final stage sexond value forreproduction 
    mat[1, nStages]<- r1[2]
  }
       
    for(n in seq_len(nStages-1))                        # loop for as many lifestages before final (growth)
    {         
  mat[n+1, n]<- g1[n]   # inputs growth rates into respective matrix location
  }   
  mat[(nStages),(nStages)]<- S[1]                      # final stage remaining is adult survival
  
  # this section only runs for 2 sex matrix
   if(nSex==2){   # only runs if 2 sex model
         for(m in (nStages+1):(2*(nStages)+1))                   # loop for m= 5-9
  {         
    mat[m, m-1]<- g2[m-1]  # inputs growth rates into respective matrix location  
  }
     } 
  if(nSex==2) 
    { 
     mat[2*nStages, 2*nStages]<- S[2]
  }
  if(nSex==2){
    mat[(nStages+1), 2*nStages]<- r2[2]    # 2nd sex reproductive rate
  }
  return(mat)  } 


# Function name = 
# Arguments:
# names= vector of lifestage names
# nSex= whether matrix includes 2 sexes, deafults to 1 if no input given. Female assumed as first sex, male second
# r= reproduction, if 2 sex then vector length 2 with (f,m) rates
# r2= reproduction for second sex
# g1 = vector of growth 
# g2 = growth 
# S= Adult survival vector
# Output = mat, a matrix of vital rates for each lifestage

# trial
stageNames<- c("Cub", "Yearling", "Adult", "Elder")
mat.test<- create.matrix(names= stageNames,
                         nSex=1,
                         r1=c(0.48, 0.35), 
                         g1=c(0.29, 0.49, 0.81), 
                         S=0.2)
mat.test    # SUCCESS

# 2 SEX TEST
mat2.test<- create.matrix(names= stageNames,
                          nSex=2,
                          r1=c(0.48, 0.35), 
                          r2=c(0.6, 0.58),
                          g1=c(0.22, 0.3, 0.6), 
                          S=c(0.2, 0.1))
mat2.test

# NOTES----
# Assumes reproductive phases are the oldest stages, and that all others are immature
# Assumes equal number of stages in both sexes if 2 sexes included
# WORKS UP TO 2 REPRODUCTIVE STAGES