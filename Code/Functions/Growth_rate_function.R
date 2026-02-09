# Function- Population growth rate 
# 23/10/25

# function to calculate growth rate (+pop size)
growth<- function(projected) 
{popN<- apply(projected, 2, sum) # intermediate pop size calculation. 2= by col
gr<- popN[2:length(popN)]/popN[1:(length(popN)-1)]  # 
return(gr)  
}   

# Function name: growth
# Arguments    : projected matrix (made with mat.proj)
# Purpose      : calculates intermediate population size (popN) for future use
#                generates growth rate per year (gr) by dividing pop size by previous year 

growth(project1)  # SUCCESS!



# creating similar function that works with my output format in DDproj
growth.rate <- function(out){
  N <- out$pop    # isolating pop size vector
  lambda <- N[2:length(N)]/N[1:(length(N)-1)]
  
  return(lambda)
}
