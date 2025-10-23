# Population growth rate - Function
# 23/10/25

# function to calculate growth rate (+pop size)
growth<- function(projected) 
{popN<- apply(projected, 2, sum) # intermediate pop size calculation
gr<- popN[2:length(popN)]/popN[1:(length(popN)-1)]
return(gr)
}   

# Function name: growth
# Arguments    : projected matrix (made with mat.proj)
# Purpose      : caculates intermediate population size (popN) for future use
#                generates growth rate per year (gr) by dividing pop size by previous year 