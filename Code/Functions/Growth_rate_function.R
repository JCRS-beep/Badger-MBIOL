# Function- Population growth rate 
# 23/10/25

# growth rate calculate with output format in DDproj
growth.rate <- function(out, vis = FALSE, rem_year = NULL){
  N <- out$pop    # isolating pop size vector
  lambda <- N[2:length(N)]/N[1:(length(N)-1)]
  # output proj
  lambda_out <- list(lambda = vector(), lamb = NA)
  lambda_out$lambda <- lambda
  
  if(is.numeric(rem_year)){
  ry <- rem_year
  }
  
  if(!isFALSE(vis)){
  lamb <- plot(x= 1:length(lambda), y= lambda, 
               xlab = "Year", 
               ylab = "lambda") 
    lines(lambda, col= "#94D673")     # keeping cols consistent - red = removals, green = growth rates
    
    title("Pop Growth rate per year", )
    
    lambda_out$lamb <- lamb   # issue here - not loading plot as object within list
  }
  
  return(lambda_out)
}
# Function name: growth.rate
# Arguments    : dd_proj or rem_proj output, vis = TRUE or FALSE for whether you want graph plotted
# Purpose      : generates growth rate per year (gr) by dividing pop size by previous year 
