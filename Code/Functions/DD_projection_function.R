# DD projection function
# 17/11/2025

# WARNING: MUST RUN MATING SYSTEM FUNCTION BEFORE THIS FUNCTION
# Needs completing

dd.proj<-function(# Fmat,  # MAX FERTILITY MATRIX - needed or calculated later?
                  Umat,   # MAX SURVIVAL
                  initial, 
                  params, 
                  stagenames, 
                  time, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="matrix", 
                  Mfunction= "Min",
                  return.vec= TRUE) 
{
  # Time 
  if (time<=1) stop("Time must be a positive integer")
   else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector for each stage abundance")
   else(n0 <- as.numeric(initial))
  
  
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("You must provide the parameters that match your selected function for density-dependent reproduction.")
  
  
  nStages<-length(stagenames)/2      # how many stages
  
  # population size
  if (is.null(memberN)){
    memberN <- 1:length(stagenames)   # all members contribute to pop size
    
  } else if(length(memberN) >= 1){
    memberN <- c(memberN) 
  }
  
  # Calculating Unions with mating function
    Nf <- n0[nStages]       # pulls adult female entry from initial vector
    Nm <- n0[2*nStages]      # adult male in vector (final entry)
    
    # mating func gives initial Fmat for first year
    mating.out <- mating.func(params, stagenames, Nf, Nm, Mfunction,  return.mat= TRUE)  # Npairs and fmat given 
    U <- mating.out$U 
    Fmat <- mating.out$Fmat
    
  # creating Amat - too early?
  Amat <- Fmat+Umat
  
  # Initialize the output for population vector:
  out<- list(pop=vector(), vec=matrix(), mat=matrix())
  
  Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- rep(NA, (time + 1))       # vector to fill with total pop size each year
  colnames(Vec) <- stagenames   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- n0                   
  Pop[1] <- sum(n0)    
  
  # Loop = density dependent matrix application for each year
  for (i in 1:time) {  
    
    # Mating Fmat creation 
      thisNf <- Vec[i,nStages]    # Nf is mid col in Vec matrix
      thisNm <- Vec[i,2*nStages]  # Nm 
       
    # apply mating func to calculate pairs
    thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
    
    thisFmat <-thisMating$Fmat
    thisU <- thisMating$U
      
    # ricker density dependence each year
    thisN <- sum(Vec[i,memberN])  # pop sizes sums row i for cols included in N
    thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
    Pop[i + 1] <- sum(Vec[(i + 1), ])
    # add thisU ?
  }
  
  # out objects
  out$pop <- Pop
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 
  
  out$mat <- thisAmat
  
  return(out)
}

# Notes----
# Function name- DD.proj
# Inputs
# Fmat: blank matrix that is filled using mating system function for vital rates based on pop structure
# Umat: max predicted survival rates
# initial: population vector for abundance of each stage class
# params: defined in other funcs, importantly b, 
# stagenames: vec of stage classes and sex (2 stage classes should have stagenames length 6)
# time: projection interval
# memberN= which stages included in total pop size count? Defautls NULL, all individuals included. Opts- c(a,b,c) will select those entries of initial vec
# AdultNx= are adults/ oldest class the only individuals included in Nf and Nm? If elderly incl and yearlings not, tweaks needed. TRUE default includes only oldest age class, FALSE incl all stages except first 
# DDapply= how is density dep applied (matrix, fertility, survival, recruitment)
# return.vec= abundance stage class matrix returned? defualts FALSE





# Function to plot projection output
dd_plot <- function(out, 
                    y_val= "N", 
                    ylab = "abundance", 
                    xlab = "time (t)",
                    col= "black",
                    legend.pos = "topright",
                    cex.legend = 0.8) # plot pop size over time, structure (or U?)
  {
  Vec <- out$vec
  Pop <- out$pop
  time_vec <- c(0: nrow(Vec)-1)    # time = 0:time
 
  if(y_val %in% c("N", "Pop Size", "pop size", "Pop", "pop")) {
  matplot(Pop,      # transposing the matrix to plot col on x and abundance on y
          type="l",  # type l for line
          lty= 1,
          col=col, 
          ylab = ylab, 
          xlab = xlab)

  }else if(y_val %in% c("Vec", "Pop Structure", "Stages", "vec")){
    # creating line and colour vectors based on sex and age
    lty_vec<- rep(seq_len(ncol(Vec)/2), times=2)   # lty = age classes
    col_vec<- c(rep(col[1], ncol(Vec)/2), rep(col[2], ncol(Vec)/2))  # repeats colour for stages so each colour = each sex
    
    # 
    matplot(Vec,      # transposing the matrix to plot col on x and abundance on y
            type="l",  # type l for line
            lty= lty_vec,
            col=col_vec, 
            ylab = ylab, 
            xlab = xlab)
    
    legend(legend.pos, 
           legend= colnames(Vec),  
           cex=0.65,
           lty= lty_vec,  # line type matches graph
           col=col_vec)            # colours match graph
    
  }

}
# Notes----
# Function name = dd_plot
# inputs out, 
# y_val= "N", method to plot, pop size vec or age structured matrix (vec)
# ylab = "abundance", labs
# xlab = "time (t)",   
# col= "black",    colours used
# legend.pos = "topright",   legend position on graph
# cex.legend = 0.8) # plot pop size over time, structure (or U?) legend size
# 
# improvements = should use points plus line to join- see dots on ggplot but not function plot. 
# Add theme for plot
# Title? 