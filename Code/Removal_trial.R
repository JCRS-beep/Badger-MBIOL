# Removal methods
# 31.12.2025

# initial structure = yf, af, ym, am (25, 10, 25, 10)
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
n0 <- c(25, 10, 25, 10)
# Creating my control - population projected over 100 years
Umat <- matrix(0, nrow=4, ncol=4)
Umat[2,1]<- 0.67   # yearling f survival
Umat[2,2]<- 0.78   # adult f survival
Umat[4,3]<- 0.64   # yearling m survival
Umat[4,4]<- 0.70

params<- data.frame(      # dataframe 
  fmax= 3.2,     # F fecundity max (max cubs per adult female) 
  Sc_f_max=0.65,   # max cub survival (equal for sexes)
  Sc_m_max=0.65,
  b=0.003,       # temp value- must be calculated from provided datasets
  rep_K= 4,          #litter size (K)
  h= 6   # harem size per male
)

# Setting up with projection function
proj1 <- dd.proj(Umat,      # seems to reach stability quickly - some kind of stocasticity needed?
                 initial = n0, 
                 params, 
                 stagenames = stages, 
                 time = 100, 
                 memberN=NULL,  # which individuals contribute to pop size? (as vec)
                 DDapply="Fmat", 
                 Mfunction= "Min",
                 return.vec= TRUE) 

col_vec <- c("#FF6A6A", "#87CEEB")

(proj_plot <- dd_plot(proj1, 
                     y_val= "Vec", 
                     ylab = "Abundance", 
                     xlab = "Time (t)",
                     theme = theme_classic(), 
                     cols= col_vec,    # can be vector of cols
                     legend.pos = "topright",
                     cex.legend = 0.8))
# reaches plateau around 25 years - removal at year = 25?

# First removal strat = calculate X% of initial pop (r), in single year remove Nr? ----
# Must write removal funcs similar to DDproj
rem.proj <- function(Umat,   # MAX SURVIVAL
                     initial, 
                     params, 
                     stagenames, 
                     time, 
                     memberN=NULL,  # which individuals contribute to pop size? (as vec)
                     DDapply="matrix", 
                     Mfunction= "Min",
                     intensity= 50,  # percentage you want REMOVED from pop at time T=ry
                     remyear = 25,  # removal year = decrease from following year
                     return.vec= TRUE) 
{
  # Time 
  if (time<=1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  
  nStages<-length(stagenames)/2      # how many stages
  ry <- as.numeric(remyear)           # shortening name for future use
  
  # population size
  if (is.null(memberN)){
    memberN <- 1:length(stagenames)   # NULL = all members contribute to pop size
    
  } else if(length(memberN) >= 1){
    memberN <- c(memberN) 
  }
  
  # Calculating Unions with mating function
  Nf <- n0[nStages]     # pulls yearling and adult female entry from initial vector
  Nm <- n0[2*nStages]      # adult male in vector (final entry)
  
  # mating func gives initial Fmat for first year
  mating.out <- mating.func(params, stagenames, Nf, Nm, Mfunction,  return.mat= TRUE)  # Npairs and fmat given 
  U <- mating.out$U 
  Fmat <- mating.out$Fmat

  # Initialize the output for population vector:
  out<- list(pop=vector(), vec=matrix(), mat=matrix())
  
  Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- rep(NA, (time + 1))       # vector to fill with total pop size each year
  colnames(Vec) <- stagenames   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- n0                   
  Pop[1] <- sum(n0)    
  
  # Loop = density dependent matrix application for each year
  for (i in 1:ry) {   # normal projection until remyear (if removal at t=5, project until t=5, final entry inserted to row 6) remember vec[5,] holds entries for year=4 (0:time)
    
    # Mating Fmat creation 
    thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
    thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
    
    # apply mating func to calculate pairs
    thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
    
    thisFmat <-thisMating$Fmat
    thisU <- thisMating$U
    
    # ricker density dependence each year
    thisN <- sum(Vec[i,memberN])  # pop sizes sums row i for cols included in N
    thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
    Pop[i + 1] <- sum(Vec[(i + 1), ])
  }
 # ------------------------------------------------
  # removal year 
  Rem <- Pop[ry + 1] * intensity/100   # ry +1 IS ROW REPRESENTING YEAR remyear! Nremoved by calculating latest N, multiplying by percentage 
  
 # rem_vec <- rep(NA, length(stagenames))  # generate number vec length= stages # sum(rem_vec) == Rem
  
 # sample_vec <- c(1:Rem)
  # for (r in 1:Rem){    # sampling loop until numbers generated
  #  rem_vec <- sample(sample_vec, size = length(stagenames), replace = TRUE)
  # }   
  
  # for now, simply divide rem/4 and remove from each class
  rem_vec <- rep(Rem/4, length(stagenames))
  
  Vec[ry + 2,] <- Vec[ry + 1,]  - rem_vec  # year after remyear = 2 rows later filled with new stage vec
  Pop[ry + 2] <- Pop[ry+1] - Rem  # remyear starting pop = prev N - Nremoved. *edit after*
# -------------------------------------------------

# repeat loop after removal year
  for (j in (ry + 2):time) {    # starts filling from 27th row (26th year)
    
    # Mating Fmat creation 
    thisNf <- sum(Vec[j,nStages-1], Vec[j,nStages-1])    # Nf is mid col in Vec matrix
    thisNm <- sum(Vec[j,2*nStages-1],Vec[j,2*nStages])   # Nm 
    
    # apply mating func to calculate pairs
    thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
    
    thisFmat <-thisMating$Fmat

    # ricker density dependence each year
    thisN <- sum(Vec[j,memberN])  # pop sizes sums row i for cols included in N
    thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
    
    Vec[(j + 1), ] <- thisAmat %*% Vec[j, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
    Pop[j + 1] <- sum(Vec[(j + 1), ])
  }
  
  # out objects
  out$pop <- Pop
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 
  
  out$mat <- thisAmat
  
  return(out)
}
# SUCESS !!

proj2 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stocasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages, 
                  time = 100, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="Fmat", 
                  Mfunction= "Min",
                  intensity= 50,  # percentage you want REMOVED from pop at time T=ry
                  remyear = 25, 
                  return.vec= TRUE) 

(proj2_plot <- dd_plot(proj2, 
                     y_val= "Vec", 
                     ylab = "Abundance", 
                     xlab = "Time (t)",
                     theme = theme_classic(), 
                     cols= col_vec,    # can be vector of cols
                     legend.pos = "topright",
                     cex.legend = 0.8))
# NOTES - population very quickly returns to optima - some environmental stochasticity needed?

proj3 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stocasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages, 
                  time = 100, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="Fmat", 
                  Mfunction= "Min",
                  intensity= 90,  # percentage you want REMOVED from pop at time T=ry
                  remyear = 25, 
                  return.vec= TRUE) 

(proj3_plot <- dd_plot(proj3, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       cex.legend = 0.8))


# Second removal strat = calculate X% of initial pop (Nr), over 3 years remove Nr/3?

# Third removal strat = calculate X% of initial pop (Nr), over y years remove r until Nt+y = (1-X)* N0


