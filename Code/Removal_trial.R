# Removal methods
# 31.12.2025

# vec structure = yf, af, ym, am (25, 10, 25, 10)
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
n0 <- c(25, 10, 25, 10)

# first time using params from extraction 
# Creating my control - population projected over 100 years
Umat <- matrix(0, nrow=4, ncol=4)
Umat[2,1]<- 0.851   # yearling f survival
Umat[2,2]<- 0.851   # adult f survival
Umat[4,3]<- 0.809   # yearling m survival
Umat[4,4]<- 0.749   # adult m survival

params<- data.frame(fmax= 0.8436,   # F fecundity max (max cubs per adult female) 
                    Sc_f_max=0.65,   # max cub survival (equal for sexes)
                    Sc_m_max=0.65,
                    b=0.004,       # temp value- must be calculated from provided datasets
                    rep_K= 1.8,          #litter size (K)
                    h= 6   # harem size per male
                    )

# Setting up with projection function
proj1 <- dd.proj(Umat,      # seems to reach stability quickly - some kind of stocasticity needed?
                 initial = n0, 
                 params, 
                 stagenames = stages, 
                 time = 50, 
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


(proj_Nplot <- dd_plot(proj1, 
                      y_val= "Pop", 
                      ylab = "Pop size", 
                      xlab = "Time (t)",
                      theme = theme_classic(), 
                      cols= ,    # can be vector of cols
                      legend.pos = "topright",
                      cex.legend = 0.8))   # why doesn't this run? 

# reaches plateau around 20 years - removal at year = 25?

# First removal strat = generate normal dist of 0 - intensity. Mean = intensity. sample 4 numbers and remove that proportion of stage
# generating rnorm
intensity = 50
prop <- rnorm(4, mean = intensity/100, sd= 0.05)  # 4 samples from dist mean 0.5, sd 0.05
removed <- proj1$vec* prop
remaining <- proj1$vec - removed 



#removal funcs similar to DDproj - replacement? Need to stop simulation once pop reaches 0 - in some cases becomes negative 
rem.proj <- function(Umat,   # MAX SURVIVAL
                     initial, 
                     params, 
                     stagenames, 
                     time, 
                     memberN=NULL,  # which individuals contribute to pop size? (as vec)
                     DDapply="matrix", 
                     Mfunction= "Min",
                     #   removal = FALSE,  # allows function versitility, proj without removals
                     intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                     remyear = NULL,  # removal year = decrease from following year
                     return.vec= TRUE, 
                     return.remvec = FALSE) 
{
  # Time 
  if (time<=1) stop("Time must be a positive integer")
  else(time <- as.integer(time))
  
  if (is.null(initial)) stop("You must provide initial population vector with each stage abundance")
  else(n0 <- as.numeric(initial))
  
  
  if (length(initial) != length(stagenames)) stop("initial pop vector must equal length(stagenames)")
  if(is.null(stagenames) || length(stagenames)== 0) stop("stagenames must be provided for correct matrix dimensions")
  
  if (is.null(params)) stop("You must provide parameters for selected density-dependent function")
  
  
  nStages <- length(stagenames)/2      # how many stages
  
  
  # population size
  if (is.null(memberN)){
    memberN <- 1:length(stagenames)   # NULL = all members contribute to pop size
    
  } else if(length(memberN) >= 1){
    memberN <- c(memberN) 
  }
  
  # Calculating Unions with mating function
  Nf <- n0[nStages]     # pulls ONLY adult female entry from initial vector
  Nm <- n0[2*nStages]      # adult male in vector (final entry)
  
  # mating func gives initial Fmat for first year
  mating.out <- mating.func(params, stagenames, Nf, Nm, Mfunction,  return.mat= TRUE)  # Npairs and fmat given 
  U <- mating.out$U 
  Fmat <- mating.out$Fmat
  
  # Set up the output
  out <- list(pop = vector(), 
              vec = matrix(), 
              mat = matrix(), 
              Nremoved = numeric(), # how to leave blank if no removals?
              remvec = vector())  # including number removed and removals from each stage
  
  Vec <- matrix(0, ncol = length(stagenames), nrow = time + 1)  # matrix to fill with stage abundance.  row= time, col= stage
  Pop <- rep(NA, (time + 1))       # vector to fill with total pop size each year
  Nremoved <- 0
  remvec <- rep(NA, length(stagenames))
  
  colnames(Vec) <- stagenames   # naming cols matrix as stages 
  rownames(Vec) <- 0:(time)   # rows correspond to each year of projection. Row 0 = initial or n0
  Vec[1, ] <- n0                   
  Pop[1] <- sum(n0)    
  
  
  if(is.numeric(intensity)){   # removal scenarios only
    ry <- as.numeric(remyear)           # shortening name for future use
  } else if(is.null(intensity)){
    ry <- time                   # no removals = run until time
  }
  # Loop = density dependent matrix application for each year
  for (i in 1:ry) {   # normal projection until remyear (if removal at t=5, project until t=5, final entry inserted to row 6) remember vec[5,] holds entries for year=4 (0:time)
    
    # Mating Fmat creation 
    thisNf <- sum(Vec[i,nStages-1], Vec[i,nStages])    # Nf sums yearling and adult fems in Vec matrix
    thisNm <- sum(Vec[i,2*nStages-1],Vec[i,2*nStages])   # Nm 
    
    # apply mating func to calculate pairs
    thisMating<- mating.func(params, stagenames, thisNf, thisNm, Mfunction, return.mat=TRUE)     # how to use?
    
    thisFmat <- thisMating$Fmat
    thisU <- thisMating$U
    
    # ricker density dependence each year
    thisN <- sum(Vec[i,memberN])  # pop sizes sums row i for cols included in N
    thisAmat <- apply.DD(params, thisFmat, Umat, thisN, DDapply)  # entire pop size used to calculate ricker, apply to recruitment - how to limit births based on U?
    
    # If the projection matrix has any negative values in it, stop iterating but
    # return the projection up until this point.
    if (sum(thisAmat<0)>0){
      warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
      break
    }
    # if pop size <= 0, stop and return
    if(Pop[i]<= 0){
      warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
      break
    }
    # if any stage becomes negative, set to zero and continue
    if (any(Vec < 0)) {
      warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
      Vec[Vec < 0] <- 0
    }
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
    Pop[i + 1] <- sum(Vec[(i + 1), ])
  }
  
  if (is.numeric(intensity)){
    # ------------------------------------------------
    # removal year 
    Rem <- Pop[ry + 1] * intensity/100   # ry +1 IS ROW REPRESENTING YEAR remyear! Nremoved by calculating latest N, multiplying by percentage 
    
    #generating the distribution
    prop <- rnorm(length(stagenames), mean = intensity/100, sd= intensity/1000)  # 4 samples from dist mean 0.5, sd 0.05
    rem_vec <- Vec[ry+1,] * prop    # remaining <- Vec[ry+1] - rem_vec
    
    Vec[ry + 2,] <- Vec[ry + 1,]  - rem_vec  # year after remyear = 2 rows later filled with new stage vec
    Pop[ry + 2] <- sum(Vec[ry + 2]) # filling in total pop size
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
      
      # If the projection matrix has any negative values in it, stop iterating but
      # return the projection up until this point.
      if (sum(thisAmat<0)>0){
        warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
        break
      }
      # if pop size <= 0, stop and return
      if(Pop[i]<= 0){
        warning(paste("Projection stopped at time step", i, "because pop size reached 0 or below"))
        break
      }
      # if any stage becomes negative, set to zero and continue
      if (any(Vec < 0)) {
        warning(paste("Negative abundances produced at time step", i, "setting negatives to 0 and continuing."))
        Vec[Vec < 0] <- 0
      }
      
      Vec[(j + 1), ] <- thisAmat %*% Vec[j, ]  # following year stage vector is this Amat* this year pop structure - incorporate U here for max no. births?
      Pop[j + 1] <- sum(Vec[(j + 1), ])
    }
  } 
  # out objects
  out$pop <- Pop
  if(is.numeric(intensity)){   
    out$Nremoved <- Rem    
    if (isTRUE(return.remvec)) {
      out$remvec <- rem_vec        # returns blank?
    } }
  
  if (isTRUE(return.vec)) {
    out$vec <- Vec        
  } 
  out$mat <- thisAmat
  
  return(out)
}
# testing no removal scenario
proj0 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stocasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages, 
                  time = 30, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="Fmat", 
                  Mfunction= "Min",
                  intensity= NULL,  # percentage you want REMOVED from pop at time T=ry
                  remyear = NULL, 
                  return.vec= TRUE, 
                  return.remvec = FALSE) 

(proj0_plot <- dd_plot(proj0, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       cex.legend = 0.8))

proj2 <- rem.proj(Umat,   # MAX SURVIVAL
                  initial = n0, 
                  params, 
                  stagenames = stages, 
                  time = 20, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="matrix", 
                  Mfunction= "Min",
                  intensity= 50,  # percentage you want REMOVED from pop at time T=ry
                  remyear = 5,  # removal year = decrease from following year
                  return.vec= TRUE, 
                  return.remvec = TRUE) 

(proj2_plot <- dd_plot(proj2, 
                     y_val= "Vec", 
                     ylab = "Abundance", 
                     xlab = "Time (t)",
                     theme = theme_classic(), 
                     cols= col_vec,    # can be vector of cols
                     legend.pos = "topright",
                     cex.legend = 0.8))
# NOTES - population very quickly returns to plateu - some environmental stochasticity needed?
#         Update plot function to draw red dotted line at remyear

proj3 <- rem.proj(Umat,      # seems to reach stability quickly - some kind of stochasticity needed?
                  initial = n0, 
                  params, 
                  stagenames = stages, 
                  time = 40, 
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


proj4 <- rem.proj(Umat,   # MAX SURVIVAL
                  initial = n0, 
                  params, 
                  stagenames= stages, 
                  time = 25, 
                  memberN=NULL,  # which individuals contribute to pop size? (as vec)
                  DDapply="matrix", 
                  Mfunction= "Min",
                  removal = TRUE,  
                  intensity= 98,  
                  remyear = 5,  
                  return.vec= TRUE, # again has negative pop sizes - needs solving
                  return.remvec = FALSE)

(proj4_plot <- dd_plot(proj4, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       theme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       cex.legend = 0.8)) 


















# issues - negative values for abundance should not be possible
##removal funcs similar to DDproj - combine with ddproj instead of replicates? Need to stop simulation once pop reaches 0 - in some cases becomes negative 


proj.test <- rem.proj(Umat,   # MAX SURVIVAL
                      initial = n0, 
                      params, 
                      stagenames = stages, 
                      time = 20, 
                      memberN=NULL,  # which individuals contribute to pop size? (as vec)
                      DDapply="matrix", 
                      Mfunction= "Min",
                      intensity= 95,  # percentage you want REMOVED from pop at time T=ry
                      remyear = 5,  # removal year = decrease from following year
                      return.vec= TRUE, 
                      return.remvec = TRUE) 

# almost works fine - ' Error in if (sum(thisAmat < 0) > 0) { : 
# missing value where TRUE/FALSE needed ' when values below 0?