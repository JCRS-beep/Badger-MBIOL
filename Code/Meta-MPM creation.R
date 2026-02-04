# meta-population matrix creation
# 13/01/26

# loadong library
library(tidyverse)


#loading params----
params<- data.frame(fmax= 0.8436,   # F fecundity max (max cubs per adult female) 
                    Sc_f_max=0.65,   # max cub survival (equal for sexes)
                    Sc_m_max=0.65,
                    b=0.004,       # temp value- must be calculated from provided datasets
                    rep_K= 1.8,          #litter size (K)
                    h= 6 )  # harem size per male

# stating with 2 patches, dispersal between at different rates ----
# n(t+1) = A * D * n0 
# patch 1 n0 = (10, 5, 3, 2)  
# patch 2 n0 = (3, 2, 2, 2)  
n0 <- list(c(10, 5, 3, 2), c(3, 2, 2, 2))

# create blank base matrix
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
Umat <- matrix(0, nrow=4, ncol=4)
rownames(Umat) <- stages
colnames(Umat) <- stages
# input vital rates
Umat[2,1]<- 0.851   # yearling f survival
Umat[2,2]<- 0.851   # adult f survival
Umat[4,3]<- 0.809   # yearling m survival
Umat[4,4]<- 0.749   # adult m survival


# creating Fmat for site 1 based on abundance ----
Fmat1_out <- mating.func(params,     # density dependent parameters
                     stages,   # Stages in life cycle graph (single sex)
                     Nf = 5,        # Adult and yearling females  
                     Nm = 2,             # Adult and yearling males
                     Mfunction= "min", 
                     return.mat= TRUE)    

Fmat2_out <- mating.func(params,     
                     stages,   
                     Nf = 2,        
                     Nm = 2,           
                     Mfunction= "min",  
                     return.mat= TRUE)    # both Fmats identical!

Fmat1 <- Fmat1_out$Fmat # same lol
Fmat2 <- Fmat2_out$Fmat
# no mating in dispersing individuals: count -> breed -> move


Dmat <- matrix(0, ncol= length(stages), nrow= 2) # row = number patches
colnames(Dmat) <- stages

stay1 <- c(0.7, 0.9, 0.5, 0.6)  # prob individ in patch1 remains in patch 1 for each age class
Dmat1 <- Dmat
Dmat1[1,] <- stay1 # stay (and disperse) for each class in patch 1 
Dmat1[2,] <- 1- stay1

stay2 <- c(0.3, 0.6, 0.4, 0.5)
Dmat2 <- Dmat
Dmat2[1,] <- stay2 # stay (and disperse) for each class in patch 1 
Dmat2[2,] <- 1 - stay2 

# Shortening this matrix creation - 
# Steps - create Umat - same for all patches
# n0 for all patches (also needed for projection)
# Create Fmat (matingfunc w/ Nf and Nm)
# Dmat creation using movement probs for each class from each patch calculation (distance and sex ratio and group size)
# 


# function trials
meta.mat <- function(Umat,   # vector of stage names
                     Fmat1,
                     Dmat1,
                     Fmat2, 
                     Dmat2,
                     patches = 2){    # i
  nStages<- ncol(Umat)/2  # number of stages
  stagenames <- c(colnames(Umat))
  
  # setting base for Amat
    Amat <- matrix(0, ncol= ncol(Umat), nrow= nrow(Umat))  

    # filling off diags
    # Patch 1 to 2
     Amat_off1 <- Amat   # creating patch 1 -> 2 amat
      for (i in 1:ncol(Amat)) {   # looping to fill with multiplied values
        Amat_off1[,i] <- Umat[,i] * Dmat1[2,i]  # second row = leaving patch 
        } 
    # patch 2 to 1  
     Amat_off2 <- Amat  # creating patch 2 -> 1 amat
    for (i in 1:ncol(Amat)) {   # looping to fill with multiplied values
      Amat_off2[,i] <- Umat[,i] * Dmat2[2,i] 
      }
  
     
    # for diagonals
      Amat1 <- Umat  # creating original amat combining fert and survival
      Amat1[1,2] <- Fmat1[1,2]
      Amat1[3,2] <- Fmat1[3,2]
      
      for (i in 1:ncol(Amat)) {   # looping to fill with multiplied values
        Amat1[,i] <- Amat1[,i] * Dmat1[1,i]  # first row of Dmat = stay, each col for stage
      }
      
      Amat2 <- Umat  # creating original amat combining fert and survival
      Amat2[1,2] <- Fmat2[1,2]
      Amat2[3,2] <- Fmat2[3,2]
      for (i in 1:ncol(Amat)) {   # looping to fill with multiplied values
        Amat2[,i] <- Amat2[,i] * Dmat2[1,i]  # first row of Dmat = stay, each col for stage
      }
      
      # combining Amats - patch 1 stay, next to it patch 2 move
      meta_mat <- cbind(Amat1, Amat_off2)   # leaving fertility blank in second mat
      meta_mat_low <- cbind(Amat_off1, Amat2) 
      meta_mat <- rbind(meta_mat, meta_mat_low) 
      colnames(meta_mat) <- rep(stagenames, 2)
    return(meta_mat)
  }
  
mattest <- meta.mat(Umat,   
                    Fmat1,
                    Dmat1,
                    Fmat2, 
                    Dmat2,
                    patches = 2)
# SUCCESS!

# next steps = adapting this function for multiple patch inputs = only once i create Dmat automatically from distance and Nf and Nm (logit link)

# diserpsal prob and dmat creation
# using distance to link p(move) = 0 or 1
# predictor variables = sex, stage, distance, density and sex ratio
# step 1 : load rogers movement data for movement and group size av distance for moves?
 
