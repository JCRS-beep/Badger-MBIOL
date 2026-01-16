# meta-population matrix creation
# 13/01/26

#loading params----
params<- data.frame(      # dataframe 
  fmax= 3.2,     # F fecundity max (max cubs per adult female) 
  Sc_f_max=0.65,   # max cub survival (equal for sexes)
  Sc_m_max=0.65,
  b=0.002,       # temp value- must be calculated from provided datasets
  rep_K= 4,          #litter size (K)
  h= 6   # harem size per male
)

# stating with 2 patches, dispersal between at different rates ----
# n(t+1) = A * D * n0 
# patch 1 n0 = (10, 5, 3, 2)  p(stay) = 0.8
# patch 2 n0 = (3, 2, 2, 2)   p(stay) = 0.5
n0 <- list(c(10, 5, 3, 2), c(3, 2, 2, 2))

# create blank base matrix
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
Umat <- matrix(0, nrow=4, ncol=4)
rownames(Umat) <- stages
colnames(Umat) <- stages
# input vital rates
Umat[2,1]<- 0.67   # yearling f survival
Umat[2,2]<- 0.78   # adult f survival
Umat[4,3]<- 0.65   # yearling m survival
Umat[4,4]<- 0.72  # adult m survival


# creating dispersal prob matrix (Dmat) - 2 sites 2X2, p(stay) vs p(leave)----
Dmat <- matrix(0, ncol=2, nrow=2)
rownames(Dmat) <- c("p1", "p2")
colnames(Dmat) <- c("p1", "p2")
Dmat[,1] <- c(0.8, 0.2)
Dmat[,2] <- c(0.5, 0.5)


# clone Umat x times, where x= sites^2. Multiply each by relevant Dmat entry ----
Umat1 <- Umat * Dmat[1,1]
Umat2 <- Umat * Dmat[2,2]
Umat21 <- Umat * Dmat[2,1]    # from 1 to 2 (2,1)
Umat12 <- Umat * Dmat[1,2]    # from 2 to 1

# merge matrices together, bind y1 y12
Umat <- cbind(Umat1, Umat12)
Umat_low <- cbind(Umat21, Umat2)
Umat <- rbind(Umat, Umat_low)



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

Fmat1 <- Fmat1_out$Fmat
Fmat2 <- Fmat2_out$Fmat
# how to work out mating function for Fmat when considering dispersal? Use prev year N


# combining into Amat
Amat <- Umat 
Amat[1,2] <- Fmat1[1,2]
Amat[3,2] <- Fmat1[3,2]
Amat[5,6] <- Fmat2[1,2]
Amat[7,6] <- Fmat2[3,2]

# what to enter for dispersal?

# Shortening this matrix creation
meta.mat <- function(stagenames,   # vector of stage names
                     g1,   # vector of all growth rates/ transition probabilities
                     g2= NULL, # second sex values for growth
                     s, 
                     patches = 2,
                     full_mat = TRUE){    # is Umat replicated across whole matrix or just diagonal?
  nStages<- length(stagenames)/2                            # number of stages
  # create out obj
  out <- list(U <- matrix(), mat <- matrix())
  
  # Umat creation
   Umat <- matrix(0, nrow=(2*nStages), ncol=(2*nStages))     # empty matrix with equal rows and columns as stages
   colnames(Umat) <- stagenames                                       # naming columns after stages
   rownames(Umat) <-  stagenames                             # naming rows after stages
         
   Umat[nStages, nStages - 1 ]<- g1   # inputs growth rates into respective matrix location
   Umat[2*nStages,  nStages + 1] <- g2
   Umat[nStages,nStages] <- s[1]                      # final stage remaining is adult survival
   Umat[2*nStages, 2*nStages] <- s[2]
  
   if(isTRUE(full_mat)){
  # Replicate Umat (pacthes^2)
   trans <- matrix(rep(t(Umat), times = patches^2), ncol= ncol(Umat), byrow= TRUE)       # transposed mat nrow = patches * length(stagenames)
   meta <- trans[1:(patches*nrow(Umat)),]   # forming a stack
   Metamat <- cbind(meta, meta, meta) }     # NEEDS GENERALISING
 # Metamat <- matrix(0, ncol= (patches * 2 *nStages), nrow=(patches * 2 *nStages))        # meta mat nrow = patches * length(stagenames)
 # Metamat[1:2*nStages,1:2*nStages] <- Umat
  
   else if(full_mat == FALSE) {
     Metamat <- matrix(0, ncol= (patches * 2 *nStages), nrow=(patches * 2 *nStages))        # meta mat nrow = patches * length(stagenames)
      Metamat[1:2*nStages,1:2*nStages] <- Umat  # Umat in first quadrant
   }
  # filling in outputs
  out$U <- Umat
  out$mat <- Metamat
      return(out)  
  } 

mat_test<- meta.mat(stagenames = stages, 
                     g1= 0.67, 
                     g2= 0.65,
                     s= c(0.78, 0.72), 
                    patches = 3)
outU <- mat_test$U  # Umat correct
mat_test$mat # done


trans <- matrix(rep(t(outU), times = 3^2), ncol=ncol(outU), byrow= TRUE) # 9 vertical stacks to matrix?

meta1 <- trans[1:(3*nrow(outU)),]   # getting Umat stack!
meta <- cbind(c(rep(meta1,times = 3)))  # DONE!


vec <- noquote(rep("meta1,", times = 3))
meta <- cbind(vec)

patches <- 3
trans <- matrix(rep(t(outU), times = patches^2), ncol= ncol(outU), byrow= TRUE)       # transposed mat nrow = patches * length(stagenames)
meta <- trans[1:(patches*nrow(Umat)),]   # forming a stack
chain <- noquote(rep("meta,", times = patches -1))
Metamat <- cbind(chain, meta)     
