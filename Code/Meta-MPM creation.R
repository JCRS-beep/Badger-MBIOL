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
Fmat1 <- mating.func(params,     # density dependent parameters
                     stages,   # Stages in life cycle graph (single sex)
                     Nf = 15,        # Adult and yearling females
                     Nm = 5,             # Adult and yearling males
                     Mfunction= "min", 
                     return.mat= TRUE)    

Fmat2 <- mating.func(params,     
                     stages,   
                     Nf = 5,        
                     Nm = 4,           
                     Mfunction= "min",  
                     return.mat= TRUE)    # both Fmats identical!

# how to work out mating function for Fmat when considering dispersal?


# combining into Amat
Amat <- Umat 

# what to enter for dipsersal?










