################################################################################
## Example functions for density-dependent population projections
## 
## Chrissy Hernandez
## November 2025

repro_DD<- function(N_t, params, model='Exponential'){
  if (is.null(model)){
    warning("No density-dependent model specified. Assuming density-independent reproduction.")
    percapita_repro<- params$fmax
  }
  
  if (model=="Logistic" | model=="logistic"){
    percapita_repro<- params$fmax * (1- (1-(1-params$sa_max)/(params$fmax*params$sj_max))*(N_t/params$repro_K) )
  }
  
  if (model=="Beverton-Holt" | model=="beverton-holt"){
    percapita_repro<- params$fmax / (1+ (params$fmax*params$sj_max/(1-params$sa_max) -1)*N_t/params$repro_K)
  }
  
  if (model %in% c("Exponential", "exponential", "Ricker", "ricker")){
    percapita_repro<- params$fmax * exp(-params$b*N_t)
  }
  
  return(percapita_repro)
  
}

## Function for the D factor that we can apply to the full matrix or a particular part of the matrix:
calc_Dfactor<- function(N, params=NULL, DDfunc="Ricker"){
  if (is.null(params)){
    stop("No parameters provided for density-dependent scaling.")
  }
  
  # Exponential:
  if (DDfunc %in% c("Exponential", "exponential", "Ricker", "ricker")){
    if (is.null(params$b)){
      warning("No b parameter provided for Ricker function - assuming D_factor=exp(-N).")
      params$b<-1
    }
    D_factor<- exp(-params$b*N)
  }
  
  # Beverton-Holt
  if (DDfunc %in% c("Beverton-Holt", "beverton-holt")){
    if (is.null(params$b)){
      warning("No b parameter provided for Beverton-Holt function - assuming D_factor=1/(1+N).")
      params$b<-1
    }
    D_factor<- 1/(1+params$b*N)
  }
  
  # Logistic
  if (DDfunc %in% c("Logistic", "logistic")){
    if (is.null(params$K)){
      warning("No K parameter provided for logistic equation - assuming D_factor=1-N/100.")
      params$K<- 100
    }
    D_factor<- 1-N/params$K
  }
  
  # Theta-logistic
  if (DDfunc %in% c("Theta-logistic", "theta-logistic")){
    if (is.null(params$K)){
      warning("No K parameter provided for logistic equation - assuming K=100.")
      params$K<- 100
    }
    if (is.null(params$theta)){
      warning("No theta parameter provided for logistic equation - assuming theta=1.")
      params$theta<- 1
    }
    D_factor<- 1-(N/params$K)^params$theta
  }
  
  # Logistic function for the probability of transitioning into the reproductive class
  # (generalization of the peregrine falcon example)
  if (DDfunc == "log_pb"){
    D_factor<- exp(params$x0 + params$x1*N)/(1+exp(params$x0 + params$x1*N))
  }
  
  return(D_factor)
}

## Function to apply the D_factor to a separable MPM. This function takes in
## both the Fmat and Umat, and returns Amat(N).
apply_Dfactor<- function(Fmat, Umat, N, DDfunc, DDparams, DDmethod='matrix', 
                         matelem=NULL, env.params=NULL, env.method=NULL){
  
  Amat<- Fmat + Umat
  Dfactor<- calc_Dfactor(N, DDparams, DDfunc)
  # if (Dfactor<0){
  #   Dfactor<- 0
  # }
  
  if (is.null(env.params) & is.null(env.method)){
    # Apply the DD parameter to the whole matrix:
    if (DDmethod=="matrix"){ 
      eye<- diag(1, nrow=nrow(Amat), ncol=ncol(Amat))
      #Amat_N<- Dfactor*(Amat-eye)+eye
      Amat_N<- Dfactor*Amat
    } else if (DDmethod=='survival'){
      Umat_N<- Dfactor*(Umat)
      Amat_N<- Fmat + Umat_N
    } else if (DDmethod=='fertility'){
      Fmat_N<- Dfactor*Fmat
      Amat_N<- Fmat_N + Umat
    } else if (DDmethod=='falcon'){
      Amat_N<- Amat
      Amat_N[c(1,4),3]<- Dfactor*Amat[c(1,4),3]
      Amat_N[3,3]<- (1-Dfactor)*Amat[3,3]
    } else if (DDmethod=='reprotrans'){ # for now, this is identical to the falcons. Need to generalize!
      Amat_N<- Amat
      Dfactor_0<- calc_Dfactor(0, DDparams, DDfunc)
      Amat_N[c(1,4),3]<- (Dfactor/Dfactor_0)*Amat[c(1,4),3]
      Amat_N[3,3]<- (1-Dfactor)/(1-Dfactor_0)*Amat[3,3]
    } else if (DDmethod=='element'){
      if (is.null(matelem)){
        stop("You have specificed that density dependence should apply only to a 
           specific matrix element, but you have not specified which one. 
           Provide the [i,j] row and column indices to scale particular matrix 
           elements, or choose a different method.")
      }
      Amat_N<- Amat
      Amat_N[matelem]<- Dfactor*Amat_N[matelem]
    }
  } else {
    if (is.null(env.params) | is.null(env.method)){
      stop("You have specified only environmental parameters or an environmental 
           method, but both are needed.")
    }
    if (is.null(env.params$sd)){
      warning("No standard deviation specified for the environmental effect, assuming sd=0.1.")
      env.params$sd<- 0
    }
    if (is.null(env.params$mean)){
      warning("No mean value specified for the environmental effect, assuming mean=0.")
      env.params$mean<- 0
    }
    if (env.method=='matrix'){
      eye<- diag(1, nrow=nrow(Amat), ncol=ncol(Amat))
      this_env<- rnorm(1, mean=env.params$mean, sd=env.params$sd)
      # Apply the DD parameter and environmental stochasticity to the whole matrix:
      if (DDmethod=="matrix"){ 
        Amat_N<- Dfactor*(Amat-eye)+eye
      } else if (DDmethod=='survival'){ 
        # Apply DD only to survival, and then environmental stochasticity to the
        # whole matrix.
        Umat_N<- Dfactor*(Umat)
        Amat_N<- Fmat + Umat_N 
      } else if (DDmethod=='fertility'){
        Fmat_N<- Dfactor*Fmat
        Amat_N<- Fmat_N + Umat
      } else if (DDmethod=='element'){
        if (is.null(matelem)){
          stop("You have specificed that density dependence should apply only to a 
           specific matrix element, but you have not specified which one. 
           Provide the [i,j] row and column indices to scale particular matrix 
           elements, or choose a different method.")
        }
        Amat_N<- Amat
        Amat_N[matelem]<- Dfactor*Amat_N[matelem]
      }
      
      # Add the environmental stochasticity:
      Amat_N<- Amat_N + this_env*eye
      
    } else {
      stop("The only method currently defined for applying environmental stochasticity is 'matrix'.")
    }
  }
  
  return(Amat_N)
  
}

## Function to project a population using our very general approach to density
## dependence using the Dfactor
project_Dfactor<- function(Fmat, Umat, vector=NULL, time=100, return.vec=TRUE, 
                           DDfunc="Ricker", DDparams=NULL, DDmethod='matrix', 
                           matelem=NULL, ismemberN=NULL){
  # This function takes a single set of low-density vital rates in the form of a
  # reproduction matrix (Fmatrix) and transition matrix (Umatrix), and then
  # projects using the chosen DD function and DD method
  
  # Check that the Fmatrix and Umatrix are the same dimensions
  if (sum(dim(Fmat)==dim(Umat))!=2){
    stop("Reproductive (Fmat) and survival/growth transition (Umat) matrices 
         must have the same dimensions.")
  }
  
  # Make the A matrix:
  Amat<- Fmat+Umat
  order <- dim(Amat)[1]
  
  # For the DD part of this function:
  if (is.null(DDparams)){
    stop("You must provide the parameters that match your selected function for density-dependent reproduction.")
  }
  
  # Stage names:
  stagenames <- dimnames(Fmat)[[2]]
  if (is.null(stagenames)){
    stagenames <- paste("S", as.character(1:order), sep = "")
  }
  
  # Which classes of the population contribute to the count of N?
  if (DDfunc=="log_pb" & is.null(ismemberN)){
    ismemberN<- 4
  }
  else if (is.null(ismemberN)){
    ismemberN<- 1:length(stagenames)
  }
  
  # Initialize an output object:
  out<- list(pop=vector(), vec=matrix(), mat=matrix())
  
  # Initial vector:
  if (is.null(vector)){
    n0<- rep(1/order, order)
  } else{
    n0 <- vector
  }
  
  # Initialize the output for population vector:
  Vec <- matrix(0, ncol = order, nrow = time + 1)
  Pop <- rep(NA, (time + 1))
  dimnames(Vec)[[2]] <- stagenames
  dimnames(Vec)[[1]] <- 1:(time+1)
  Vec[1, ] <- n0
  Pop[1] <- sum(n0)
  for (i in 1:time) {
    thisN<- sum(Vec[i,ismemberN])
    thisAmat<- apply_Dfactor(Fmat, Umat, thisN, DDfunc, DDparams, DDmethod, 
                               matelem)
    
    # If the projection matrix has any negative values in it, stop iterating but
    # return the projection up until this point.
    if (sum(thisAmat<0)>0){
      warning(paste("Projection stopped at time step", i, "because the density-dependent projection matrix has negative values."))
      break
    }
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]
    Pop[i + 1] <- sum(Vec[(i + 1), ])
    
  }
    
  if (is.null(dim(Pop))) 
    dim(Pop) <- time + 1
  out$pop <- Pop
  if (return.vec) 
    out$vec <- Vec
  out$mat <- thisAmat
  
  return(out)
  
}

## Function to project a population using our very general approach to density
## dependence using the Dfactor
project_Dfactor_stochenv<- function(Fmat, Umat, vector=NULL, time=100, return.vec=TRUE,
                                    DDfunc="Ricker", DDparams=NULL, DDmethod='matrix',
                                    matelem=NULL, env.params=NULL, env.method='matrix'){
  # This function takes a single set of low-density vital rates in the form of a
  # reproduction matrix (Fmatrix) and transition matrix (Umatrix), and then
  # projects using the chosen DD function and DD method, and including simple
  # environmental stochasticity.
  
  # Check that the Fmatrix and Umatrix are the same dimensions
  if (sum(dim(Fmat)==dim(Umat))!=2){
    stop("Reproductive (Fmat) and survival/growth transition (Umat) matrices 
         must have the same dimensions.")
  }
  
  # Make the A matrix:
  Amat<- Fmat+Umat
  order <- dim(Amat)[1]
  
  # For the DD part of this function:
  if (is.null(DDparams)){
    stop("You must provide the parameters that match your selected function for density-dependent reproduction.")
  }
  
  # Stage names:
  stagenames <- dimnames(Fmat)[[2]]
  if (is.null(stagenames)){
    stagenames <- paste("S", as.character(1:order), sep = "")
  }
  
  # Initialize an output object:
  out<- list(pop=vector(), vec=matrix(), mat=matrix())
  
  # Initial vector:
  if (is.null(vector)){
    n0<- rep(1/order, order)
  } else{
    n0 <- vector
  }
  
  # Initialize the output for population vector:
  Vec <- matrix(0, ncol = order, nrow = time + 1)
  Pop <- numeric(time + 1)
  dimnames(Vec)[[2]] <- stagenames
  Vec[1, ] <- vector
  Pop[1] <- sum(vector)
  for (i in 1:time) {
    thisAmat<- apply_Dfactor(Fmat, Umat, Pop[i], DDfunc, DDparams, DDmethod, 
                             env.params = env.params, env.method=env.method)
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]
    Pop[i + 1] <- sum(Vec[(i + 1), ])
  }
  if (is.null(dim(Pop))) 
    dim(Pop) <- time + 1
  out$pop <- Pop
  if (return.vec) 
    out$vec <- Vec
  out$mat <- thisAmat
  
  return(out)
  
}


project_DD<- function(A, vector=NULL, time=100, return.vec=TRUE, 
                      DD="Exponential", DD_params=NULL){
  # This function can take a list of population matrices (where they differ in
  # adult survival), but now reproductive output is a density-dependent
  # function.
  
  if (is.list(A) & length(A) == 1) {
    A <- A[[1]] }
  if (is.matrix(A)) {
    M1 <- A
    dim(A) <- c(dim(A), 1)
    dimnames(A)[[1]] <- dimnames(M1)[[1]]
    dimnames(A)[[2]] <- dimnames(M1)[[2]]
    dimnames(A)[[3]] <- NULL
  }
  if (is.list(A) & length(A) > 1) {
    numA <- length(A)
    alldim <- sapply(A, dim)
    if (!diff(range(alldim)) == 0) {
      stop("all matrices in A must be square and have the same dimension as each other")
    }
    dimA <- mean(alldim)
    L <- A
    A <- numeric(dimA * dimA * numA)
    dim(A) <- c(dimA, dimA, numA)
    for (i in 1:numA) {
      A[, , i] <- L[[i]]
    }
    dimnames(A)[[1]] <- dimnames(L[[1]])[[1]]
    dimnames(A)[[2]] <- dimnames(L[[1]])[[2]]
    dimnames(A)[[3]] <- names(L)
  }
  
  order <- dim(A)[1]
  stagenames <- dimnames(A)[[2]]
  if (is.null(stagenames)) {
    stagenames <- paste("S", as.character(1:order), sep = "") }
  nmat <- dim(A)[3]
  
  # For the DD part of this function:
  if (is.null(DD_params)){
    stop("You must provide the parameters that match your selected function for density-dependent reproduction.")
  }
  if (DD=="Exponential" | DD=="exponential"){
    if (is.null(DD_params$fmax)){
      stop("You must provide DD_params$fmax for an exponential decay DD function.") }
    if (is.null(DD_params$b)){
      b<- 1
    }
  }
  ## Note: We can change this to be Beverton-Holt, etc. See function above.
  ## Exponential has the fewest number of parameters, so I'm starting with that
  ## one.
  
  # Initialize an output object:
  out<- list(pop=vector(), vec=matrix(), mat=matrix())
  
  # Initial vector:
  if (is.null(vector)){
    n0<- rep(1/order, order)
  } else{
    n0 <- vector
  }
  
  # Set the order that the matrices will be applied:
  if (dim(A)[3]==1){
    MC<- rep(1, time)
  } else if (dim(A)[3]>=time){
    MC<- 1:time
  } else { # if A has fewer matrices than the requested number of timesteps, 
    # then loop through them until you reach the number of requested timesteps
    MC<- rep(1:nmat, ceiling(time/nmat))[1:time]
  }
  
  # Initialize the output for population vector:
  Vec <- matrix(0, ncol = order, nrow = time + 1)
  Pop <- numeric(time + 1)
  dimnames(Vec)[[2]] <- stagenames
  Vec[1, ] <- vector
  Pop[1] <- sum(vector)
  for (i in 1:time) {
    thisAmat<- A[, , MC[i]]
    thisrepro<- repro_DD(Pop[i], DD_params)
    thisAmat[1,order]<- thisrepro
    
    Vec[(i + 1), ] <- thisAmat %*% Vec[i, ]
    Pop[i + 1] <- sum(Vec[(i + 1), ])
  }
  if (is.null(dim(Pop))) 
    dim(Pop) <- time + 1
  out$pop <- Pop
  if (return.vec) 
    out$vec <- Vec
  out$mat <- A
  
  return(out)
  
}