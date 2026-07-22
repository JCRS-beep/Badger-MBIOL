# custom functions for calculating output metrics:
# Written by Jay Creese
# April 2026
# Last updated: 14/07


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
# Function name : growth.rate
# Arguments : 
# dd_proj or rem_proj output
# vis = TRUE or FALSE for whether you want graph plotted
# Use : generates growth rate per year (gr) by dividing pop size by previous year 




# calculating average lambda, av ssd in single function
# Splitting calculations into their own functions? Av lambda, relative pop sizes?

# calculating pop size per rep function
N.extract <- function(proj_list) {
  reps <- length(proj_list)
  N_list <- list()   # list of N vecs for each repeat
  
  for (t in 1:reps){
    N_list[[t]] <- proj_list[[t]]$pop    # isolating pop size vector
  }
  
  return(N_list)
}


lamb.av <- function(proj_list, return.Lambda = FALSE){
  reps <- length(proj_list)
  time <- length(proj_list[[1]]$pop)- 1 # minus 1 for the initial population at t = 0
  
  N_list <- N.extract(proj_list)
  
  # calculating lambda for each rep lambda calculation as a function applied to pop size vec
  lambda <- function(N){
    N[2:length(N)]/N[1:(length(N)-1)]
  } 
  
  lambda_list <- lapply(N_list , lambda)
  av_lambda <- sapply(lambda_list, mean)
  
  # setting up out obj
  lamb_out <- list(lambda_list = vector(), 
                   av_lambda = vector() )
  
  if(return.Lambda == TRUE){
    lamb_out$lambda <- lambda_list
  }
  lamb_out$av_lambda <- av_lambda
}
# lambda : lambda per year per rep (list of 100 vecs length 20)
# av_lambda : average lambda per rep (vec length 100)


# ssd function  - 
ssd.av <- function(proj_list, return.Mats = FALSE){
  
  reps <- length(proj_list)
  time <- length(proj_list[[1]]$pop)- 1
  
  # separating a matrix from out obj
  abun_mat <- vector("list", reps)   # list of abundances for each repeat
  
  for (t in 1:reps){
    abun_mat[[t]] <- proj_list[[t]]$vec    # isolating pop size vector
  }
  
  mat <- matrix(0, ncol = ncol(abun_mat[[1]]), nrow = nrow(abun_mat[[1]]))    # list of empty matrix to fill with stage props
  
  
  stageMat <- vector("list", reps)   # fill each stage mats list with mat?
  
  # out obj is a list with matrix and av prop 1 row matrix
  
  for (t in 1:reps){
    for(i in 1:nrow(abun_mat[[1]])) {   # loop for each column 
      mat[i,] <- abun_mat[[t]][i,]/sum(abun_mat[[t]][i,])          # column i of matrix filled with row i divided by col sum
    } 
    stageMat[[t]] <- mat # ssd for each year per rep
  }
  
  list_prop <- lapply(stageMat, colMeans) # list per rep, each single vec of av abundances
  
  # combining lists into a single matrix
  av_prop <- matrix(unlist(list_prop), byrow = TRUE, nrow = reps, ncol = 4) # filling each row with list vec
  colnames(av_prop) <- colnames(abun_mat[[1]])
  
  # calculating sex ratio (proportion female in pop)
  fems <- function(x) {
    x[2]/(x[2] + x[4])   #  entries 2/ sum adult entries = proportion adult fems in adult pop
  }
  sr <- sapply(list_prop, fems)  # vector of adult fem proportion 
  
  
  # set up out obj
  ssd_out <- list(ssdMat = matrix(), av_prop = matrix(), sex_ratio = vector())
  
  if(return.Mats == TRUE){
    ssd_out$ssdMat <- stageMat # stage proportions by year
  }
  
  ssd_out$av_prop <- av_prop  # single stage proportion for each rep 
  
  ssd_out$sex_ratio <- sr  
  
  return(ssd_out)
}
# Outputs 
# ssdMat : stage dist mat
# av_prop : average proportions of each stage in matrix   




relative.pop <- function(proj_list,   
                         baseline_list = NULL) # baseline for comparison. if empty, no relative population sizes calculated
{
  # basic checks
  if(length(baseline_list)!= length(proj_list)){
    stop(paste("length of baseline and projection differ - comparison not possible. Ensure repetitions are equal before continuing."))
  }
  
  reps <- length(proj_list)
  time <- length(proj_list[[1]]$pop)- 1  # minus 1 for the initial population at t = 0
  
  N_list <- N.extract(proj_list)   # list of N vecs for each repeat
  base_list <- N.extract(baseline_list)
  
  fin_N <- sapply(N_list, function(x) x[20])
  
  # average pop size proportions
  av_N <- sapply(N_list, mean)   # vector of mean pop sizes per rep
  av_base_N <- mean(sapply(base_list, mean))   # SINGLE value for mean pop size across all reps of baseline
  
  if(!is.null(baseline_list)){   # only if there is comparison proj provided
    # calculating final pop size in each scenario  (av?)
    fin_props <- vector(length = reps)   # vector of proportions
    
    for (t in 1:reps){
      fin_props[t] <- N_list[[t]][time]/base_list[[t]][time]     #  list calculate proportions
    }
    
    av_base_N <- sapply(base_list, mean)  
    relative_meanN <- av_N / av_base_N   # issues = if baseline proj goes to 0 for some reason, causes
  }
  
  # output proj
  pop_out <- list(fin_N = vector(), av_fin_N = numeric(), 
                  fin_props = vector(), av_fin = vector(), 
                  relative_meanN= vector()) # number of matrices = rep, rows = time
  
  
  pop_out$fin_N <- fin_N
  pop_out$av_fin_N <- mean(fin_N)
  pop_out$fin_props <- fin_props
  pop_out$av_fin <- mean(fin_props, na.rm = TRUE)
  pop_out$relative_meanN <- relative_meanN
  
  
  return(pop_out)
}
# NOTES
# Name =  av_pop
# proj_list,   
# baseline_list : the baseline projection to compare projections
# return.Lambda : return lambda per year for each rep? T/F
# return.Mats : return stage distribution per year for each rep? T/F
# 
# Outputs
# fin_props : final year proportions of projection compared to baseline (vec length 100)
# av_fin : average final pop size of proj compared to baseline (numeric val).
# relative_meanN : mean pop size across proj compared to mean pop size of baseline (vec length 100)


### ISSUES WITH COMPARISONS - if there is a year when baseline goes to extinction, this skews that years comparison. 
# better to have an average pop size for baseline and compare all repeats to this value?


# Vulnerability function - how many of our repetitions fall below pop size 10?
extinction.risk <- function(proj_list){
  reps <- length(proj_list)
  time <- length(proj_list[[1]]$pop)- 1 # minus 1 for the initial population at t = 0
  
  # isolating pop size from list
  pop_sizes <- N.extract(proj_list)
  
  
  # searching for any projection where: pop size falls below 10
  vul.idx <- function(pop) {    # population size vector 
    if(any(pop <= 20)){
      which(20 >= pop)[[1]]    # returns first index in vector less than 10
      
    } else {return(FALSE)}    # returns FALSE
  }
  
  vul_it <- sapply(pop_sizes, vul.idx)  # apply across whole list to count 
  # vector with 0 for normal projection, number when low pop size occured
  
  # number of numeric values is number of vulnerable projs
  nVul <- length(vul_it[vul_it >=1])
  prop_vul <- nVul/ reps
  
  return(prop_vul)
}



# Stage distribution calculation
# Updating function to match output of ddproj obj
ssd <- function(out, vis = FALSE, cols) {     # input projected matrix
  # separating a matrix from out obj
  mat <- out$vec   
  
  stageMat<- matrix(0, ncol=ncol(mat), nrow=nrow(mat))    # empty matrix to fill
  rownames(stageMat) <- rownames(mat)
  colnames(stageMat) <- colnames(mat)
  
  # out obj is a list with matrix and plot
  ssd_out <- list(stageMat = matrix(), plot = NA)
  
  for(i in 1:nrow(mat)) {   # loop for each column 
    stageMat[i,]<- mat[i,]/sum(mat[i,])          # column i of matrix filled with row i divided by col sum
  } 
  ssd_out$stageMat <- stageMat
  
  if(!isFALSE(vis)){
    # stageMat as df
    nStage <- ncol(stageMat)   # number of classes and sexes (if nStages = 2 and sex =2, x =4)
    # turning into dataframe
    df <- as.data.frame(stageMat)
    df$Year <- as.numeric(rownames(stageMat))    # year column from 0 to t years
    
    # tidy data - converting to long format so each row is a single observation 
    df_long <- gather(df, key= "Stage", value = "Proportion", 1:nStage)   # creating a stage col in df with abundance
    df_long <- separate(df_long, col= "Stage", into= c("Stage", "Sex"), sep='_')   # splliting by sex, seperated by _
    
    plot <- ggplot(data = df_long, 
                   aes(x = Year, y = Proportion, colour = Sex,  # qhy has year ordered so weird?
                       linetype = Stage, shape = Stage)) +  # sexes diff cols, shapes and lines diff for stages
      geom_line(data = df_long, position = "jitter") + # why is this not joining as other 
      scale_colour_manual(values = cols,
                          labels = c("Female", "Male")) +
      labs(title = "Stage Proportions over Time", 
           x = "Year", y = "Proportion of population") 
    
    ssd_out$plot <- plot   # issue here - not loading plot as object within list
  }
  
  return(ssd_out)   # returns matrix of each stage as proportion of total pop
  
}
# Function name= prop.stage
# Arguments: 
# out = stage abundance matrix over a time interval, produced by proj function
# vis = whether to print graph
# cols = vector of 2 colours for each sex
# Purpose:  Calculates the proportion of each stage class out of the total pop size in a given year





# Used in analysis script in dataframe construction

# Turning relative pop size outputs into comparison dataframe for plotting
rel.df <- function(rel_projs = "list")  # list of all projections to compare, length = n projs
{
  nProj <- length(rel_projs)  # how many projections do we have?
  relative_df <- data.frame()  # to store other dfs in
  
  # set up individual dfs for eac proj
  for (p in 1:nProj){   # repeat for each projection
    # assigning strategy 1,2,3 own names. Worry = if names do not match actual projection
    if(p == 1) strat = "Random" 
    if(p == 2) strat = "Adult males" 
    if(p == 3) strat = "Adult females"
    
    Strategy <- as.character(rep(strat, 100))   # better with an informative name (projection1)
    relative_mean_N <- rel_projs[[p]]$relative_meanN
    relative_final_N <- rel_projs[[p]]$fin_props
    final_N <- rel_projs[[p]]$fin_N
    av_final_N <-  rel_projs[[p]]$av_fin_N
    df <- data.frame(Strategy,final_N, av_final_N, relative_mean_N, relative_final_N)  # what to do with this df? Store in list?
    
    # storing df in our relative df - intial 
    relative_df <- rbind(relative_df, df)
  }
  return(relative_df)
}

# Turning list outputs into dataframe with av stage proportion per rep and sex ratio for comparison 
sex.df <- function(sex_projs = "list")  # list of all projections to compare, length = n projs
{
  nProj <- length(sex_projs)  # how many projections do we have?
  sex_df <- data.frame()  # to store other dfs in
  
  # set up individual dfs for eac proj
  for (p in 1:nProj){   # repeat for each projection
    n <- p-1   # incl proj0, need to name 0
    if(n == 0) strat = "Baseline" 
    if(n == 1) strat = "Random" 
    if(n == 2) strat = "Adult males" 
    if(n == 3) strat = "Adult females"
    
    Strategy <- as.character(rep(strat, 100))   # better with an informative name (projection1)
    av_prop <- sex_projs[[p]]$av
    sex_ratio <- sex_projs[[p]]$sex_ratio
    
    df <- data.frame(Strategy, av_prop, sex_ratio)  # what to do with this df? Store in list?
    
    # storing df in our relative df - intial 
    sex_df <- rbind(sex_df, df)
  }
  return(sex_df)
}