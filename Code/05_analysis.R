# Analysis script - metrics and rate calculation
# 20.04

# custom functions for this script -----
# # Turning relative pop size outputs into comparison dataframe for plotting
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
    df <- data.frame(Strategy, relative_mean_N, relative_final_N)  # what to do with this df? Store in list?
    
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



# baseline checks ----
popN <- N.extract(rep_proj0)
Nfin <- sapply(popN, function(x) x[20])
summary(Nfin)

av_lamb <- summary(lamb.av(rep_proj0))

av_ssd0 <- ssd.av(rep_proj0)
base_prop <- colMeans(av_ssd0$av_prop)


# single removal scenarios comparison -------
rel_proj1 <- relative.pop(rep_proj1,   
                          baseline_list = rep_proj0) 
rel_proj2 <- relative.pop(rep_proj2,   
                          baseline_list = rep_proj0) 
rel_proj3 <- relative.pop(rep_proj3,   
                          baseline_list = rep_proj0) 

rel_projs <- list(rel_proj1, rel_proj2, rel_proj3)  # will this always list in order? 

rel_df <- rel.df(rel_projs)  # using prev defined function to turn into comparison df


# comparing proportion female across projections
rep_list <- list(rep_proj0, rep_proj1, rep_proj2, rep_proj3)
sex_list <- lapply(rep_list, ssd.av) 
sex_df <- sex.df(sex_list)  # custom function to turn lists into a dataframe


# Double removals -----
rel_du1 <- relative.pop(du_proj1,   
baseline_list = rep_proj0) 
rel_du2 <- relative.pop(du_proj2,   
                        baseline_list = rep_proj0) 
rel_du3 <- relative.pop(du_proj2,   
                        baseline_list = rep_proj0) 
du_list <- list(rel_du1, rel_du2, rel_du3) 

du_rel_df <- rel.df(du_list)


# sex ratio 
du_rep_list <- list(du_proj1, du_proj2, du_proj3)

# need to add baseline for comparison in plot
av_ssd0 <- sex_list[[1]]
av_du_ssd1 <- ssd.av(du_proj1, return.Mats = FALSE)
av_du_ssd2 <- ssd.av(du_proj2, return.Mats = FALSE)
av_du_ssd3 <- ssd.av(du_proj3, return.Mats = FALSE)

du_sex_list <- list(av_ssd0, av_du_ssd1, av_du_ssd2, av_du_ssd3)
du_sex_df <- sex.df(du_sex_list)  # success!


# Multi removals over 5 years
multi_rep_list <- list(multi_proj1, multi_proj2, multi_proj3)
rel_multi1 <- relative.pop(multi_proj1,   
                           baseline_list = rep_proj0) 
rel_multi2 <- relative.pop(multi_proj2,   
                           baseline_list = rep_proj0) 
rel_multi3 <- relative.pop(multi_proj2,   
                           baseline_list = rep_proj0) 
multi_list <- list(rel_multi1, rel_multi2, rel_multi3) 

multi_rel_df <- rel.df(multi_list)

# sex ratio 
av_multi_ssd1 <- ssd.av(multi_proj1, return.Mats = FALSE)
av_multi_ssd2 <- ssd.av(multi_proj2, return.Mats = FALSE)
av_multi_ssd3 <- ssd.av(multi_proj3, return.Mats = FALSE)

# combine in list
multi_sex_list <- list(av_ssd0, av_multi_ssd1, av_multi_ssd2, av_multi_ssd3)
multi_sex_df <- sex.df(multi_sex_list)  


# combined df analyses -----
sing <- rel_df
sing$rem_freq <- rep(as.character("Single", nrow(du_rel_df)))
duplo <- du_rel_df 
duplo$rem_freq <- rep(as.character("Double", nrow(du_rel_df)))
multi <- multi_rel_df 
multi$rem_freq <- rep(as.character("Multiple", nrow(multi)))


comb_rel_df <- rbind(sing, duplo, multi)


# combined sex ratio plot
# baseline sex ratio
base_ratio <- mean(sex_df[1:100,]$sex_ratio)
sing_sr <- sex_df[101:nrow(sex_df),] # exclude baseline sex ratio? 
sing_sr$rem_freq <- rep(as.character("Single", nrow(sing_sr)))

du_sr <- du_sex_df[101:nrow(du_sex_df),]   # specify rows 100 onward?
du_sr$rem_freq <- rep(as.character("Double", nrow(du_sr)))

multi_sr <- multi_sex_df[101:nrow(multi_sex_df),]
multi_sr$rem_freq <- rep(as.character("Multiple", nrow(multi_sr)))

comb_sex_df <- rbind(sing_sr, du_sr, multi_sr)  


# Extinction risk and vulnerabilities (low pop sizes) ----
single_list <- list(rep_proj1, rep_proj2, rep_proj3)
single_vul <- sapply(single_list, extinction.risk)  # no extinctions in baseline or single rem scenarios
du_vul <- sapply(du_rep_list, extinction.risk)  # some extinctions in consecutive removals
multi_vul <- sapply(multi_rep_list, extinction.risk)  # freq extinctions in multi removals

# how to visualise - in bar plot?
vul_df <- data.frame(cbind(single_vul, du_vul, multi_vul))   # proj, remfreq and resulting risk val
# setting 
vul_df <- gather(vul_df, key = "Frequency", value = "Vulnerability")   # rename risk later?
vul_df$Strategy <- as.character(rep(c(1,2,3), 3))
vul_df$Frequency <- factor(vul_df$Frequency,
                           levels = c("single_vul", "du_vul", "multi_vul"),
                           labels = c("Single", "Double", "Multi"))


