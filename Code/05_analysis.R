# Analysis script - metrics and rate calculation
# 20.04

# custom functions for this script -----
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


# baseline checks  -------
# average pop size and lambda
popN <- N.extract(rep_proj0)    # pop sizes per year per repetition
Nfin <- sapply(popN, function(x) x[20])   # extracting final entry per rep
summary(Nfin)
av_fin_N <- round(mean(Nfin), 1)    # mean final N value used in results and comparison

av_lamb <- round(summary(lamb.av(rep_proj0)), 2)  # lambda value - mean and range

av_ssd0 <- ssd.av(rep_proj0)
base_prop <- round(colMeans(av_ssd0$av_prop), 2)  # proportion adult used in results



# single removal scenarios comparison -------
# Pop size
rel_proj1 <- relative.pop(rep_proj1,   
                          baseline_list = rep_proj0)    # extracting pop size metrics for single rems
rel_proj2 <- relative.pop(rep_proj2,   
                          baseline_list = rep_proj0) 
rel_proj3 <- relative.pop(rep_proj3,   
                          baseline_list = rep_proj0) 

rel_projs <- list(rel_proj1, rel_proj2, rel_proj3)  
rel_df <- rel.df(rel_projs)  # using prev defined function to turn into comparison df

# sex ratio
rep_list <- list(rep_proj0, rep_proj1, rep_proj2, rep_proj3)
sex_list <- lapply(rep_list, ssd.av) 
sex_df <- sex.df(sex_list)  # custom function to turn lists into a dataframe



# Double rems --------
# Pop size 
rel_du1 <- relative.pop(du_proj1,   
                        baseline_list = rep_proj0) 
rel_du2 <- relative.pop(du_proj2,   
                        baseline_list = rep_proj0) 
rel_du3 <- relative.pop(du_proj3,   
                        baseline_list = rep_proj0) 
du_list <- list(rel_du1, rel_du2, rel_du3) 
du_rel_df <- rel.df(du_list)

# sex ratio 
av_ssd0 <- sex_list[[1]]   # extracting baseline, stored in single list already
av_du_ssd1 <- ssd.av(du_proj1, return.Mats = FALSE)
av_du_ssd2 <- ssd.av(du_proj2, return.Mats = FALSE)
av_du_ssd3 <- ssd.av(du_proj3, return.Mats = FALSE)

# combine in list with baseline 
du_sex_list <- list(av_ssd0, av_du_ssd1, av_du_ssd2, av_du_ssd3)
du_sex_df <- sex.df(du_sex_list)  # turning list into df



# Multiple removals
# Pop size
rel_multi1 <- relative.pop(multi_proj1,   
                           baseline_list = rep_proj0) 
rel_multi2 <- relative.pop(multi_proj2,   
                           baseline_list = rep_proj0) 
rel_multi3 <- relative.pop(multi_proj3,   
                           baseline_list = rep_proj0) 
multi_list <- list(rel_multi1, rel_multi2, rel_multi3) 
multi_rel_df <- rel.df(multi_list)

# sex ratio 
av_multi_ssd1 <- ssd.av(multi_proj1, return.Mats = FALSE)
av_multi_ssd2 <- ssd.av(multi_proj2, return.Mats = FALSE)
av_multi_ssd3 <- ssd.av(multi_proj3, return.Mats = FALSE)

# combine in list with baseline
multi_sex_list <- list(av_ssd0, av_multi_ssd1, av_multi_ssd2, av_multi_ssd3)
multi_sex_df <- sex.df(multi_sex_list)  



# Continuos removals ------
rel_cont1 <- relative.pop(cont_proj1,   
                          baseline_list = rep_proj0) 
rel_cont2 <- relative.pop(cont_proj2,   
                          baseline_list = rep_proj0) 
rel_cont3 <- relative.pop(cont_proj3,   
                          baseline_list = rep_proj0) 

cont_list <- list(rel_cont1, rel_cont2, rel_cont3) 
cont_rel_df <- rel.df(cont_list)  # comparison df

# sex ratio 
av_cont_ssd1 <- ssd.av(cont_proj1, return.Mats = FALSE)
av_cont_ssd2 <- ssd.av(cont_proj2, return.Mats = FALSE)
av_cont_ssd3 <- ssd.av(cont_proj3, return.Mats = FALSE)

# combine in list with baseline
cont_sex_list <- list(av_ssd0, av_cont_ssd1, av_cont_ssd2, av_cont_ssd3)
cont_sex_df <- sex.df(cont_sex_list)  



# combined analyses -----
# creating new dfs with pop size info to bind together
sing <- rel_df    
sing$Frequency <- rep(as.character("Single", nrow(du_rel_df)))

duplo <- du_rel_df 
duplo$Frequency <- rep(as.character("Double", nrow(du_rel_df)))

multi <- multi_rel_df 
multi$Frequency <- rep(as.character("Multiple", nrow(multi)))

cont <- cont_rel_df
cont$Frequency <- rep(as.character("Continuous", nrow(cont)))

comb_rel_df <- rbind(sing, duplo, multi, cont)

# combined final pop N
# fit the two-way ANOVA model
finN_model <- aov(relative_mean_N ~ Strategy * Frequency, data = comb_rel_df)

#view the model output
summary(finN_model)

anova_N_tab <- anova(finN_model)
N_var <- round((anova_N_tab$`Sum Sq` / sum(anova_N_tab$`Sum Sq`)), 2)    # values to use in paper - how much variance explained by each factor?


# combined sex ratio 
base_ratio <- mean(av_ssd0$sex_ratio)     # baseline sex ratio as proportion female - calculated in beginning

sing_sr <- sex_df[101:nrow(sex_df),] # exclude baseline sex ratio for sinlge rem df
sing_sr$Frequency <- rep(as.character("Single", nrow(sing_sr)))

du_sr <- du_sex_df[101:nrow(du_sex_df),]   # specify rows 100 onward
du_sr$Frequency <- rep(as.character("Double", nrow(du_sr)))

multi_sr <- multi_sex_df[101:nrow(multi_sex_df),]
multi_sr$Frequency <- rep(as.character("Multiple", nrow(multi_sr)))

cont_sr <- cont_sex_df[101:nrow(cont_sex_df),]
cont_sr$Frequency <- rep(as.character("Continuous", nrow(cont_sr)))

comb_sex_df <- rbind(sing_sr, du_sr, multi_sr, cont_sr)  

# fit the two-way ANOVA model
sr_model <- aov(sex_ratio ~ Strategy * Frequency, data = comb_sex_df)

#view the model output
summary(sr_model)

anova_sr_tab <- anova(sr_model)
sr_var <- round((anova_sr_tab$`Sum Sq` / sum(anova_sr_tab$`Sum Sq`)), 2) # use in paper - how much variance explained by each factor?


# Combined Extinction risk and vulnerabilities (low pop sizes) ----
single_rep_list <- list(rep_proj1, rep_proj2, rep_proj3)  
single_vul <- sapply(single_rep_list, extinction.risk)  # no extinctions in baseline or single rem scenarios

du_rep_list <- list(du_proj1, du_proj2, du_proj3)   # listing repeated projection
du_vul <- sapply(du_rep_list, extinction.risk)  # some extinctions in consecutive removals

multi_rep_list <- list(multi_proj1, multi_proj2, multi_proj3)   # SAVE FOR COMB
multi_vul <- sapply(multi_rep_list, extinction.risk)  # freq extinctions in multi removals

cont_rep_list <- list(cont_proj1, cont_proj2, cont_proj3)
cont_vul <- sapply(cont_rep_list, extinction.risk) 

vul_df <- data.frame(cbind(single_vul, du_vul, multi_vul, cont_vul))   # proj, remfreq and resulting risk val

# setting 
vul_df <- gather(vul_df, key = "Frequency", value = "Vulnerability")   # rename risk later?
vul_df$Strategy <- as.character(rep(c("Random", "Adult male", "Adult female"), 4))
vul_df$Frequency <- factor(vul_df$Frequency,
                           levels = c("single_vul", "du_vul", "multi_vul", "cont_vul"),
                           labels = c("Single", "Double", "Multiple", "Continuous"))

# performing 2 way ANOVA
vul_model <- aov(Vulnerability ~ Strategy * Frequency, data = vul_df)
summary(vul_model)   #view the model output

anova_vul_tab <- anova(vul_model)
vul_var <- round((anova_vul_tab$`Sum Sq` / sum(anova_vul_tab$`Sum Sq`)), 2) 
# warning In anova.lm(vul_model) :
 # ANOVA F-tests on an essentially perfect fit are unreliable






