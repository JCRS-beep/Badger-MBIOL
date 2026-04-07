# Script for results draft 
# 15.03.26

library(ggplot2)
library(gridExtra)
library(tidyr)
library(tidyverse)
library(dplyr)
library(readr)
library(here)

# sourcing from model projection scripts
source(here("Code/03_model_projections.R"))  # loads all params, model projections

# future improvements - final pop size plot, sex ratio plot, extinction plot with all strats and freqs 

# Functions used in this script -----
# Turning relative pop size outputs into comparison dataframe for plotting
rel.df <- function(rel_projs = "list")  # list of all projections to compare, length = n projs
{
  nProj <- length(rel_projs)  # how many projections do we have?
  relative_df <- data.frame()  # to store other dfs in
  
  # set up individual dfs for eac proj
  for (p in 1:nProj){   # repeat for each projection
    # assigning strategy 1,2,3 own names. Worry = names do not match actual projection
    if(p == 1) strat = "random" 
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
    if(n == 0) strat = "baseline" 
    if(n == 1) strat = "random" 
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

# final N plot creation
# rel plot function - WHY not wokring?
rel.plot <- function(rel_df, yval, save_name = FALSE){
  
  baseplot <- ggplot(data = rel_df, aes(x = Strategy, y = yval)) +
    geom_boxplot()
  
  plot <- baseplot + 
    geom_hline(yintercept = 1, aes(colour = "grey20") ) +
    labs(y = "Final population size relative to baseline average") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme_minimal() +
    theme(
      text = element_text(size = 16),
      plot.title = element_text(size = 16 + 2, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16 - 2),
    ) 
  
  if(save_name != FALSE){ # saving as ggsave
    ggsave(filename = save_name,
           plot = plot,
           device = "png",
           path = here("Figs"), 
           bg = "white")
  }
  
  return(plot)
}
#  Error in `geom_boxplot()`:
# ! Problem while computing aesthetics.
# ℹ Error occurred in the 1st layer.
# Caused by error:
#  ! object 'relative_final_N' not found


# average pop size and lambda of baseline -------
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

rel_df <- rel.df(rel_projs)  # using prev defined function to turn inot comparison df

# visualising in a boxplot
finN_box <- rel.plot(rel_df, yval = relative_final_N)  # why this err?


  
finN_box <-  ggplot(rel_df, aes(x = Strategy, y = relative_final_N)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, aes(colour = "grey20") ) +
  labs(y = "Final population size relative to baseline average") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 
    # ggsave these to figs folder
ggsave(filename = "finalN_boxplot.png",
       plot = finN_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")


meanN_box <- ggplot(rel_df, aes(x = Strategy, y = relative_mean_N)) + 
  geom_boxplot(outlier.colour="red") + # need to remove outliers as these skew axis too much
  geom_hline(yintercept = 1, aes(colour = "grey20")) +
  labs(title = "Average population size relative to baseline average",
       y = "Relative mean pop size") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 


# comparing ssd and sex ratio across scenarios
# combine in list
rep_list <- list(rep_proj0, rep_proj1, rep_proj2, rep_proj3)

sex_list <- lapply(rep_list, ssd.av) 

sex_df <- sex.df(sex_list)  # custom function to turn lists into a dataframe


# how to visualise this?
# boxplot for sex ratio?
sr_box <- ggplot(sex_df, aes(x = Strategy, y = sex_ratio))+
  geom_boxplot(outlier.colour="red") +
  labs( y = "Sex ratio (proportion female)") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 

ggsave(filename = "sex_ratio_boxplot.png",
       plot = sr_box,
       device = "png",
       path = here("Figs"), 
       bg = "white") 
      


# looking at repeated rems --------
# I want a single box plot, colour coded by temporal frequency and seperated by strategy ()
# 
rel_du1 <- relative.pop(du_proj1,   
                          baseline_list = rep_proj0) 
rel_du2 <- relative.pop(du_proj2,   
                          baseline_list = rep_proj0) 
rel_du3 <- relative.pop(du_proj2,   
                          baseline_list = rep_proj0) 
du_list <- list(rel_du1, rel_du2, rel_du3) 

du_rel_df <- rel.df(du_list)

du_finN_box <- ggplot(du_rel_df, aes(x = Strategy, y = relative_final_N)) +
  geom_boxplot(outlier.colour="red") +
  geom_hline(yintercept = 1, aes(colour = "grey20") ) +
  labs(title = "Final population size relative to baseline average",
       y = "Relative final pop size") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 

ggsave(filename = "double_rem_finalN_boxplot",
       plot = du_finN_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")

du_meanN_box <- ggplot(du_rel_df, aes(x = Strategy, y = relative_mean_N)) + 
  geom_boxplot(outlier.colour="red") + # need to remove outliers as these skew axis too much
  geom_hline(yintercept = 1, aes(colour = "grey20")) +
  labs(title = "Average population size relative to baseline average",
       y = "Relative mean pop size") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 


# sex ratio 
du_rep_list <- list(du_proj1, du_proj2, du_proj3)

# du_sex_list <- lapply(du_list, ssd.av) 
# need to add baseline for comparison in plot
av_ssd0 <- sex_list[[1]]
av_du_ssd1 <- ssd.av(du_proj1, return.Mats = FALSE)
av_du_ssd2 <- ssd.av(du_proj2, return.Mats = FALSE)
av_du_ssd3 <- ssd.av(du_proj3, return.Mats = FALSE)

# combine in list with baseline 
# du_sex_box_list <- list(av_ssd0, lapply(du_list, ssd.av))
# du_sex_box_list[[1]] <- av_ssd0  # baseline has to be in first for projection labelling in correct order
# du_sex_box_list[[2:4]] <-  du_sex_list   #  ?

du_sex_list <- list(av_ssd0, av_du_ssd1, av_du_ssd2, av_du_ssd3)
du_sex_df <- sex.df(du_sex_list)  # success!


# boxplot for sex ratio?
du_sr_box <- ggplot(du_sex_df, aes(x = Strategy, y = sex_ratio))+
  geom_boxplot(outlier.colour="red") +
  labs(title = "Average Sex ratio across years",
       y = "Sex ratio averaged across years") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 

ggsave(filename = "double_rem_sex_ratio_boxplot",
       plot = du_sr_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")

# removals over 5 years
multi_rep_list <- list(multi_proj1, multi_proj2, multi_proj3)
rel_multi1 <- relative.pop(multi_proj1,   
                        baseline_list = rep_proj0) 
rel_multi2 <- relative.pop(multi_proj2,   
                        baseline_list = rep_proj0) 
rel_multi3 <- relative.pop(multi_proj2,   
                        baseline_list = rep_proj0) 
multi_list <- list(rel_multi1, rel_multi2, rel_multi3) 

multi_rel_df <- rel.df(multi_list)

multi_finN_box <- ggplot(multi_rel_df, aes(x = Strategy, y = relative_final_N)) +
  geom_boxplot(outlier.colour="red") +
  geom_hline(yintercept = 1, aes(colour = "grey20") ) +
  labs(title = "Final population size relative to baseline average",
       y = "Relative final pop size") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 

ggsave(filename = "multi_rem_finalN_boxplot",
       plot = multi_finN_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")

multi_meanN_box <- ggplot(multi_rel_df, aes(x = Strategy, y = relative_mean_N)) + 
  geom_boxplot(outlier.colour="red") + # need to remove outliers as these skew axis too much
  geom_hline(yintercept = 1, aes(colour = "grey20")) +
  labs(title = "Average population size relative to baseline average",
       y = "Relative mean pop size") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 

# sex ratio 
av_multi_ssd1 <- ssd.av(multi_proj1, return.Mats = FALSE)
av_multi_ssd2 <- ssd.av(multi_proj2, return.Mats = FALSE)
av_multi_ssd3 <- ssd.av(multi_proj3, return.Mats = FALSE)

# combine in list
multi_sex_list <- list(av_ssd0, av_multi_ssd1, av_multi_ssd2, av_multi_ssd3)
multi_sex_df <- sex.df(multi_sex_list)  


# boxplot for sex ratio?
multi_sr_box <- ggplot(sex_df, aes(x = Strategy, y = sex_ratio))+
  geom_boxplot(outlier.colour="red") +
  labs(title = "Average Sex ratio across years",
       y = "Sex ratio averaged across years") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 
ggsave(filename = "multi_rem_sex_ratio_boxplot",
       plot = multi_sr_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")




# combined df - 
sing <- rel_df
sing$rem_freq <- rep(as.character("Single", nrow(du_rel_df)))
duplo <- du_rel_df 
duplo$rem_freq <- rep(as.character("Double", nrow(du_rel_df)))
multi <- multi_rel_df 
multi$rem_freq <- rep(as.character("Multi", nrow(multi)))

comb_rel_df <- rbind(sing, duplo, multi)

# combined plots


# plot - how to make x axis wider
comb_finN_box <- ggplot(comb_rel_df, aes(x = Strategy, y = relative_final_N, fill = rem_freq)) +
  geom_boxplot(outlier.colour="red") +
  geom_hline(yintercept = 1, aes(colour = "grey20") ) +
  labs(y = "Final population size relative to baseline average") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 

ggsave(filename = "comb_finalN_boxplot.png",
       plot = comb_finN_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")



# checking for vulnerabilities (low pop sizes)
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

cols <- c("#1F78B4", "#EEAD0E", "#E31A1C")
vul_plot <- ggplot(vul_df, aes(x = Strategy , y =  Vulnerability, fill = Frequency)) +
  geom_bar(position = 'dodge', stat = "identity") +
  labs(title = "Extinction risk by strategy and frequency",
       y = "Extinction probability") +
  scale_fill_manual(values = cols) +   # colour assigned as single, dbl, multi
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 

ggsave(filename = "extinction_bar.png",
       plot = vul_plot,
       device = "png",
       path = here("Figs"), 
       bg = "white")


# combining final N in a table





## for appendix
# # lambda comparisons 
proj0 <- av_proj0$av_lambda
df <-  as.data.frame(proj0)   # val 1 = 1.036787
df$proj1 <- av_proj1$av_lambda
df$proj2 <- av_proj2$av_lambda
df$proj3 <- av_proj3$av_lambda   # 27, 30, 71 are extremely large/ low - what went wrong?

# merging into single col of lambda values split by Strategy
df_long <- gather(df, key = "Strategy" , value = "av_lambda") 

# boxlpot - comparing rem scenario lambda values
lamb_box <- ggplot(df_long, aes(x = Strategy, y = av_lambda)) +
  geom_boxplot(outlier.colour="red") 
# why are there sm outliers from all Strategies?
