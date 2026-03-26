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
source(here("Code/03_model_projections.R"))  # should load all params, model projections
# issues - if the projection reaches 0 and stops, this means objs are not loaded

# need to combine summaries into a dataframe to be plotted together 
# df1 = av lambda per year for all projections (incl proj0)
# df2 = comparison of proportion final N and av N compared to baseline?  (plotted with hline = 1, anything above has increased, anything below has decreased)


# adding relative mean N and final N in a df
# easier = set up 3 df and rbind?  combining in projection df first, then merging
rel_proj1 <- relative.pop(rep_proj1,   
                          baseline_list = rep_proj0) 
rel_proj2 <- relative.pop(rep_proj2,   
                          baseline_list = rep_proj0) 
rel_proj3 <- relative.pop(rep_proj3,   
                          baseline_list = rep_proj0) 


# function to generalise this process
rel.df <- function(rel_projs = "list")  # list of all projections to compare, length = n projs
{
  nProj <- length(rel_projs)  # how many projections do we have?
  relative_df <- data.frame()  # to store other dfs in
  
  # set up individual dfs for eac proj
  for (p in 1:nProj){   # repeat for each projection
    Projection <- as.character(rep(p, 100))   # better with an informative name (projection1)
    relative_mean_N <- rel_projs[[p]]$relative_meanN
    relative_final_N <- rel_projs[[p]]$fin_props
    df <- data.frame(Projection, relative_mean_N, relative_final_N)  # what to do with this df? Store in list?
    
    # storing df in our relative df - intial 
      relative_df <- rbind(relative_df, df)
    }
  return(relative_df)
}

rel_projs <- list(rel_proj1, rel_proj2, rel_proj3)
rel_df <- rel.df(rel_projs)

# visualising in a boxplot
finN_box <- ggplot(rel_df, aes(x = Projection, y = relative_final_N)) +
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
  

meanN_box <- ggplot(rel_df, aes(x = Projection, y = relative_mean_N)) + 
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


# comparing ssd and sex ratio across scenarios?
av_ssd0 <- ssd.av(rep_proj1, return.Mats = TRUE)
av_ssd1 <- ssd.av(rep_proj1, return.Mats = TRUE)
av_ssd2 <- ssd.av(rep_proj2, return.Mats = TRUE)
av_ssd3 <- ssd.av(rep_proj3, return.Mats = TRUE)

# combine in list
sex_list <- list(av_ssd0, av_ssd1, av_ssd2, av_ssd3)

# function to turn lists into a dataframe
sex.df <- function(sex_projs = "list")  # list of all projections to compare, length = n projs
  {
  nProj <- length(sex_projs)  # how many projections do we have?
  sex_df <- data.frame()  # to store other dfs in
  
  # set up individual dfs for eac proj
  for (p in 1:nProj){   # repeat for each projection
    n <- p-1   # incl proj0, need to name 0
    Projection <- as.character(rep(n, 100))   # better with an informative name (projection1)
    av_prop <- sex_projs[[p]]$av
    sex_ratio <- sex_projs[[p]]$sex_ratio
    
    df <- data.frame(Projection, av_prop, sex_ratio)  # what to do with this df? Store in list?
    
    # storing df in our relative df - intial 
    sex_df <- rbind(sex_df, df)
  }
  return(sex_df)
}

sex_df <- sex.df(sex_list)  # success!


# how to visualise this?
# boxplot for sex ratio?
summary(sex_df$sex_ratio)   # seems fairly constant?
sr_box <- ggplot(sex_df, aes(x = Projection, y = sex_ratio))+
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




# lambda comparisons - not so important, but leave for appendix?
proj0 <- av_proj0$av_lambda
df <-  as.data.frame(proj0)   # val 1 = 1.036787
df$proj1 <- av_proj1$av_lambda
df$proj2 <- av_proj2$av_lambda
df$proj3 <- av_proj3$av_lambda   # 27, 30, 71 are extremely large/ low - what went wrong?

# merging into single col of lambda values split by projection
df_long <- gather(df, key = "Projection" , value = "av_lambda") 

# boxlpot - comparing rem scenario lambda values
lamb_box <- ggplot(df_long, aes(x = Projection, y = av_lambda))+
  geom_boxplot(outlier.colour="red") 
# why are there sm outliers from all projections?


