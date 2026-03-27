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
# issues - du_proj3 not found?

# need to combine summaries into a dataframe to be plotted together 
# Functions used in this script -----
# Turning relative pop size outputs into comparison dataframe for plotting
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

# Turning list outputs into dataframe with av stage proportion per rep and sex ratio for comparison 
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

# average pop size and lambda of baseline -------
popN <- N.extract(rep_proj0)
Nfin <- sapply(popN, function(x) x[20])
summary(Nfin)

base_lamb <- lamb.av(rep_proj0)
av_lamb <- summary(base_lamb)

av_ssd0 <- ssd.av(rep_proj0)
base_prop <- colMeans(av_ssd0$av_prop)

# adding relative mean N and final N in a df
# easier = set up 3 df and rbind?  combining in projection df first, then merging
rel_proj1 <- relative.pop(rep_proj1,   
                          baseline_list = rep_proj0) 
rel_proj2 <- relative.pop(rep_proj2,   
                          baseline_list = rep_proj0) 
rel_proj3 <- relative.pop(rep_proj3,   
                          baseline_list = rep_proj0) 

rel_projs <- list(rel_proj1, rel_proj2, rel_proj3)

rel_df <- rel.df(rel_projs)  # using prev defined function to turn inot comparison df

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
 
   # ggsave these to figs folder
ggsave(filename = "finalN_boxplot",
       plot = finN_box,
       device = "png",
       path = here("Figs")
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


# comparing ssd and sex ratio across scenarios
# combine in list
rep_list <- list(rep_proj0, rep_proj1, rep_proj2, rep_proj3)

sex_list <- lapply(rep_list, ssd.av) 

sex_df <- sex.df(sex_list)  # custom function to turn lists into a dataframe


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

ggsave(filename = "sex_ratio_boxplot",
       plot = finN_box,
       device = "png",
       path = here("Figs")
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

du_finN_box <- ggplot(du_rel_df, aes(x = Projection, y = relative_final_N)) +
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
       path = here("Figs")
)

du_meanN_box <- ggplot(du_rel_df, aes(x = Projection, y = relative_mean_N)) + 
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
# du_list <- list(du_proj1, du_proj2, du_proj3)

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
du_sr_box <- ggplot(du_sex_df, aes(x = Projection, y = sex_ratio))+
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
       path = here("Figs")
)

# removals over 5 years
rel_multi1 <- relative.pop(multi_proj1,   
                        baseline_list = rep_proj0) 
rel_multi2 <- relative.pop(multi_proj2,   
                        baseline_list = rep_proj0) 
rel_multi3 <- relative.pop(multi_proj2,   
                        baseline_list = rep_proj0) 
multi_list <- list(rel_multi1, rel_multi2, rel_multi3) 

multi_rel_df <- rel.df(multi_list)

multi_finN_box <- ggplot(multi_rel_df, aes(x = Projection, y = relative_final_N)) +
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
       path = here("Figs")
)

multi_meanN_box <- ggplot(multi_rel_df, aes(x = Projection, y = relative_mean_N)) + 
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
multi_sex_df <- sex.df(multi_sex_list)  # success!


# boxplot for sex ratio?
multi_sr_box <- ggplot(sex_df, aes(x = Projection, y = sex_ratio))+
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
       path = here("Figs")
)

# combined df
duplo <- du_rel_df 
duplo$rem_freq <- rep(as.character(2, nrow(du_rel_df)))
multi <- multi_rel_df 
multi$rem_freq <- rep(as.character(5, nrow(multi)))

comb_rel_df <- rbind(duplo, multi)

# combined plots
comb_finN_box <- ggplot(comb_rel_df, aes(x = Projection, y = relative_final_N, fill = rem_freq)) +
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

ggsave(filename = "comb_finalN_boxplot",
       plot = comb_finN_box,
       device = "png",
       path = here("Figs")
)

# checking for vulnerabilities (low pop sizes)
