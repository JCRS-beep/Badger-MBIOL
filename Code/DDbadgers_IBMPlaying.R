# Function Conversion to IBM
# Written by Jay C
# 13/11/26

# testing conversion of inital function from IBM file
out <- proj_meta_DDdisp_norem_IBM (threepatch_Umat,    
                            n0_3patch,
                            params,
                            stages,
                            npatches = 3,
                            time = 30,
                            DDapply="fertility",
                            dispersal_stages = c("Adult_f", "Adult_m"),
                            return.vec= TRUE)

(test_plot  <- meta_plot(out,   # output obj of dd.proj
                           y_val= "N",   # plot type - N or Vec 
                           ylab = "Population size", 
                           xlab = "time (t)",
                           rem_year = NULL,
                           mytheme = theme_classic(), 
                           cols= "black",    # can be vector of 2 cols
                           legend.pos = "top",
                           base_size = 16)) 


time <- as.numeric(c(0:(nrow(out$vec) - 1)) )      # vector 0:time

# for group size over time graph - must be as df
df <- as.data.frame(out$group)
colnames(df) <- paste0("Group", seq_len(ncol(df)))
df$Year <- time
group_df <- pivot_longer(df, cols = starts_with("Group"),
                         names_to = "Group", values_to = "Value") 
  
  # creating group projection plot
(base.plot <- ggplot(group_df, aes(Year, Value, colour = Group)) +
    geom_line(size = 1) +
    geom_point(size = 1.5) +
    labs(
      x = "xlab",
      y = "ylab")) 


# test with tighter params  - issues in groups overshooting 30 individuals, needs further tightening?
out_np <- proj_meta_DDdisp_norem_IBM (threepatch_Umat,    
                                   n0_3patch,
                                   meta_params,
                                   stages,
                                   npatches = 3,
                                   time = 30,
                                   DDapply="fertility",
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec= TRUE)

(np_plot  <- meta_plot(out_np,   # output obj of dd.proj
                         y_val= "N",   # plot type - N or Vec 
                         ylab = "Population size", 
                         xlab = "time (t)",
                         rem_year = NULL,
                         mytheme = theme_classic(), 
                         cols= "black",    # can be vector of 2 cols
                         legend.pos = "top",
                         base_size = 16)) 


df <- as.data.frame(out_np$group)
colnames(df) <- paste0("Group", seq_len(ncol(df)))
df$Year <- time
group_df <- pivot_longer(df, cols = starts_with("Group"),
                         names_to = "Group", values_to = "Value") 

# creating group projection plot
(base.plot <- ggplot(group_df, aes(Year, Value, colour = Group)) +
    geom_line(size = 1) +
    geom_point(size = 1.5) +
    labs(
      x = "xlab",
      y = "ylab")) 



# testing and improving claude written code
rem_test <- proj_meta_DDdisp_rem_IBM(threepatch_Umat,    
                                     n0_3patch,
                                     stagedist,
                                     params,
                                     stages,
                                     npatches = 3,
                                     time = 30,
                                     DDapply="fertility",
                                     dispersal_stages = c("Adult_f", "Adult_m"),
                                     remyear = c(10:13),  # removal year = vector of years 
                                     rem_strat = "random",  # not used
                                     bias = NULL,           # not used
                                     return.vec=TRUE, 
                                     return.remvec = TRUE)

(rem_plot  <- meta_plot(rem_test,   # output obj of dd.proj
                       y_val= "N",   # plot type - N or Vec 
                       ylab = "Population size", 
                       xlab = "time (t)",
                       rem_year = NULL,
                       mytheme = theme_classic(), 
                       cols= "black",    # can be vector of 2 cols
                       legend.pos = "top",
                       base_size = 16)) 


df <- as.data.frame(rem_test$group)
colnames(df) <- paste0("Group", seq_len(ncol(df)))
df$Year <- time
group_df <- pivot_longer(df, cols = starts_with("Group"),
                         names_to = "Group", values_to = "Value") 

# creating group projection plot
(baseRem_op.plot <- ggplot(group_df, aes(Year, Value, colour = Group)) +
    geom_line(size = 1) +
    geom_point(size = 1.5) +
    labs(
      x = "xlab",
      y = "ylab")) 

# saving plot to show
ggsave(filename = "metaIBM_DDdisp_rem_oldparams.png",
       plot = baseRem_op.plot,
       device = "png",
       path = here("Figs"), 
       bg = "white")

# removals with new parameters
rem_mp_test <- proj_meta_DDdisp_rem_IBM(threepatch_Umat,    
                                     n0_3patch,
                                     stagedist,
                                     meta_params,
                                     stages,
                                     npatches = 3,
                                     time = 30,
                                     DDapply="fertility",
                                     dispersal_stages = c("Adult_f", "Adult_m"),
                                     remyear = c(10:13),  # removal year = vector of years 
                                     rem_strat = "random",  # not used
                                     bias = NULL,           # not used
                                     return.vec=TRUE, 
                                     return.remvec = TRUE)

(rem_plot  <- meta_plot(rem_mp_test,   # output obj of dd.proj
                        y_val= "N",   # plot type - N or Vec 
                        ylab = "Population size", 
                        xlab = "time (t)",
                        rem_year = NULL,
                        mytheme = theme_classic(), 
                        cols= "black",    # can be vector of 2 cols
                        legend.pos = "top",
                        base_size = 16)) 


df <- as.data.frame(rem_mp_test$group)
colnames(df) <- paste0("Group", seq_len(ncol(df)))
df$Year <- time
group_df <- pivot_longer(df, cols = starts_with("Group"),
                         names_to = "Group", values_to = "Value") 

# creating group projection plot
(baseRem_np.plot <- ggplot(group_df, aes(Year, Value, colour = Group)) +
    geom_line(size = 1) +
    geom_point(size = 1.5) +
    labs(
      x = "xlab",
      y = "ylab")) 

# saving plot to share
ggsave(filename = "metaIBM_DDdisp_rem_newparams.png",
       plot = baseRem_np.plot,
       device = "png",
       path = here("Figs"), 
       bg = "white")