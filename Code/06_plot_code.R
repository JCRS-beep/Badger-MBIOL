# all plots from my project
# 
# basic visualisations = the trajectory of each scenario investigated --------

# setting colours for male and female plots
col_vec <- c("#FF6A6A", "#87CEEB")

# baseline projection = no removals, 20 years
# stage abundance over time - remove titles
proj0_plot <- dd.plot(proj0, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = NULL,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16)
proj0_Nplot <- dd.plot(proj0, 
                        y_val= "N", 
                        ylab = "Population size",
                       xlab = "Time (t)",
                       rem_year = NULL,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16)
multipanel <- proj0_Nplot + proj0_plot
# saving tihs plot
ggsave(filename = "baseline_projection.png",
       plot = multipanel,
       device = "png",
       path = here("Figs"), 
       bg = "white")


    

    
scale_factor <- max(pop_df$Pop, na.rm = TRUE) / max(df_long$Abundance, na.rm = TRUE)
(comb_plot <- ggplot() +    # improve by making y axis range bigger - gap between pop size and abundance. Change shape of pop size
  # Population size (primary axis)
  geom_line(data = pop_df, aes(x = Year, y = Pop),
            size = 1.2, colour = "black") +
  geom_point(data = pop_df, aes(x = Year, y = Pop),
             size = 3, 
             shape = 3,
             colour = "black") +
  # Abundance (rescaled)
  geom_line(data = df_long, aes(x = Year, 
                                y = Abundance * scale_factor,
                                colour = Sex, 
                                linetype = Stage),
            size = 1.2, alpha = 0.7) +
  geom_point(data = df_long, 
             aes(x = Year, 
                 y = Abundance * scale_factor,
                 colour = Sex,
                 shape = Stage),
             size = 2, alpha = 0.8,
             position = position_jitter(width = 0.2)) +
  
  # Dual axis
# scale_y_continuous(name = "Population size", sec.axis = sec_axis(~ . / scale_factor, name = "Abundance")) +
  scale_y_continuous(
    name = "Population size",
    expand = expansion(mult = c(0.05, 0.3)),  # 5% below, 30% above
    sec.axis = sec_axis(~ . / scale_factor, name = "Abundance")
  ) +
  scale_colour_manual(values = col_vec,
                      labels = c("Female", "Male")) +
  labs(x = "time (t)",
       colour = "Sex",
       linetype = "Stage",
       shape = "Stage") +
  theme_classic() +
  theme(text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top"))





# analysis plots
# rel final N combined plot
(comb_rel_finN_box <- ggplot(comb_rel_df, aes(x = Strategy, y = relative_final_N, 
                                              fill = factor(Frequency,
                                                            levels = c("Single", "Double", "Multiple", "Continuous"),
                                                            labels = c("Single", "Double", "Multiple", "Continuous"))))+
    geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single", padding = 0.2),  # changing width changes 
                 width = 0.75,   # change whole plot width 
                 outlier.size = 1.2) +
    geom_hline(yintercept = 1, colour = "grey20",     # int = av final_N baseline
               linetype = "dashed",
               linewidth = 0.7)  +
    labs(y = "Final population size relative to baseline",
         fill = "Removal frequency") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 16) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.grid.major.x = element_blank(),
    ) )
ggsave(filename = "comb_relative_finalN_boxplot.png",
       plot = comb_rel_finN_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")


# sex ratio combined plor
(comb_sex_box <- ggplot(comb_sex_df, aes(x = Strategy, y = sex_ratio,     # why no adult fem cont plot?
                                         fill = factor(Frequency,
                                                       levels = c("Single", "Double", "Multiple", "Continuous"),
                                                       labels =c("Single", "Double", "Multiple", "Continuous")))) +
    geom_boxplot(position = position_dodge2(width = 1, preserve = "single", padding = 0.2),
                 width = 0.75,
                 outlier.size = 1.2) +
    geom_hline(yintercept = base_ratio, 
               linetype = "dashed",
               linewidth = 0.7) +  # line for baseline av sex ratio
    labs(y = "Proportion of females in the population", 
         fill = "Removal frequency") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 16) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.grid.major.x = element_blank(),
    ) )
ggsave(filename = "comb_sex_ratio_boxplot.png",
       plot = comb_sex_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")



# extinction risk plot

(vul_plot <- ggplot(vul_df, aes(x = Strategy , y =  Vulnerability, fill = Frequency)) +
    geom_bar(position = 'dodge', stat = "identity") +
    labs(y = "Extinction probability") +
    scale_fill_brewer(palette = "Set2") +   # colour assigned to removal frequency
    theme_minimal() +
    theme(
      text = element_text(size = 16),
      plot.title = element_text(size = 16 + 2, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16 - 2),
    )) 
ggsave(filename = "extinction_bar.png",
       plot = vul_plot,
       device = "png",
       path = here("Figs"), 
       bg = "white")
