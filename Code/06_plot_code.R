# all plots from my project
# 
# basic visualisations = the trajectory of each scenario investigated --------
av_baseplot <-  ggplot() +
  geom_line(data = pop_df, 
            aes(mapping = rep, 
               x = Year, 
               y = N), 
            size = 1.2, 
            colour = "grey80", 
            alpha = 0.2) +
  labs(x = "time",
       y = "Population size per repetition")  +      # add legend showing dots = mean?
  geom_point(data = N_dist, aes(x = Year, y = av_N), 
             size = 2, colour = "black", alpha = 0.7) +
  geom_line(data = N_dist, aes(x = Year, y = av_N),
            size = 1.2, colour = "black", alpha = 0.8) +
  geom_errorbar(data = N_dist, 
                aes(ymin = (av_N - sd), ymax = (av_N + sd),
                    y = av_N, 
                    x = Year), 
                width = 0.2,
                alpha = 0.6,
                position=position_dodge(0.05)
  )+               # include sd as err bars
  labs(
    x = "Year",
    y = "Population size") +
  scale_y_continuous(
    name = "Population size",
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.2)),
    limits = c(0,NA))+  
  theme_classic() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16 + 2, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  ) 



v_baseplot <-  ggplot() +   #  increase point size, remove jitter, increase point size
  geom_line(data= df_long, 
            aes(x =Year, 
                y =Abundance, 
                colour= Sex, # pink or blue
                mapping = Stage,
                mapping = rep),
            size = 1.1, 
            alpha=0.1) +
  geom_point(data = av_df, 
             aes(x = Year, 
                 y = Av_abundance,
                 shape = Stage), 
             size = 2.5, 
             alpha=0.7
  ) + 
  geom_line(data= av_df,
            aes(x=Year, 
                y=Av_abundance, 
                mapping= Stage, 
                mapping= Sex
            ),
            size = 1,
            alpha=0.6) +
  geom_errorbar(data = av_df,      # err bars for each sex and stage
                aes( x= Year,
                     ymin = Av_abundance - sd, ymax = Av_abundance + sd), 
                width =.2,
                alpha = 0.6,
                position = position_dodge(0.05)) +
  labs(
    x = "Year",
    y = "Abundance",
    shape = "Stage", 
    colour = "Sex"
  ) +
  scale_y_continuous(
    name = "Abundance",
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.05)),  # 5% below, 30% above
    limits = c(0, NA)) +
  theme_classic() +
  theme(
    text = element_text(size = base_size),
    plot.title = element_text(size = base_size + 2, face = "bold"),
    axis.title = element_text(size = base_size),
    axis.text = element_text(size = base_size - 2),
    legend.position = "top",
    # journal-style compact legend:
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(face = "bold", size = base_size),
    legend.text = element_text(size = base_size - 2),
    legend.spacing.x = unit(0.2, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  # make legend compact and show titles above keys (journal style)
  guides(
    colour = guide_legend(title.position = "top",
                          title.hjust = 0.5,
                          nrow = 1,
                          byrow = TRUE,
                          override.aes = list(size = 3, linetype = 1, shape = 16, alpha = 1)),
    shape = guide_legend(title.position = "top",
                         title.hjust = 0.5,
                         nrow = 1,
                         byrow = TRUE,
                         override.aes = list(size = 3, alpha = 1))
    
  )


base_joint_plot <- av_baseplot + v_baseplot         # improve = gap between plots



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
      axis.title.y = element_text(vjust = 2),
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
    labs(y = "Proportion of Adult females in the population", 
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
