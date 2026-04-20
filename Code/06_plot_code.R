# all plots from my project
# 
# basic visualisations = the trajectory of each scenario investigated --------

# setting colours for male and female plots
col_vec <- c("#FF6A6A", "#87CEEB")

# baseline projection = no removals, 20 years
# stage abundance over time
(proj0_plot <- dd_plot(proj0, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = NULL,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16))






# analysis plots
# final N boxplots

# final N combined plots
comb_finN_box <- ggplot(comb_rel_df, aes(x = Strategy, y = relative_final_N, fill = rem_freq)) +
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
ggsave(filename = "comb_finalN_boxplot.png",
       plot = comb_finN_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")


# sex ratio
comb_sex_box <- ggplot(comb_sex_df, aes(x = Strategy, y = sex_ratio, 
                                        fill = factor(rem_freq,
                                                      levels = c("Single", "Double", "Multiple"),
                                                      labels = c("Single", "Double", "Multiple")))) +
  geom_boxplot(position = position_dodge2(width = 1, preserve = "single", padding = 0.2),
               width = 0.75,
               outlier.size = 1.2) +
  geom_hline(yintercept = base_ratio, 
             linetype = "dashed",
             linewidth = 0.7) +  # line for baseline av sex ratio
  labs(y = "Female proportion of the population", 
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
  ) 
ggsave(filename = "comb_sex_ratio_boxplot.png",
       plot = comb_sex_box,
       device = "png",
       path = here("Figs"), 
       bg = "white")



# extinction risk plot
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
