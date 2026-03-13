# all plots from my project
# 
# basic visualisations = the trajectory of each scenario investigated

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

# population size over time
(proj0_Nplot <- dd_plot(proj0, 
                        y_val= "N", 
                        ylab = "Abundance", 
                        xlab = "Time (t)",
                        rem_year = NULL,
                        mytheme = theme_classic(), 
                        cols= col_vec,    # can be vector of cols
                        legend.pos = "top" ,
                        base_size = 16))


# scenario 1 = 70% random removals 
# Stage Abundance
(proj1_plot <- dd_plot(proj1, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 5,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16))
# Population size
(proj1_Nplot <- dd_plot(proj1, 
                       y_val= "N", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 5,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "top",
                       base_size = 16))

# Scenario 2 = biased male removals (free shooting to target moving individuals)
(proj2_plot <- dd_plot(proj2, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       base_size = 16))
(proj2_plot <- dd_plot(proj2, 
                       y_val= "N", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       base_size = 16))

# Scenario 3 = Sett killings (adult females)
(proj3_plot <- dd_plot(proj3, 
                       y_val= "Vec", 
                       ylab = "Abundance", 
                       xlab = "Time (t)",
                       rem_year = 10,
                       mytheme = theme_classic(), 
                       cols= col_vec,    # can be vector of cols
                       legend.pos = "topright",
                       base_size = 16))
