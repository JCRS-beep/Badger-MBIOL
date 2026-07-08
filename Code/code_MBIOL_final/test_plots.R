# example graphs


# pop removals
# Visualising how removals may impact disease prevalence----
rem <- dd.plot(rep_proj1[[1]],
        y_val= "N",   # plot type - N or Vec 
        ylab = "abundance", 
        xlab = "time (t)",
        rem_year = 5,
        mytheme = theme_classic(), 
        cols= "black",    # can be vector of 2 cols
        legend.pos = "top",
        base_size = 16)



t <- nrow(out$vec) - 1                 # t= n years (0-t = t+1 entries)
time <- as.numeric(c(0:t))       # vector 0:t

# for pop size over time graph - must be as df
  pop_df <- data.frame(Year = time,
                       Pop  = out$pop )# converting pop to a df
  
d = vector()
d[c(1:5)] <- 0.7 *(pop_df$Pop[c(1:5)])
d[c(6:11)] <- seq(from = 30, to = 20, by = -2)
d[c(11:21)] <- seq(from = 19, to = 8, by = -1.1)
  # adding hypothetical disease line
  rem +
    geom_line(aes(d) ) # follow pop size until year 5, then slowly decrease
  
  rem+
    geom_line(aes(x = pop_df$Year,y =  d), 
              colour = "red3", 
              linetype = "dashed", linewidth = 1.5)
  
  # creating pop projection plot
  base.plot <- ggplot(pop_df, aes(x = Year, y = Pop)) +
    geom_line(size = 1.2, colour = "grey30") +
    geom_point(size = 2, colour = "black") +
    labs(
      x = xlab,
      y = ylab) +
    scale_y_continuous(
      name = "Population size",
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0.1, 0.2)))+  # 5% below, 30% above
    mytheme +
    theme(
      text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2, face = "bold"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 2),
      legend.position = "top"
    ) 
  
  if (!is.null(rem_year)){
    plot <- base.plot +
      geom_vline(xintercept = rem_year,                       # Adding a line to show removal year
                 colour = "red3", 
                 linetype = "dashed", linewidth = 1, 
                 alpha = 0.5) 
    
  } else if (is.null(rem_year)){
    plot <- base.plot 
  }
  
rem_plot <-  ggplot() +    # improve by making y axis range bigger - gap between pop size and abundance. Change shape of pop size
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
             position = position_jitter(width = 0.2))









# plotting baseline projection averages ----
# pop projection line with err, vec with err
popN <- N.extract(rep_proj0)
# turn into df - rows = rep (list number) =  c("rep", "Year", "pop_N")
t <- length(popN[[1]]) - 1                 # t= n years (0-t = t+1 entries)
time <- as.numeric(c(0:t))       # vector 0:t

# for pop size over time graph - must be as df
pop_df <- data.frame(Year = rep(time, 100))   # empty df

# fill each row with rep number (r), year and N
for (r in 1:length(popN)){   # for each rep
  r_lims <- c((21*r-20):(21*r))     # rows to fill with each rep
  pop_df$rep[r_lims] <- as.character(rep(r, 21))
  pop_df$N[r_lims] <- popN[[r]]
}

base_size <- 16

# plotting N across years

base_line_plot <- ggplot(pop_df, aes(x = Year, y = N)) +
  geom_line(data = pop_df, aes(mapping = rep), 
            size = 1.2, colour = "grey80", alpha = 0.2) +
  labs(x = "time",
       y = "Population size per repetition") +
  scale_y_continuous(
    name = "Population size",
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.05)),  # 5% below, 30% above
    limits = c(0, 250)) +  
  theme_classic() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16 - 2),
  )


# adding line for average pop size each year, include error
av_N <- sapply(popN, mean) # av pop size PER REP


year_av_N <- function(n_list){   # list of pop sizes per year, calculate  mean per year with some measure of spread
 av_n <- vector()     # yearly mean n 
 sd_n <- vector()
 
   for(i in 1:(t+1)){
    year_N <- sapply(n_list, function(x) x[i])  # extract all pop sizes this year
   av_n[i] <- mean(year_N)
   sd_n[i] <- sd(year_N)    # sd of  yearly pop sizes
   }
 
 t <- length(n_list[[1]]) - 1                 # t= n years (0-t = t+1 entries)

 df <- data.frame(av_N = av_n, sd = sd_n)
 df$Year <- as.numeric(c(0:t))    
 
  return(df)
} 
  
N_dist <- year_av_N(popN)    



# av pop size per year plotted with stadard dev
av_plot <- ggplot(data = N_dist, aes(x = Year, y = av_N)) +
  geom_point() +
  geom_line(size = 1.2, colour = "black", alpha = 0.7) +
  geom_errorbar(data = N_dist, 
                aes(ymin= av_N - sd, ymax= av_N + sd), 
                width= 0.2,
                alpha = 0.6,
                position=position_dodge(0.05))  +                  # include sd as err bars
  labs(
    x = "time",
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

  
av_baseplot <-   base_line_plot +      # add legend showing dots = mean?
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




# doing same for stage abundances ---
v.extract <- function(x){     # like N extract but separates stage abundances from each lsit
  reps <- length(x)
  v_list <- list()   # list of N vecs for each repeat
  
  for (r in 1:reps){
    v_list[[r]] <- x[[r]]$vec    # isolating pop size vector
  }
  
  return(v_list)
}

v_list <- v.extract(rep_proj0)

x_val <- ncol(v_list[[1]])   # number of classes and sexes (if nStages = 2 and sex =2, x =4)

# turning into dataframe - rep, year and abundance for each stage 
vec_df <- data.frame(Year = rep(time, 100)) # years repeated per rep

for (r in 1:length(v_list)){   # for each rep
  r_lims <- c((21*r-20):(21*r))     # rows to fill with each rep
  
  vec_df$rep[r_lims] <- as.character(rep(r, 21))  # filling in rep column with this loop iteration
  vec_df$Yearling_Female[r_lims] <- v_list[[r]][,1]
  vec_df$Adult_Female[r_lims] <- v_list[[r]][,2]
  vec_df$Yearling_Male[r_lims] <- v_list[[r]][,3]
  vec_df$Adult_Male[r_lims] <- v_list[[r]][,4] 
}

# long format so each row is a single observation 
df_long <- gather(vec_df, key = "Stage", value = "Abundance", c(3,4,5,6))   # creating a stage col in df with abundance
df_long <- separate(df_long, 
                    col= "Stage", 
                    into= c("Stage", "Sex"), 
                    sep='_')   # splitting by sex, separated by _

# plotting graph 
vec_plot <- ggplot(data= df_long) +  # sexes diff cols, shapes and lines diff for stages
 # geom_point(position= "jitter", size = 2, alpha=0.3) +  # jitter to avoid overlap of yearlings
  geom_line(data= df_long, aes(x=Year, 
                               y=Abundance, 
                               colour= Sex, 
                               linetype= Stage,
                               mapping = rep),
            size = 1.2, alpha=0.1) +
  labs(
    x = "Year",
    y = "Abundance",
    colour = "Sex",
    linetype = "Stage") +
  
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
    linetype = guide_legend(title.position = "top",
                            title.hjust = 0.5,
                            nrow = 1,
                            byrow = TRUE,
                            override.aes = list(size = 1.2, alpha = 1))
   
  )
 

# using mean per year and err bars --------
# calculating mean per year per stage - storing in df ready to plot
names <- c("Yearling_Female", "Adult_Female",    "Yearling_Male", "Adult_Male"  )


# function to store in df
year_av_V <- function(v_list, names){   # list of pop sizes per year, calculate  mean per year with some measure of spread
  t <- nrow(v_list[[1]]) - 1                 # t= n years (0-t = t+1 entries)
  time <- c(0:t)
  
  av_v <- matrix(0, ncol = 4, nrow = length(time))     # yearly mean n 
  sd_v <- matrix(0, ncol = 4, nrow = length(time))     
  
  # naming cols and rows
  rownames(av_v) <- time
  rownames(sd_v) <- time
  colnames(av_v) <- names
  colnames(sd_v) <- names
  
  for(i in 1:(t+1)){  # for each row representing year 0 to t
    year_v <- sapply(v_list, function(x) x[i,])  # extract stage abundance this year across 100 reps - matrix with 4 rows, 100 cols
   
    av_v[i,] <- round(rowMeans(year_v), 1)         # each row = stage, row means = mean per stage per year
   
    for (n in 1:nrow(year_v)){
      sd_v[i,n] <- round(sd(year_v[n,]), 2)   # sd of each row of
    }
  }

 
  # creating df with av abundance and sd for each stage and sex
  v_df <- data.frame(av_v = av_v)    # matrices into df - how to stop naming as av.stage?
  colnames(v_df) <- names
  v_df$Year <- time
  
  # sd df
  v_df2 <- data.frame(sd = sd_v)
  colnames(v_df2) <- names
  v_df2$Year <- time
   
  # long format so each row is a single observation 
  v_long <- gather(v_df, key = "Stage", value = "Av_abundance", c(1:4))     # creating a stage col in df with abundance
  # why is col name for year now NA?
  v_long2 <- gather(v_df2, key = "Stage", value = "sd", c(1:4))
  
  v_big <- merge(v_long, v_long2, by = c("Stage", "Year")) 
  
  v_big <- separate(v_big, 
                     col= "Stage", 
                     into= c("Stage", "Sex"), 
                     sep='_')   # splitting by sex, separated by _
  
  return(v_big)
} 

av_df <-  year_av_V(v_list, names)



# plotting av points and line with err bars - linetype = sex ?
av_vplot <- ggplot(data= av_df) +  # sexes diff cols, shapes and lines diff for stages
  geom_point(aes(x = Year, y = Av_abundance, shape = Stage,), 
             position= "jitter", 
             size = 3, alpha=0.7) +  # jitter to avoid overlap of yearlings
  geom_line(data= av_df, aes(x=Year, 
                             y=Av_abundance, 
                             colour= Stage, 
                             linetype= Sex),
            size = 1.2, alpha=0.8) +
  labs(
    x = "Year",
    y = "Abundance",
    colour = "Stage",
    shape = "Stage",
    linetype = "Sex") +
  scale_colour_manual(values = c("black", "grey1")) +
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
    shape = guide_legend(title.position = "top",
                         title.hjust = 0.5,
                          nrow = 1,
                          byrow = TRUE,
                          override.aes = list(size = 3, alpha = 1)),
    linetype = guide_legend(title.position = "top",
                            title.hjust = 0.5,
                            nrow = 1,
                            byrow = TRUE,
                            override.aes = list(size = 1.2, alpha = 1))
  )



v_baseplot <-  ggplot() +   # keep shape same, yearlings dotted line?
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
                 shape = Stage
                 ), 
             size = 2, 
             alpha=0.7
             ) + 
  geom_line(data= av_df,
            aes(x=Year, 
                y=Av_abundance, 
                linetype= Stage, 
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
ggsave(filename = "base_ranges.png",
       plot = base_joint_plot,
       device = "png",
       path = here("Figs"), 
       bg = "white")



