# just visualising a normal dist



library(ggplot2)

x <- rnorm(1000, mean = 0.1, sd= 0.05) # sample of 1000 
# trim 0 - 1
x[x < 0] <- 0
x[x > 1] <- 1
df <- as.data.frame(x)


n_obs <-  sum(!is.na(df$x))
bw <- 0.05



base.plot <- ggplot(df, aes(x)) + 
  geom_histogram(aes(y = after_stat(density)),  
                 binwidth = bw,
                 fill = "grey50", alpha= 0.3) + 
  stat_function(fun = dnorm, args = list(mean = mean(df$x), sd = sd(df$x)), 
                size = 1.5, colour = "black")+ 
  geom_vline(aes(xintercept = 0.1),                       # Adding a line to show removal year
             colour = "red", linetype = "dashed", linewidth= 1.5, alpha = 0.8) +
  labs(x = "Removal probability") +
  theme_classic(base_size = 22) +
  theme(
    axis.text = element_text(size = 22, color = "black"), 
    
  )

# rescaling y
ybreaks = seq(0,1000,200) 
# On primary axis
base.plot + 
  scale_y_continuous("", breaks = round(ybreaks / (bw * n_obs),3), labels = ybreaks )


# remove y axis?


# biased trials = multiple plots needed with varying means
df$unbias <- rnorm(1000, mean = 0.05, sd= 0.05)
df$bias <- rnorm(1000, mean = 0.2, sd= 0.05)

# graph limits
means <- sapply(df[c("x", "unbias", "bias")], mean, na.rm = TRUE)
sds   <- sapply(df[c("x", "unbias", "bias")], sd,   na.rm = TRUE)
min_x <- min(means - 4 * sds)
max_x <- max(means + 4 * sds)

# want to add vis of these extra means on top of base plot
bias.plot <- ggplot(df, aes(x)) + 
  geom_histogram(aes(y = after_stat(density)),  
                 binwidth = bw,
                 fill = "grey50", alpha= 0.15) + 
 # stat_function(fun = dnorm, args = list(mean = mean(df$x), sd = sd(df$x)), 
  #              size = 1.2, colour = "black", alpha = 0.5) + 
  #geom_vline(aes(xintercept = mean(df$x)),                       # Adding a line to show removal year
   #          colour = "black", linetype = "dashed", size= 1, alpha = 0.5) +
  
  stat_function(fun = dnorm, args = list(mean = mean(df$unbias), sd = sd(df$unbias)), 
                size = 1.2, colour = "blue") + 
  geom_vline(aes(xintercept = mean(df$unbias)),                       # Adding a line to show removal year
             colour = "blue", linetype = "dashed", size= 1, alpha = 0.5) +
  
  stat_function(fun = dnorm, args = list(mean = mean(df$bias), sd = sd(df$bias)), 
                size = 1.2, colour = "red") + 
  geom_vline(aes(xintercept = mean(df$bias)),                       # Adding a line to show removal year
             colour = "red", linetype = "dashed", size= 1, alpha = 0.5) +
  
  labs(x = "Removal probability") +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(color = "black")
  )+
 # coord_cartesian(xlim = c(min_x, max_x))  + # makes sure full range of x axis included
  scale_y_continuous("", breaks = round(ybreaks / (bw * n_obs),3), labels = ybreaks )
# almost works, why are my curves not complete?




# plot example for trajectories
# removal plot
x <- rep_proj1[[1]]
base_size = 16
rem_year = 5

# 3 alt dfs with varying trajectories after removals
t <- nrow(x$vec) - 1                 # t= n years (0-t = t+1 entries)
time <- as.numeric(c(0:t))       # vector 0:t

# for pop size over time graph - must be as df
  pop_df <- data.frame(Year = time,
                       Pop  = x$pop )# converting pop to a df
  
  # creating pop projection plot
p <- ggplot(pop_df, aes(x = Year, y = Pop)) +
  geom_line(size = 1.2, colour = "black") +
  geom_point(size = 2, colour = "grey7") +
  labs(
    x = "Year",
    y = "Population size") +
  scale_y_continuous(
    name = "Population size",
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.1, 0.2)))+  # 5% below, 30% above
  theme_classic() +
  theme(
    text = element_text(size = base_size),
    plot.title = element_text(size = base_size + 2, face = "bold"),
    axis.title = element_text(size = base_size),
    axis.text = element_text(size = base_size - 2),
    legend.position = "top"
  ) +
    geom_vline(xintercept = 5,                       # Adding a line to show removal year
               colour = "red3", 
               linetype = "dashed", linewidth = 1, 
               alpha = 0.5) 

head(pop_df)
pop_df$Pop_low <- rep(0, 21)
pop_df$Pop_low[1:6] <- pop_df$Pop[1:6]

pop_df$Pop_low[7:14] <- floor(seq(from = 25, to = 0, length.out = 8))
# creating a sequence as a curve
t = c(1:6)
y = 35 - 15 * log(t, base = exp(1))
pop_df$Pop_low[7:12] <- y

#plot(x = pop_df$Year , y = exp(-pop_df$Year))
#plot(x = pop_df$Year , y = 35/(0.5 *pop_df$Year))
#plot(x = pop_df$Year , y = 120 - 40 *log(pop_df$Year + 6))


pop_df$Pop_mid <- pop_df$Pop
pop_df$Pop_mid[7:21] <- seq(from = 40, to = 75, length.out = 15)

p + 
  geom_point(data = pop_df, aes(x = Year, y = Pop_low), 
             size = 2, colour = "red", alpha = 0.6) +
  geom_line(data = pop_df, aes(x = Year, y = Pop_low),  
            size = 1.2, colour= "red", alpha = 0.7) +
  geom_point(data = pop_df, aes(x = Year, y = Pop_mid), 
             size = 2, colour = "orange", alpha = 0.6) +
  geom_line(data = pop_df, aes(x = Year, y = Pop_mid),  
            size = 1.2, colour= "orange", alpha = 0.7) 
  
  