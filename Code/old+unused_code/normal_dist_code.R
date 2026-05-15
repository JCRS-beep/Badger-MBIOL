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
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(color = "black")
  )

# rescaling y
ybreaks = seq(0,1000,200) 
# On primary axis
base.plot + 
  scale_y_continuous("", breaks = round(ybreaks / (bw * n_obs),3), labels = ybreaks )


# remove y axis?


# biased trials = multiple plots needed with varying means
df$unbias <- rnorm(1000, mean = 0.35, sd= 0.1)
df$bias <- rnorm(1000, mean = 0.65, sd= 0.1)

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
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(color = "black")
  )+
 # coord_cartesian(xlim = c(min_x, max_x))  + # makes sure full range of x axis included
  scale_y_continuous("", breaks = round(ybreaks / (bw * n_obs),3), labels = ybreaks )
# almost works, why are my curves not complete?