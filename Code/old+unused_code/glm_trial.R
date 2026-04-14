# testing conditions and movement link
install.packages("msm")
library(msm) # truncated normal distribution - could just clip regular
library(ggplot2)

# for now, movement prob sampled from random dist 0-1
# p(move) for each stage or patch? stage and patch?
prob <- rtnorm(1000, mean = 0.5, sd= 0.3) # sample of 1000 
prob_df <- as.data.frame(prob)

n_obs <-  sum(!is.na(prob_df$prob))

norm.plot <- ggplot(prob_df, aes(prob))  + 
  geom_histogram(aes(y = after_stat(density)),  binwidth = 0.05) + 
  stat_function(fun = dnorm, args = list(mean = mean(prob_df$prob), sd = sd(prob_df$prob)))

# rescaling y
ybreaks = seq(0,50,5) 
# On primary axis
norm.plot + scale_y_continuous("Counts", breaks = round(ybreaks / (bw * n_obs),3), labels = ybreaks)

# randomly generating movement probs for each class and patch
dmat <- matrix(0, ncol= 4, nrow= 2) # row 1 = pstay, row 2 = move
colnames(dmat) <-  c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
rownames(dmat) <- c(1:2)
dmat1 <- dmat            # stay and move for patch 1 - assuming better than patch 2
dmat2 <- dmat

# generating probs by sampling truncated norm. Assume patch 1 better than 2, so move low for 1, high for 2
move1 <- rtnorm(4, mean = 0.2, sd= 0.2, lower = 0, upper = 1) # vector of movement probs -  shouldnt fems disperse less?
dmat1[2,] <- move1
dmat1[1,] <- 1-move1
print(dmat1)
move2<- rtnorm(4, mean = 0.6, sd= 0.2, lower = 0, upper = 1) 
dmat2[2,] <- move2
dmat2[1,] <- 1-move2
print(dmat2)


# logit link - for sex ratio of current and movement prob, (not distance!)
# example sr = just over 1 as relative fem abundance (Nf/Nm)
sr <- rnorm(n=50, mean = 1.5, sd= 0.5)  # can get up to 3:1!
if(any(sr <0) ){
  sr[sr < 0] <- 0
}

# visualising ratio
sr_df <- as.data.frame(sr)
obs_sr = sum(!is.na(sr_df$sr))

logit.plot <- ggplot(sr_df, aes(sr))  + 
  geom_histogram(aes(y = ..density..),  binwidth = 0.05) + 
  stat_function(fun = dnorm, args = list(mean = mean(sr_df$sr), sd = sd(sr_df$sr)))

# rescaling y
ybreaks = seq(0,50,5) 
## On primary axis
logit.plot + scale_y_continuous("Counts", breaks = round(ybreaks / (0.05 * obs_sr),3), labels = ybreaks)

# defining binary movement based off this
y <- ifelse(sr_df > 1, 0, 1)    # males only move when fem abundance greater than male
sr_df$Y <- y 


# visualising distance against Y (move or not)
ggplot(data = sr_df, aes(x= sr, y=Y))+
  geom_point() +
  labs(title = "Male movement probability distribution",
       y = "Movement prob", 
       x = "Sex ratio (number of females per male)")

mylogit <- glm(Y ~ sr, data = sr_df, family = "binomial")  # creating gen linear model to predict effect of sex ratio on movement
summary(mylogit)

new_df <- data.frame(sample_ratio = seq(min(sr_df$sr), max(sr_df$sr), length.out=100))   # creating df with possible sex ratios (between min and max generated)
new_df$prob <- predict(mylogit, data = new_df, type = "response")   # applying my model to these sample ratios

# visualising model
ggplot(sr_df, aes(x = sr, y = Y)) +
  geom_jitter(height = 0.03, width = 0, alpha = 0.6) +
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE, color = "blue", size = 1.2) +
  labs(title = "Male movement probability distribution",
       y = "Movement prob", 
       x = "Sex ratio (number of females per male)")   # represents a single male sample - what does he do at diff sex ratios?


# increasing replicates 
# treat each y as movement data for an individ MALE badger
y1 <- ifelse(sr > 1.2, 0, 1)   # setting various different movement thresholds, represents larger sample size
y2 <- ifelse(sr > 1, 0, 1)
y3 <- ifelse(sr > 1.25, 0, 1)
y4 <- ifelse(sr > 0.9, 0, 1)
y5 <- ifelse(sr > 0.7, 0, 1)

y_list <- vector("list")  # blank vecor for p move FOR EACH INDIVIDUAL - 50 responses (length sr) for each individual

for (i in 1:10) {
min_f <- rnorm(1, mean = 1.5, sd = 0.2)  # random number to set as threshold

y_list[[i]] <- ifelse(sr > min_f, 0, 1)    # for each invidual list, stay if greater than this random number, move if not
}

# blank df with 10 cols
move_df <- data.frame(matrix(0, nrow = 50, ncol = 10))
colnames(move_df) <- c("Male_1","Male_2","Male_3", "Male_4", "Male_5", "Male_6", "Male_7", "Male_8", "Male_9", "Male_10")

for ( i in 1:10){
  move_df[,i] <- y_list[[i]]  # first col = first list from y_list
}

move_df$Y <- rowMeans(move_df) # mean response for each sr (from sr_df), creates continuous var from binary

move_df$sr <- sr # sr must add in later so we don't interfer with rowmeans

# new logit function based of these 
newlogit <- glm(Y ~ sr, data = move_df, family = "binomial")  # creating gen linear model to predict effect of sex ratio on movement
summary(newlogit)

new <- data.frame(sr = seq(min(move_df$sr), max(move_df$sr), length.out=100))
new$prob <- predict(newlogit, data = new, type = "response") # prediction of movement given sex ratio based on new glm  
                                                 # - doesn't seem to correlate movement with sr?

# visualising model
ggplot(move_df, aes(x = sr, y = Y)) +
  geom_jitter(height = 0.03, width = 0, alpha = 0.6) +
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE, color = "blue", size = 1.2) +
  labs(title = "Male movement probability distribution",
       y = "Movement prob", 
       x = "Sex ratio (number of females per male)")    # more samples from individuals = curvier


# need to extract equation of curve 
# then, can input Nf. Nm and auto calculate p(move) for all males



# Defining movement function for Dmat
# base logistic equation = L / (1 - e^(-k(x-x0)))
# where L = 1 (max value, prob bound between 0 and 1)
# K = steepness of curve (vis using desmos). We want pop size as some density val?
# x0 = 0.5 (midpoint) 




# how to set up more patches 
set.seed(123)  # setting our number 
reps <- 10
sizes <- runif(reps, min= 3, max=30) 
# multiply each of these pop sizes with stage dist

init_mat <- matrix(0, ncol = 4, nrow = reps) # fill each row with an inital vec per patch 
for(i in 1:reps){
  init_mat[i, ] <- floor(sizes[i]* stagedist)
}

initial <- c(t(init_mat))



# FUTURE - adding density dependent dispersal
ddDmat.create <- function(colnames, nPatches, 
                          group_size, sex_ratio = NULL, 
                          k) {
  Dmat <- matrix(0, ncol= length(colnames), nrow = nPatches) # row = number patches
  colnames(Dmat) <- colnames
  
  max_group_size <- 28 # based on papers recording 26 and 27 
  
  # logistic equation - total group size
  for (p in 1:nPatches){
    x <- group_size/ max_group_size    # proportion of max size
    move <- 1 / (1 - exp((-k(x-0.5)))) 
    Dmat[p,] <- move   # equal for all indivs or variable between sexes? vary K value for males and females?
  }
  
  if(!is.null(sex_ratio) == FALSE){
    for (p in 1:nPatches){
      x <- sex_ratio    # proportion of max size
      move <- 1 / (1 - exp((-k[2](x-0.5)))) # second value for k needed
      Dmat[p,2] <- move   # male movement prob 
    }
  }
  return(Dmat)
}

# NOTES a more versitile way to construct in case we are deciding if yearlings are also moving is to have 
# 4X4 Dmat with 0 in non-moving stages. Multiply relevant entries and if not moving, will be multiplying 0. 





# creating ndist plots for removal visualisation -----
dis <- rnorm(1000, mean = 0.5, sd= 0.05) # sample of 1000 between 0 and 1
df <- as.data.frame(dis)
bw <- 0.005
n_obs <-  sum(!is.na(df$dis))

base.plot <- ggplot(df, aes(dis)) + 
  geom_histogram(aes(y = after_stat(density)),  
                 binwidth = bw,
                 fill = "grey50", alpha= 0.3) + 
  stat_function(fun = dnorm, args = list(mean = mean(df$dis), sd = sd(df$dis)), 
                size = 1.5, colour = "black")+ 
  geom_vline(aes(xintercept = 0.5),                       # Adding a line to show removal year
             colour = "red", linetype = "dashed", size= 1, alpha = 0.8) +
  labs(x = "Removal probability") +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(color = "black")
  )

# rescaling y
ybreaks = seq(0,50,5) 
# On primary axis
base.plot + scale_y_continuous("Counts", breaks = round(ybreaks / (bw * n_obs),3), labels = ybreaks)


# biased trials = multiple plots needed with varying means
df$unbias <- rnorm(1000, mean = 0.45, sd= 0.05)
df$bias <- rnorm(1000, mean = 0.55, sd= 0.05)

# graph limits
means <- sapply(df[c("dis", "unbias", "bias")], mean, na.rm = TRUE)
sds   <- sapply(df[c("dis", "unbias", "bias")], sd,   na.rm = TRUE)
min_x <- min(means - 4 * sds)
max_x <- max(means + 4 * sds)

# want to add vis of these extra means on top of base plot
bias.plot <- ggplot(df, aes(dis)) + 
  geom_histogram(aes(y = after_stat(density)),  
                 binwidth = bw,
                 fill = "grey50", alpha= 0.15) + 
  stat_function(fun = dnorm, args = list(mean = mean(df$dis), sd = sd(df$dis)), 
                size = 1.2, colour = "black", alpha = 0.5) + 
  geom_vline(aes(xintercept = mean(df$dis)),                       # Adding a line to show removal year
             colour = "black", linetype = "dashed", size= 1, alpha = 0.5) +
  
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
  coord_cartesian(xlim = c(min_x, max_x))   # makes sure full range of x axis included
# almost works, why are my curves not complete?