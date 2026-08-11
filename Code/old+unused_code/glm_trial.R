# testing conditions and movement link
library(ggplot2)

# for now, movement prob sampled from random dist 0-1
# p(move) for each stage or patch? stage and patch?
prob <- rnorm(1000, mean = 0.5, sd= 0.3) # sample of 1000 
prob_df <- as.data.frame(prob)

n_obs <-  sum(!is.na(prob_df$prob))

norm.plot <- ggplot(prob_df, aes(prob))  + 
  geom_histogram(aes(y = after_stat(density)),  binwidth = 0.05) + 
  stat_function(fun = dnorm, args = list(mean = mean(prob_df$prob), sd = sd(prob_df$prob)))

# rescaling y
ybreaks = seq(0,50,5) 
# On primary axis
norm.plot + scale_y_continuous("Counts", breaks = round(ybreaks / (0.05 * n_obs),3), labels = ybreaks)

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



# FUTURE - density dependent dispersal (p leave affected by crowdedness)
ddDmat.create <- function(colnames, npatches, 
                          group_size, sex_ratio = NULL, 
                          k) {  # K value strength of 
  Dmat <- matrix(0, ncol= length(colnames), nrow = npatches) # row = number patches
  colnames(Dmat) <- colnames
  
  max_group_size <- 28 # based on papers recording 26 and 27 
  
  # total group size
  for (p in 1:npatches){
    x <- group_size/ max_group_size    # proportion each group of max size
    move <- 1 / (1 - exp((-k*(x-0.5))))   # logistic equation 
    move[move<0] <- 0   # setting any negatives to 0
    Dmat[p,] <- move   # equal for all indivs or variable between sexes? vary K value for males and females?
  }
  
  if(!is.null(sex_ratio) == FALSE){
    for (p in 1:npatches){
      x <- sex_ratio    # proportion of max size
      move <- 1 / (1 - exp((-k[2]*(x-0.5)))) # second value for k needed
      Dmat[p,2] <- move   # male movement prob 
    }
  }
  return(Dmat)
}

# NOTES a more versitile way to construct in case we are deciding if yearlings are also moving is to have 
# 4X4 Dmat with 0 in non-moving stages. Multiply relevant entries and if not moving, will be multiplying 0. 



ddDmat.create(colnames = c("Adult_f", "Adult_m"), npatches = 3, 
               group_size = c(21, 3, 11), sex_ratio = NULL, 
               k = 0.1)

