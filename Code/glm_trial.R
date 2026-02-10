# testing distance and movement link
install.packages("msm")
library(msm)
library(ggplot2)

# for now, movement prob sampled from random dist 0-1
# p(move) for each stage or patch? stage and patch?
prob <- rtnorm(1000, mean = 0.5, sd= 0.3, lower = 0, upper = 1) # sample of 1000 between 0 and 1
prob_df <- as.data.frame(prob)
bw <- 0.05
n_obs <-  sum(!is.na(prob_df$prob))

norm.plot <- ggplot(prob_df, aes(prob))  + 
  geom_histogram(aes(y = after_stat(density)),  binwidth = bw) + 
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
# example sr = should be around 0.5 as proportion fem Nf/(Nf+Nm)
sr <- rtnorm(n=1000, mean = 1.2, sd=0.4, lower = 0)  # can get up to 3:1!

# visualising distance points
sr_df <- as.data.frame(sr)
bw_sr = 0.05
obs_sr = sum(!is.na(sr_df$sr))

logit.plot <- ggplot(sr_df, aes(sr))  + 
  geom_histogram(aes(y = ..density..),  binwidth = bw_sr) + 
  stat_function(fun = dnorm, args = list(mean = mean(sr_df$sr), sd = sd(sr_df$sr)))

# rescaling y
ybreaks = seq(0,50,5) 
## On primary axis
logit.plot + scale_y_continuous("Counts", breaks = round(ybreaks / (bw_sr * obs_sr),3), labels = ybreaks)

## Or on secondary axis
#logit.plot + scale_y_continuous("Density", sec.axis = sec_axis(
#  transform = ~ . * bw_sr * obs_sr, name = "Counts", breaks = ybreaks))



# taking smaller sample of normal (100 points)
sample_ratio <- rtnorm(n=100, mean = 1.2, sd=0.4, lower = 0)
df <- as.data.frame(sample_ratio)

# treat each y as movement data for an individ MALE badger
y1 <- ifelse(sample_ratio > 1.2, 0, 1)
y2 <- ifelse(sample_ratio > 1, 0, 1)
y3 <- ifelse(sample_ratio > 1.25, 0, 1)
y4 <- ifelse(sample_ratio > 0.9, 0, 1)
y5 <- ifelse(sample_ratio > 0.7, 0, 1)

move_data <-as.data.frame(y1) # 5 rows, each 
move_data$y2 <- y2
move_data$y3 <- y3
move_data$y4 <- y4
move_data$y5 <- y5
Y <- rowMeans(move_data)

df$Y <- Y # av per row across y vecs


# visualising distance against Y (move or not)
ggplot(data = df, aes(x= sample_ratio, y=Y))+
  geom_point() 
  
  
# creating logit model
mylogit <- glm(Y ~ sample_ratio, data = df, family = "binomial")
summary(mylogit)

new <- data.frame(sample_ratio=seq(min(df$sample_ratio), max(df$sample_ratio), length.out=100))
new$prob <- predict(mylogit, newdata= new, type = "response")

# visualising model
ggplot(df, aes(x = sample_ratio, y = Y)) +
  geom_jitter(height = 0.03, width = 0, alpha = 0.6) +
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE, color = "blue", size = 1.2) +
  labs(x = "Proportion fem", y = "Probability of Y = 1 (males)")   # more samples from individuals = curvier




# test syntax ----
n <- list(c(10,11,12,13), c(1,2,3,4))
length(n)
# how to split into objects? 
n[[1]][2]
lapply(n, sum) 

Reduce(`+`, n)  # sums entries within lists
sum(n[[1]])  # must have double brackets to sum

n1 <- c(10, 10, 13, 12)
length(n1)

mat.test <- list(matrix(0, ncol = 4, nrow = 11))
Vec.test <- rep(mat.test, 2)  # why no rep this time?
size.test <- list(rep(NA, (11)))       # vector to fill with total pop size each year
Pop.test <- rep(size.test, 2)

for (m in 1:2){
  Vec.test[[m]][1,] <- n[[m]]          # for each list in Vec, first row is initial         
  Pop.test[[m]][1] <- lapply(n, sum)   
  
}

Vec.test
Pop.test[1:2] <- lapply(n0, sum)  
