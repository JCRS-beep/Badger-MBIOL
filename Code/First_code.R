# Beginning of script
# 13/10/25

# install.packages("popdemo") (do once)
library(popdemo)

# loading my data from wytham badger pops (Macdonald and Newman, 2018)
badger_raw<- read.csv("Data/wytham.data.csv")   # more info on this df in "Docs" folder. Survival rates of badgers by ages

# Creating my matrix based off this (info in doc )
badger_matrix<- matrix(c(0,    0,	   3, 
                         0.12, 0,    0, 
                         0,    0.66, 0.32), 
                       byrow=TRUE, 3,3)



# creating a population projection
initial<- c(20, 20, 20) # vector of starting pops for each stage
project<- popdemo::project(badger_matrix, vector = initial, time = 10) 


project.plot<- vec(project)
matplot(project.plot, type="l", log="y")
legend("topright", legend= c("Cub", "Yearling", "Adult"), 
       col= 1:ncol(project.plot), lty= 1:ncol(project.plot))
# projection shows decline to extinction over 10 time intervals

# Asking for eigen values and vectors
eigs(badger_matrix)
#lamdba= 0.7464394 stable stage distribution= 0.7094330 0.1140507 0.1765163
#reproductive rates = 0.375847 2.337892 2.644082 
############
## Comment from Chrissy: The left eigenvector is reproductive *value*, not
#reproductive rate. There's a subtle difference - check that you understand (and
#we can discuss in a future meeting as needed).
############

# Elasticity analysis
popdemo::elas(badger_matrix)
#           [,1]      [,2]      [,3]
# [1,]         0         0 0.2666382
# [2,] 0.2666382         0         0
# [3,]         0 0.2666382 0.2000853
# suggests 3 transitions are equally important: birth, and cub and yearling development
############
## Comment from Chrissy: Equally important to/for what? What does elasticity tell us?
############


# Life table response experiments- USEFUL IN FINAL REMOVAL EXPERIMENTS
# How to model removal
removal<- matrix(c(0,      0,	     3*1.5,     # if reproduction increases by 50%
                   0.12*2, 0,        0,       # if cub survival doubles
                   0,      0.66*0.5, 0.32*0.5),     # if adult and yearling survival falls 50%
                byrow=TRUE, 3,3)

popdemo::eigs(removal, what = "lambda")
# 0.7665447
elas(removal)
#         [,1]      [,2]       [,3]
# [1,] 0         0         0.30639233
# [2,] 0.3063923 0 ,       0
# [3,] 0         0.3063923 0.08082302

############ 
# Comment from Chrissy: I like the direction that this is going! It's not quite
# a Life Table Response Experiment (LTRE) yet as written here. An LTRE would
# take two (or more) matrices and compare them, decomposing the difference or
# variance in population growth rate (or another metric of interest) into
# contributions from the different matrix elements (or underlying vital rates).
# It's a pretty flexible tool that I am an expert in, so we can discuss more in
# a future meeting!
############


