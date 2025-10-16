# Beginning of script
# 13/10/25

# install.packages("popdemo") (do once)
library(popdemo)

# loading my data from wytham badger pops (Macdonald and Newman, 2018)
badger_raw<- read.csv("Data/wytham.data.csv")   # more info on this df in "Docs" folder. Survival rates of badgers by ages

# Creating my matrix based off this (info in doc )
badger_matrix<- matrix(c(0,    0,	   3, 
                         0.66, 0,    0, 
                         0,    0.66, 0.29), 
                       byrow=TRUE, 3,3)



# creating a population projection
initial<- c(20, 20, 20) # vector of starting pops for each stage
project<- popdemo::project(badger_matrix, vector = initial, time = 10) 


project.plot<- vec(project)
matplot(project.plot, type="l", log="y")
legend("topright", legend= c("Cub", "Yearling", "Adult"), 
       col= 1:ncol(project.plot), lty= 1:ncol(project.plot))
# projection shows rise over 10 time intervals

# Asking for eigen values and vectors
eigs(badger_matrix)
# lamdba= 1.199006
# stable stage distribution= 0.5127877 0.2822671 0.2049451
# reproductive value = 0.5875586 1.0674033 1.9391252  

# Elasticity analysis
popdemo::elas(badger_matrix)
#           [,1]      [,2]      [,3]
# [1,]         0         0 0.30129286
# [2,]  0.3012929         0         0
# [3,]         0  0.3012929 0.09612143
# equal values for birth, and cub and yearling development, changes in these vital rates affect lambda 



# Life table response experiments- USEFUL IN FINAL REMOVAL EXPERIMENTS
# How to model removal
removal<- matrix(c(0,      0,	     3*1.5,     # if reproduction increases by 50%
                   0.66*2, 0,        0,       # if cub survival doubles
                   0,      0.66*0.5, 0.29*0.5),     # if adult and yearling survival falls 50%
                byrow=TRUE, 3,3)

popdemo::eigs(removal, what = "lambda")
#1.301755
elas(removal)
#         [,1]      [,2]       [,3]
# [1,] 0         0          0.31996410
# [2,] 0.3199641 0 ,        0
# [3,] 0         0.3199641  0.04010769