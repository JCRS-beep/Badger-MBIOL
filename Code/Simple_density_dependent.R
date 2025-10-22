# Adding density dependence with lefko
# 17/10/25


# Installing lefko package for MPM analysis
install.packages("lefko3", dependencies = TRUE)
library(lefko3)


# recreating my simple MPM with lefko
# creating the matrix
A<- matrix(c(0,    0,	   3, # pasting prev matrix
             0.66, 0,    0, 
             0,    0.66, 0.29), 
          nrow=3, byrow=TRUE)
stages<- c("Cub", "Yealing", "Adult")
dimnames(A)<- list(stages, stages)

# creating stageframe, dataframe needed for MPM analysis
sframe<- sf_create( 
  sizes=c(1,2,3),   # place holder for required size col
  stagenames= stages , # lifestages for Leskovitch matrix
  repstatus=c(0, 0, 1),   # only adults reproductive
  matstatus=c(0, 0, 1),   # adults mature 
  immstatus=c(1, 0, 0 ), # cubs immature
  propstatus = c(0, 0, 0),
)

# Creating our MPM
MPM<- create_lM(mats=list(A),   # must be converted to list even if length=1
                stageframe = sframe) 
summary(MPM)

# recreating our 10 year projection
project<- projection3(MPM, times=10)  # projecting MPM over 10 timesteps
summary(project) # assumes starting vector of 1 per stage- very small pop!
plot(project)

project2<- projection3(MPM, times=10, integeronly = TRUE)  # allowing only integer values
summary(project2) 
plot(project2)  # quickly falls to extinction- only use if interested in small pops, abolsute population size and impacts of demographic stochasticity


# creating start vector 
project3<- projection3(MPM, times=10, 
                      start_vec=c(20, 20, 20), 
                      ) # used in simple first analysis
summary(project3) 
plot(project3)

project4<- projection3(MPM, times=10, 
                       start_vec=c(20, 20, 20), 
                       integeronly = TRUE)   # this time, no extinction!
summary(project4) 
plot(project4)

# creating a starting pop more similar to expected 
initial<- c(5, 10, 0)
project_final<- projection3(MPM, times=10, 
                       start_vec=initial,   # using prev defined vec
                       integeronly = TRUE) 
summary(project_final) 
plot(project_final)