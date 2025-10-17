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
class(A)

# creating stageframe, dataframe needed for MPM analysis
sframe<- sf_create( 
  sizes=c(1,2,3),   # place holder for required size col
  stagenames= , # lifestages for Leskovitch matrix
  repstatus=c(0, 0, 1),   # only adults reproductive
  matstatus=c(0, 0, 1),   # adults mature 
  immstatus=c(1, 0, 0 ), # cubs immature
  propstatus = c(0, 0, 0),
)


MPM<- create_lM(mats=list(A), # must be converted to list even if length=1
                stageframe = sframe) 
summary(MPM)
