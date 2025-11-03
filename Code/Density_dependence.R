# Density dependence and mating systems
# 29/10/2025

# creating a 2 sex matrix with create.matrix
stageNames<- c("Cub", "Yearling", "Adult") # vector of stages

mat1<-  create.matrix(names= stageNames,   # vector of stage names    
                      nSex=2,    # 1 or 2 sex matrix (default 1)             
                      r1= 0.27 ,             # reproduction sex 1
                      r2= 0.32,              # reproduction sex 2
                      g1= c(0.65, 0.65),   # vector of all growth rates/ transition probabilities
                      g2= c(0.67, 0.67), # second sex values for growth
                      S= c(0.86, 0.82)
                      )





#  vital rates must vary across years as a function of pop size 
#  Rates affected by desnity dependence= reproduction and recruitment. Adult and yearling survival stable

