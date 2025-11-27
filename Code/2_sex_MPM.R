# 2 Sex Matrix population Models
# 22.10.25

# NOTE- incorrect matrix used, should be 4X4 final not 6X6.


# creating seperate male and female leslie matrix using survival table data----
Fmat<- matrix(0, nrow=2, ncol=2)# blank female (F) matrix of zeros to fill
Fmat[1,2]<- 0.27   # female reproduction- from matrix so no updates needed
Fmat[2,1]<- 0.65    # cub to yearling 
Fmat[2,2]<-  0.86  # yearling to adult, adult survival as(n-1/n)

Mmat<-  matrix(0, nrow=2, ncol=2)   # blank male matrix 
Mmat[1,2]<- 0.32  #  reproduction
Mmat[2,1]<- 0.67   # cub to yearling 
Mmat[2,2]<-  0.82  # yearling to adult, male adult survival lower than females

# combining into 2 sex matrix for stage specific vital rates

mat2<- matrix(0, nrow=4, ncol=4)
mat2[1:2, 1:2]<- Fmat   # adding females
mat2[3:4, 3:4]<- Mmat    # adding males
names<- c("Yearling(f)", "Adult(f)", "Yearling(m)", "Adult(m)")
colnames(mat2)<- names
rownames(mat2)<- names
mat2   # female produce only female cubs, and males produce only male cubs




# Projections----
initial<- matrix(c(10, 5, 10, 4), ncol=1)  # initial pop structure matrix (no cubs)
N1<- mat2 %*% initial  # 1 year projection. %*% for matrix multiplication

years<- 10  # over 10 years
project1<- matrix(0, nrow=nrow(mat2), ncol=years+1)  # blank matrix to fill without our function
rownames(project1)<- rownames(mat2)   # rows show sex and stage
colnames(project1)<- c(0:10)          # year of projection (initial =0)
project1[,1]<- initial    # first row of blank matrix is initial vector

# loop calculating pop structure - (later as function mat.proj)
for(t in 2:(years+1)){      # t representing row, each row a year. first col is initial so starts from 2
  project1[,t]<- mat2 %*% project1[,t-1]    # fill rows of blank projection matrix by multiplying prev col by vital rate matrix
}
project1  


# plotting this 10 year projection and colour coding life stages and sex
cols<- c( "red", "red", "blue", "blue")   # custom colour vector
matplot(t(project1),   # transposing matrix- swaps cols for rows
        type="l",  # type l for line
        col=cols, lty= 1:2,           # red female, blue male, line type
        ylab= "abundance", xlab="time (t)")
legend("topright", legend= c("Yearling", "Adult"),
       cex=0.65,
       lty= 1:2)   # line type matches graph


library(popdemo)
eigs(mat2)
# lambda 1.030333
# $ss
# 0.2076391 0.7923609 0.0000000 0.0000000
# $rv
# 0.6832316 1.0830095 0.0000000 0.0000000


#manually calculating pop growth rate
popN<- apply(project1, 2, sum)   # apply sum func over matrix array (2= by col). Pop size per year
plot(1:11, popN, 
    ylab= "Population size", xlab="year")  # population size falls over the years (lambda<1)

lambda<- popN[2:11]/popN[1:10]    # each year divided by prev year abundance
lambda  # calculates growth rate for each year
#        1         2         3         4         5         6         7         8         9        10 
# 0.8072414 1.0814310 1.0202555 1.0311666 1.0291067 1.0294971 1.0294235 1.0294387 1.0294368 1.0294381 

# stage distribution
propStage<- project1/popN # divide each value by total pop size to get proportion each year
matplot(1:11, propStage, type= "l")

# elasticity and sensitivity analysis
popdemo:: elas(mat2)
#             Cub(f) Yearling(f) Adult(f)    Cub(m) Yearling(m)  Adult(m)
# Cub(f)           0           0        0 0.0000000   0.0000000 0.0000000
# Yearling(f)      0           0        0 0.0000000   0.0000000 0.0000000
# Adult(f)         0           0        0 0.0000000   0.0000000 0.0000000
# Cub(m)           0           0        0 0.0000000   0.0000000 0.1108541
# Yearling(m)      0           0        0 0.1108541   0.0000000 0.0000000
# Adult(m)         0           0        0 0.0000000   0.1108541 0.6674378

# How to do sensitivity analysis?