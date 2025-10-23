# 2 Sex Matrix population Models
# 22.10.25

# creating seperate male and female leslie matrix using survival table data----
Fmat<- matrix(0, nrow=3, ncol=3)# blank female (F) matrix of zeros to fill
Fmat[1,3]<- 0.27   # female reproduction
Fmat[2,1]<- 0.65    # cub to yearling 
Fmat[3,]<- c(0, 0.65, 0.33)  # yearling to adult, adult survival

Mmat<-  matrix(0, nrow=3, ncol=3)   # blank male matrix 
Mmat[1,3]<- 0.32  #  reproduction
Mmat[2,1]<- 0.67   # cub to yearling 
Mmat[3,]<- c(0, 0.67, 0.25)  # yearling to adult, adult survival

# combining into 2 sex matrix for stage specific vital rates
mat2<- matrix(0, nrow=6, ncol=6)
mat2[1:3, 1:3]<- Fmat   # adding females
mat2[4:6, 4:6]<- Mmat    # adding males
cnames<- c("Cub(f)", "Yearling(f)", "Adult(f)", "Cub(m)", "Yearling(m)", "Adult(m)")
rnames<- c("Cub(f)", "Yearling(f)", "Adult(f)", "Cub(m)", "Yearling(m)", "Adult(m)")
colnames(mat2)<- cnames
rownames(mat2)<- rnames
mat2   # female produce only female cubs, and males produce only male cubs




# Projections----
initial<- matrix(c(0,20, 5, 0, 20, 4), ncol=1)  # initial pop structure matrix (no cubs)
N1<- mat2 %*% initial  # 1 year projection. %*% for matrix multiplication

years<- 10  # over 10 years
project1<- matrix(0, nrow=nrow(mat2), ncol=years+1)  # blank matrix to fill without our function
rownames(project1)<- rownames(mat2)   # rows show sex and stage
colnames(project1)<- c(0:10)          # year of projection (initial =0)
project1[,1]<- initial    # first row of blank matrix is initial vector

# loop interactive calculating pop structure - create function
for(t in 2:(years+1)){      # t representing row, each row a year, first col is initial
  project1[,t]<- mat2 %*% project1[,t-1]    # fill rows of blank projection matrix by multiplying prev col by vital rate matrix
}
project1

# plotting this 10 year projection and colour coding life stages and sex
cols<- c("red", "red", "red", "blue", "blue", "blue")   # custom colour vector
matplot(project1, type="l",  # type l for line
        col=cols, lty= 1:3,           # red female, blue male, line type
        ylab= "abundance", xlab="time (t)")
legend("topright", legend= c("Cub", "Yearling", "Adult"),
       cex=0.65,
       lty= 1:3)   # line type matches graph


library(popdemo)
eigs(mat2)
# lambda= 0.6234685 
# ss 0.2297954 0.2395743 0.5306303 0.0000000 0.0000000 0.0000000
# rv 1.0550894 1.0120231 0.9707146 0.0000000 0.0000000 0.000000

#manually calculating pop growth rate
popN<- apply(project1, 2, sum)   # apply sum func over matrix array (2= by col). Pop size per year
plot(1:11, popN, 
    ylab= "Population size", xlab="year")  # population size falls over the years (lambda<1)

lambda<- popN[2:11]/popN[1:10]    # each year divided by prev year abundance
lambda  # calculates growth rate for each year


# stage distribution - incorrect?
propStage<- project1/popN    # divide each value by total pop size to get proportion each year
matplot(1:11, propStage, type= "l")

# function to calculate stage distribution per year
propStage<- function(projected, popN)     # inputs= projected matrix, calculated pop size for each year
    {
   for(i in 1:ncol(projected))  # repeat for each column
   
    }
