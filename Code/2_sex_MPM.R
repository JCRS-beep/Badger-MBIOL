# 2 Sex Matrix population Models
# 22.10.25

# creating seperate male and female leslie matrix using survival table data----
Fmat<- matrix(0, nrow=3, ncol=3)# blank female (F) matrix of zeros to fill
Fmat[1,3]<- 0.27   # female reproduction
Fmat[2,1]<- 0.65    # cub to yearling 
Fmat[3,]<- c(0, 0.65, 0.86)  # yearling to adult, adult survival as(n-1/n)

Mmat<-  matrix(0, nrow=3, ncol=3)   # blank male matrix 
Mmat[1,3]<- 0.32  #  reproduction
Mmat[2,1]<- 0.67   # cub to yearling 
Mmat[3,]<- c(0, 0.67, 0.82)  # yearling to adult, male adult survival lower than females

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

# loop calculating pop structure - (later as function mat.proj)
for(t in 2:(years+1)){      # t representing row, each row a year. first col is initial so starts from 2
  project1[,t]<- mat2 %*% project1[,t-1]    # fill rows of blank projection matrix by multiplying prev col by vital rate matrix
}
project1  


# plotting this 10 year projection and colour coding life stages and sex
cols<- c("red", "red", "red", "blue", "blue", "blue")   # custom colour vector
matplot(t(project1),   # transposing matrix- swaps cols for rows
        type="l",  # type l for line
        col=cols, lty= 1:3,           # red female, blue male, line type
        ylab= "abundance", xlab="time (t)")
legend("topright", legend= c("Cub", "Yearling", "Adult"),
       cex=0.65,
       lty= 1:3)   # line type matches graph


library(popdemo)
eigs(mat2)
# lambda 1.002837
# ss 0.0000000 0.0000000 0.0000000 0.2082479 0.1391314 0.6526206
# rv 0.0000000 0.0000000 0.0000000 0.5323178 0.7967578 1.1925640

#manually calculating pop growth rate
popN<- apply(project1, 2, sum)   # apply sum func over matrix array (2= by col). Pop size per year
plot(1:11, popN, 
    ylab= "Population size", xlab="year")  # population size falls over the years (lambda<1)

lambda<- popN[2:11]/popN[1:10]    # each year divided by prev year abundance
lambda  # calculates growth rate for each year


# stage distribution - incorrect?
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

# How to do sensitivity analysis



# creating a 2 sex matrix with create.matrix----
stageNames<- c("Cub", "Yearling", "Adult") # vector of stages

mat1<-  create.matrix(names= stageNames,   # vector of stage names    
                      nSex=2,    # 1 or 2 sex matrix (default 1)             
                      r1= 0.32 ,             # reproduction sex 1
                      r2= 0.4,              # reproduction sex 2
                      g1= c(0.67, 0.67),   # vector of all growth rates/ transition probabilities
                      g2= c(0.65, 0.65), # second sex values for growth
                      S= c(0.86, 0.82)
)
mat1


# projecting matrix for population size estimates----
n0<- c(10, 10, 10, 10, 10, 10)  # initial pop structure
cols<- c("red", "blue")
proj1<- mat.proj(n0, mat1,20)  # 10 year projection
plot.proj(proj1, cols, nSex=2)   
growth(proj1)  # varies below 1
 