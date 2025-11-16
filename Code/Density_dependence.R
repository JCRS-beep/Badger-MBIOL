# Density dependent recruitment
# 12/11/25

# creating blank base matrix
stages<- c("Cub_f", "Yearling_f", "Adult_f", "Cub_m", "Yearling_m", "Adult_m")
mat.blank<- matrix(0, ncol=6, nrow=6)
rownames(mat.blank)<- stages
colnames(mat.blank)<- stages

# creating Umat, martix with only survival and transitions
ex_Umat<- mat.blank
ex_Umat[2,1]<- "e^-bN*G1f"    # only cub survival density dependent
ex_Umat[3,2]<- "G2f"
ex_Umat[3,3]<- "Sf"
ex_Umat[5,4]<- "e^-bN*G1m"   
ex_Umat[6,5]<- "G2m"
ex_Umat[6,6]<- "Sm"



Umat<-mat.blank  # creating my working matrix. Values from 
Umat[2,1]<- 0.32   # only cub survival density dependent
Umat[3,2]<- 0.67  
Umat[3,3]<- 0.86
Umat[5,4]<- 0.27 
Umat[6,5]<- 0.65
Umat[6,6]<- 0.77

# loading Parameters for density and Fmat creation 
params<- data.frame(      # dataframe 
  fmax= 3.2,     # F fecundity max (max cubs per adult female)
  Sc_f_max=0.65,   # max female  survival
  Sc_m_max=0.65, 
  b=0.002,       # temp value- must be calculated from provided datasets
  rep_K= 3,          #litter size (K)
  h= 6   # harem size per male
)

# creating Fmat basic, with max reproductive output
Fmat<- mat.blank # Fmat will be blank, each year of projection reproductive rates are calculated
f<- params$fmax
Fmat[1,3]<- 0.5*f
Fmat[4,3]<- 0.5*f

Fmat2<- Mating.Fmat(params, stages, Nf=20, Nm=20)  # expected based on N=40
Fmat2

# Ricker function----
ricker<- function(params, N){   # Using params (b) and population size input
  dd_fun<-exp(-params$b*N)
  return(dd_fun)      # returns value of multiplier
}

# Applying to matrix with function
Amat<-apply.DD(params, Fmat, Umat, N=40, DDapply="recruitment") # N includes only adults and yearlings, same as Nf and Nm
Amat
# Issues with recruitment application- Amat full of NA?




