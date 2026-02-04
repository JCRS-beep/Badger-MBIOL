# Density dependent recruitment
# 12/11/25

# creating blank base matrix
stages<- c("Yearling_f", "Adult_f", "Yearling_m", "Adult_m")
mat.blank<- matrix(0, ncol=4, nrow=4)
rownames(mat.blank)<- stages
colnames(mat.blank)<- stages

# creating Umat, martix with only survival and transitions
ex_Umat<- mat.blank
ex_Umat[2,1]<- "G2f"
ex_Umat[2,2]<- "Sf"
ex_Umat[4,3]<- "G2m"
ex_Umat[4,4]<- "Sm"

Umat<-mat.blank  # creating my working matrix. Values from 
Umat[2,1]<- 0.67   # yearling f survival
Umat[2,2]<- 0.78   # adult f survival
Umat[4,3]<- 0.65   # yearling m survival
Umat[4,4]<- 0.72

# loading Parameters for density and Fmat creation 
params<- data.frame(      # dataframe 
  fmax= 1.2,     # F fecundity max (max cubs per adult female) 
  Sc_f_max=0.65,   # max cub survival (equal for sexes)
  Sc_m_max=0.65,
  b=0.004,       # temp value- must be calculated from provided datasets
  rep_K= 3.5,          #litter size (K)
  h= 6   # harem size per male
)


# Fmat example for density dependence
ex_Fmat<- mat.blank
ex_Umat[1,2]<- "0.5*e^-2bN*f"    # female reproduction AND cub survival density dependent. Must apply dependence twice
ex_Umat[3,2]<- "0.5*e^-2bN*f"  


# creating Fmat basic, with max reproductive output
Fmat<- mat.blank # Fmat will be blank, each year of projection reproductive rates are calculated
f<- params$fmax
Sf<- params$Sc_f_max
Sm<- params$Sc_m_max
Fmat[1,2]<- 0.5*f*Sf   # this Fmat has max values of reproduction (freq and density independent)
Fmat[3,2]<- 0.5*f * Sm

# Mating function -----
# can be loaded from script in Function folder

# Creating Fmat for max rates WITH MATING FUNCTIONS
Fmat2 <- mating.func(params, stages, Nf= 20, Nm= 20, Mfunction= "min", return.mat=TRUE)



# Applying to matrix with function applyDD for a static population - 40 adults, 20 females, 20 males -----
Amat<-apply.DD(params, Fmat2$Fmat, Umat, N=40, DDapply="Fmat") # N includes only adults and yearlings, so that N= Nf + Nm


# testing projection function----
initial <- c(10, 10, 10, 10)
test_proj <-dd.proj(Umat, 
                    initial, 
                    params, 
                    stagenames = stages, 
                    time = 20, 
                    memberN=NULL,  # which individuals contribute to pop size? (as vec)
                    DDapply= "Fmat", 
                    Mfunction= "min",
                    return.vec= TRUE) 
# SUCCESS 

# Plotting pop strucutre over time
library(ggplot2)
time_vec <- c(0:20)
 
(Nplot <- ggplot(NULL, aes(x=time_vec, y=test_proj$pop)) + 
                           xlab("time(years)") +
                           ylab("Population size") +
  geom_point(alpha=0.8)+
    geom_line() +
    theme_classic())      # nice plateau, but 300 individuals appears high? Other studied have sizes above 300, could be reasonable



# dataframe creation for plot function use

proj_df <- as.data.frame(test_proj$vec)
proj_df$Year <- c(0:20)   

# df must be in longer form - year vs stage TIDY DATA!
library(tidyr)
#  gather to create new columns stage and abundance by year?
proj_df_long <- gather(proj_df, key= "Stage", value= "Abundance", 1:4)  # key is new column you are gathering by, value is the values you want in new cols
proj_df_long2 <- separate(proj_df_long, col= "Stage", into=c("Stage", "Sex"), sep='_')    # separate stage into new col, stage and sex by _. 


# lty_v <- c(rep(1,2), rep(2,2))   # line type = stage
cols<- c("#FF6A6A", "#87CEEB", "#5CACEE") 
col_v1 <- rep(c("#FF6A6A", "#87CEEB"), 2)
col_v2 <- rep(c("#FF6A6A", "#5CACEE"), 2)

theme_custom <- function(){
  theme_classic()+
    theme(
      axis.text.x= element_text(size=10), 
      axis.text.y= element_text(size= 10), 
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, vjust = 1, hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 10, face = "italic"),
      legend.position.inside = c(0.9, 0.7))
}


(abundance_plot <- ggplot(data= proj_df_long2, aes(x=Year, y=Abundance, colour= Sex, linetype=Stage, shape=Stage))+
    geom_point(position= "jitter", alpha=0.8)+  # jitter to avoid overlap pf yearlings
    geom_line(data= proj_df_long2, alpha=0.7) +
    scale_colour_manual(values=col_v2,
                        labels=c("Female", "Male")) +
    labs(title = "Stage Abundance over Time", 
         x = "Year", y = "Abundance") + 
    theme_custom()   # added correctly here?
    )


(test_plot <- dd_plot(test_proj, 
                                  y_val= "Stages", 
                                  ylab = "Abundance", 
                                  xlab = "Time (t)",
                                  theme = theme_custom(),   
                                  cols= c("#FF6A6A", "#87CEEB"),
                                  legend.pos = "topright",
                                  cex.legend = 0.8))   
  

# using growth.rate() function to plot lambda over time
lambda <- growth.rate(test_proj)
gr_graph <- plot(lambda, xlab= "Year")  # steady decline over time - is this the pattern we want? 
SSD <- prop.stage(test_proj)

SSD_df <- as.data.frame(SSD)
colnames(SSD_df) <- c(stages)
SSD_df$Year <- c(0:20)


SSD_long <- gather(SSD_df, key= "Stage", value= "Proportion", 1:4)  # key is new column you are gathering by, value is the values you want in new cols
SSD_long <- separate(SSD_long, col= "Stage", into=c("Stage", "Sex"), sep='_')    # separate stage into new col, stage and sex by _

SSD_plot <- ggplot(data= SSD_long, aes(x= Year, y=Proportion, colour= Sex, shape = Stage))+
  geom_point()
