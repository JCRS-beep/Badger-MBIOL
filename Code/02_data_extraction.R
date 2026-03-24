# Parameter estimation for final usage

# installing digitise package
# install.packages("devtools")
# devtools::install_github("daniel1noble/metaDigitise")

# packages ----
library(metaDigitise)
library(dplyr)
library(tidyr)
library(tidyverse) # NEEDED FOR PIPE!
library(ggplot2)
library(gridExtra)


# Functions ----
# function - for scatter, rename, round, remove cols
scat.clean <- function(df,y_name = NULL, year= TRUE, remove.id= TRUE, dec = 0) {
  if(is.null(y_name)){
    yname <- unique(df$y_variable)}
  else{
    yname <- y_name
  }
  names(df)[names(df) == "y"] <- yname
  
  
  df <- dplyr::select(df, -c(col, pch, group,  y_variable, x_variable))
  
  if(year == TRUE){   
    df <- rename(df, Year = x)
    
    if(is.null(dec)){    # put dec = NULL for rounding down
      df$Year <- as.numeric(floor(df$Year)) # 1979.6 should round to 1979 not 1980
      
    }else if(is.numeric(dec)){  
      df$Year <-  as.numeric(round(df$Year), dec)  # rounding the year to specified dec places
    }
  } else if(year== FALSE){
    xname <- unique(df$x_variable)
    names(df)[names(df) == "x"] <- xname
  }
  
  if(remove.id== TRUE){
    df <- dplyr::select(df, -id)
  }
  
  names(df) <- gsub("\\s+", "_", names(df))
  df <- as.data.frame(df)
}

# converting df to wide  with id split by sex and age (use in Bright-Ross dfs)
widen <- function(df, y_val) {
  df_wide <- spread(df, "id", y_val, "Year")   # third col is survival or percentage producing offspring
  # turning second col into numeric vals
  
  cols <- names(df_wide)
  df_wide[cols] <- lapply(df_wide[cols], function(x) as.numeric(as.character(x)))
  
  names(df_wide) <- names(df_wide) |>
    gsub("\\s+", "_", x = _) |>          # replace spaces with _
    gsub("[()]", "", x = _)  |>      # remove ()  - hasn't removed successfully
    gsub("-", "_", x = _)
  
  return(df_wide)    # year, sex+age in cols
}


# manually digitising figs ----
# figs <- metaDigitise(dir ="C:/Users/jaycr/OneDrive - Nexus365/UNI/Masters/Project resources/Badger-MBIOL/Data/FigExtraction",
#                      summary = FALSE)    # summary = false imports raw data
#Version for Chrissy's computer (but I think it might work for Jay as well, if the working directory is Badger-MBIOL/)
figs <- metaDigitise(dir ="Data/FigExtraction", summary = FALSE)

# splitting into plot types is easier to work with
scatter <- figs$scatterplot
mean_err <- figs$mean_error


# organising paper 1  - Rogers 1997 ----
rogers_1997_fig1 <- figs$scatterplot$`242_Rogers_1997_fig1.png`
rogers_1997_popsize <- scat.clean(rogers_1997_fig1, "MNA", remove.id = FALSE, dec= 0)  
# order years correctly
# rogers_1997_popsize <- arrange(rogers_1997_popsize, by= Year, group_by(id)) 

# reorganising so fewer rows
rogers_1997_popsize <- spread(rogers_1997_popsize, id, MNA)
rogers_1997_popsize$Year <- as.numeric(rogers_1997_popsize$Year)
# renaming rows
names(rogers_1997_popsize) <- gsub("\\s+", "_", names(rogers_1997_popsize))

# reproductive rates
rogers_1997_rep3 <- scatter$`242_Rogers_1997_fig13.png` # Cubs / breeding fem
rogers_1997_rep3_clean <- scat.clean(rogers_1997_rep3, "Cubs_per_breeding_fem", dec = 0)
rogers_1997_rep3_clean$Cubs_per_breeding_fem[2] <- NA # removing outlier point


rogers_1997_data <- full_join(rogers_1997_rep3_clean, rogers_1997_popsize, by= "Year")
rogers_1997_data <- rename(rogers_1997_data, 
                           "Adult_density" = "Total_Adult_Badgers", "Adult_fem_density" = "Adult_Female_Badgers", "Adult_male_density" ="Adult_Male_Badgers", "Cub_density" = "Badger_Cubs")

# estimating average litter size
litter.lm <- lm(Cubs_per_breeding_fem~Adult_density, data= rogers_1997_data)  # does density affect litter size
summary(litter.lm)  # Int = 2.206774   
summary(rogers_1997_data$Cubs_per_breeding_fem, na.rm = TRUE) # mean = 2.299102, max = 3.144 ,similar to intercept

# choosing litter size value to use - mean is most representative?
rogers_k <- mean(rogers_1997_data$Cubs_per_breeding_fem, na.rm = TRUE)

# cub survival
cub_mort <- 0.24   # mortality rate from paper
rogers_cub_survival <- 1-cub_mort

# Next paper- Bright Ross 2020-----
bright_fig3ai <- scatter$`3307_Bright_2020_fig3ai.png`
bright_fig3ai_clean <- scat.clean(bright_fig3ai, remove.id= FALSE)    # survival rates

bright_fig3aii <- scatter$`3307_Bright_2020_fig3aii.png`
bright_fig3aii_clean <- scat.clean(bright_fig3aii, remove.id= FALSE)

bright_fig3aiii <- scatter$`3307_Bright_2020_fig3aiii.png`
bright_fig3aiii_clean <- scat.clean(bright_fig3aiii, remove.id= FALSE)

bright_fig3aiv <- scatter$`3307_Bright_2020_fig3aiv.png`
bright_fig3aiv_clean <- scat.clean(bright_fig3aiv, remove.id= FALSE)

# converting to wide format
bright_fig3ai <- widen(bright_fig3ai_clean, "Survival") #correct - numeric and 

bright_fig3aii <- widen(bright_fig3aii_clean, "Survival")

bright_fig3aiii <- widen(bright_fig3aiii_clean, "Survival")

bright_fig3aiv <- widen(bright_fig3aiv_clean, "Survival")

bright_fig3a.1 <- full_join(bright_fig3ai, bright_fig3aii, by= "Year")
bright_fig3a.2 <- full_join(bright_fig3aiii, bright_fig3aiv, by= "Year")
bright_survival <- full_join(bright_fig3a.1, bright_fig3a.2, by= "Year")


# bright survival rate estimation (average) ----
bright_survival <- rename(bright_survival, "Female_8+_yrs" = "Female_8+")
yf_s <- mean(bright_survival$Female_1_2yrs) # av of 1-2 year survival rates - currently a character, need as numeric
ym_s <- mean(bright_survival$Male_1_2yrs)   
af_s <- mean(c(bright_survival$Female_3_4yrs, bright_survival$Female_5_7yrs, bright_survival$`Female_8+`)) # mean of all adult survival rates
am_s <- mean(c(bright_survival$Male_3_4yrs, bright_survival$Male_5_7yrs, bright_survival$`Male_8+yrs`))

bright_survival_vec <- (c(yf_s, af_s, ym_s, am_s)) # order as vital rates - yf, af, ym ,am

# Macdonald survival----
mac_2002_survival <- as.data.frame(read.csv("Data/macdonald_2002_survival.csv"))
macs1 <- mean(mac_2002_survival$f_adult_survival)
macs2 <- mean(mac_2002_survival$m_adult_survival)
macs3 <- mean(mac_2002_survival$total_cub_survival)
mac_survival_vec <- c(macs1, macs2, macs3)  
# could use adult survival from this paper, and yearling from other?


# beta estimation
mcdonald_demo <- as.data.frame(read.csv("Data/mcdonald_2016_supinfo.csv"))  # posterior estimates from IPM

dens_posterior <- mcdonald_demo$pop_size   # posterior estimate from model
dens_posterior_mean <- mean(dens_posterior)
dens_posterior_sd <- sd(dens_posterior)

beta_reported <- 0.239  # from paper 

beta <- beta_reported/dens_posterior_sd  # explain maths in appendix!

