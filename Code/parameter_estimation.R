# Data digitisation and parameter estimation
# Focus on what you need!
# 13.12.2025

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
figs <- metaDigitise(dir ="C:/Users/jaycr/OneDrive - Nexus365/UNI/Masters/Project resources/Badger-MBIOL/Data/FigExtraction",
                     summary = FALSE)# summary = false imports raw data

# splitting into plot types is easier to work with
scatter <- figs$scatterplot
mean_err <- figs$mean_error


# organising paper 1  - Rogers 1997 ----
rogers_1997_fig1 <- figs$scatterplot$`242_Rogers_1997_fig1.png`
rogers_1997_popsize <- scat.clean(rogers_1997_fig1, remove.id = FALSE, dec= 0)  # YEARS STILL INCORRECT - ,3 onwards should all be 1 year later for 
rogers_1997_popsize <- rename(rogers_1997_popsize, MNA=  "Minimum_Number_Alive")
# order years correctly
# rogers_1997_popsize <- arrange(rogers_1997_popsize, by= Year, group_by(id)) 

# reorganising so fewer rows
rogers_1997_popsize_wide <- spread(rogers_1997_popsize, id, MNA)
rogers_1997_popsize_wide$Year <- as.numeric(rogers_1997_popsize_wide$Year)
# renaming rows
names(rogers_1997_popsize_wide) <- gsub("\\s+", "_", names(rogers_1997_popsize_wide))


# survival and population size
rogers_1997_s <- mean_err$`242_Rogers_1997_fig3.png` # prob survival
# renaming cols, joining with MNA estimates
rogers_1997_s <- rename(rogers_1997_s, Year= id, "prob of survival" = mean)
rogers_1997_s <- select(rogers_1997_s, -variable, -n)
names(rogers_1997_s) <- gsub("\\s+", "_", names(rogers_1997_s))

# reproductive rates
rogers_1997_rep1 <- scatter$`242_Rogers_1997_fig10.png` # % Af breeding
rogers_1997_rep1_clean <- scat.clean(rogers_1997_rep1, dec = 0)

rogers_1997_rep2 <- scatter$`242_Rogers_1997_fig12.png` # rep rate (cubs per adult?)
rogers_1997_rep2_clean <- scat.clean(rogers_1997_rep2, dec = 0)

rogers_1997_rep3<- scatter$`242_Rogers_1997_fig13.png` # Cubs / breeding fem
rogers_1997_rep3_clean <- scat.clean(rogers_1997_rep3, dec = 0)
rogers_1997_rep3_clean <- rename(rogers_1997_rep3_clean, "Cubs_per_breeding_fem"= "Number_of_cubs_per_breeding_female")
rogers_1997_rep3_clean$Cubs_per_breeding_fem[2] <- NA # removing outlier point


rogers_1997_rep4 <- mean_err$`242_Rogers_1997_fig11.png` # number of reproducing fems per social group 
rogers_1997_rep4_clean <- rename(rogers_1997_rep4, Year= id, "Av_reproducing_fems" = mean)
rogers_1997_rep4_clean <- dplyr::select(rogers_1997_rep4_clean, -c(n, variable))

# joining into single df 
rogers_1997_rep <- full_join(rogers_1997_rep1_clean, rogers_1997_rep3_clean,
                         by= "Year")
rogers_1997_rep <- full_join(rogers_1997_rep, rogers_1997_rep2_clean,
                                 by= "Year")
rogers_1997_rep$Year <- as.numeric(rogers_1997_rep$Year)
rogers_1997_rep4_clean$Year <- as.numeric(rogers_1997_rep4_clean$Year)  # must both be numeric 
rogers_1997_rep <- full_join(rogers_1997_rep, rogers_1997_rep4_clean,
                             by= "Year")

rogers_1997_s$Year <- as.numeric(rogers_1997_s$Year)
rogers_1997_data <- full_join(rogers_1997_s, rogers_1997_rep, by= "Year")
rogers_1997_data <- full_join(rogers_1997_data, rogers_1997_popsize_wide, by= "Year")

# visualising 
rogers_N_plot <- plot(x=rogers_1997_data$Year, y=rogers_1997_data$Total_Adult_Badgers)  # pop size over time increases linearly

# variables with density 
rogers_breeding_plot <- ggplot(rogers_1997_data, aes(x=Total_Adults, y= Percentage_of_Adult_Females_Breeding)) +
                                 geom_point()+
                                 geom_smooth(method = 'lm', formula = y~x, se= FALSE) +
    theme_classic() 
                           
rogers_reprate_plot <- ggplot(rogers_1997_data, aes(x=Adult_density, y=Reproductive_rate)) +
  geom_point()+
  geom_smooth(method = 'lm', formula = y~x, se= FALSE) +
    theme_classic()  # decline

rogers_fec_plot <- ggplot(rogers_1997_data, aes(x=Adult_density, y=Cubs_per_breeding_fem)) +
                             geom_point()+ 
    geom_smooth(method = 'lm', formula = y~x, se= FALSE) +
    theme_classic() 
                             
                       
rogers_fems_plot <- ggplot(rogers_1997_data, aes(x=Adult_density, y=Av_reproducing_fems)) +
    geom_point()+ 
    geom_smooth(method = 'lm', formula = y~x, se= FALSE) +
    theme_classic() # breeding fems per social group increases over time

grid.arrange(rogers_breeding_plot, rogers_reprate_plot, rogers_fec_plot, rogers_fems_plot, nrow=2, ncol=2)
# fairly constant % fems breeding in pop, and cohort size
# Rep rate decreases with density, reproducing fems per social group increases

breeding.lm <- lm(Percentage_of_Adult_Females_Breeding~Adult_density, data= rogers_1997_data)  # does density affect percentage adult fems reproducing
summary(breeding.lm)  # intercept sig, slope not. Low r squared, not important

cohort.lm <- lm(Number_of_cubs_per_breeding_female~Adult_density, data= rogers_1997_data)  # does density affect litter size
summary(cohort.lm)  # intercept sig, slope not. Low r squared. Int = 2.206774   

fecundity.lm <- lm(Reproductive_rate~Adult_density, data= rogers_1997_data)  # does density affect recruitment
summary(fecundity.lm)  # non significant

fems.lm <- lm(Av_reproducing_fems~Adult_density, data= rogers_1997_data)  # does pop size affect reproducing females per group?
summary(fems.lm) # sign positive slope - grr!


# Next paper - Macdonald 2013----

# extracting pop size and lambda from sup info instead of fig
library(readr)
mcdonald_demo <- as.data.frame(read_csv("Data/mcdonald_2016_supinfo.csv"))
mcdonald_demo$pop_size

# fig 2a covariate effect estimate
mcdonald_fig2a <- mean_err$`443_Mcdonald_2016_fig2a.jpg`  # this is on log scale, convert back to value?
mcdonald_fig2a_clean <- select(mcdonald_fig2a, -"n")
mcdonald_fig2a_clean <-  rename(mcdonald_fig2a_clean, "log_mean" = "mean")
  
mcdonald_fig2 <-scatter$`443_Mcdonald_2016_figS2.png`
mcdonald_fig2_clean <- select(mcdonald_fig2, -c("id", "group", "col",  "pch", "x_variable", "y_variable"))  # effect of density on recruitment
mcdonald_fig2_clean <-  rename(mcdonald_fig2_clean, "Standardised_density" = "x", "Recruitment" = "y") 


plot(x= mcdonald_demo$Year, y= mcdonald_demo$pop_size)  # decline at 170 indivs = K?
plot(x= mcdonald_fig2_clean$`Standardised Density`, y= mcdonald_fig2_clean$Recruitment, ylab= "Recruitment", xlab= "Standardised density") # issue - standardised density, how to convert back?

# how to convert standardised density to absolute? can't merge df, recruitment recordings not ordered by year
# rescaling pop size to mean 0 and sd 1. 
mcdonald_fig2_clean$Standardised_density
scaled_density <- scale(mcdonald_demo$pop_size) # values do not match - paper values have lower negative range  
sort(scaled_density, decreasing = FALSE) 



# Next = Tuyttens 2000----
# link density each year (av?) to changes in reproduction and social group size
tuyttens_fig1a <- scatter$`567_Tuyttens_2000_fig1a.png`  # density per yer in WW
tuyttens_fig1a_clean <- scat.clean(tuyttens_fig1a, y_name= "WW", remove.id= FALSE, dec= NULL)
# taking the mean of year measurements for each class
tuyttens_fig1a_mean <- aggregate(WW ~ id + Year, data = tuyttens_fig1a_clean, FUN = mean)


tuyttens_fig1b <- scatter$`567_Tuyttens_2000_fig1b.png`
tuyttens_fig1b_clean <- scat.clean(tuyttens_fig1b, y_name= "WP", remove.id= FALSE, dec = NULL)  # density per yer in WW
tuyttens_fig1b_mean <- aggregate(WP ~ id + Year, data = tuyttens_fig1b_clean, FUN = mean)

tuyttens_fig1c <- scatter$`567_Tuyttens_2000_fig1c.png`
tuyttens_fig1c_clean <- scat.clean(tuyttens_fig1c, y_name= "NN", remove.id= FALSE, dec= NULL)
tuyttens_fig1c_mean <- aggregate(NN ~ id + Year, data = tuyttens_fig1c_clean, FUN = mean)

# need estimates of mean density in each area to link rep values each year at sites
tuyttens1 <- full_join(tuyttens_fig1a_mean, tuyttens_fig1b_mean, by= c("id","Year"))
tuyttens1 <- full_join(tuyttens1, tuyttens_fig1c_mean, by= c("id","Year")) #combined dataset - density at each site per year, split by age class (wide)

# instead = year col, id col, site cols, density col
tuyttens1_long <- gather(tuyttens1, Site, Density, c("WW",  "WP",  "NN"))
tuyttens1_long <- spread(tuyttens1_long, id, Density) %>% 
  rename("Cub_density" = "Cubs",  "Adult_density" = "Adults", "Total_density" = "Total")
  
# options = mean of each year to summarise abundance each year, or keep each observation with vec naming? 
# mean better as each df after uses single recordings per site

tuyttens_fig4 <- mean_err$`567_Tuyttens_2000_fig4.png`
tuyttens_fig4_clean <- rename(tuyttens_fig4, Year= id, "Av_Badgers_per_social_group" = mean)
tuyttens_fig4_clean <- dplyr::select(tuyttens_fig4_clean, -c(n, variable))  # why is site not a colunm? need to separate out

tuyttens_fig5 <- scatter$`567_Tuyttens_2000_fig5.png` # also need to keep groupings, cant remove id
tuyttens_fig5_clean <- scat.clean(tuyttens_fig5, remove.id= FALSE, dec = NULL) # no rounding years for now
tuyttens_fig5_clean$Percentage_adult_fems_lactating[2] <- 100  # rounding one data point to 100 (down)

tuyttens_fig6b <- scatter$`567_Tuyttens_2000_fig6b.png`
tuyttens_fig6b_clean <- scat.clean(tuyttens_fig6b, remove.id= FALSE, dec = NULL)
tuyttens_fig6b_clean <- rename(tuyttens_fig6b_clean, "Fem_cubs_per_Fem_adults"= "Fem_cubs/_Fem_adults") 
tuyttens_fig6b_clean$Cubs_per_Fem_adults <- 2*(tuyttens_fig6b_clean$`Fem_cubs_per_Fem_adults`) # assuming equal sex ratio

tuyttens_fig6c <- scatter$`567_Tuyttens_2000_fig6c.png`
tuyttens_fig6c_clean <- scat.clean(tuyttens_fig6c, remove.id= FALSE, dec = NULL) 

# joining into df
tuyttens_rep <- full_join(tuyttens_fig5_clean, tuyttens_fig6b_clean, by = c("id", "Year"))
tuyttens_rep <- full_join(tuyttens_rep, tuyttens_fig6c_clean, by = c("id", "Year")) %>% 
  rename("Site" = "id")

# joining densities to rep df
tuyttens_data <- full_join(tuyttens1_long, tuyttens_rep, by = c("Year", "Site"))

# Linking changes in total or adult density to rep params - visualising
par(mfrow = c(2,2))
(breedingfems_plot <- ggplot(tuyttens_data, aes(x = Total_density, y = Percentage_adult_fems_lactating, col = Site)) +
                             geom_point() +  # correlating adult density and 
                             geom_smooth(method = 'lm'))
(fec_plot <- ggplot(tuyttens_data, aes(x = Total_density, y = Cubs_per_Fem_adults, col = Site)) +
  geom_point() +
  geom_smooth(method= 'lm'))

(litter_plot <- ggplot(tuyttens_data, aes(x = Total_density, y = Cubs_per_lactating_fem, col = Site)) +
    geom_point() +
    geom_smooth(method= 'lm'))

# linear models to get parameter estimates
fec.lm <- lm(Cubs_per_Fem_adults~Total_density + Site, data = tuyttens_data)  # linear model for fecundity (cubs / adult fem) and density, 
                                                                              # site as a covariate
summary(fec.lm) # Intercept (0.84), slope 0.07

litter.lm <- lm(Cubs_per_lactating_fem~Total_density + Site, data = tuyttens_data)  # linear model for cohort size (cubs / breeding fem) and density, 
                                                                                    # site as a covariate
summary(litter.lm)  # intercept 1.8, slope = 0.187. using max cohort size as cohort size when density = 0

breeding.lm <- lm(Percentage_adult_fems_lactating~Total_density + Site, data = tuyttens_data)  # linear model for percemtage reproducing and density, 
                                                                                               # site as a covariate
summary(breeding.lm)  # intercept = 68.7, slope = 3.31  


# Next paper- Bright Ross 2020-----
bright_fig2 <- scatter$`3307_Bright_2020_fig2.png`
bright_fig2_clean <- scat.clean(bright_fig2, remove.id = FALSE, dec=0)  # years are not in order :(
bright_fig2_wide <- spread(bright_fig2_clean, "id", "Density")  # MNA for Adult males and fems, and total 
# can't spread?
# ----

bright_fig3ai <- scatter$`3307_Bright_2020_fig3ai.png`
bright_fig3ai_clean <- scat.clean(bright_fig3ai, remove.id= FALSE)    # survival rates

bright_fig3aii <- scatter$`3307_Bright_2020_fig3aii.png`
bright_fig3aii_clean <- scat.clean(bright_fig3aii, remove.id= FALSE)

bright_fig3aiii <- scatter$`3307_Bright_2020_fig3aiii.png`
bright_fig3aiii_clean <- scat.clean(bright_fig3aiii, remove.id= FALSE)

bright_fig3aiv <- scatter$`3307_Bright_2020_fig3aiv.png`
bright_fig3aiv_clean <- scat.clean(bright_fig3aiv, remove.id= FALSE)

# converting to wide format
bright_fig3ai_wide <- widen(bright_fig3ai_clean, "Survival") #correct - numeric and 

bright_fig3aii_wide <- widen(bright_fig3aii_clean, "Survival")

bright_fig3aiii_wide <- widen(bright_fig3aiii_clean, "Survival")

bright_fig3aiv_wide <- widen(bright_fig3aiv_clean, "Survival")

bright_fig3a.1 <- full_join(bright_fig3ai_wide, bright_fig3aii_wide, by= "Year")
bright_fig3a.2 <- full_join(bright_fig3aiii_wide, bright_fig3aiv_wide, by= "Year")
bright_survival <- full_join(bright_fig3a.1, bright_fig3a.2, by= "Year")



# Rep by year
bright_fig3bi <- scatter$`3307_Bright_2020_fig3bi.png`
bright_fig3bi_clean <- scat.clean(bright_fig3bi, remove.id= FALSE)  # percentage producing offspring - need to keep ids seperated

bright_fig3bii <- scatter$`3307_Bright_2020_fig3bii.png`
bright_fig3bii_clean <- scat.clean(bright_fig3bii, remove.id= FALSE)

bright_fig3biii <- scatter$`3307_Bright_2020_fig3biii.png`
bright_fig3biii_clean <- scat.clean(bright_fig3biii, remove.id= FALSE)

bright_fig3biv <- scatter$`3307_Bright_2020_fig3biv.png`
bright_fig3biv_clean <- scat.clean(bright_fig3biv, remove.id= FALSE)


# widen to seperate sex in col, age in col to join df
bright_fig3bi_wide <- widen(bright_fig3bi_clean, "Percentage_producing_offspring")

bright_fig3bii_wide <- widen(bright_fig3bii_clean, "Percentage_producing_offspring")

bright_fig3biii_wide <- widen(bright_fig3biii_clean, "Percentage_producing_offspring")

bright_fig3biv_wide <- widen(bright_fig3biv_clean, "Percentage_producing_offspring")

# joining into reproduction df (percentage producing offspring)
bright_fig3.1 <- full_join(bright_fig3bi_wide, bright_fig3bii_wide, by= "Year")
bright_fig3.2 <- full_join(bright_fig3biii_wide, bright_fig3biv_wide, by= "Year")
bright_rep <- full_join(bright_fig3.1, bright_fig3.2, by= "Year")

names(bright_survival)
names(bright_rep)  # cols dont match first in s is 1-2, rep has only 2+

# survival rate estimation (average)
yf_s <- mean(bright_survival$Female_1_2yrs) # av of 1-2 year survival rates - currently a character, need as numeric
ym_s <- mean(bright_survival$Male_1_2yrs)   
af_s <- mean(c(bright_survival$Female_3_4yrs, bright_survival$Female_5_7yrs, bright_survival$`Female_8+`)) # mean of all adult survival rates
am_s <- mean(c(bright_survival$Male_3_4yrs, bright_survival$Male_5_7yrs, bright_survival$`Male_8+yrs`))


# Mcdonald 2002 Wytham paper ----
macdonald_fig2 <- scatter$Macdonald_2002_fig2.png
macdonald_fig2_clean <- scat.clean(macdonald_fig2)
macdonald_fig2_clean <- rename(macdonald_fig2_clean, "Total_MNA" = "Number_of_badgers") 

macdonald_fig7 <- scatter$Macdonald_2002_fig7.png
macdonald_fig7_clean <- scat.clean(macdonald_fig7)

macdonald_fig8 <- scatter$Macdonald_2002_fig8.png
macdonald_fig8_clean <- scat.clean(macdonald_fig8)
macdonald_fig8_clean <- rename(macdonald_fig8_clean, "Net_reproductive_rate" = "Reproductive_rates") 

macdonald_fig8b <- scatter$Macdonald_2002_fig8b.png
macdonald_fig8b_clean <- scat.clean(macdonald_fig8b)

# combining into large dataframe
macdonald_data <- full_join(macdonald_fig2_clean, macdonald_fig7_clean, by= "Year")
macdonald_data <- full_join(macdonald_data, macdonald_fig8_clean, by= "Year")
macdonald_data <- full_join(macdonald_data, macdonald_fig8b_clean, by= "Year")
macdonald_data <- arrange(macdonald_data, Year) #reorder by ascending year 

# visualising
plot(macdonald_data$Year, y=macdonald_data$Total_MNA)
plot(macdonald_data$Total_MNA, y=macdonald_data$Mean_group_size_MNA) # increasing group size with higher pop size
plot(macdonald_data$Total_MNA, y=macdonald_data$Net_reproductive_rate)
plot(macdonald_data$Total_MNA, y=macdonald_data$Number_of_cubs_per_reproductive_fem) 
# increasing reproduction with group size??