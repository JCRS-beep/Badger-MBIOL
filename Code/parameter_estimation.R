# Calculating parameters for density dependence
# 13.12.2025

# installing digitise package
# install.packages("devtools")
# devtools::install_github("daniel1noble/metaDigitise")

# packages ----
library(metaDigitise)
library(dplyr)
library(tidyr)

# function - for scatter, rename, round, remove cols.----
scat.clean <- function(df,y_name = NULL, year= TRUE, remove.id= TRUE) {
  if(is.null(y_name)){
  yname <- unique(df$y_variable)}
  else{
  yname <- y_name
  }
  names(df)[names(df) == "y"] <- yname
  
  
  df <- dplyr::select(df, -c(col, pch, group,  y_variable, x_variable))
  
  if(year == TRUE){   
    df <- rename(df, Year = x)
    df$Year <- as.numeric(round(df$Year, 0))
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

# manually digitising figs ----
figs <- metaDigitise(dir ="C:/Users/jaycr/OneDrive - Nexus365/UNI/Masters/Project resources/Badger-MBIOL/Data/FigExtraction",
                     summary = FALSE)# summary = false imports raw data

# splitting into plot types is easier to work with
scatter <- figs$scatterplot
mean_err <- figs$mean_error


# organising paper 1  - Rogers 1997 ----
rogers_1997_fig1 <- figs$scatterplot$`242_Rogers_1997_fig1.png`
rogers_1997_popsize <- rename(rogers_1997_fig1, Year= x, MNA= y) # renaming to clearer cols

rogers_1997_popsize$Year <- as.numeric(round(rogers_1997_popsize$Year, 0)) # removing decimal places from Year

# DESIRED = id, Year, MNA
rogers_1997_popsize <- dplyr::select(rogers_1997_popsize, id, Year, MNA) 

# reorganising so fewer rows
rogers_1997_popsize_wide <- spread(rogers_1997_popsize, id, MNA)
rogers_1997_popsize_wide$Year <- as.numeric(rogers_1997_popsize_wide$Year)
# renaming rows
names(rogers_1997_popsize_wide) <- gsub("\\s+", "_", names(rogers_1997_popsize_wide))


# survival and population size ----
rogers_1997_s <- mean_err$`242_Rogers_1997_fig3.png` # prob survival
# renaming cols, joining with MNA estimates
rogers_1997_s <- rename(rogers_1997_s, Year= id, "prob of survival" = mean)
rogers_1997_s <- select(rogers_1997_s, -variable, -n)
names(rogers_1997_s) <- gsub("\\s+", "_", names(rogers_1997_s))

# reproductive rates ----
rogers_1997_rep1 <-scatter$`242_Rogers_1997_fig10.png` # % Af breeding
rogers_1997_rep1_clean <- scat.clean(rogers_1997_rep1)


rogers_1997_rep2 <- scatter$`242_Rogers_1997_fig12.png` # rep rate 
rogers_1997_rep2_clean <- scat.clean(rogers_1997_rep2)

rogers_1997_rep3<- scatter$`242_Rogers_1997_fig13.png` # Cubs / breeding fem
rogers_1997_rep3_clean <- scat.clean(rogers_1997_rep3)

rogers_1997_rep4 <- mean_err$`242_Rogers_1997_fig11.png` # rep fems
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

rogers_plot <- plot(x=rogers_1997_data$Percentage_of_Adult_Females_Breeding, y=rogers_1997_data$Total_Adult_Badgers) 
# no clear relationship shown in these plots
# Analysis, test which variable is greatest predictor for cubs per female, percentage breeding and

# Next paper - Macdonald 2013----
mcdonald_fig1 <- scatter$`443_Mcdonald_2016_fig1.png`
mcdonald1_clean <- scat.clean(mcdonald_fig1)

mcdonald_fig2e <- mean_err$`443_Mcdonald_2016_fig2e.png` # Regression slopes (b) for relationship between demographic rates;recruitment (R) and survival (S), and covariate effect density (N). The posterior mean is displayed alongside the corresponding 95%credible intervals.
mcdonald_fig2e_clean <- select(mcdonald_fig2e, -variable, -n)

mcdonald_fig2 <-scatter$`443_Mcdonald_2016_figS2.png`
mcdonald_fig2_clean <- select(mcdonald_fig2, -c("group", "col",  "pch"))  # mb shouldnt be a df? effect of density on recruitment



# Tuyttens 2000 - must link density each year (av?) to changes in reproduction and social group size
tuyttens_fig1a <- scatter$`567_Tuyttens_2000_fig1a.png`  # density per yer in WW
tuyttens_fig1a_clean <- scat.clean(tuyttens_fig1a, y_name= "WW_badgers/km^2", remove.id= FALSE)  # rename x= Year, y= yvar, remove cols
tuyttens_1a_wide <- spread(tuyttens_fig1a_clean, id, "WW_badgers/km^2") # doesn't work if multiple readings per year


tuyttens_fig1b <- scatter$`567_Tuyttens_2000_fig1b.png`
tuyttens_fig1b_clean <- scat.clean(tuyttens_fig1b, y_name= "WP_badgers/km^2", remove.id= FALSE)  # density per yer in WW

tuyttens_fig1c <- scatter$`567_Tuyttens_2000_fig1c.png`
tuyttens_fig1c_clean <- scat.clean(tuyttens_fig1c, y_name= "NN_badgers/km^2", remove.id= FALSE)
# need estimates of density in each area to link rep values each year at sites

tuyttens_fig4 <- mean_err$`567_Tuyttens_2000_fig4.png`
tuyttens_fig4_clean <- rename(tuyttens_fig4, Year= id, "Av_Badgers_per_social_group" = mean)
tuyttens_fig4_clean <- dplyr::select(tuyttens_fig4_clean, -c(n, variable))

tuyttens_fig5 <- scatter$`567_Tuyttens_2000_fig5.png` # also need to keep groupings, cant remove id
tuyttens_fig5_clean <- scat.clean(tuyttens_fig5, remove.id= FALSE)

tuyttens_fig6b <- scatter$`567_Tuyttens_2000_fig6b.png`
tuyttens_fig6b_clean <- scat.clean(tuyttens_fig6b, remove.id= FALSE)

tuyttens_fig6c <- scatter$`567_Tuyttens_2000_fig6c.png`
tuyttens_fig6c_clean <- scat.clean(tuyttens_fig6c, remove.id= FALSE)

tuyttens_rep <- full_join(tuyttens_fig5_clean, tuyttens_fig6b_clean, by = c("id", "Year"))
tuyttens_rep <- full_join(tuyttens_rep, tuyttens_fig6c_clean, by = c("id", "Year"))

# Final paper
bright_fig2 <- scatter$`3307_Bright_2020_fig2.png`
bright_fig2_clean <- scat.clean(bright_fig2, remove.id = FALSE) # years are not in order

bright_fig3ai <- scatter$`3307_Bright_2020_fig3ai.png`
bright_fig3ai_clean <- scat.clean(bright_fig3ai, remove.id= FALSE)

bright_fig3aii <- scatter$`3307_Bright_2020_fig3aii.png`
bright_fig3aii_clean <- scat.clean(bright_fig3aii, remove.id= FALSE)

bright_fig3aiii <- scatter$`3307_Bright_2020_fig3aiii.png`
bright_fig3aiii_clean <- scat.clean(bright_fig3aiii, remove.id= FALSE)

bright_fig3aiv <- scatter$`3307_Bright_2020_fig3aiv.png`
bright_fig3aiv_clean <- scat.clean(bright_fig3aiv, remove.id= FALSE)

bright_fig3bi <- scatter$`3307_Bright_2020_fig3bi.png`
bright_fig3bi_clean <- scat.clean(bright_fig3bi, remove.id= FALSE)

bright_fig3bii <- scatter$`3307_Bright_2020_fig3bii.png`
bright_fig3bii_clean <- scat.clean(bright_fig3bii, remove.id= FALSE)

bright_fig3biii <- scatter$`3307_Bright_2020_fig3biii.png`
bright_fig3biii_clean <- scat.clean(bright_fig3biii, remove.id= FALSE)

bright_fig3biv <- scatter$`3307_Bright_2020_fig3biv.png`
bright_fig3biv_clean <- scat.clean(bright_fig3biv, remove.id= FALSE)