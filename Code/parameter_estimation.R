# Calculating parameters for density dependence
# 11.11.2025


# load the data
DD_data<- read.csv("Data/DD_vital_rates.csv")   # data from Tuyttens 2000
head(DD_data)


byear<- plot(x=DD_data$Year, y= DD_data$Births_f) # changes in births per year- link? 


# linking cub production AND SURVIVAL to population density 
bd_plot<- plot(x= DD_data$Density, y=DD_data$Births_f)    # births against density- no clear pattern seen

# loading ggplot for nicer graphs
library(ggplot2)

(birth_plot<- ggplot(data= DD_data, aes(x=Density, y=Births_f, colour= Site))+  # seeing if site explains some variation
  geom_point())     
# 
