# run script 
# temporay solution to issues in report .rmd

## section 1 -  loading libraries

library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)
library(dplyr)
library(readr)
library(here)


# Sourcing other scripts - 02 not working? no digitised papers?

# # sourcing functions used throughout 
source(here("Code/Functions/01_all_functions.R")) 

# data extraction
source(here("Code/code_MBIOL_final/02_data_extraction.R"))  # select options 2 (existing data) and 1 (all data)

# parameter definition and df set ups
source(here("Code/code_MBIOL_final//03_parameter_rate_setup.R"))


# following sections only relevant for MBIOL project - will be updating the removal scenarios with new functions and different removal biases.
# setting up model scenarios
source(here("Code/code_MBIOL_final//04_model_projections.R"))  

# analysis script
source(here("Code/code_MBIOL_final//05_analysis.R"))  


