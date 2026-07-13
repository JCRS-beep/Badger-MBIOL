# Run script for density-dependent badger population modeling with removals
# Written by Jay Creese with input from Chrissy Hernandez
# July 2026
# Last edited: July 2026

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
source(here("Code/02_data_extraction.R"))  # select options 2 (existing data) and 1 (all data)

# parameter definition and df set ups
source(here("Code/03_parameter_rate_setup.R"))


# following sections only relevant for MBIOL project - will be updating the removal scenarios with new functions and different removal biases.
# setting up model scenarios
source(here("Code/code_MBIOL_final//04_model_projections.R"))  

# analysis script
source(here("Code/code_MBIOL_final//05_analysis.R"))  


