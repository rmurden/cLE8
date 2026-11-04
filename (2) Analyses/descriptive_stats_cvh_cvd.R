##########################################################################
######  Descriptive Statistics Table for cLE8 Scores and CVD        ######
######    Outcomes, Stratified by State                             ######
#################  Author: Raphiel J. Murden;           ##################
##########################################################################
rm(list = ls(all = T))
set.seed(217)

## Data Prep
# Load required packages
library(tidyverse); library(readxl); library(arsenal);

# 1. Load your dataset
setwd( "C:/Users/rmurden/OneDrive - Emory/Documents/GitHub/cLE8")
data.1 <- read.csv("Data/cle8_three_scores_txgaca_within_state.csv")
names(data.1)
dim(data.1) # 471 obs. (3=471 counties)
data.1$GEOID <- sprintf("%05d", as.numeric(data.1$CountyFIPS))

## Merge with outcome data
outcome.dat.cvd <- read.csv("Data/cvd_mortality_2021-2023.csv") %>%
  dplyr::select(County.Code, Crude.Rate, Deaths, Population) %>%
  rename(c("County.Code" = "County.Code", "CVD_mortality" = "Crude.Rate")) %>%
  transform(CVD_mortality = as.numeric(CVD_mortality), Deaths = as.numeric(Deaths), Population = as.numeric(Population))
outcome.dat.all.cause <- read.csv("Data/all_cause_mortality_2021-2023.csv") %>%
  dplyr::select(County.Code, Crude.Rate, Deaths, Population) %>%
  rename(c("County.Code" = "County.Code", "AllCause_mortality" = "Crude.Rate")) %>%
  transform(AllCause_mortality = as.numeric(AllCause_mortality), Deaths = as.numeric(Deaths), Population = as.numeric(Population))

outcome.dat.morbidity <- read.csv("Data/PLACES__County_Data_(GIS_Friendly_Format),_2024_release_20260622.csv") %>%
  dplyr::select(CHD_AdjPrev, CHD_CrudePrev, STROKE_CrudePrev, STROKE_AdjPrev, CountyFIPS)

## Merge into analysis data set
data <- data.1 %>%
  left_join(outcome.dat.cvd, by = c("CountyFIPS" = "County.Code")) %>%
  filter(!is.na(CVD_mortality))

# 1a. Descriptive Statistics Table
descriptive_stats <- data %>%
  left_join(outcome.dat.all.cause, by = c("CountyFIPS" = "County.Code")) %>%
  left_join(outcome.dat.morbidity, by = c("CountyFIPS")) %>%
  tableby(StateAbbr ~ cle8_point + CVD_mortality + AllCause_mortality + CHD_AdjPrev + STROKE_AdjPrev, data = .,
  control = tableby.control(test = FALSE, total = TRUE, numeric.stats = c("meansd", "range"), digits = 2))
summary(descriptive_stats)
# write2word(summary(descriptive_stats), file = "Output/Results/Descriptive_Stats_CVH_CVD.doc")
