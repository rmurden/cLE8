##########################################################################
######  Fit models for MDPI Paper with CVD Outcomes and report fit  ######
######    statistics                                                ######
#################  Author: Raphiel J. Murden;           ##################
##########################################################################
rm(list = ls(all = T))
set.seed(217)

## Data Prep
# Load required packages
library(tidyverse); library(sf); library(tigris); library(spdep); library(ggplot2); 
library(ggspatial); library(readxl); library(RColorBrewer); library(classInt);
require(INLA)
library(gridExtra)
options(tigris_use_cache = TRUE)

# 1. Load your dataset
setwd( "C:/Users/rmurden/OneDrive - Emory/Documents/GitHub/cLE8/")
data.1 <- read.csv("Data/cle8_three_scores_txgaca_within_state.csv")
names(data.1)
dim(data.1) # 471 obs. (3=471 counties)
data.1$GEOID <- sprintf("%05d", as.numeric(data.1$CountyFIPS))

## Create 65+ dataset
chrr_2024 <- read_xlsx("Data/analytic_data2024.xlsx")

pct_65plus <- chrr_2024 %>%
  select(County.Code = `5-digit FIPS Code`,pct_65_older = `% 65 and Older raw value`) %>%
  mutate(County.Code = as.numeric(County.Code),
         pct_65_older = as.numeric(str_replace(as.character(pct_65_older), "%", "")))

## Merge with outcome data
outcome.dat.all.cause <- read.csv("Data/all_cause_mortality_2021-2023.csv") %>%
  select(County.Code, Crude.Rate, Deaths, Population) %>%
  rename(c("County.Code" = "County.Code", "AllCause_mortality" = "Crude.Rate")) %>%
  transform(AllCause_mortality = as.numeric(AllCause_mortality), Deaths = as.numeric(Deaths), Population = as.numeric(Population))

outcome.covar.dat.all.cause<-outcome.dat.all.cause %>%
  left_join(pct_65plus, by = "County.Code")

## Merge both into analysis data set
data <- data.1 %>%
  left_join(outcome.covar.dat.all.cause, by = c("CountyFIPS" = "County.Code"))  %>%
  filter(!is.na(AllCause_mortality)) 
data_CA <- data %>% filter(StateAbbr == "CA") %>%
            mutate(county_int = as.numeric(factor(CountyFIPS)),
                    county_slope = county_int, 
                    ExpectedDeaths = sum(Deaths)/sum(Population)*Population)
data_GA <- data %>% filter(StateAbbr == "GA") %>%
            mutate(county_int = as.numeric(factor(CountyFIPS)),
                    county_slope = county_int, 
                    ExpectedDeaths = sum(Deaths)/sum(Population)*Population)
data_TX <- data %>% filter(StateAbbr == "TX") %>%
            mutate(county_int = as.numeric(factor(CountyFIPS)),
                    county_slope = county_int, 
                    ExpectedDeaths = sum(Deaths)/sum(Population)*Population)


# 2. Load and transform shapefiles
counties_sf <- counties(state = c("CA", "GA", "TX"), cb = TRUE, resolution = "5m", year = 2022) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  select(GEOID, geometry)

counties_sf_CA <- counties(state = c("CA"), cb = TRUE, resolution = "5m", year = 2022) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  select(GEOID, geometry)

counties_sf_GA <- counties(state = c("GA"), cb = TRUE, resolution = "5m", year = 2022) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  select(GEOID, geometry)

counties_sf_TX <- counties(state = c("TX"), cb = TRUE, resolution = "5m", year = 2022) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  select(GEOID, geometry)

# 3. Merge data with counties
merged_sf <- counties_sf %>%
  left_join(data, by = c("GEOID")) %>%
  filter(!is.na(AllCause_mortality))

merged_sf_CA <- counties_sf_CA %>%
  left_join(data_CA, by = c("GEOID")) %>%
  filter(!is.na(AllCause_mortality))

merged_sf_GA <- counties_sf_GA %>%
  left_join(data_GA, by = c("GEOID")) %>%
  filter(!is.na(AllCause_mortality))

merged_sf_TX <- counties_sf_TX %>%
  left_join(data_TX, by = c("GEOID")) %>%
  filter(!is.na(AllCause_mortality))

# 4. Choropleth map with color-blind friendly red-to-green palette
cb_palette <- c("#D73027", "#FC8D59", "#FEE08B", "#D9EF8B", "#91CF60", "#1A9850")

choropleth <- ggplot() +
  geom_sf(data = merged_sf, aes(fill = AllCause_mortality), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette, name = "cLE8 Score", na.value = "grey90") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  labs(
    title = "County-Level CVH Score (Repositioned AK & HI)",
    subtitle = "Color-blind friendly red-to-green scale with state borders",
    caption = "Data source: all_merged_2025-07-07 subset.xlsx"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(hjust = 0.5),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8)
  )

# ggsave("CVH_Score_Map_with_States.pdf", plot = choropleth, width = 11, height = 8.5, units = "in", dpi = 300)
# ggsave("CVH_Score_Map_with_States.png", plot = choropleth, width = 11, height = 8.5, units = "in", dpi = 300)

# 6. Project for spatial analysis
merged_proj <- st_transform(merged_sf, crs = 5070)
merged_proj_CA <- st_transform(merged_sf_CA, crs = 5070)
merged_proj_GA <- st_transform(merged_sf_GA, crs = 5070)
merged_proj_TX <- st_transform(merged_sf_TX, crs = 5070)

# 7. Create spatial weights
nb <- poly2nb(merged_proj)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

nb_CA <- poly2nb(merged_proj_CA)
nb_GA <- poly2nb(merged_proj_GA)
nb_TX <- poly2nb(merged_proj_TX)
lw_CA <- nb2listw(nb_CA, style = "W", zero.policy = TRUE)
lw_GA <- nb2listw(nb_GA, style = "W", zero.policy = TRUE)
lw_TX <- nb2listw(nb_TX, style = "W", zero.policy = TRUE)

# 8. Global Moran's I
moran_global <- moran.test(merged_proj$AllCause_mortality, lw, zero.policy = TRUE)
print(moran_global)

moran_global_CA <- moran.test(merged_proj_CA$AllCause_mortality, lw_CA, zero.policy = TRUE)
print(moran_global_CA)
moran_global_GA <- moran.test(merged_proj_GA$AllCause_mortality, lw_GA, zero.policy = TRUE)
print(moran_global_GA)
moran_global_TX <- moran.test(merged_proj_TX$AllCause_mortality, lw_TX, zero.policy = TRUE)
print(moran_global_TX)

# 9. Moran scatterplot
moran.plot(merged_proj$AllCause_mortality, lw, zero.policy = TRUE)

# 10. Local Moran's I
local_moran <- localmoran(merged_proj$AllCause_mortality, lw, zero.policy = TRUE)
merged_proj$local_I <- local_moran[, 1]
merged_proj$p_value <- local_moran[, 5]

local_moran_CA <- localmoran(merged_proj_CA$AllCause_mortality, lw_CA, zero.policy = TRUE)
merged_proj_CA$local_I <- local_moran_CA[, 1]
merged_proj_CA$p_value <- local_moran_CA[, 5]

local_moran_GA <- localmoran(merged_proj_GA$AllCause_mortality, lw_GA, zero.policy = TRUE)
merged_proj_GA$local_I <- local_moran_GA[, 1]
merged_proj_GA$p_value <- local_moran_GA[, 5]

local_moran_TX <- localmoran(merged_proj_TX$AllCause_mortality, lw_TX, zero.policy = TRUE)
merged_proj_TX$local_I <- local_moran_TX[, 1]
merged_proj_TX$p_value <- local_moran_TX[, 5]


# 11. Map of Local Moran's I (hotspots)
hotspot_map <- ggplot(merged_proj) +
  geom_sf(aes(fill = local_I), color = NA) +
  scale_fill_viridis_c(option = "plasma", name = "Local Moran's I") +
  labs(title = "Local Moran's I (Hotspot Detection)") +
  theme_minimal()

# ggsave("Local_Morans_I_Hotspots.png", plot = hotspot_map, width = 11, height = 8.5, units = "in", dpi = 300)

# Prepare data for modeling
nb2INLA("Data/map_CA.adj", nb_CA)
g_CA <- inla.read.graph("Data/map_CA.adj")

nb2INLA("Data/map_GA.adj", nb_GA)
g_GA <- inla.read.graph("Data/map_GA.adj")

nb2INLA("Data/map_TX.adj", nb_TX)
g_TX <- inla.read.graph("Data/map_TX.adj")

################          Models          ################
## Same GLM formula for all states
formula.0.glm <- Deaths ~ 1 + scale(pct_65_older)
formula.1.glm <- Deaths ~ 1 + scale(pct_65_older) + scale(cle8_point) 

## CA Models
formula.0.icar.CA <- Deaths ~ 1 + scale(pct_65_older) + f(county_int, model = "besag", graph = g_CA)
formula.1.icar.CA <- Deaths ~ 1 + scale(pct_65_older) + scale(cle8_point) + f(county_int, model = "besag", graph = g_CA)
formula.0.bym.CA <- Deaths ~ 1 + scale(pct_65_older) +f(county_int, model = "bym2", graph = g_CA)
formula.1.bym.CA <- Deaths ~ 1 + scale(pct_65_older) +scale(cle8_point) + f(county_int, model = "bym2", graph = g_CA)

## GA Models
formula.0.icar.GA <- Deaths ~ 1 + scale(pct_65_older) +f(county_int, model = "besag", graph = g_GA)
formula.1.icar.GA <- Deaths ~ 1 + scale(pct_65_older) +scale(cle8_point) + f(county_int, model = "besag", graph = g_GA)
formula.0.bym.GA <- Deaths ~ 1 + scale(pct_65_older) +f(county_int, model = "bym2", graph = g_GA)
formula.1.bym.GA <- Deaths ~ 1 + scale(pct_65_older) +scale(cle8_point) + f(county_int, model = "bym2", graph = g_GA)

## TX Models
formula.0.icar.TX <- Deaths ~ 1 + scale(pct_65_older) +f(county_int, model = "besag", graph = g_TX)
formula.1.icar.TX <- Deaths ~ 1 + scale(pct_65_older) +scale(cle8_point) + f(county_int, model = "besag", graph = g_TX)
formula.0.bym.TX <- Deaths ~ 1 + scale(pct_65_older) +f(county_int, model = "bym2", graph = g_TX)
formula.1.bym.TX <- Deaths ~ 1 + scale(pct_65_older) +scale(cle8_point) + f(county_int, model = "bym2", graph = g_TX)

##########   CA   ##########    
##### GLMs for All-cause mortality in CA
### Model 0: Intercept-only model
mod.0.CA.glm <- inla(formula.0.glm, family = "nbinomial", data = merged_proj_CA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.0.CA.glm)
sum(mod.0.CA.glm$cpo$cpo) 
hist(mod.0.CA.glm$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

### Model 1: CVH score as a covariate
mod.1.CA.glm <- inla(formula.1.glm, family = "nbinomial", data = merged_proj_CA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.1.CA.glm)
sum(mod.1.CA.glm$cpo$cpo) 
hist(mod.1.CA.glm$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

##### GLMMs with ICAR for  mortality in CA
### Model 0: Intercept-only model
mod.0.CA.icar <- inla(formula.0.icar.CA, family = "poisson", data = merged_proj_CA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.0.CA.icar)
sum(mod.0.CA.icar$cpo$cpo) 
hist(mod.0.CA.icar$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

### Model 1: CVH score as a covariate
mod.1.CA.icar <- inla(formula.1.icar.CA, family = "poisson", data = merged_proj_CA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.1.CA.icar)
sum(mod.1.CA.icar$cpo$cpo) 
hist(mod.1.CA.icar$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

##### GLMMs with Besag-York-Mollie for  mortality in CA
### Model 0: Intercept-only model
mod.0.CA.bym <- inla(formula.0.bym.CA, family = "poisson", data = merged_proj_CA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.0.CA.bym)
sum(mod.0.CA.bym$cpo$cpo) 
hist(mod.0.CA.bym$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

### Model 1: CVH score as a covariate
mod.1.CA.bym <- inla(formula.1.bym.CA, family = "poisson", data = merged_proj_CA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.1.CA.bym)
sum(mod.1.CA.bym$cpo$cpo) 
hist(mod.1.CA.bym$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

CA.mod.results = list("GLM 0 -CA" = mod.0.CA.glm, "GLM 1 -CA" = mod.1.CA.glm,
                     "ICAR 0 -CA" = mod.0.CA.icar, "ICAR 1 -CA" = mod.1.CA.icar,
                     "BYM 0 -CA" = mod.0.CA.bym, "BYM 1 -CA" = mod.1.CA.bym)

extract.model.results <- function(model) {
  return(c("DIC" = model$dic$dic, "pD" = model$dic$p.eff, "WAIC" = model$waic$waic,
  "pWAIC" = model$waic$p.eff,  "Log Marginal Likelihood" = model$mlik[2,1], "Sum CPO" = sum(model$cpo$cpo)))
}
extract.fixed.effects <- function(model) {
  return(model$summary.fixed)
}

CA.out = t(sapply(CA.mod.results, extract.model.results))
CA.out

CA.fixed = do.call(rbind, lapply(CA.mod.results, extract.fixed.effects))
##########   GA   ##########    
##### GLMs for All-cause mortality in GA
### Model 0: Intercept-only model
mod.0.GA.glm <- inla(formula.0.glm, family = "nbinomial", data = merged_proj_GA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.0.GA.glm)
sum(mod.0.GA.glm$cpo$cpo) 
hist(mod.0.GA.glm$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

### Model 1: CVH score as a covariate
mod.1.GA.glm <- inla(formula.1.glm, family = "nbinomial", data = merged_proj_GA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.1.GA.glm)
sum(mod.1.GA.glm$cpo$cpo) 
hist(mod.1.GA.glm$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

##### GLMMs with ICAR for All-cause mortality in GA
### Model 0: Intercept-only model
mod.0.GA.icar <- inla(formula.0.icar.GA, family = "poisson", data = merged_proj_GA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.0.GA.icar)
sum(mod.0.GA.icar$cpo$cpo) 
hist(mod.0.GA.icar$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

### Model 1: CVH score as a covariate
mod.1.GA.icar <- inla(formula.1.icar.GA, family = "poisson", data = merged_proj_GA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.1.GA.icar)
sum(mod.1.GA.icar$cpo$cpo) 
hist(mod.1.GA.icar$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

##### GLMMs with Besag-York-Mollie for All-cause mortality in GA
### Model 0: Intercept-only model
mod.0.GA.bym <- inla(formula.0.bym.GA, family = "poisson", data = merged_proj_GA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.0.GA.bym)
sum(mod.0.GA.bym$cpo$cpo) 
hist(mod.0.GA.bym$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

### Model 1: CVH score as a covariate
mod.1.GA.bym <- inla(formula.1.bym.GA, family = "poisson", data = merged_proj_GA, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.1.GA.bym)
sum(mod.1.GA.bym$cpo$cpo) 
hist(mod.1.GA.bym$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

GA.mod.results = list("GLM 0 -GA" = mod.0.GA.glm, "GLM 1 -GA" = mod.1.GA.glm,
                     "ICAR 0 -GA" = mod.0.GA.icar, "ICAR 1 -GA" = mod.1.GA.icar,
                     "BYM 0 -GA" = mod.0.GA.bym, "BYM 1 -GA" = mod.1.GA.bym)

GA.out = t(sapply(GA.mod.results, extract.model.results))
GA.out

GA.fixed = do.call(rbind, lapply(GA.mod.results, extract.fixed.effects))
##########   TX   ##########    
##### GLMs for All-cause mortality in TX
### Model 0: Intercept-only model
mod.0.TX.glm <- inla(formula.0.glm, family = "nbinomial", data = merged_proj_TX , E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.0.TX.glm)
sum(mod.0.TX.glm$cpo$cpo) 
hist(mod.0.TX.glm$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

### Model 1: CVH score as a covariate
mod.1.TX.glm <- inla(formula.1.glm, family = "nbinomial", data = merged_proj_TX, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.1.TX.glm)
sum(mod.1.TX.glm$cpo$cpo) 
hist(mod.1.TX.glm$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

##### GLMMs with ICAR for All-cause mortality in TX
### Model 0: Intercept-only model
mod.0.TX.icar <- inla(formula.0.icar.TX, family = "poisson", data = merged_proj_TX, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.0.TX.icar)
sum(mod.0.TX.icar$cpo$cpo) 
hist(mod.0.TX.icar$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

### Model 1: CVH score as a covariate
mod.1.TX.icar <- inla(formula.1.icar.TX, family = "poisson", data = merged_proj_TX, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.1.TX.icar)
sum(mod.1.TX.icar$cpo$cpo) 
hist(mod.1.TX.icar$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

##### GLMMs with Besag-York-Mollie for All-cause mortality in TX
### Model 0: Intercept-only model
mod.0.TX.bym <- inla(formula.0.bym.TX, family = "poisson", data = merged_proj_TX, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.0.TX.bym)
sum(mod.0.TX.bym$cpo$cpo) 
hist(mod.0.TX.bym$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

### Model 1: CVH score as a covariate
mod.1.TX.bym <- inla(formula.1.bym.TX, family = "poisson", data = merged_proj_TX, E = ExpectedDeaths,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(return.marginals.predictor = TRUE, 
                                         dic = TRUE, cpo = TRUE, waic = TRUE))
summary(mod.1.TX.bym)
sum(mod.1.TX.bym$cpo$cpo) 
hist(mod.1.TX.bym$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)
tau.bym.0.GA <- mod.0.GA.bym$summary.hyperpar$mean[1]
phi.bym.0.GA <- mod.0.GA.bym$summary.hyperpar$mean[2]
tau.bym.1.GA <- mod.1.GA.bym$summary.hyperpar$mean[1]
phi.bym.1.GA <- mod.1.GA.bym$summary.hyperpar$mean[2]
tau.bym.0.CA <- mod.0.CA.bym$summary.hyperpar$mean[1]
phi.bym.0.CA <- mod.0.CA.bym$summary.hyperpar$mean[2]
tau.bym.1.CA <- mod.1.CA.bym$summary.hyperpar$mean[1]
phi.bym.1.CA <- mod.1.CA.bym$summary.hyperpar$mean[2]
tau.bym.0.TX <- mod.0.TX.bym$summary.hyperpar$mean[1]
phi.bym.0.TX <- mod.0.TX.bym$summary.hyperpar$mean[2]
tau.bym.1.TX <- mod.1.TX.bym$summary.hyperpar$mean[1]
phi.bym.1.TX <- mod.1.TX.bym$summary.hyperpar$mean[2]

tau.m0 <- round(c(tau.bym.0.CA, tau.bym.0.GA, tau.bym.0.TX), 2)
phi.m0 <- round(c(phi.bym.0.CA, phi.bym.0.GA, phi.bym.0.TX), 2)
marg.var.m0 <- 1/tau.m0
spat.var.m0 <- marg.var.m0*phi.m0


tau.m1 <- round(c(tau.bym.1.CA, tau.bym.1.GA, tau.bym.1.TX), 2)
phi.m1 <- round(c(phi.bym.1.CA, phi.bym.1.GA, phi.bym.1.TX), 2)
marg.var.m1 <- 1/tau.m1
spat.var.m1 <- marg.var.m1*phi.m1

tau.change <- round((tau.m1 - tau.m0)/tau.m0, 3)*100
phi.change <- round((phi.m1 - phi.m0)/phi.m0, 3)*100
marg.var.change <- round((marg.var.m1 - marg.var.m0)/marg.var.m0, 3)*100
spat.var.change <- round((spat.var.m1 - spat.var.m0)/spat.var.m0, 3)*100

Var.res <- cbind("Tau (M0)" = tau.m0,
                  "Tau (M1)" = tau.m1,
                  "%-Change in Tau" = tau.change,
                  "Phi (M0)" = phi.m0,
                  "Phi (M1)" = phi.m1,
                  "%-Change in Phi" =  phi.change,
                  "1/Tau (M0)" = round(marg.var.m0,3),
                  "1/Tau (M1)" = round(marg.var.m1,3),
                  "%-Change in 1/Tau" = marg.var.change,
                  "Spatial Var (M0)" = round(spat.var.m0,3),
                  "Spatial Var (M1)" = round(spat.var.m1,3),
                  "%-Change in Spatial Var" =  spat.var.change
                 )
rownames(Var.res) <- c("CA", "GA", "TX" )
write.csv(Var.res, "Output/Results/AllCause_Mortality_variance_components.csv", row.names = TRUE)

TX.mod.results = list("GLM 0 -TX" = mod.0.TX.glm, "GLM 1 -TX" = mod.1.TX.glm,
                     "ICAR 0 -TX" = mod.0.TX.icar, "ICAR 1 -TX" = mod.1.TX.icar,
                     "BYM 0 -TX" = mod.0.TX.bym, "BYM 1 -TX" = mod.1.TX.bym)

TX.out = t(sapply(TX.mod.results, extract.model.results))
TX.out

TX.fixed = do.call(rbind, lapply(TX.mod.results, extract.fixed.effects))

write.csv(rbind(CA.out %>% as.data.frame() %>% arrange(DIC),
                GA.out %>% as.data.frame() %>% arrange(DIC),
                TX.out %>% as.data.frame() %>% arrange(DIC)),
          "Output/Results/AllCause_Mortality_NB-model_results.csv",
          row.names = TRUE)

write.csv(round(rbind(CA.fixed, GA.fixed, TX.fixed),5),
          "Output/Results/AllCause_Mortality_NB-model_fixed_effects.csv",
          row.names = TRUE)


merged_proj_CA$RE_bym0 <- mod.0.CA.bym$summary.random$county_int$mean[1:58]
merged_proj_GA$RE_bym0 <- mod.0.GA.bym$summary.random$county_int$mean[1:159]
merged_proj_TX$RE_bym0 <- mod.0.TX.bym$summary.random$county_int$mean[1:251]
merged_proj_CA$RE_bym1 <- mod.1.CA.bym$summary.random$county_int$mean[1:58]
merged_proj_GA$RE_bym1 <- mod.1.GA.bym$summary.random$county_int$mean[1:159]
merged_proj_TX$RE_bym1 <- mod.1.TX.bym$summary.random$county_int$mean[1:251]

mod.1.CA.bym$summary.hyperpar

(var(merged_proj_CA$RE_bym1) - var(merged_proj_CA$RE_bym0))/var(merged_proj_CA$RE_bym0)
# = -0.1345998
(var(merged_proj_GA$RE_bym1) - var(merged_proj_GA$RE_bym0))/var(merged_proj_GA$RE_bym0)
# = -0.4458711
(var(merged_proj_TX$RE_bym1) - var(merged_proj_TX$RE_bym0))/var(merged_proj_TX$RE_bym0)
# = 0.05285844

##################      Maps of Random Effects      ##################
#####   NULL Models    #####
merged_proj_CA$RE_bym0 <- mod.0.CA.bym$summary.random$county_int$mean[1:58]
merged_proj_GA$RE_bym0 <- mod.0.GA.bym$summary.random$county_int$mean[1:159]
merged_proj_TX$RE_bym0 <- mod.0.TX.bym$summary.random$county_int$mean[1:251]

# merged_proj_all$RE_bym0 <- 
# names(mod.0.CA.bym$summary.random$county_int$ID)

mod.0.CA.map <- ggplot() +
  geom_sf(data = merged_proj_CA, aes(fill = RE_bym0), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette[6:1], name = "Random Effect", na.value = "grey90", 
                        limits = c(-0.4, 0.4)) +
  annotation_scale(location = "bl", width_hint = 0.3) #+
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering())
# ggsave("Output/Figures/CA_BYM_AllCauseMort_0.jpeg", plot = mod.0.CA.map, width = 4, height = 5, units = "in", dpi = 300)

mod.0.TX.map <- ggplot() +
  geom_sf(data = merged_proj_TX, aes(fill = RE_bym0), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette[6:1], name = "Random Effect", na.value = "grey90") +
  annotation_scale(location = "bl", width_hint = 0.3) #+
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering()) 
# ggsave("Output/Figures/TX_BYM_AllCauseMort_0.jpeg", plot = mod.0.TX.map, width = 4, height = 5, units = "in", dpi = 300)

mod.0.GA.map <- ggplot() +
  geom_sf(data = merged_proj_GA, aes(fill = RE_bym0), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette[6:1], name = "Random Effect", na.value = "grey90", 
                        limits = c(-0.075, 0.05)) +
  annotation_scale(location = "bl", width_hint = 0.3) #+
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering()) 
# ggsave("Output/Figures/GA_BYM_AllCauseMort_0.jpeg", plot = mod.0.GA.map, width = 4, height = 5, units = "in", dpi = 300)

combined.0 <- grid.arrange(mod.0.CA.map, mod.0.TX.map, mod.0.GA.map,  ncol = 2) 
ggsave("Output/Figures/Combined_BYM_AllCauseMort_0.jpeg", plot = combined.0, width = 6.67, height = 5, units = "in", dpi = 300)

#####   CVH Models    #####
merged_proj_CA$RE_bym1 <- mod.1.CA.bym$summary.random$county_int$mean[1:58]
merged_proj_GA$RE_bym1 <- mod.1.GA.bym$summary.random$county_int$mean[1:159]
merged_proj_TX$RE_bym1 <- mod.1.TX.bym$summary.random$county_int$mean[1:251]

# merged_proj_all$RE_bym1 <- 
# names(mod.1.CA.bym$summary.random$county_int$ID)

mod.1.CA.map <- ggplot() +
  geom_sf(data = merged_proj_CA, aes(fill = RE_bym1), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette[6:1], name = "Random Effect", na.value = "grey90", 
                        limits = c(-0.4, 0.4)) +
  annotation_scale(location = "bl", width_hint = 0.3) #+
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering())
# ggsave("Output/Figures/CA_BYM_AllCauseMort_1.jpeg", plot = mod.1.CA.map, width = 4, height = 5, units = "in", dpi = 300)

mod.1.TX.map <- ggplot() +
  geom_sf(data = merged_proj_TX, aes(fill = RE_bym1), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette[6:1], name = "Random Effect", na.value = "grey90") +
  annotation_scale(location = "bl", width_hint = 0.3) #+
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering()) 
# ggsave("Output/Figures/TX_BYM_AllCauseMort_1.jpeg", plot = mod.1.TX.map, width = 4, height = 5, units = "in", dpi = 300)

mod.1.GA.map <- ggplot() +
  geom_sf(data = merged_proj_GA, aes(fill = RE_bym1), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette[6:1], name = "Random Effect", na.value = "grey90", 
                        limits = c(-0.075, 0.05)) +
  annotation_scale(location = "bl", width_hint = 0.3) #+
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering()) 
# ggsave("Output/Figures/GA_BYM_AllCauseMort_1.jpeg", plot = mod.1.GA.map, width = 4, height = 5, units = "in", dpi = 300)

combined.1<-grid.arrange(mod.1.CA.map, mod.1.TX.map, mod.1.GA.map,  ncol = 2) 
ggsave("Output/Figures/Combined_BYM_AllCauseMort_1.jpeg", plot = combined.1, width = 6.67, height = 5, units = "in", dpi = 300)

#####   Change from null to CVH models    #####
merged_proj_CA$RE_change <- (merged_proj_CA$RE_bym1 - merged_proj_CA$RE_bym0)/abs(merged_proj_CA$RE_bym0)
merged_proj_GA$RE_change <- (merged_proj_GA$RE_bym1 - merged_proj_GA$RE_bym0)/abs(merged_proj_GA$RE_bym0)
merged_proj_TX$RE_change <- (merged_proj_TX$RE_bym1 - merged_proj_TX$RE_bym0)/abs(merged_proj_TX$RE_bym0)
summary(merged_proj_CA$RE_change)
hist(merged_proj_CA$RE_change, breaks = 30)
summary(merged_proj_GA$RE_change)
summary(merged_proj_TX$RE_change)

var(merged_proj_CA$RE_bym1) - var(merged_proj_CA$RE_bym0)
var(merged_proj_GA$RE_bym1) - var(merged_proj_GA$RE_bym0)
var(merged_proj_TX$RE_bym1) - var(merged_proj_TX$RE_bym0)

merged_proj_CA$RE_abs_magchange <- (abs(merged_proj_CA$RE_bym1) - abs(merged_proj_CA$RE_bym0))/abs(merged_proj_CA$RE_bym0)
merged_proj_GA$RE_abs_magchange <- (abs(merged_proj_GA$RE_bym1) - abs(merged_proj_GA$RE_bym0))/abs(merged_proj_GA$RE_bym0)
merged_proj_TX$RE_abs_magchange <- (abs(merged_proj_TX$RE_bym1) - abs(merged_proj_TX$RE_bym0))/abs(merged_proj_TX$RE_bym0)
summary(merged_proj_CA$RE_abs_magchange)
summary(merged_proj_GA$RE_abs_magchange)
summary(merged_proj_TX$RE_abs_magchange)

# merged_proj_all$RE_change <- 
# names(mod.1.CA.bym$summary.random$county_int$ID)

mod.1.CA.map <- ggplot() +
  geom_sf(data = merged_proj_CA, aes(fill = RE_change), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette[6:1], name = "%-Change in Random Effect", na.value = "grey90", 
                        limits = c(-70, 70)) #+
  # annotation_scale(location = "bl", width_hint = 0.3) #+
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering())
# ggsave("Output/Figures/CA_BYM_AllCauseMort_Change.jpeg", plot = mod.1.CA.map, width = 4, height = 5, units = "in", dpi = 300)

mod.1.TX.map <- ggplot() +
  geom_sf(data = merged_proj_TX, aes(fill = RE_change), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette[6:1], name = "%-Change in Random Effect", na.value = "grey90", 
                        limits = c(-70, 70)) #+
  # annotation_scale(location = "bl", width_hint = 0.3) #+
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering()) 
# ggsave("Output/Figures/TX_BYM_AllCauseMort_Change.jpeg", plot = mod.1.TX.map, width = 4, height = 5, units = "in", dpi = 300)

mod.1.GA.map <- ggplot() +
  geom_sf(data = merged_proj_GA, aes(fill = RE_change), color = "white", size = 0.1) +
 #  geom_sf(data = states_sf, fill = NA, color = "grey30", size = 0.3) +
  scale_fill_gradientn(colors = cb_palette[6:1], name = "%-Change in Random Effect", na.value = "grey90", 
                        limits = c(-70, 70)) #+
  # annotation_scale(location = "bl", width_hint = 0.3) #+
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering()) 
# ggsave("Output/Figures/GA_BYM_AllCauseMort_Change.jpeg", plot = mod.1.GA.map, width = 4, height = 5, units = "in", dpi = 300)

combined.1<-grid.arrange(mod.1.CA.map + theme(legend.title = element_text(size = 8), axis.text = element_text(size = 8, angle = 45)),
                         mod.1.TX.map + theme(legend.title = element_text(size = 8), axis.text = element_text(size = 8, angle = 45)),
                         mod.1.GA.map + theme(legend.title = element_text(size = 8), axis.text = element_text(size = 8, angle = 45)),  ncol = 2) 

ggsave("Output/Figures/Combined_RE_AllCauseMort_Change_Maps.jpeg", plot = combined.1, width = 6.67, height = 5, units = "in", dpi = 300)



# ################          AllCause Mortality - Model 0          ################
# ##########   CA   ##########    
# ##### Exchangeable random effect model for AllCause mortality in CA
# formula.0.exc.cvh <- AllCause_mortality ~ 1 +
#   f(county_int, model = "iid", graph = g_CA)
# formula.1.exc.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "iid", graph = g_CA)
# formula.2.exc.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "iid", graph = g_CA) + f(county_slope, scale(cle8_point), model = "iid", graph = g_CA)

# ### Model 0: Exchangeable random effect model
# mod.0.CA.exc <- inla(formula.0.exc.cvh, family = "gaussian", data = merged_proj_CA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.0.CA.exc)
# # Marginal likelihood = 2.75e+10
# # DIC = -1585.51 | Saturated DIC = 115.91 | Eff. # Parms = 57.95
# sum(mod.0.CA.exc$cpo$cpo) #=22027236
# hist(mod.0.CA.exc$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 1: Exchangeable random effect model
# mod.1.CA.exc <- inla(formula.1.exc.cvh, family = "gaussian", data = merged_proj_CA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.1.CA.exc)
# # Marginal likelihood = -357.94
# # DIC = 672.23 | Saturated DIC = 115.91 | Eff. # Parms = 57.95
# sum(mod.1.CA.exc$cpo$cpo) #=22027236
# hist(mod.1.CA.exc$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 2: Exchangeable random effect model
# mod.2.CA.exc <- inla(formula.2.exc.cvh, family = "gaussian", data = merged_proj_CA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.2.CA.exc)
# # Marginal likelihood = -353.11
# # DIC = 671.53 | Saturated DIC = -47.32 | Eff. # Parms = 2.54
# sum(mod.2.CA.exc$cpo$cpo) #=22027236
# hist(mod.2.CA.exc$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ##### ICAR random effect model for AllCause mortality in CA
# formula.0.besag.cvh <- AllCause_mortality ~ 1 +
#   f(county_int, model = "besag", graph = g_CA)
# formula.1.besag.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "besag", graph = g_CA)
# formula.2.besag.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "besag", graph = g_CA) + f(county_slope, scale(cle8_point), model = "iid", graph = g_CA)

# ### Model 0: Besag-ICAR
# mod.0.CA.besag <- inla(formula.0.besag.cvh, family = "gaussian", data = merged_proj_CA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.0.CA.besag)
# # Marginal likelihood = -387.55
# # DIC = 670.88 | Saturated DIC = -50.05 | Eff. # Parms = 1.89
# sum(mod.0.CA.besag$cpo$cpo) #=0.21
# hist(mod.0.CA.besag$cpo$pit) # closer to uniform 

# ### Model 1: Besag-ICAR
# mod.1.CA.besag <- inla(formula.1.besag.cvh, family = "gaussian", data = merged_proj_CA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.1.CA.besag)
# # Marginal likelihood = -389.49
# # DIC = 672.00 | Saturated DIC = -47.32 | Eff. # Parms = 2.78
# sum(mod.1.CA.besag$cpo$cpo) #=0.21
# hist(mod.1.CA.besag$cpo$pit) # closer to uniform 

# ### Model 2: Besag-ICAR
# mod.2.CA.besag <- inla(formula.2.besag.cvh, family = "gaussian", data = merged_proj_CA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.2.CA.besag)
# # Marginal likelihood = -390.97
# # DIC = 672.00 | Saturated DIC = -47.32 | Eff. # Parms = 2.78
# sum(mod.2.CA.besag$cpo$cpo) #=0.21
# hist(mod.2.CA.besag$cpo$pit) # closer to uniform 

# ##### Besag-York-Mollie random effect model for AllCause mortality in CA
# formula.0.bym.cvh <- AllCause_mortality ~ 1 +
#   f(county_int, model = "bym2", graph = g_CA, scale.model = TRUE)
# formula.1.bym.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "bym2", graph = g_CA, scale.model = TRUE)
# formula.2.bym.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "bym2", graph = g_CA, scale.model = TRUE) +
#   f(county_slope, scale(cle8_point), model = "iid", graph = g_CA)

# ### Model 0: BYM2
# mod.0.CA.bym <- inla(formula.0.bym.cvh, family = "gaussian", data = merged_proj_CA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.0.CA.bym)
# # Marginal likelihood = -323.89
# # DIC = 670.19 | Saturated DIC = -48.48 | Eff. # Parms = 1.48
# sum(mod.0.CA.bym$cpo$cpo) #=0.21
# hist(mod.0.CA.bym$cpo$pit) # closer to uniform distrution, but still not a good fit

# ### Model 1: BYM2
# mod.1.CA.bym <- inla(formula.1.bym.cvh, family = "gaussian", data = merged_proj_CA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.1.CA.bym)
# # Marginal likelihood = -326.23
# # DIC = 671.41 | Saturated DIC = -47.84 | Eff. # Parms = 2.48
# sum(mod.1.CA.bym$cpo$cpo) #=0.21
# hist(mod.1.CA.bym$cpo$pit) # closer to uniform distrution, but still not a good fit

# ### Model 2: BYM2
# mod.2.CA.bym <- inla(formula.2.bym.cvh, family = "gaussian", data = merged_proj_CA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.2.CA.bym)
# # Marginal likelihood = -324.36
# # DIC = 671.19 | Saturated DIC = -50.48 | Eff. # Parms = 2.40
# sum(mod.2.CA.bym$cpo$cpo) #=0.21
# hist(mod.2.CA.bym$cpo$pit) # closer to uniform distrution, but still not a good fit
# CA.mod.results = list(mod.0.CA.exc, mod.1.CA.exc, mod.2.CA.exc,
#                      mod.0.CA.besag, mod.1.CA.besag, mod.2.CA.besag,
#                      mod.0.CA.bym, mod.1.CA.bym, mod.2.CA.bym)

# CA.out = t(sapply(CA.mod.results, function(x){c(x$dic$dic, x$dic$p.eff,x$waic$waic, x$waic$p.eff, x$mlik[2,1])}))
# colnames(CA.out) = c("DIC", "pD", "WAIC", "pWAIC", "Log Marginal Likelihood")
# rownames(CA.out) = c("Model 0: Exchangeable", "Model 1: Exchangeable", "Model 2: Exchangeable",
#                      "Model 0: Besag-ICAR", "Model 1: Besag-ICAR", "Model 2: Besag-ICAR",
#                      "Model 0: BYM2", "Model 1: BYM2", "Model 2: BYM2")
# CA.out

# ##########   TX   ##########    
# formula.0.exc.cvh <- AllCause_mortality ~ 1 +
#   f(county_int, model = "iid", graph = g_TX)
# formula.1.exc.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "iid", graph = g_TX)
# formula.2.exc.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "iid", graph = g_TX) + f(county_slope, scale(cle8_point), model = "iid", graph = g_TX)

# ##### Exchangeable random effect model for AllCause mortality in TX
# ### Model 0: Exchangeable random effect model
# mod.0.TX.exc <- inla(formula.0.exc.cvh, family = "gaussian", data = merged_proj_TX,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.0.TX.exc)
# # Marginal likelihood = -1558.89
# # DIC = 3082.03 | Saturated DIC = 115.91 | Eff. # Parms = 1.86
# sum(mod.0.TX.exc$cpo$cpo) #=0.55
# hist(mod.0.TX.exc$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 1: Exchangeable random effect model
# mod.1.TX.exc <- inla(formula.1.exc.cvh, family = "gaussian", data = merged_proj_TX,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.1.TX.exc)
# # Marginal likelihood = -1560.81
# # DIC = 3083.48 | Saturated DIC = -616.76 | Eff. # Parms = 2.80
# sum(mod.1.TX.exc$cpo$cpo) #=22027236
# hist(mod.1.TX.exc$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 2: Exchangeable random effect model
# mod.2.TX.exc <- inla(formula.2.exc.cvh, family = "gaussian", data = merged_proj_TX,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.2.TX.exc)
# # Marginal likelihood = -1565.46
# # DIC = -1585.51 | Saturated DIC = 115.91 | Eff. # Parms = 57.95
# sum(mod.2.TX.exc$cpo$cpo) #=22027236
# hist(mod.2.TX.exc$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ##### Besag-ICAR random effect model for AllCause mortality in TX
# formula.0.besag.cvh <- AllCause_mortality ~ 1 +
#   f(county_int, model = "besag", graph = g_TX)
# formula.1.besag.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "besag", graph = g_TX)
# formula.2.besag.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "besag", graph = g_TX) + f(county_slope, scale(cle8_point), model = "besag", graph = g_TX)

# ### Model 0: Besag-ICAR
# mod.0.TX.besag <- inla(formula.0.besag.cvh, family = "gaussian", data = merged_proj_TX,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.0.TX.besag)
# # Marginal likelihood = ______
# # DIC = ______ | Saturated DIC = ______ | Eff. # Parms = ______
# sum(mod.0.TX.besag$cpo$cpo) #=0.21
# hist(mod.0.TX.besag$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 1: Besag-ICAR
# mod.1.TX.besag <- inla(formula.1.besag.cvh, family = "gaussian", data = merged_proj_TX,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.1.TX.besag)
# # Marginal likelihood = ______
# # DIC = ______ | Saturated DIC = ______ | Eff. # Parms = ______
# sum(mod.1.TX.besag$cpo$cpo) #=0.21
# hist(mod.1.TX.besag$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 2: Besag-ICAR
# mod.2.TX.besag <- inla(formula.2.besag.cvh, family = "gaussian", data = merged_proj_TX,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.2.TX.besag)
# # Marginal likelihood = ______
# # DIC = ______ | Saturated DIC = ______ | Eff. # Parms = ______
# sum(mod.2.TX.besag$cpo$cpo) #=0.21
# hist(mod.2.TX.besag$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ##### Mesag-York_Mollie random effect model for AllCause mortality in TX
# formula.0.bym.cvh <- AllCause_mortality ~ 1 +
#   f(county_int, model = "bym2", graph = g_TX, scale.model = TRUE)
# formula.1.bym.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "bym2", graph = g_TX, scale.model = TRUE)
# formula.2.bym.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "bym2", graph = g_TX, scale.model = TRUE) + f(county_slope, scale(cle8_point), model = "bym2", graph = g_TX, scale.model = TRUE)

# ### Model 0: BYM2
# mod.0.TX.bym <- inla(formula.0.bym.cvh, family = "gaussian", data = merged_proj_TX,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.0.TX.bym)
# # Marginal likelihood = -323.89
# # DIC = 670.19 | Saturated DIC = -48.48 | Eff. # Parms = 1.48
# sum(mod.0.TX.bym$cpo$cpo) #=0.21
# hist(mod.0.TX.bym$cpo$pit) # closer to uniform distrution, but still not a good fit

# ### Model 1: BYM2
# mod.1.TX.bym <- inla(formula.1.bym.cvh, family = "gaussian", data = merged_proj_TX,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.1.TX.bym)
# # Marginal likelihood = -326.23
# # DIC = 671.41 | Saturated DIC = -47.84 | Eff. # Parms = 2.48
# sum(mod.1.TX.bym$cpo$cpo) #=0.21
# hist(mod.1.TX.bym$cpo$pit) # closer to uniform distrution, but still not a good fit

# ### Model 2: BYM2
# mod.2.TX.bym <- inla(formula.2.bym.cvh, family = "gaussian", data = merged_proj_TX,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.2.TX.bym)
# # Marginal likelihood = -323.89
# # DIC = 670.19 | Saturated DIC = -48.48 | Eff. # Parms = 1.48
# sum(mod.2.TX.bym$cpo$cpo) #=0.21
# hist(mod.2.TX.bym$cpo$pit) # closer to uniform distrution, but still not a good fit

# mod.2.TX.bym = list(NA)
# TX.mod.results = list(mod.0.TX.exc, mod.1.TX.exc, 
#                      mod.0.TX.besag, mod.1.TX.besag, 
#                      mod.0.TX.bym, mod.1.TX.bym)

# TX.out = t(sapply(TX.mod.results, function(x){c(x$dic$dic, x$dic$p.eff,x$waic$waic, x$waic$p.eff, x$mlik[2,1])}))
# colnames(TX.out) = c("DIC", "pD", "WAIC", "pWAIC", "Log Marginal Likelihood")
# rownames(TX.out) = c("Model 0: Exchangeable", "Model 1: Exchangeable",
#                      "Model 0: Besag-ICAR", "Model 1: Besag-ICAR",
#                      "Model 0: BYM2", "Model 1: BYM2")

# ##########   GA   ##########    
# formula.0.exc.cvh <- AllCause_mortality ~ 1 +
#   f(county_int, model = "iid", graph = g_GA)
# formula.1.exc.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "iid", graph = g_GA)
# formula.2.exc.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "iid", graph = g_GA) + f(county_slope, scale(cle8_point), model = "iid", graph = g_GA)

# ##### Exchangeable random effect model for AllCause mortality in GA
# ### Model 0: Exchangeable random effect model
# mod.0.GA.exc <- inla(formula.0.exc.cvh, family = "gaussian", data = merged_proj_GA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.0.GA.exc)
# # Marginal likelihood = -1558.89
# # DIC = 3082.03 | Saturated DIC = 115.91 | Eff. # Parms = 1.86
# sum(mod.0.GA.exc$cpo$cpo) #=0.55
# hist(mod.0.GA.exc$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 1: Exchangeable random effect model
# mod.1.GA.exc <- inla(formula.1.exc.cvh, family = "gaussian", data = merged_proj_GA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.1.GA.exc)
# # Marginal likelihood = -1560.81
# # DIC = 3083.48 | Saturated DIC = -616.76 | Eff. # Parms = 2.80
# sum(mod.1.GA.exc$cpo$cpo) #=22027236
# hist(mod.1.GA.exc$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 2: Exchangeable random effect model
# mod.2.GA.exc <- inla(formula.2.exc.cvh, family = "gaussian", data = merged_proj_GA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.2.GA.exc)
# # Marginal likelihood = -1565.46
# # DIC = -1585.51 | Saturated DIC = 115.91 | Eff. # Parms = 57.95
# sum(mod.2.GA.exc$cpo$cpo) #=22027236
# hist(mod.2.GA.exc$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ##### Besag-ICAR random effect model for AllCause mortality in GA
# formula.0.besag.cvh <- AllCause_mortality ~ 1 +
#   f(county_int, model = "besag", graph = g_GA)
# formula.1.besag.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "besag", graph = g_GA)
# formula.2.besag.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "besag", graph = g_GA) + f(county_slope, scale(cle8_point), model = "besag", graph = g_GA)

# ### Model 0: Besag-ICAR
# mod.0.GA.besag <- inla(formula.0.besag.cvh, family = "gaussian", data = merged_proj_GA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.0.GA.besag)
# # Marginal likelihood = ______
# # DIC = ______ | Saturated DIC = ______ | Eff. # Parms = ______
# sum(mod.0.GA.besag$cpo$cpo) #=0.21
# hist(mod.0.GA.besag$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 1: Besag-ICAR
# mod.1.GA.besag <- inla(formula.1.besag.cvh, family = "gaussian", data = merged_proj_GA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.1.GA.besag)
# # Marginal likelihood = ______
# # DIC = ______ | Saturated DIC = ______ | Eff. # Parms = ______
# sum(mod.1.GA.besag$cpo$cpo) #=0.21
# hist(mod.1.GA.besag$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ### Model 2: Besag-ICAR
# mod.2.GA.besag <- inla(formula.2.besag.cvh, family = "gaussian", data = merged_proj_GA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.2.GA.besag)
# # Marginal likelihood = ______
# # DIC = ______ | Saturated DIC = ______ | Eff. # Parms = ______
# sum(mod.2.GA.besag$cpo$cpo) #=0.21
# hist(mod.2.GA.besag$cpo$pit) # data do not fit well (https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf)

# ##### Mesag-York_Mollie random effect model for AllCause mortality in GA
# formula.0.bym.cvh <- AllCause_mortality ~ 1 +
#   f(county_int, model = "bym2", graph = g_GA, scale.model = TRUE)
# formula.1.bym.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "bym2", graph = g_GA, scale.model = TRUE)
# formula.2.bym.cvh <- AllCause_mortality ~ 1 + scale(cle8_point) +
#   f(county_int, model = "bym2", graph = g_GA, scale.model = TRUE) + f(county_slope, scale(cle8_point), model = "bym2", graph = g_GA, scale.model = TRUE)

# ### Model 0: BYM2
# mod.0.GA.bym <- inla(formula.0.bym.cvh, family = "gaussian", data = merged_proj_GA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.0.GA.bym)
# # Marginal likelihood = -323.89
# # DIC = 670.19 | Saturated DIC = -48.48 | Eff. # Parms = 1.48
# sum(mod.0.GA.bym$cpo$cpo) #=0.21
# hist(mod.0.GA.bym$cpo$pit) # closer to uniform distrution, but still not a good fit

# ### Model 1: BYM2
# mod.1.GA.bym <- inla(formula.1.bym.cvh, family = "gaussian", data = merged_proj_GA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.1.GA.bym)
# # Marginal likelihood = -326.23
# # DIC = 671.41 | Saturated DIC = -47.84 | Eff. # Parms = 2.48
# sum(mod.1.GA.bym$cpo$cpo) #=0.21
# hist(mod.1.GA.bym$cpo$pit) # closer to uniform distrution, but still not a good fit

# ### Model 2: BYM2
# mod.2.GA.bym <- inla(formula.2.bym.cvh, family = "gaussian", data = merged_proj_GA,
#                   control.predictor = list(compute = TRUE),
#                   control.compute = list(return.marginals.predictor = TRUE, 
#                                          dic = TRUE, cpo = TRUE, waic = TRUE))
# summary(mod.2.GA.bym)
# # Marginal likelihood = -323.89
# # DIC = 670.19 | Saturated DIC = -48.48 | Eff. # Parms = 1.48
# sum(mod.2.GA.bym$cpo$cpo) #=0.21
# hist(mod.2.GA.bym$cpo$pit) # closer to uniform distrution, but still not a good fit

# GA.mod.results = list(mod.0.GA.exc, mod.1.GA.exc, mod.2.GA.exc,
#                      mod.0.GA.besag, mod.1.GA.besag, mod.2.GA.besag,
#                      mod.0.GA.bym, mod.1.GA.bym, mod.2.GA.bym)

# GA.out = t(sapply(GA.mod.results, function(x){c(x$dic$dic, x$dic$p.eff,x$waic$waic, x$waic$p.eff, x$mlik[2,1])}))
# colnames(GA.out) = c("DIC", "pD", "WAIC", "pWAIC", "Log Marginal Likelihood")
# rownames(GA.out) = c("Model 0: Exchangeable", "Model 1: Exchangeable", "Model 2: Exchangeable",
#                      "Model 0: Besag-ICAR", "Model 1: Besag-ICAR", "Model 2: Besag-ICAR",
#                      "Model 0: BYM2", "Model 1: BYM2", "Model 2: BYM2")

# GA.out
# write.csv(CA.out, "Output/Results/CA_model_AllCause_results.csv", row.names = TRUE)
# write.csv(TX.out, "Output/Results/TX_model_AllCause_results.csv", row.names = TRUE)  
# write.csv(GA.out, "Output/Results/GA_model_AllCause_results.csv", row.names = TRUE)  
