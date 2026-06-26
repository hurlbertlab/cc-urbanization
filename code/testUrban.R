
# Clear R's memory
rm(list = ls())

library(rjags)
library(coda)
library(tidyverse)
library(tidybayes)
library(viridis)
library(raster)
library(geodata)
library(sf)
library(terra)
library(rvest)
library(lubridate)
library(terra)
library(ggplot2)
library(ggpubr)
library(png)
library(maps)
library(viridisLite)
library(vioplot)
library(tibble)
library(jsonlite)
library(magick)
require(vegan)
require(ggimage) 
library(gt) # to create tables



# (1) Read in latest Caterpillars Count! raw dataset from the caterpillars-analysis-public repo----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))  %>% filter(Year <= 2025)


# urban500m = read.csv("data/urbanization_USA_Canada_500m.csv") %>% 
#   rename(dev = urban_percent) %>% dplyr::select(-ID, -X, -Longitude)


urban1km = read.csv("data/urbanization_USA_Canada_1km.csv") %>% 
  rename(dev = urban_percent) %>% dplyr::select(-ID, -X, -Longitude)

# (2) Filter dataset ----
# to sites east of 100W with a minimum of 50 branch surveys during June and July (juliandays 152-213)
minSurveys = 50
julianWindow = 152:213



# to sites east of 100W with a minimum of 50 branch surveys during June and July (juliandays 152-213)
minSurveys = 50
julianWindow = 152:213

goodSites = fullDataset %>%
  filter(julianday %in% julianWindow,
         Longitude > -100,
         WetLeaves == 0) %>%
  group_by(Name, ObservationMethod) %>%
  summarize(nSurvs = n_distinct(ID)) %>%
  filter(nSurvs >= minSurveys) 


goodData = goodSites %>% 
  left_join(fullDataset %>%
              filter(
                julianday %in% julianWindow,
                WetLeaves == 0,
                !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>%
              group_by(Name, ID, ObservationMethod) %>%
              summarize(caterpillar = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
                        spider = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
                        beetle = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
                        truebug = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
                        hopper = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
                        ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0),
                        grasshopper = ifelse(sum(Group == "grasshopper", na.rm = TRUE) > 0, 1, 0),
                        fly = ifelse(sum(Group == "fly", na.rm = TRUE) > 0, 1, 0),
                        daddylonglegs = ifelse(sum(Group == "daddylonglegs", na.rm = TRUE) > 0, 1, 0)),
            by = c("Name", "ObservationMethod")
  )
sites = distinct(fullDataset, Name, Region, Longitude, Latitude) %>% 
  inner_join(urban1km, by = c("Name", "Latitude"))
  

dataset = inner_join(goodData, sites, by = 'Name') %>% 
  filter(!is.na(ID))



dataset %>% 
  filter(nSurvs < minSurveys) # All good!






prop_fullDataset = fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>% 
  group_by(Name, ID, ObservationMethod) %>%
  summarize(caterpillar = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
            spider = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
            beetle = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
            truebug = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
            hopper = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
            ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0),
            grasshopper = ifelse(sum(Group == "grasshopper", na.rm = TRUE) > 0, 1, 0),
            fly = ifelse(sum(Group == "fly", na.rm = TRUE) > 0, 1, 0),
            daddylonglegs = ifelse(sum(Group == "daddylonglegs", na.rm = TRUE) > 0, 1, 0),
            nSurv = n_distinct(ID)) %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(caterpillar = mean(caterpillar),
            spider = mean(spider),
            beetle = mean(beetle),
            truebug = mean(truebug),
            hopper  = mean(hopper),
            ant = mean(ant),
            grasshopper = mean(grasshopper),
            fly = mean(fly),
            daddylonglegs = mean(daddylonglegs),
            nSurv = sum(nSurv))  




# Bayesian analysis----

dataset$dev_c = as.numeric(scale(dataset$dev))
dataset$lat_c  = as.numeric(scale(dataset$Latitude))
dataset$dev_lat = dataset$dev_c * dataset$lat_c
dataset$method = as.numeric(dataset$ObservationMethod == "Visual")




cat("
model {

  # --------------------
  # Likelihood
  # --------------------
 for (i in 1:n) {
  y[i] ~ dbern(p[i])
  logit(p[i]) <- 
    beta_0 +                # global intercept
    beta_dev * dev[i] +     # urban development slope
    beta_lat * lat[i] +     # Latitide slope
    beta_int * dev_lat[i] + # Interaction effect of dev and latitude
    beta_method * method[i] # Observation method ( visual or beeting sheet; visual is the reference)
}

  # --------------------
  # Priors
  # --------------------
  beta_0      ~ dnorm(0.0, 0.01)
  beta_dev    ~ dnorm(0.0, 0.01)
  beta_lat    ~ dnorm(0.0, 0.01)
  beta_int    ~ dnorm(0.0, 0.01)
  beta_method ~ dnorm(0.0, 0.01)

}
", file = "urbanDevLat.jags")


# this function takes each arthropod group as reponse_var and ensures that it is a numeric 
# variable. The JAGs data list is created and we use a specified burn in.
fit_arthropod_model <- function(response_var) {
  
  dataset$y <- as.numeric(dataset[[response_var]])
  
  jags_data <- list(
    y       = dataset$y,
    dev     = as.numeric(dataset$dev_c),
    lat     = as.numeric(dataset$lat_c),
    dev_lat = as.numeric(dataset$dev_lat),
    method  = as.numeric(dataset$method),
    n       = nrow(dataset)
  )
  
  model <- jags.model(
    file = "urbanDevLat.jags",
    data = jags_data,
    n.chains = 3
  )
  
  update(model, 1000) # I think this burn-in takes to much computation time
  
  fit <- coda.samples(
    model,
    variable.names = c("beta_0","beta_dev",
                       "beta_lat","beta_int",
                       "beta_method"),
    n.iter = 2000
  )
  
  return(fit)
}

caterpillarFit  = fit_arthropod_model("caterpillar")
spiderFit       = fit_arthropod_model("spider")
beetleFit       = fit_arthropod_model("beetle")
truebugFit      = fit_arthropod_model("truebug")
hopperFit       = fit_arthropod_model("hopper")
antFit          = fit_arthropod_model("ant")
grasshopperFit  = fit_arthropod_model("grasshopper")
# flyFit          = fit_arthropod_model("fly")
# daddylonglegsFit = fit_arthropod_model("daddylonglegs")

all_fits = list(
  caterpillar = caterpillarFit,
  spider = spiderFit,
  beetle = beetleFit,
  truebug = truebugFit,
  hopper = hopperFit,
  ant = antFit,
  grasshopper = grasshopperFit
  # fly = flyFit,
  # daddylonglegs = daddylonglegsFit
)

saveRDS(all_fits, "arthropod_fits_all.rds")


lapply(all_fits, gelman.diag)
lapply(all_fits, effectiveSize)


pdf("all_traceplots.pdf", width = 10, height = 8)

for (name in names(all_fits)) {
  traceplot(all_fits[[name]],
            main = paste("Traceplot:", name))
}

dev.off()


# Save each fit as an .rds file
saveRDS(caterpillarFit, file = "caterpillarFit1000.rds")
saveRDS(spiderFit, file = "spiderFit1000.rds")
saveRDS(beetleFit, file = "beetleFit1000.rds")
saveRDS(truebugFit, file = "truebugFit1000.rds")
saveRDS(hopperFit, file = "hopperFit1000.rds")
saveRDS(antFit, file = "antFit1000.rds")
saveRDS(grasshopperFit, file = "grasshopperFit1000.rds")
# saveRDS(flyFit, file = "flyFit500.rds") # so, fly and daddylongleg might have the old and wrong urbanization analysis saved in it.
# saveRDS(daddylonglegsFit, file = "daddylonglegsFit500.rds") # right now, I do not care about it, so to save time I will not run them.

