
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

fullDataset <- read.csv(paste0(github_raw, latest_file))




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
sites = distinct(fullDataset, Name, Region, Longitude, Latitude)

dataset = left_join(goodData, sites, by = 'Name')

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



# (4) Read in maps and Landcover data----
#     Sources: NLCD
#              Canada Land Cover 2020

# This raster is in geographic coverage, but at 3 km resolution
us_devo3km_geog = rast('data/devo_geographic_3km.tif')

# Rasters at 1 km resolution
canada_devo_1km = rast('data/devo_1km_toronto.tif')
us_devo_1km = rast('data/devo_1km_east.tif')
canada_forest_1km = rast('data/forest_1km_toronto.tif')
us_forest_1km = rast('data/forest_1km_east.tif')

us_devo_geog = project(us_devo_1km, crs(us_devo3km_geog))
can_devo_geog = project(canada_devo_1km, crs(us_devo_geog))
can_devo_albers = project(canada_devo_1km, crs(us_devo_1km))

can_forest_geog = project(canada_forest_1km, crs(us_devo_geog))
us_forest_geog = project(us_forest_1km, crs(us_devo_geog))



northam = st_read("data/na_base_Lambert_Azimuthal.shp")

us_devo_la = project(us_devo_1km, crs(northam))

us_devo_na <- clamp(us_devo_1km, lower=0.0001, value=FALSE)

northam_geom = st_cast(northam$geometry,"POLYGON")

northam_albers = st_transform(northam_geom, crs(us_devo_1km))

northam_albers_vect = vect(northam_albers)

# Crop out ocean
us_devo_crop <- crop(us_devo_1km, northam_albers_vect, mask= T)
can_devo_crop <- crop(can_devo_albers, northam_albers_vect, mask= T)

# Crop to similar bounding box
e <- ext(-1e+06, 2.5e+06, 170000, 3.1e+06)

us_devo_ext = extend(us_devo_crop, can_devo_crop, snap = "out") %>%
  crop(e, snap = "out")

can_devo_ext = extend(can_devo_crop, us_devo_crop, snap = "in") %>%
  crop(e, snap = "out")

# Still not working because even after projecting to same projection,
# resolution and extents end up differing

# nam_devo = us_devo_ext + can_devo_ext


# Mosaic not working at the moment, as projected rasters don't have the same resolution. Need to fix...
#devo1km = terra::mosaic(can_devo_geog, us_devo_geog)


sites = distinct(fullDataset, Name, Region, Longitude, Latitude)

## Add urban cover and forest covrer information to sites
CAN_sites = sites[sites$Region %in% c('AB', 'ON'),]
US_sites = sites[!sites$Region %in% c('AB', 'ON'),]

us_devo = terra::extract(us_devo_geog, US_sites[, c('Longitude', 'Latitude')])
can_devo = terra::extract(can_devo_geog, CAN_sites[, c('Longitude', 'Latitude')])

us_forest = terra::extract(us_forest_geog, US_sites[, c('Longitude', 'Latitude')])
can_forest = terra::extract(can_forest_geog, CAN_sites[, c('Longitude', 'Latitude')])


CAN_sites$dev = can_devo$Canada2020                     
CAN_sites$forest = can_forest$Canada2020
US_sites$dev = us_devo$`NLCD Land Cover Class`
US_sites$forest = us_forest$`NLCD Land Cover Class`

sites = rbind(US_sites, CAN_sites)

write.csv(sites, file = "data/sites.csv") # keep a local copy.

sites_sf <- st_as_sf(sites[sites$Name %in% goodSites$Name, ], 
                     coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs(us_devo_1km))

us_devo_df = as.data.frame(us_devo_crop, xy = TRUE)


# Doesn't currently include Canada landcover but good enough for now.
ggplot() +
  geom_raster(data = us_devo_df, aes(x = x, y = y, fill = `NLCD Land Cover Class`)) + 
  scale_fill_gradientn(colors = viridis(10)) + 
  geom_sf(data = sites_sf, col = "magenta", size = 2.5)


dataset = left_join(goodData, sites, by = 'Name')

head(dataset)


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
    n.chains = 4
  )
  
  update(model, 2000) # I think this burn-in takes to much computation time
  
  fit <- coda.samples(
    model,
    variable.names = c("beta_0","beta_dev",
                       "beta_lat","beta_int",
                       "beta_method"),
    n.iter = 5000
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
flyFit          = fit_arthropod_model("fly")
daddylonglegsFit = fit_arthropod_model("daddylonglegs")

all_fits = list(
  caterpillar = caterpillarFit,
  spider = spiderFit,
  beetle = beetleFit,
  truebug = truebugFit,
  hopper = hopperFit,
  ant = antFit,
  grasshopper = grasshopperFit,
  fly = flyFit,
  daddylonglegs = daddylonglegsFit
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
saveRDS(caterpillarFit, file = "caterpillarFit.rds")
saveRDS(spiderFit, file = "spiderFit.rds")
saveRDS(beetleFit, file = "beetleFit.rds")
saveRDS(truebugFit, file = "truebugFit.rds")
saveRDS(hopperFit, file = "hopperFit.rds")
saveRDS(antFit, file = "antFit.rds")
saveRDS(grasshopperFit, file = "grasshopperFit.rds")
saveRDS(flyFit, file = "flyFit.rds")
saveRDS(daddylonglegsFit, file = "daddylonglegsFit.rds")

