
# Clear R's memory
rm(list = ls())


library(geodata)
library(sf)
library(terra)
library(rvest)
library(stringr)
library(lubridate)
library(terra)
library(ggplot2)
library(interactions)
library(ggpubr)
library(png)
library(maps)
library(viridisLite)
library(vioplot)
library(tibble)
library(tidyverse)
library(jsonlite)
library(magick)
require(vegan)
require(ggimage) 



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

dataset %>% 
  filter(nSurvs < minSurveys) # All good!





# Compare survey method across urbanization gradient and Latitude


prop_dataset = left_join(prop_fullDataset, sites, by = 'Name') %>% 
  filter(nSurv >= minSurveys)


prop_dataset %>%
  ggplot(aes(x = ObservationMethod, y = dev, fill = ObservationMethod)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 3, color = "black") +
  geom_violin(trim = TRUE, alpha = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
    guides(fill = "none")+
  theme_minimal() +
  labs(x = "Observation Method", y = "% Development")

prop_dataset %>%
  ggplot(aes(x = ObservationMethod, y = Latitude, fill = ObservationMethod)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 3, color = "black") +
  geom_violin(trim = TRUE, alpha = 0.3) +
  scale_fill_manual(values = c("red", "yellow"))+
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  guides(fill = "none")+
  theme_minimal() +
  labs(x = "Observation Method", y = "Latitude",
      # subtitle = "Site observers preference for observation method dooes not differ with latitude"
       )







prop.cat.Dev = glm(caterpillar ~ dev + ObservationMethod, 
                   data = prop_dataset, weights = nSurv,  family = "binomial")

prop.spi.Dev = glm(spider ~ dev + ObservationMethod,
                   data = prop_dataset, weights = nSurv,  family = "binomial")

prop.bet.Dev= glm(beetle ~ dev + ObservationMethod, 
                  data = prop_dataset, weights = nSurv,  family = "binomial")

prop.bug.Dev = glm(truebug ~ dev + ObservationMethod, 
                   data = prop_dataset, weights = nSurv,  family = "binomial")


prop.hop.Dev = glm(hopper ~ dev + ObservationMethod, 
                   data = prop_dataset, weights = nSurv,  family = "binomial")


prop.ant.Dev = glm(ant ~ dev + ObservationMethod, 
                   data = prop_dataset, weights = nSurv,  family = "binomial")

prop.fly.Dev = glm(fly ~ dev + ObservationMethod, 
                   data = prop_dataset, weights = nSurv,  family = "binomial")

prop.grasshopper.Dev = glm(grasshopper ~ dev + ObservationMethod, 
                           data = prop_dataset, weights = nSurv,  family = "binomial")

prop.daddylonglegs.Dev = glm(daddylonglegs ~ dev + ObservationMethod, 
                             data = prop_dataset, weights = nSurv,  family = "binomial")



prop.devOnlyOutput = data.frame(rbind(summary(prop.cat.Dev)$coefficients, 
                                      summary(prop.spi.Dev)$coefficients, 
                                      summary(prop.bet.Dev)$coefficients,
                                      summary(prop.bug.Dev)$coefficients, 
                                      summary(prop.hop.Dev)$coefficients,
                                      summary(prop.ant.Dev)$coefficients,
                                      summary(prop.fly.Dev)$coefficients,
                                      summary(prop.grasshopper.Dev)$coefficients,
                                      summary(prop.daddylonglegs.Dev)$coefficients))
prop.devOnlyOutput$term = rep(c('Intercept', 'dev', 'Method'), times = 9)
prop.devOnlyOutput$Group = rep(c('Caterpillar', 'Spider', 'Beetle', 
                                 'Leafhopper', 'Truebugs', 'Ant',
                                 'Fly', 'Grasshopper', 'Daddylonglegs'), each = 3)

rownames(prop.devOnlyOutput) = NULL


# Interaction effect of dev and lat on arthropod ----

prop.cat.DevLat = glm(caterpillar ~ dev * Latitude  + ObservationMethod, 
                      data = prop_dataset, weights = nSurv,  family = "binomial")

prop.spi.DevLat = glm(spider ~ dev * Latitude + ObservationMethod,
                      data = prop_dataset, weights = nSurv,  family = "binomial")

prop.bet.DevLat= glm(beetle ~ dev * Latitude + ObservationMethod, 
                     data = prop_dataset, weights = nSurv,  family = "binomial")

prop.bug.DevLat = glm(truebug ~ dev * Latitude + ObservationMethod, 
                      data = prop_dataset, weights = nSurv,  family = "binomial")


prop.hop.DevLat = glm(hopper ~ dev * Latitude + ObservationMethod, 
                      data = prop_dataset, weights = nSurv,  family = "binomial")


prop.ant.DevLat = glm(ant ~ dev * Latitude + ObservationMethod, 
                      data = prop_dataset, weights = nSurv,  family = "binomial")

prop.fly.DevLat = glm(fly ~ dev * Latitude + ObservationMethod, 
                      data = prop_dataset, weights = nSurv,  family = "binomial")

prop.grasshopper.DevLat = glm(grasshopper ~ dev * Latitude + ObservationMethod, 
                              data = prop_dataset, weights = nSurv,  family = "binomial")

prop.daddylonglegs.DevLat = glm(daddylonglegs ~ dev * Latitude + ObservationMethod, 
                                data = prop_dataset, weights = nSurv,  family = "binomial")






prop.devLatOutput = data.frame(rbind(summary(prop.cat.DevLat)$coefficients, 
                                     summary(prop.spi.DevLat)$coefficients, 
                                     summary(prop.bet.DevLat)$coefficients,
                                     summary(prop.bug.DevLat)$coefficients, 
                                     summary(prop.hop.DevLat)$coefficients,
                                     summary(prop.ant.DevLat)$coefficients,
                                     summary(prop.fly.DevLat)$coefficients,
                                     summary(prop.grasshopper.DevLat)$coefficients,
                                     summary(prop.daddylonglegs.DevLat)$coefficients))
prop.devLatOutput$term = rep(c('Intercept', 'dev', 'Latitude', 'dev*Latitude', 'Method'),
                             times = 9)
prop.devLatOutput$Group = rep(c('Caterpillar', 'Spider', 'Beetle', 
                                'Leafhopper', 'Truebugs', 'Ant',
                                'Fly', 'Grasshopper', 'Daddylonglegs'), each = 5)

rownames(prop.devLatOutput) = NULL



### urbanization effect on occurrence by observationMethod

prop_dataset %>% 
  pivot_longer(cols = -c("Name", "ObservationMethod", "Region", 
                         "Longitude", "Latitude", "dev", "forest", "nSurv"),
               names_to = "Group",
               values_to = "Occurence") %>% 
  ggplot(aes(x = dev, y = Occurence)) +
  geom_point(aes(colour = ObservationMethod, size = nSurv, alpha = 0.2))+
  geom_smooth(
    method = "glm",
    method.args = list(family = binomial),
    aes(weight = nSurv, colour = ObservationMethod ),
    se = TRUE
  )+
  facet_wrap(~ Group, scales = "free_y") +
  labs(
    colour = "Observation method",
    size   = "Number of surveys",
    x      = "% Urban development",
    y      = "Proportion of cccurrence"
  ) +
  theme_bw()


 

# Plots

### Arthropod images ----
catImage = readPNG('images/caterpillar.png')
antImage = readPNG('images/ant.png')
beetleImage = readPNG('images/beetle.png')
spiderImage = readPNG('images/spider.png')
hopperImage = readPNG('images/leafhopper.png')
truebugImage = readPNG('images/truebugs.png')

# https://www.phylopic.org/images/7174527d-3060-4c78-ab59-c3ccf76075de/opilio-saxatilis
daddylonglegImage = image_read('images/daddylongleg.png') %>% 
  image_convert(type = "truecolor") %>% 
  as.raster()

# https://www.phylopic.org/images/2bec28d1-3d13-4512-8b4d-00cd0fcef6e4/clogmia-albipunctata
flyImage = image_read('images/fly.png') %>% 
  image_convert(type = "truecolor") %>% 
  as.raster()

# https://www.phylopic.org/images/0de35750-d9ba-472a-b4f2-3c630777fcc3/acrididae
grasshopperImage = image_read('images/grasshopper.png') %>% 
  image_convert(type = "truecolor") %>% 
  as.raster()

catDevPlot = interact_plot(prop.cat.DevLat, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with caterpillars",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(catImage, ymin = .085, ymax = .1, xmin = 60, xmax = 100)

betDevPlot = interact_plot(prop.bet.DevLat, pred = dev, modx = Latitude,
                            y.label = "Prop. of surveys with beetles",
                            x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                            colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  #ylim(0.12, 0.26) +
  annotation_raster(beetleImage, ymin = .28, ymax = .30, xmin = 60, xmax = 100) 

bugDevPlot = interact_plot(prop.bug.DevLat, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with true bugs",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(truebugImage, ymin = .09, ymax = .11, xmin = 0, xmax = 40)

spiDevPlot = interact_plot(prop.spi.DevLat, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with spiders",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(spiderImage, ymin = .29, ymax = .35, xmin = 60, xmax = 100)

hopDevPlot = interact_plot(prop.hop.DevLat, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with hoppers",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(hopperImage, ymin = .07, ymax = .095, xmin = 60, xmax = 100)

antDevPlot = interact_plot(prop.ant.DevLat, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with ants",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(antImage, ymin = .14, ymax = .16, xmin = 60, xmax = 99)


flyDevPlot = interact_plot(prop.fly.DevLat, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with flys",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(flyImage, ymin = .14, ymax = .16, xmin = 20, xmax = 60)



grasshoppersDevPlot = interact_plot(prop.grasshopper.DevLat, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with grasshoppers",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(grasshopperImage, ymin = .12, ymax = .16, xmin = 15, xmax = 60)


daddylonglegsDevPlot = interact_plot(prop.daddylonglegs.DevLat, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with daddylonglegs",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(daddylonglegImage, ymin = .04, ymax = .06, xmin = 60, xmax = 90)


ggarrange(catDevPlot, spiDevPlot, betDevPlot, bugDevPlot, hopDevPlot, antDevPlot,
          flyDevPlot, grasshoppersDevPlot, daddylonglegsDevPlot,
          ncol=3, nrow=3, common.legend = TRUE, legend="bottom")

## Supplementary analysis----

groups <- seq(min(dataset$Latitude), max(dataset$Latitude), length.out = 9)
mod.x = groups[-c(1,2, length(groups), length(groups)-1)]
 


dataset %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(nSurvs = sum(nSurvs)) %>% 
  filter(nSurvs < minSurveys)

# Developed cover models
cat.Dev.Latitude = glm(caterpillar ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")

spi.Dev.Latitude = glm(spider ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")

beet.Dev.Latitude = glm(beetle ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                        data = dataset, family = "binomial")

hop.Dev.Latitude = glm(hopper ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")

bug.Dev.Latitude = glm(truebug ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")

ant.Dev.Latitude = glm(ant ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")


fly.Dev.Latitude = glm(fly ~ dev * Latitude + ObservationMethod, 
                      data = dataset, family = "binomial")

grasshopper.Dev.Latitude = glm(grasshopper ~ dev * Latitude + ObservationMethod, 
                              data = dataset, family = "binomial")

daddylonglegs.Dev.Latitude = glm(daddylonglegs ~ dev * Latitude + ObservationMethod, 
                                data = dataset, family = "binomial")







simDev_caterpillar = sim_slopes(cat.Dev.Latitude, pred = dev, modx = Latitude,
                                 johnson_neyman = FALSE, digits = 4)

simDev_beetle = sim_slopes(beet.Dev.Latitude, pred = dev, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simDev_spider = sim_slopes(spi.Dev.Latitude, pred = dev, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simDev_hopper = sim_slopes(hop.Dev.Latitude, pred = dev, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simDev_ant = sim_slopes(ant.Dev.Latitude, pred = dev, modx = Latitude,
                        johnson_neyman = FALSE, digits = 4)

simDev_truebug = sim_slopes(bug.Dev.Latitude, pred = dev, modx = Latitude,
                            johnson_neyman = FALSE, digits = 4)

simDev_fly = sim_slopes(fly.Dev.Latitude, pred = dev, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simDev_grasshopper = sim_slopes(grasshopper.Dev.Latitude, pred = dev, modx = Latitude,
                        johnson_neyman = FALSE, digits = 4)

simDev_daddylonglegs = sim_slopes(daddylonglegs.Dev.Latitude, pred = dev, modx = Latitude,
                            johnson_neyman = FALSE, digits = 4)

Dev_sims =
  bind_rows(
    simDev_caterpillar$slopes     %>% mutate(Group = "caterpillar"),
    simDev_beetle$slopes          %>% mutate(Group = "beetle"),
    simDev_ant$slopes             %>% mutate(Group = "ant"),
    simDev_hopper$slopes          %>% mutate(Group = "hopper"),
    simDev_truebug$slopes         %>% mutate(Group = "truebug"),
    simDev_spider$slopes          %>% mutate(Group = "spider"),
    simDev_fly$slopes             %>% mutate(Group = "fly"),
    simDev_grasshopper$slopes     %>% mutate(Group = "grasshopper"),
    simDev_daddylonglegs$slopes   %>% mutate(Group = "daddylonglegs")
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "caterpillar",
        "beetle",
        "ant",
        "hopper",
        "truebug",
        "spider",
        "fly",
        "grasshopper",
        "daddylonglegs"
      )
    )
  )



 
par(mfrow = c(3, 3))  # 9 panels

johnson_neyman(cat.Dev.Latitude, pred = dev, modx = Latitude, control.fdr = TRUE)
johnson_neyman(spi.Dev.Latitude, pred = dev, modx = Latitude, control.fdr = TRUE)
johnson_neyman(beet.Dev.Latitude, pred = dev, modx = Latitude, control.fdr = TRUE)
johnson_neyman(hop.Dev.Latitude, pred = dev, modx = Latitude, control.fdr = TRUE)
johnson_neyman(bug.Dev.Latitude, pred = dev, modx = Latitude, control.fdr = TRUE)
johnson_neyman(ant.Dev.Latitude, pred = dev, modx = Latitude, control.fdr = TRUE)
johnson_neyman(fly.Dev.Latitude, pred = dev, modx = Latitude, control.fdr = TRUE)
johnson_neyman(grasshopper.Dev.Latitude, pred = dev, modx = Latitude, control.fdr = TRUE)
johnson_neyman(daddylonglegs.Dev.Latitude, pred = dev, modx = Latitude, control.fdr = TRUE)

par(mfrow = c(1, 1))  # reset

#######################################################################################################

# Community level analysis----


# Note that theoretically, the P-value from using lm() on dissimilarity computation, in this case, cannot be valid. This is because pair-wise correlations are not independent. That is what a mantel test corrects for by doing permutations.
########################################################################################################################

prop.num <- prop_dataset[,3:11]
site.info <- prop_dataset[,c(1,2,12:17)]

prop.num <- prop.num %>% as.data.frame()

prop.num_ilogit <- plogis(as.matrix(prop.num)) %>% as.data.frame
prop.num_ilogit


site.z <- scale(prop_dataset[,c( "Latitude", "dev", "forest")]) %>% 
  as.data.frame()

Observ <- prop_dataset [, c("ObservationMethod","Longitude")]


part.prop.rda.LD <-rda(
  prop.num_ilogit ~ Latitude * dev +  Condition(ObservationMethod),
  data = cbind(site.z, Observ))

summary(part.prop.rda.LD)
RsquareAdj(part.prop.rda.LD)$adj.r.squared
# Test whether the model is statistically significant


anova.cca(part.prop.rda.LD, step = 999) # good!
anova.cca(part.prop.rda.LD, step = 999, by = "axis")
anova.cca(part.prop.rda.LD, step = 999, by= 'terms')

anova.cca(part.prop.rda.LD, step = 999, by= 'margin')

vif.cca(part.prop.rda.LD) # Below 1.1


vpart.prop.LD <- varpart(prop.num_ilogit, site.z %>%select(-forest),  Observ)
vpart.prop.LD


species_score.ld <- scores(part.prop.rda.LD, 
                           display = "species",
                           choices = 1:2) %>% 
  as.data.frame()
site_score.ld <- scores(part.prop.rda.LD, 
                        display = "sites",
                        choices = 1:2) %>% 
  as.data.frame()

site_score.sites.ld  <- cbind(prop_dataset [,c("Name", "Region", "Longitude", 
                                               "Latitude", "dev", "forest")], site_score.ld)

rda_axis.ld <- scores(part.prop.rda.LD, display = "bp", choices = 1:3) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  filter(!grepl(":", term))  %>% # Remove all interaction terms
  mutate(
    term = recode(term,
                  dev = "Development")) 




taxa_scores <- part.prop.rda.LD$CCA$v %>%
  as.data.frame() %>%
  tibble::rownames_to_column("taxon")

taxa_scores <- taxa_scores %>%
  mutate(image = case_when(
    taxon == "caterpillar"   ~ "images/caterpillar.png",
    taxon == "spider"        ~ "images/spider.png",
    taxon == "beetle"        ~ "images/beetle.png",
    taxon == "truebug"       ~ "images/truebugs.png",
    taxon == "hopper"        ~ "images/leafhopper.png",
    taxon == "ant"           ~ "images/ant.png",
    taxon == "grasshopper"   ~ "images/grasshopper.png",
    taxon == "fly"           ~ "images/fly.png",
    taxon == "daddylonglegs" ~ "images/daddylongleg.png"
  ))


ggplot() +
  geom_segment(
    data = rda_axis.ld,
    aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
    arrow = arrow(length = unit(0.3, "cm")), 
    color = "black",
    size = 1, alpha = .5
  ) +
  geom_point(
    data = site_score.sites.ld, 
    aes(x = RDA1, y = RDA2, fill = dev),
    shape = 21, size = 4, alpha = 0.5
  ) +
  geom_image(
    data = taxa_scores,
    aes(x = RDA1, y = RDA2, image = image),
    size = 0.065   # adjust this for figure scale
  ) +
  geom_text(data = rda_axis.ld, 
            aes(x = RDA1, y = RDA2, 
                label = term), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("green", "lightyellow", "blue"),
    name = "Development"
  ) +
  labs(x = "RDA1", y = "RDA2") +
  guides(color = "none", shape = "none", alpha = "none",
         fill = guide_colorbar(title = "% Development")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal()





ggplot() +
  geom_segment(data = rda_axis.ld,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black",
               size = 1, alpha = .5) +
  geom_point(data = site_score.sites.ld, 
             aes(x = RDA1, y = RDA2, 
                 color = dev, fill = dev, alpha = 0.5), 
             shape = 21, size = 4) +
  
  geom_text(data = part.prop.rda.LD$CCA$v, 
            aes(x = RDA1, y = RDA2, label = rownames(part.prop.rda.LD$CCA$v)), 
            color = "black") +
  geom_text(data = rda_axis.ld, 
            aes(x = RDA1, y = RDA2, label = term), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("green", "lightyellow", "blue"),
    name = "Development"
  ) +
  labs(x = "RDA1", y = "RDA2") +
  guides(color = "none", shape = "none", alpha = "none",
         fill = guide_colorbar(title = "% Development")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_minimal()
