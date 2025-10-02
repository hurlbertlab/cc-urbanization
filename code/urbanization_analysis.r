library(dplyr)
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


# (1) Read in latest Caterpillars Count! raw dataset from the caterpillars-analysis-public repo
data_repo <- "https://github.com/hurlbertlab/caterpillars-analysis-public/tree/master/data"
webpage <- read_html(data_repo)
repo_links <- unique(html_attr(html_nodes(webpage, "a"), "href"))
dataset_url <- repo_links[grepl("fullDataset", repo_links)]

latest_file <- word(dataset_url, -1, sep = "/")

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/refs/heads/master/data/"

fullDataset = read.csv(paste0(github_raw, latest_file), header = TRUE, quote = '\"', fill = TRUE)


# (2) Filter dataset to sites east of 100W with a minimum of 40 branch surveys during June and July (juliandays 152-213)
minSurveys = 50
julianWindow = 152:213

goodSites = fullDataset %>%
  filter(julianday %in% julianWindow,
         Longitude > -100,
         WetLeaves == 0) %>%
  group_by(Name, ObservationMethod) %>%
  summarize(nSurvs = n_distinct(ID)) %>%
  filter(nSurvs >= minSurveys) %>%
  arrange(desc(nSurvs))

# Dataset thru 2024 includes 154 sites with at least 50 surveys in the seasonal window


# (3) Reorganize data to presence-absence of arthropod groups on each survey for these good sites
goodData = fullDataset %>%
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
            ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0))
    
goodDataBS = fullDataset %>%  # This would not have up to 50 in each set
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         ObservationMethod == 'Beat sheet',
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>%
  group_by(Name, ID) %>%
  summarize(caterpillar = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
            spider = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
            beetle = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
            truebug = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
            hopper = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
            ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0))

goodDataVis = fullDataset %>% # This would not have up to 50 in each set
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         ObservationMethod == 'Visual',
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>%
  group_by(Name, ID) %>%
  summarize(caterpillar = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
            spider = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
            beetle = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
            truebug = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
            hopper = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
            ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0))

# Data only on Acer rubrum
goodDataACRU = fullDataset %>%
  filter(julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         sciName == "Acer rubrum") %>%
  group_by(Name, ID, ObservationMethod) %>%
  summarize(caterpillar = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
            spider = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
            beetle = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
            truebug = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
            hopper = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
            ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0))

# (4) Read in maps and Landcover data
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
nam_devo = us_devo_ext + can_devo_ext


# Mosaic not working at the moment, as projected rasters don't have the same resolution. Need to fix...
#devo1km = terra::mosaic(can_devo_geog, us_devo_geog)


sites = distinct(fullDataset, Name, Region, Longitude, Latitude)

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


sites_sf <- st_as_sf(sites[sites$Name %in% goodSites$Name, ], 
                     coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs(us_devo_1km))

us_devo_df = as.data.frame(us_devo_crop, xy = TRUE)


# Doesn't currently include Canada landcover but good enough for now.
ggplot() +
  geom_raster(data = us_devo_df, aes(x = x, y = y, fill = `NLCD Land Cover Class`)) + 
  scale_fill_gradientn(colors = viridis(10)) + 
  geom_sf(data = sites_sf, col = "magenta", size = 2.5)
  

# (5) Join landcover data to arthropod data

dataset = left_join(goodData, sites, by = 'Name')

datasetBS = left_join(goodDataBS, sites, by = 'Name')

datasetVis = left_join(goodDataVis, sites, by = 'Name')

datasetACRU = left_join(goodDataACRU, sites, by = 'Name')


# Note that 3 historical sites in western NC from Coweeta only surveyed for caterpillars. These sites should be excluded from analyses of other arthropod groups. Perhaps it makes sense even to include them from caterpillar analyses so that results for caterpillars are geographically comparable to other groups

siteSummary = dataset %>%
  group_by(Name, Region, Longitude, Latitude, dev, forest) %>%
  summarize(nSurvs = n_distinct(ID),
            nSurvsBS = n_distinct(ID[ObservationMethod == 'Beat sheet']),
            nSurvsVis = n_distinct(ID[ObservationMethod == 'Visual']),
            propCat = n_distinct(ID[caterpillar == 1])/nSurvs,
            propBeet = n_distinct(ID[beetle == 1])/nSurvs,
            propTruebug = n_distinct(ID[truebug == 1])/nSurvs,
            propspider = n_distinct(ID[spider == 1])/nSurvs,
            propHopper = n_distinct(ID[hopper == 1])/nSurvs,
            propAnt = n_distinct(ID[ant == 1])/nSurvs,
            propCatBS = n_distinct(ID[caterpillar == 1 & ObservationMethod == 'Beat sheet'])/nSurvsBS,
            propBeetBS = n_distinct(ID[beetle == 1 & ObservationMethod == 'Beat sheet'])/nSurvsBS,
            propTruebugBS = n_distinct(ID[truebug == 1 & ObservationMethod == 'Beat sheet'])/nSurvsBS,
            propspiderBS = n_distinct(ID[spider == 1 & ObservationMethod == 'Beat sheet'])/nSurvsBS,
            propHopperBS = n_distinct(ID[hopper == 1 & ObservationMethod == 'Beat sheet'])/nSurvsBS,
            propAntBS = n_distinct(ID[ant == 1 & ObservationMethod == 'Beat sheet'])/nSurvsBS,
            propCatVis = n_distinct(ID[caterpillar == 1 & ObservationMethod == 'Visual'])/nSurvsVis,
            propBeetVis = n_distinct(ID[beetle == 1 & ObservationMethod == 'Visual'])/nSurvsVis,
            propTruebugVis = n_distinct(ID[truebug == 1 & ObservationMethod == 'Visual'])/nSurvsVis,
            propspiderVis = n_distinct(ID[spider == 1 & ObservationMethod == 'Visual'])/nSurvsVis,
            propHopperVis = n_distinct(ID[hopper == 1 & ObservationMethod == 'Visual'])/nSurvsVis,
            propAntVis = n_distinct(ID[ant == 1 & ObservationMethod == 'Visual'])/nSurvsVis)

vioplot(list(siteSummary$propHopperBS, siteSummary$propHopperVis,
             siteSummary$propTruebugBS, siteSummary$propTruebugVis,
             siteSummary$propCatBS, siteSummary$propCatVis,
             siteSummary$propspiderBS, siteSummary$propspiderVis,
             siteSummary$propAntBS, siteSummary$propAntVis,
             siteSummary$propBeetBS, siteSummary$propBeetVis
             ), xaxt = "n", 
        col = rep(c('gray90', 'gray30'), times = 6), las = 1)

mtext(c('hopper', 'true bug', 'caterpillar', 'spider', 'ant', 'beetle'), 1, 
      at = (2*1:6) -.5, cex = 1.8, line = 1)

legend('topleft', legend = c('Beat sheet', 'Visual'), pch = 15, pt.cex = 4, cex = 1.5, col = c('gray90', 'gray30'))



# (6) glms examining presence as a function of % developed or forest cover, latitude, and interaction



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

# No strong interaction, so model without interaction here:
bug.Dev.Latitude2 = glm(truebug ~ dev + Latitude + ObservationMethod, 
                        data = dataset, family = "binomial")

ant.Dev.Latitude = glm(ant ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")


# Forest cover models
cat.For.Latitude = glm(caterpillar ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")

spi.For.Latitude = glm(spider ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")

beet.For.Latitude = glm(beetle ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                        data = dataset, family = "binomial")

hop.For.Latitude = glm(hopper ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")

bug.For.Latitude = glm(truebug ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")

ant.For.Latitude = glm(ant ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = dataset, family = "binomial")



# GLM output

devOutput = data.frame(rbind(summary(cat.Dev.Latitude)$coefficients, 
                  summary(spi.Dev.Latitude)$coefficients, 
                  summary(beet.Dev.Latitude)$coefficients,
                  summary(hop.Dev.Latitude)$coefficients, 
                  summary(bug.Dev.Latitude)$coefficients,
                  summary(ant.Dev.Latitude)$coefficients))
devOutput$term = rep(c('Intercept', 'dev', 'Latitude', 'dev*Latitude', 'Method'), times = 6)
devOutput$Group = rep(c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'), each = 5)


forOutput = data.frame(rbind(summary(cat.For.Latitude)$coefficients, 
                             summary(spi.For.Latitude)$coefficients, 
                             summary(beet.For.Latitude)$coefficients,
                             summary(hop.For.Latitude)$coefficients, 
                             summary(bug.For.Latitude)$coefficients,
                             summary(ant.For.Latitude)$coefficients))
forOutput$term = rep(c('Intercept', 'forest', 'Latitude', 'forest*Latitude', 'Method'), times = 6)
forOutput$Group = rep(c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'), each = 5)


# Plots

# Arthropod images
catImage = readPNG('images/caterpillar.png')
antImage = readPNG('images/ant.png')
beetleImage = readPNG('images/beetle.png')
spiderImage = readPNG('images/spider.png')
hopperImage = readPNG('images/leafhopper.png')
truebugImage = readPNG('images/truebugs.png')


catDevPlot = interact_plot(cat.Dev.Latitude, pred = dev, modx = Latitude,
              y.label = "Prop. of surveys with caterpillars",
              x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
              colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(catImage, ymin = .085, ymax = .1, xmin = 60, xmax = 100)

beetDevPlot = interact_plot(beet.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with beetles",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  #ylim(0.12, 0.26) +
  annotation_raster(beetleImage, ymin = .3, ymax = .325, xmin = 60, xmax = 100) 

bugDevPlot = interact_plot(bug.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with true bugs",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(truebugImage, ymin = .09, ymax = .11, xmin = 0, xmax = 40)

spiDevPlot = interact_plot(spi.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with spiders",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(spiderImage, ymin = .32, ymax = .37, xmin = 60, xmax = 100)

hopDevPlot = interact_plot(hop.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with hoppers",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(hopperImage, ymin = .06, ymax = .09, xmin = 60, xmax = 100)

antDevPlot = interact_plot(ant.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with ants",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(antImage, ymin = .13, ymax = .15, xmin = 60, xmax = 100)

ggarrange(catDevPlot, spiDevPlot, beetDevPlot, bugDevPlot, hopDevPlot, antDevPlot, 
          ncol=3, nrow=2, common.legend = TRUE, legend="bottom")



catForPlot = interact_plot(cat.For.Latitude, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with caterpillars",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(catImage, ymin = .11, ymax = .13, xmin = 0, xmax = 40)

beetForPlot = interact_plot(beet.For.Latitude, pred = forest, modx = Latitude,
                            y.label = "Prop. of surveys with beetles",
                            x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                            colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  #ylim(0.15, 0.26) +
  annotation_raster(beetleImage, ymin = .295, ymax = .315, xmin = 0, xmax = 40)

bugForPlot = interact_plot(bug.For.Latitude, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with true bugs",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(truebugImage, ymin = .03, ymax = .05, xmin = 0, xmax = 40)

spiForPlot = interact_plot(spi.For.Latitude, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with spiders",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(spiderImage, ymin = .31, ymax = .38, xmin = 0, xmax = 40)

hopForPlot = interact_plot(hop.For.Latitude, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with hoppers",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(hopperImage, ymin = .06, ymax = .085, xmin = 0, xmax = 40)

antForPlot = interact_plot(ant.For.Latitude, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with ants",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(antImage, ymin = .135, ymax = .15, xmin = 5, xmax = 45)

ggarrange(catForPlot, spiForPlot, beetForPlot, bugForPlot, hopForPlot, antForPlot, 
          ncol=3, nrow=2, common.legend = TRUE, legend="bottom")






par(mar = c(4, 4, 0, 0), cex.axis = 2)
plot(us_devo_geog, ylim = c(25, 50), xlim = c(-100, -66), las = 1, xaxt = "n", background = 'white')
points(dataset$Longitude, dataset$Latitude, pch = 16, col = 'red')

# Confirming that there are a range of % developed cover values across all latitudinal bands
par(mar = c(6, 7, 1, 1), mgp = c(4, 1, 0))
plot(siteSummary$dev, siteSummary$Latitude, pch = 16, col = viridis(101)[floor(siteSummary$dev)+1], xlab = '% developed cover', ylab = 'Latitude', cex.lab = 3, cex.axis = 2, cex = 2, las = 1)
cor(siteSummary$Latitude, siteSummary$dev)



bluemono = colorRampPalette(c("#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF"))

plot(siteSummary$dev, siteSummary$propCat, 
     col = bluemono(100)[round(100*(max(siteSummary$Latitude) - siteSummary$Latitude)/(max(siteSummary$Latitude) - min(siteSummary$Latitude)))], 
     cex = log10(siteSummary$nSurvs), xlab = "% developed cover", 
     ylab = "Prop. caterpillars", pch = 16)



# (7.1) -------------- SIMPER SLOP ANALAYSIS -----------------------------------------------
# read more at: https://cran.r-project.org/web/packages/interactions/vignettes/interactions.html

# Urban cover

simDev_caterpillar <- sim_slopes(cat.Dev.Latitude, pred = dev, modx = Latitude,
                                 johnson_neyman = FALSE, digits = 4)

simDev_beetle <-sim_slopes(beet.Dev.Latitude, pred = dev, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simDev_spider <-sim_slopes(spi.Dev.Latitude, pred = dev, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simDev_hopper <-sim_slopes(hop.Dev.Latitude, pred = dev, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simDev_ant <-sim_slopes(ant.Dev.Latitude, pred = dev, modx = Latitude,
                        johnson_neyman = FALSE, digits = 4)

simDev_truebug <-sim_slopes(ant.Dev.Latitude, pred = dev, modx = Latitude,
                            johnson_neyman = FALSE, digits = 4)

# Forest
simFor_caterpillar <- sim_slopes(cat.For.Latitude, pred = forest, modx = Latitude,
                                 johnson_neyman = FALSE, digits = 4)

simFor_beetle <-sim_slopes(beet.For.Latitude, pred = forest, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simFor_spider <-sim_slopes(spi.For.Latitude, pred = forest, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simFor_hopper <-sim_slopes(hop.For.Latitude, pred = forest, modx = Latitude,
                           johnson_neyman = FALSE, digits = 4)

simFor_ant <-sim_slopes(ant.For.Latitude, pred = forest, modx = Latitude,
                        johnson_neyman = FALSE, digits = 4)

simFor_truebug <-sim_slopes(ant.For.Latitude, pred = forest, modx = Latitude,
                            johnson_neyman = FALSE, digits = 4)


For_sims = data.frame(rbind(simFor_caterpillar$slopes %>% mutate(Group = "caterpillar"), 
                            simFor_beetle$slopes %>% mutate(Group = "beetle"), 
                            simFor_spider$slopes %>% mutate(Group = "spider"),
                            simFor_ant$slopes %>% mutate(Group = "ant"), 
                            simFor_hopper$slopes %>% mutate(Group = "hopper"), 
                            simFor_truebug$slopes %>% mutate(Group = "truebug")))

Dev_sims = data.frame(rbind(simDev_caterpillar$slopes %>% mutate(Group = "caterpillar"), 
                            simDev_hopper$slopes %>% mutate(Group = "hopper"), 
                            simDev_beetle$slopes %>% mutate(Group = "beetle"), 
                            simDev_truebug$slopes %>% mutate(Group = "truebug"),
                            simDev_ant$slopes %>% mutate(Group = "ant"), 
                            simDev_spider$slopes %>% mutate(Group = "spider"))) %>% 
            mutate(Group = factor(Group, levels = c("caterpillar", "beetle", "ant",
                                             "hopper", "truebug", "spider")))




ggplot(Dev_sims, aes(x = Est., 
               y = Value.of.Latitude, 
               xmin = X2.5., 
               xmax = X97.5., 
               )) +
  geom_point(size = 3, aes(colour = Est.)) +
  geom_errorbarh(width = 0.2, size = 1, aes(colour = Est.)) +   # horizontal error bars
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2)+
  scale_color_gradientn(
    colours = c("#084594", "#2171B5", "#4292C6", "#9ECAE1"),
  ) +
  facet_wrap(~ Group, scales = "free_x") + # separate panel for each group
  labs(
    x = "Estimated effect of urban cover",
    y = "Latitude",
    title = "Effect of Urban Cover on Arthropod Groups by Latitude",
    subtitle = "Points = estimates. Estimate bars (=95% CI) touching the dashed line are not statistically significant."
  ) +
  guides(color = "none")+
  scale_y_continuous(limits = c(32, 44))+
  theme_bw()
 



# (7.2) -------------- SIMPER SLOP ANALAYSIS -----------------------------------------------
# read more at: https://cran.r-project.org/web/packages/interactions/vignettes/interactions.html

groups <- seq(min(dataset$Latitude), max(dataset$Latitude), length.out = 9)
mod.x = groups[-c(1,2, length(groups), length(groups)-1)]

# Urban cover

simDev_5caterpillar <- sim_slopes(cat.Dev.Latitude, pred = dev, modx = Latitude,modx.values = mod.x,
                                  johnson_neyman = FALSE, digits = 4)

simDev_5beetle <-sim_slopes(beet.Dev.Latitude, pred = dev, modx = Latitude,modx.values = mod.x,
                            johnson_neyman = FALSE, digits = 4)

simDev_5spider <-sim_slopes(spi.Dev.Latitude, pred = dev, modx = Latitude,modx.values = mod.x,
                            johnson_neyman = FALSE, digits = 4)

simDev_5hopper <-sim_slopes(hop.Dev.Latitude, pred = dev, modx = Latitude,modx.values = mod.x,
                            johnson_neyman = FALSE, digits = 4)

simDev_5ant <-sim_slopes(ant.Dev.Latitude, pred = dev, modx = Latitude,modx.values = mod.x,
                         johnson_neyman = FALSE, digits = 4)

simDev_5truebug <-sim_slopes(ant.Dev.Latitude, pred = dev, modx = Latitude,modx.values = mod.x,
                             johnson_neyman = FALSE, digits = 4)

# Forest
simFor_5caterpillar <- sim_slopes(cat.For.Latitude, pred = forest, modx = Latitude,modx.values = mod.x,
                                  johnson_neyman = FALSE, digits = 4)

simFor_5beetle <-sim_slopes(beet.For.Latitude, pred = forest, modx = Latitude,modx.values = mod.x,
                            johnson_neyman = FALSE, digits = 4)

simFor_5spider <-sim_slopes(spi.For.Latitude, pred = forest, modx = Latitude,modx.values = mod.x,
                            johnson_neyman = FALSE, digits = 4)

simFor_5hopper <-sim_slopes(hop.For.Latitude, pred = forest, modx = Latitude,modx.values = mod.x,
                            johnson_neyman = FALSE, digits = 4)

simFor_5ant <-sim_slopes(ant.For.Latitude, pred = forest, modx = Latitude,modx.values = mod.x,
                         johnson_neyman = FALSE, digits = 4)

simFor_5truebug <-sim_slopes(ant.For.Latitude, pred = forest, modx = Latitude,modx.values = mod.x,
                             johnson_neyman = FALSE, digits = 4)


For_sims5 = data.frame(rbind(simFor_5caterpillar$slopes %>% mutate(Group = "caterpillar"), 
                             simFor_5beetle$slopes %>% mutate(Group = "beetle"), 
                             simFor_5spider$slopes %>% mutate(Group = "spider"),
                             simFor_5ant$slopes %>% mutate(Group = "ant"), 
                             simFor_5hopper$slopes %>% mutate(Group = "hopper"), 
                             simFor_5truebug$slopes %>% mutate(Group = "truebug")))

Dev_sims5 = data.frame(rbind(simDev_5caterpillar$slopes %>% mutate(Group = "caterpillar"), 
                             simDev_5hopper$slopes %>% mutate(Group = "hopper"), 
                             simDev_5beetle$slopes %>% mutate(Group = "beetle"), 
                             simDev_5truebug$slopes %>% mutate(Group = "truebug"),
                             simDev_5ant$slopes %>% mutate(Group = "ant"), 
                             simDev_5spider$slopes %>% mutate(Group = "spider"))) %>% 
  mutate(Group = factor(Group, levels = c("caterpillar", "beetle", "ant",
                                          "hopper", "truebug", "spider")))




ggplot(Dev_sims5, aes(x = Est., 
                      y = Value.of.Latitude, 
                      xmin = X2.5., 
                      xmax = X97.5.)) +
  geom_point(size = 2, aes(colour = Est.)) +
  geom_errorbarh(width = 0.02, size = 1, aes(colour = Est.)) +   # horizontal error bars
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2)+
  scale_color_gradientn(
    colours = c("red", "#424242", "green"),
  ) +
  facet_wrap(~ Group, scales = "free_x") + # separate panel for each group
  labs(
    x = "Estimated effect of urban cover",
    y = "Latitude",
    title = "Effect of Urban Cover on Arthropod Groups by Latitude",
    subtitle = "Points = estimates. Estimate bars (=95% CI) touching the dashed line are not statistically significant."
  ) +
  guides(color = "none")+
  scale_y_continuous(limits = c(34, 45))+
  scale_x_continuous(limits = c(min(Dev_sims$X2.5.), max(Dev_sims$X97.5.)))+
  theme_bw()



# Forest cover %
ggplot(For_sims5, aes(x = Est., 
                      y = Value.of.Latitude, 
                      xmin = X2.5., 
                      xmax = X97.5.)) +
  geom_point(size = 2, aes(colour = Est.)) +
  geom_errorbarh(width = 0.02, size = 1, aes(colour = Est.)) +   # horizontal error bars
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2)+
  scale_color_gradientn(
    colours = c("red", "#424242", "green"),
  ) +
  facet_wrap(~ Group, scales = "free_x") + # separate panel for each group
  labs(
    x = "Estimated effect of Forest cover",
    y = "Latitude",
    title = "Effect of Forest Cover on Arthropod Groups by Latitude",
    subtitle = "Points = estimates. Estimate bars (=95% CI) touching the dashed line are not statistically significant."
  ) +
  guides(color = "none")+
  scale_y_continuous(limits = c(34, 45))+
  scale_x_continuous(limits = c(min(For_sims5$X2.5.), max(For_sims5$X97.5.)))+
  theme_bw()


# (8) Analysis restricted to a single tree species, Acer rubrum, the most widespread species in the dataset

minSurveysACRU = 10


siteSummaryACRU = datasetACRU %>%
  group_by(Name, Region, Longitude, Latitude, dev, forest) %>%
  summarize(nSurvs = n_distinct(ID),
            nSurvsBS = n_distinct(ID[ObservationMethod == 'Beat sheet']),
            nSurvsVis = n_distinct(ID[ObservationMethod == 'Visual']),
            propCat = n_distinct(ID[caterpillar == 1])/nSurvs,
            propBeet = n_distinct(ID[beetle == 1])/nSurvs,
            propTruebug = n_distinct(ID[truebug == 1])/nSurvs,
            propspider = n_distinct(ID[spider == 1])/nSurvs,
            propHopper = n_distinct(ID[hopper == 1])/nSurvs,
            propAnt = n_distinct(ID[ant == 1])/nSurvs) %>%
  filter(nSurvs >= minSurveysACRU)

ACRUdata = datasetACRU %>%
  filter(Name %in% siteSummaryACRU$Name)

# GLMs
cat.Dev.Lat.ACRU = glm(caterpillar ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")

spi.Dev.Lat.ACRU = glm(spider ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")

beet.Dev.Lat.ACRU = glm(beetle ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                        data = ACRUdata, family = "binomial")

hop.Dev.Lat.ACRU = glm(hopper ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")

bug.Dev.Lat.ACRU = glm(truebug ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")

# No strong interaction, so model without interaction here:
bug.Dev.Lat.ACRU2 = glm(truebug ~ dev + Latitude + ObservationMethod, 
                        data = ACRUdata, family = "binomial")

ant.Dev.Lat.ACRU = glm(ant ~ dev + Latitude + dev*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")


# Forest cover models
cat.For.Lat.ACRU = glm(caterpillar ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")

spi.For.Lat.ACRU = glm(spider ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")

beet.For.Lat.ACRU = glm(beetle ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                        data = ACRUdata, family = "binomial")

hop.For.Lat.ACRU = glm(hopper ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")

bug.For.Lat.ACRU = glm(truebug ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")

ant.For.Lat.ACRU = glm(ant ~ forest + Latitude + forest*Latitude + ObservationMethod, 
                       data = ACRUdata, family = "binomial")


# GLM output

devOutputACRU = data.frame(rbind(summary(cat.Dev.Lat.ACRU)$coefficients, 
                             summary(spi.Dev.Lat.ACRU)$coefficients, 
                             summary(beet.Dev.Lat.ACRU)$coefficients,
                             summary(hop.Dev.Lat.ACRU)$coefficients, 
                             summary(bug.Dev.Lat.ACRU)$coefficients,
                             summary(ant.Dev.Lat.ACRU)$coefficients))
devOutputACRU$term = rep(c('Intercept', 'dev', 'Latitude', 'dev*Latitude', 'Method'), times = 6)
devOutputACRU$Group = rep(c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'), each = 5)


forOutputACRU = data.frame(rbind(summary(cat.For.Lat.ACRU)$coefficients, 
                             summary(spi.For.Lat.ACRU)$coefficients, 
                             summary(beet.For.Lat.ACRU)$coefficients,
                             summary(hop.For.Lat.ACRU)$coefficients, 
                             summary(bug.For.Lat.ACRU)$coefficients,
                             summary(ant.For.Lat.ACRU)$coefficients))
forOutputACRU$term = rep(c('Intercept', 'forest', 'Latitude', 'forest*Latitude', 'Method'), times = 6)
forOutputACRU$Group = rep(c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'), each = 5)

# Interaction plots
catDevPlotACRU = interact_plot(cat.Dev.Lat.ACRU, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with caterpillars",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(catImage, ymin = .045, ymax = .06, xmin = 0, xmax = 40)

beetDevPlotACRU = interact_plot(beet.Dev.Lat.ACRU, pred = dev, modx = Latitude,
                            y.label = "Prop. of surveys with beetles",
                            x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                            colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  #ylim(0.12, 0.26) +
  annotation_raster(beetleImage, ymin = .18, ymax = .23, xmin = 33, xmax = 67) 

bugDevPlotACRU = interact_plot(bug.Dev.Lat.ACRU, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with true bugs",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(truebugImage, ymin = .16, ymax = .22, xmin = 0, xmax = 40)

spiDevPlotACRU = interact_plot(spi.Dev.Lat.ACRU, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with spiders",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(spiderImage, ymin = .32, ymax = .4, xmin = 60, xmax = 100)

hopDevPlotACRU = interact_plot(hop.Dev.Lat.ACRU, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with hoppers",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(hopperImage, ymin = .03, ymax = .07, xmin = 60, xmax = 100)

antDevPlotACRU = interact_plot(ant.Dev.Lat.ACRU, pred = dev, modx = Latitude,
                           y.label = "Prop. of surveys with ants",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(antImage, ymin = .075, ymax = .1, xmin = 60, xmax = 100)

ggarrange(catDevPlotACRU, spiDevPlotACRU, beetDevPlotACRU, bugDevPlotACRU, hopDevPlotACRU, antDevPlotACRU, 
          ncol=3, nrow=2, common.legend = TRUE, legend="bottom")



catForPlotACRU = interact_plot(cat.For.Lat.ACRU, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with caterpillars",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(catImage, ymin = .105, ymax = .125, xmin = 33, xmax = 67)

beetForPlotACRU = interact_plot(beet.For.Lat.ACRU, pred = forest, modx = Latitude,
                            y.label = "Prop. of surveys with beetles",
                            x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                            colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  #ylim(0.15, 0.26) +
  annotation_raster(beetleImage, ymin = .16, ymax = .205, xmin = 0, xmax = 35)

bugForPlotACRU = interact_plot(bug.For.Lat.ACRU, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with true bugs",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(truebugImage, ymin = .1, ymax = .13, xmin = 60, xmax = 100)

spiForPlotACRU = interact_plot(spi.For.Lat.ACRU, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with spiders",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(spiderImage, ymin = .31, ymax = .4, xmin = 0, xmax = 40)

hopForPlotACRU = interact_plot(hop.For.Lat.ACRU, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with hoppers",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(hopperImage, ymin = .01, ymax = .05, xmin = 0, xmax = 40)

antForPlotACRU = interact_plot(ant.For.Lat.ACRU, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with ants",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(antImage, ymin = .085, ymax = .11, xmin = 5, xmax = 45)

ggarrange(catForPlotACRU, spiForPlotACRU, beetForPlotACRU, bugForPlotACRU, hopForPlotACRU, antForPlotACRU, 
          ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

