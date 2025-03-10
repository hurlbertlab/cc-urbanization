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
         Longitude > -100) %>%
  group_by(Name) %>%
  summarize(nSurvs = n_distinct(ID)) %>%
  filter(nSurvs >= minSurveys)

# Dataset thru 2024 includes 156 sites with at least 50 surveys in the seasonal window


# (3) Reorganize data to presence-absence of arthropod groups on each survey for these good sites
goodData = fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow) %>%
  group_by(Name, ID) %>%
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
us_forest_geog = project(us_forest_1km, crs(us_devo3_geog))



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
  geom_sf(data = sites_sf, aes(color = "red"), size = 2)
  

# (5) Join landcover data to arthropod data

dataset = left_join(goodData, sites, by = 'Name')


# Note that 3 historical sites in western NC from Coweeta only surveyed for caterpillars. These sites should be excluded from analyses of other arthropod groups. Perhaps it makes sense even to include them from caterpillar analyses so that results for caterpillars are geographically comparable to other groups

dataset2 = dataset %>%
  filter(!Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'))

siteSummary = dataset %>%
  group_by(Name, Region, Longitude, Latitude, dev, forest) %>%
  summarize(nSurvs = n_distinct(ID),
            propCat = n_distinct(ID[caterpillar == 1])/nSurvs,
            propBeet = n_distinct(ID[beetle == 1])/nSurvs,
            propTruebug = n_distinct(ID[truebug == 1])/nSurvs,
            propspider = n_distinct(ID[spider == 1])/nSurvs,
            propHopper = n_distinct(ID[hopper == 1])/nSurvs,
            propAnt = n_distinct(ID[ant == 1])/nSurvs,)

# (6) glms examining presence as a function of % developed or forest cover, latitude, and interaction



# Developed cover models
cat.Dev.Latitude = glm(caterpillar ~ dev + Latitude + dev*Latitude, 
                          data = dataset2, family = "binomial")

spi.Dev.Latitude = glm(spider ~ dev + Latitude + dev*Latitude, 
                       data = dataset2, family = "binomial")

beet.Dev.Latitude = glm(beetle ~ dev + Latitude + dev*Latitude, 
                       data = dataset2, family = "binomial")

hop.Dev.Latitude = glm(hopper ~ dev + Latitude + dev*Latitude, 
                       data = dataset2, family = "binomial")

bug.Dev.Latitude = glm(truebug ~ dev + Latitude + dev*Latitude, 
                       data = dataset2, family = "binomial")

# No strong interaction, so model without interaction here:
bug.Dev.Latitude2 = glm(truebug ~ dev + Latitude, data = dataset2, family = "binomial")

ant.Dev.Latitude = glm(ant ~ dev + Latitude + dev*Latitude, 
                       data = dataset2, family = "binomial")

# Forest cover models
cat.For.Latitude = glm(caterpillar ~ forest + Latitude + forest*Latitude, 
                       data = dataset2, family = "binomial")

spi.For.Latitude = glm(spider ~ forest + Latitude + forest*Latitude, 
                       data = dataset2, family = "binomial")

beet.For.Latitude = glm(beetle ~ forest + Latitude + forest*Latitude, 
                        data = dataset2, family = "binomial")

hop.For.Latitude = glm(hopper ~ forest + Latitude + forest*Latitude, 
                       data = dataset2, family = "binomial")

bug.For.Latitude = glm(truebug ~ forest + Latitude + forest*Latitude, 
                       data = dataset2, family = "binomial")

ant.For.Latitude = glm(ant ~ forest + Latitude + forest*Latitude, 
                       data = dataset2, family = "binomial")


# GLM output

devOutput = data.frame(rbind(summary(cat.Dev.Latitude)$coefficients, 
                  summary(spi.Dev.Latitude)$coefficients, 
                  summary(beet.Dev.Latitude)$coefficients,
                  summary(hop.Dev.Latitude)$coefficients, 
                  summary(bug.Dev.Latitude)$coefficients,
                  summary(ant.Dev.Latitude)$coefficients))
devOutput$term = rep(c('Intercept', 'dev', 'Latitude', 'dev*Latitude'), times = 6)
devOutput$Group = rep(c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'), each = 4)


forOutput = data.frame(rbind(summary(cat.For.Latitude)$coefficients, 
                             summary(spi.For.Latitude)$coefficients, 
                             summary(beet.For.Latitude)$coefficients,
                             summary(hop.For.Latitude)$coefficients, 
                             summary(bug.For.Latitude)$coefficients,
                             summary(ant.For.Latitude)$coefficients))
forOutput$term = rep(c('Intercept', 'forest', 'Latitude', 'forest*Latitude'), times = 6)
forOutput$Group = rep(c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'), each = 4)


# Plots

# Arthropod images
catImage = readPNG('images/caterpillar.png')
antImage = readPNG('images/ant.png')
beetleImage = readPNG('images/beetle.png')
spiderImage = readPNG('images/spider.png')
hopperImage = readPNG('images/leafhopper.png')
truebugImage = readPNG('images/truebugs.png')


plot1 = intplotOriginLatitudeName + 
  theme_bw() +
  annotation_raster(caterpillar, ymin = .041, ymax= .052,xmin = 39, xmax = 43.4) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(margin = margin(t = 6)), 
        axis.title.y = element_text(margin = margin(l = 12), vjust = 2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 20))

catDevPlot = interact_plot(cat.Dev.Latitude, pred = dev, modx = Latitude,
              y.label = "Prop. of surveys with caterpillars",
              x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
              colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(catImage, ymin = .085, ymax = .1, xmin = 60, xmax = 100)

beetDevPlot = interact_plot(beet.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with beetles",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  ylim(0.12, 0.26) +
  annotation_raster(beetleImage, ymin = .12, ymax = .16, xmin = 60, xmax = 100) 

bugDevPlot = interact_plot(bug.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with true bugs",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(truebugImage, ymin = .08, ymax = .1, xmin = 0, xmax = 40)

spiDevPlot = interact_plot(spi.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with spiders",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(spiderImage, ymin = .27, ymax = .33, xmin = 60, xmax = 100)

hopDevPlot = interact_plot(hop.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with hoppers",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(hopperImage, ymin = .06, ymax = .09, xmin = 60, xmax = 100)

antDevPlot = interact_plot(ant.Dev.Latitude, pred = dev, modx = Latitude,
                        y.label = "Prop. of surveys with ants",
                        x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                        colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) +
  annotation_raster(antImage, ymin = .11, ymax = .13, xmin = 60, xmax = 100)

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
  ylim(0.15, 0.26) +
  annotation_raster(beetleImage, ymin = .15, ymax = .175, xmin = 0, xmax = 40)

bugForPlot = interact_plot(bug.For.Latitude, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with true bugs",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(truebugImage, ymin = .025, ymax = .045, xmin = 0, xmax = 40)

spiForPlot = interact_plot(spi.For.Latitude, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with spiders",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(spiderImage, ymin = .28, ymax = .35, xmin = 0, xmax = 40)

hopForPlot = interact_plot(hop.For.Latitude, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with hoppers",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(hopperImage, ymin = .05, ymax = .075, xmin = 0, xmax = 40)

antForPlot = interact_plot(ant.For.Latitude, pred = forest, modx = Latitude,
                           y.label = "Prop. of surveys with ants",
                           x.lab = "% Forest cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkgreen', 'green', 'lightgreen'), line.thickness = 2) +
  annotation_raster(antImage, ymin = .11, ymax = .125, xmin = 5, xmax = 45)

ggarrange(catForPlot, spiForPlot, beetForPlot, bugForPlot, hopForPlot, antForPlot, 
          ncol=3, nrow=2, common.legend = TRUE, legend="bottom")







par(mar = c(4, 4, 0, 0), cex.axis = 2)
plot(us_devo_geog, ylim = c(25, 50), xlim = c(-100, -66), las = 1, xaxt = "n", background = 'white')
points(dataset2$Longitude, dataset2$Latitude, pch = 16, col = 'red')

# Confirming that there are a range of % developed cover values across all latitudinal bands
par(mar = c(6, 7, 1, 1), mgp = c(4, 1, 0))
plot(siteSummary$dev, siteSummary$Latitude, pch = 16, col = viridis(101)[floor(siteSummary$dev)+1], xlab = '% developed cover', ylab = 'Latitude', cex.lab = 3, cex.axis = 2, cex = 2, las = 1)
cor(siteSummary$Latitude, siteSummary$dev)



bluemono = colorRampPalette(c("#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF"))

plot(siteSummary$dev, siteSummary$propCat, 
     col = bluemono(100)[round(100*(max(siteSummary$Latitude) - siteSummary$Latitude)/(max(siteSummary$Latitude) - min(siteSummary$Latitude)))], 
     cex = log10(siteSummary$nSurvs), xlab = "% developed cover", 
     ylab = "Prop. caterpillars", pch = 16)

