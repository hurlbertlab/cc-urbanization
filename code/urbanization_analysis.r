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
    

# (4) Read in Landcover data
#     Sources: NLCD
#              Canada Land Cover 2020

# This raster is in geographic coverage, but at 3 km resolution
us_devo3km_geog = rast('data/devo_geographic_3km.tif')

# Rasters at 1 km resolution
canada_devo_1km = rast('data/devo_1km_toronto.tif')
us_devo_1km = rast('data/devo_1km_east.tif')
canada_forest_1km = rast('data/forest_1km_toronto.tif')
us_forest_1km = rast('data/forest_1km_east.tif')

can_devo_geog = project(canada_devo_1km, crs(us_devo3km_geog))
us_devo_geog = project(us_devo_1km, crs(us_devo3km_geog))

can_forest_geog = project(canada_forest_1km, crs(us_devo3km_geog))
us_forest_geog = project(us_forest_1km, crs(us_devo3km_geog))


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


# (5) Join landcover data to arthropod data

dataset = left_join(goodData, sites, by = 'Name')


# Note that 3 historical sites in western NC from Coweeta only surveyed for caterpillars. These sites should be excluded from analyses of other arthropod groups. Perhaps it makes sense even to include them from caterpillar analyses so that results for caterpillars are geographically comparable to other groups

nonCoweetaDataset = dataset %>%
  filter(!Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'))

# (6) glms examining presence as a function of % developed or forest cover, latitude, and interaction



# Developed cover models
cat.Dev.Latitude = glm(caterpillar ~ dev + Latitude + dev*Latitude, 
                          data = nonCoweetaDataset, family = "binomial")

spi.Dev.Latitude = glm(spider ~ dev + Latitude + dev*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")

beet.Dev.Latitude = glm(beetle ~ dev + Latitude + dev*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")

hop.Dev.Latitude = glm(hopper ~ dev + Latitude + dev*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")

bug.Dev.Latitude = glm(truebug ~ dev + Latitude + dev*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")

ant.Dev.Latitude = glm(ant ~ dev + Latitude + dev*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")

# Forest cover models
cat.For.Latitude = glm(caterpillar ~ forest + Latitude + forest*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")

spi.For.Latitude = glm(spider ~ forest + Latitude + forest*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")

beet.For.Latitude = glm(beetle ~ forest + Latitude + forest*Latitude, 
                        data = nonCoweetaDataset, family = "binomial")

hop.For.Latitude = glm(hopper ~ forest + Latitude + forest*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")

bug.For.Latitude = glm(truebug ~ forest + Latitude + forest*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")

ant.For.Latitude = glm(ant ~ forest + Latitude + forest*Latitude, 
                       data = nonCoweetaDataset, family = "binomial")


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








plot(us_devo_geog, ylim = c(25, 50), xlim = c(-100, -66))
points(nonCoweetaDataset$Longitude, nonCoweetaDataset$Latitude, pch = 16, col = 'red')

# Confirming that there are a range of % developed cover values across all latitudinal bands
plot(nonCoweetaDataset$dev, nonCoweetaDataset$Latitude, pch = 16, col = 'red', xlab = '% developed cover', ylab = 'Latitude')
cor(nonCoweetaDataset$Latitude, nonCoweetaDataset$dev)




