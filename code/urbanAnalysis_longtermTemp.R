
library(mgcv)
library(gratia)
library(patchwork)
library(ENMeval)
library(dismo)
library(ENMeval)
library(terra)
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
library(gt)

# (1) Read in latest Caterpillars Count! raw dataset from the caterpillars-analysis-public repo----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file)) %>% 
  filter(Year <= 2025)

urban500m = read.csv("data/urbanization_USA_Canada_500m.csv") %>% 
   rename(dev = urban_percent) %>% dplyr::select(-ID, -X, -Longitude, -Region)

clim2021_2040 <- stack("largeFile/climate/wc2.1_10m_bioc_ACCESS-CM2_ssp245_2021-2040.tif")
clim2061_2080 <- stack("largeFile/climate/wc2.1_10m_bioc_ACCESS-CM2_ssp245_2061-2080.tif")
clim2081_2100 <- stack("largeFile/climate/wc2.1_10m_bioc_ACCESS-CM2_ssp245_2081-2100.tif")

myCRS1 = CRS("+init=epsg:4326") # WGS 84


clim_hist_list = list.files("largeFile/climate/Historical_10m",pattern=".tif$",full.names = T)
clim_hist = raster::stack(clim_hist_list) 
names(clim_hist) = gsub("wc2.1_10m_","",names(clim_hist))

plot(clim_hist[[2]])


# (2) Filter dataset ----
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

sites = distinct(fullDataset, Name, Region, Longitude, Latitude) %>% 
  inner_join(urban500m, by = c("Name", "Latitude"))

myCRS1 = CRS("+init=epsg:4326") # WGS 84

sitesLatLon = sites %>% 
  rename(lon = Longitude,
         lat = Latitude) 
# make occurrences spatial
coordinates(sitesLatLon) <- ~ lon + lat
crs(sitesLatLon) <- myCRS1

sitesLatLon$BIO10_hist <- raster::extract(clim_hist[[2]], sitesLatLon)
sitesLatLon$BIO05_hist <- raster::extract(clim_hist[[15]], sitesLatLon)

sites = sites %>%
  mutate(
    BIO10 = sitesLatLon$BIO10_hist,
    BIO05 = sitesLatLon$BIO05_hist)

sites$BIO10_2021_2040 <-
  raster::extract(clim2021_2040[["bio10"]], sitesLatLon)

sites$BIO10_2061_2080 <-
  raster::extract(clim2061_2080[["bio10"]], sitesLatLon)

sites$BIO10_2081_2100 <-
  raster::extract(clim2081_2100[["bio10"]], sitesLatLon)

sites$BIO05_2021_2040 <-
  raster::extract(clim2021_2040[["bio05"]], sitesLatLon)

sites$BIO05_2061_2080 <-
  raster::extract(clim2061_2080[["bio05"]], sitesLatLon)

sites$BIO05_2081_2100 <-
  raster::extract(clim2081_2100[["bio05"]], sitesLatLon)


dataset = inner_join(goodData, sites, by = 'Name') %>% 
  filter(!is.na(ID))

dataset %>% 
  filter(nSurvs < minSurveys) # All good!

# 
# clim2021_2040_30s <- cmip6_world(
#   var = "bio",
#   model = "ACCESS-CM2",
#   ssp = "245",
#   time = "2061-2080",
#   res = 0.5,
#   path = "largeFile/climate"
# )
# The geodata server is temporary out of service for maintenance. It should be back on 22 June


catBIO10 <- gam(caterpillar ~ dev + s(BIO10,  k = 4) + ObservationMethod,
  family = binomial(link = "logit"),
  data = dataset)

summary(catBIO10)
draw(catBIO10, select = "s(BIO10)")+
  geom_hline(yintercept = 0, linetype = "dashed")


catBIO05 <- gam(caterpillar ~ dev + s(BIO05,  k = 4) + ObservationMethod,
                family = binomial(link = "logit"),
                data = dataset)
summary(catBIO05)
draw(catBIO05, select = "s(BIO05)")+
  geom_hline(yintercept = 0, linetype = "dashed")

spiBIO10 <- gam(spider ~ dev + s(BIO10,  k = 4) + ObservationMethod,
                family = binomial(link = "logit"),
                data = dataset)
summary(spiBIO10)
draw(spiBIO10, select = "s(BIO10)")+
  geom_hline(yintercept = 0, linetype = "dashed")



spiBIO05 <- gam(spider ~ dev + s(BIO05,  k = 4) + ObservationMethod,
                family = binomial(link = "logit"),
                data = dataset)
summary(spiBIO05)
draw(spiBIO05, select = "s(BIO05)")+
  geom_hline(yintercept = 0, linetype = "dashed")


antBIO10 <- gam(ant ~ dev + s(BIO10,  k = 4) + ObservationMethod,
                family = binomial(link = "logit"),
                data = dataset)
summary(antBIO10)
draw(antBIO10, select = "s(BIO10)")+
  geom_hline(yintercept = 0, linetype = "dashed")


antBIO05 <- gam(ant ~ dev + s(BIO05,  k = 4) + ObservationMethod,
                family = binomial(link = "logit"),
                data = dataset)
summary(antBIO05)
draw(antBIO05, select = "s(BIO05)")+
  geom_hline(yintercept = 0, linetype = "dashed")




truebugBIO10 <- gam(truebug ~ dev + s(BIO10,  k = 4) + ObservationMethod,
                family = binomial(link = "logit"),
                data = dataset)
summary(truebugBIO10)
draw(truebugBIO10, select = "s(BIO10)")+
  geom_hline(yintercept = 0, linetype = "dashed")


truebugBIO05 <- gam(truebug ~ dev + s(BIO05,  k = 4) + ObservationMethod,
                family = binomial(link = "logit"),
                data = dataset)
summary(truebugBIO05)

draw(truebugBIO05, select = "s(BIO05)") +
  geom_hline(yintercept = 0, linetype = "dashed")


hopperBIO10 <- gam(hopper ~ dev + s(BIO10,  k = 4) + ObservationMethod,
                    family = binomial(link = "logit"),
                    data = dataset)
summary(hopperBIO10)
draw(hopperBIO10, select = "s(BIO10)")+
  geom_hline(yintercept = 0, linetype = "dashed")


hopperBIO05 <- gam(hopper ~ dev + s(BIO05,  k = 4) + ObservationMethod,
                    family = binomial(link = "logit"),
                    data = dataset)
summary(hopperBIO05)

draw(hopperBIO05, select = "s(BIO05)") +
  geom_hline(yintercept = 0, linetype = "dashed")

grasshopperBIO10 <- gam(grasshopper ~ dev + s(BIO10,  k = 4) + ObservationMethod,
                   family = binomial(link = "logit"),
                   data = dataset)
summary(grasshopperBIO10)
draw(grasshopperBIO10, select = "s(BIO10)")+
  geom_hline(yintercept = 0, linetype = "dashed")


grasshopperBIO05 <- gam(grasshopper ~ dev + s(BIO05,  k = 4) + ObservationMethod,
                   family = binomial(link = "logit"),
                   data = dataset)
summary(grasshopperBIO05)

draw(grasshopperBIO05, select = "s(BIO05)") +
  geom_hline(yintercept = 0, linetype = "dashed")



beetleBIO10 <- gam(beetle ~ dev + s(BIO10,  k = 4) + ObservationMethod,
                   family = binomial(link = "logit"),
                   data = dataset)
summary(beetleBIO10)
draw(beetleBIO10, select = "s(BIO10)")+
  geom_hline(yintercept = 0, linetype = "dashed")


beetleBIO05 <- gam(beetle ~ dev + s(BIO05,  k = 4) + ObservationMethod,
                   family = binomial(link = "logit"),
                   data = dataset)
summary(beetleBIO05)

draw(beetleBIO05, select = "s(BIO05)") +
  geom_hline(yintercept = 0, linetype = "dashed")






plot_gam_response <- function(model, dataset, var, title_text) {
  newdat <- data.frame(
    x = seq(min(dataset[[var]], na.rm = TRUE),
            max(dataset[[var]], na.rm = TRUE),
            length.out = 200),
    dev = mean(dataset$dev, na.rm = TRUE),
    ObservationMethod = "Visual"
  )
  # rename x back to the variable name expected by the model
  names(newdat)[1] <- var
  pred <- predict(model,
                  newdata = newdat,
                  type = "response",
                  se.fit = TRUE) # standard error
  
  newdat$fit <- pred$fit
  newdat$lower <- pmax(0, pred$fit - 1.96 * pred$se.fit) # ensure p don't get below 0
  newdat$upper <- pmin(1, pred$fit + 1.96 * pred$se.fit) # ensure p don't get above 1
  
  ggplot(newdat, aes(x = .data[[var]], y = fit)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.2) +
    ylab("Predicted probability of occurrence") +
    xlab(var) +
    ggtitle(title_text) +
    theme_bw()
}

p_cat10 <- plot_gam_response(catBIO10, dataset,"BIO10", "Caterpillar")
p_spi10 <- plot_gam_response(spiBIO10, dataset, "BIO10","Spider")
p_ant10 <- plot_gam_response(antBIO10, dataset,"BIO10", "Ant")
p_true10 <- plot_gam_response(truebugBIO10, dataset,"BIO10", "True bug")
p_hopper10 <- plot_gam_response(hopperBIO10, dataset,"BIO10", "Leafhopper")
p_beetle10 <- plot_gam_response(beetleBIO10, dataset,"BIO10", "Beetle")
p_orthopteran10 <- plot_gam_response(grasshopperBIO10, dataset, "BIO10","Orthoptera")

p_cat05 <- plot_gam_response(catBIO05, dataset, "BIO05", "Caterpillar")
p_spi05 <- plot_gam_response(spiBIO05, dataset, "BIO05","Spider")
p_ant05 <- plot_gam_response(antBIO05, dataset, "BIO05","Ant")
p_true05 <- plot_gam_response(truebugBIO05, dataset, "BIO05","True bug")
p_hopper05 <- plot_gam_response(hopperBIO05, dataset, "BIO05","Leafhopper")
p_beetle05 <- plot_gam_response(beetleBIO05, dataset, "BIO05","Beetle")
p_orthopteran05 <- plot_gam_response(grasshopperBIO05, dataset, "BIO05","Orthoptera")

BIO10_panel =
  wrap_plots(
    p_cat10,
    p_orthopteran10,
    p_spi10,
    p_ant10,
    p_beetle10,
    p_hopper10,
    p_true10,
    ncol = 2
  ) +
  plot_annotation(
    title = "Smooth effects of BIO10 on arthropod occurrence")

BIO10_panel


BIO05_panel =
  wrap_plots(
    p_cat05,
    p_orthopteran05,
    p_spi05,
    p_ant05,
    p_beetle05,
    p_hopper05,
    p_true05,
    ncol = 2
  ) +
  plot_annotation(
    title = "Smooth effects of BIO05 on arthropod occurrence")

BIO05_panel






catBIO10_2021_2040 <- gam(caterpillar ~ dev + s(BIO10_2021_2040,  k = 4) + ObservationMethod,
                          family = binomial(link = "logit"),
                          data = dataset)

summary(catBIO10_2021_2040)
draw(catBIO10_2021_2040, select = "s(BIO10_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")


catBIO05_2021_2040 <- gam(caterpillar ~ dev + s(BIO05_2021_2040,  k = 4) + ObservationMethod,
                          family = binomial(link = "logit"),
                          data = dataset)
summary(catBIO05_2021_2040)
draw(catBIO05_2021_2040, select = "s(BIO05_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")

spiBIO10_2021_2040 <- gam(spider ~ dev + s(BIO10_2021_2040,  k = 4) + ObservationMethod,
                          family = binomial(link = "logit"),
                          data = dataset)
summary(spiBIO10_2021_2040)
draw(spiBIO10_2021_2040, select = "s(BIO10_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")



spiBIO05_2021_2040 <- gam(spider ~ dev + s(BIO05_2021_2040,  k = 4) + ObservationMethod,
                          family = binomial(link = "logit"),
                          data = dataset)
summary(spiBIO05_2021_2040)
draw(spiBIO05_2021_2040, select = "s(BIO05_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")


antBIO10_2021_2040 <- gam(ant ~ dev + s(BIO10_2021_2040,  k = 4) + ObservationMethod,
                          family = binomial(link = "logit"),
                          data = dataset)
summary(antBIO10_2021_2040)
draw(antBIO10_2021_2040, select = "s(BIO10_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")


antBIO05_2021_2040 <- gam(ant ~ dev + s(BIO05_2021_2040,  k = 4) + ObservationMethod,
                          family = binomial(link = "logit"),
                          data = dataset)
summary(antBIO05_2021_2040)
draw(antBIO05_2021_2040, select = "s(BIO05_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")




truebugBIO10_2021_2040 <- gam(truebug ~ dev + s(BIO10_2021_2040,  k = 4) + ObservationMethod,
                              family = binomial(link = "logit"),
                              data = dataset)
summary(truebugBIO10_2021_2040)
draw(truebugBIO10_2021_2040, select = "s(BIO10_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")


truebugBIO05_2021_2040 <- gam(truebug ~ dev + s(BIO05_2021_2040,  k = 4) + ObservationMethod,
                              family = binomial(link = "logit"),
                              data = dataset)
summary(truebugBIO05_2021_2040)

draw(truebugBIO05_2021_2040, select = "s(BIO05_2021_2040)") +
  geom_hline(yintercept = 0, linetype = "dashed")


hopperBIO10_2021_2040 <- gam(hopper ~ dev + s(BIO10_2021_2040,  k = 4) + ObservationMethod,
                             family = binomial(link = "logit"),
                             data = dataset)
summary(hopperBIO10_2021_2040)
draw(hopperBIO10_2021_2040, select = "s(BIO10_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")


hopperBIO05_2021_2040 <- gam(hopper ~ dev + s(BIO05_2021_2040,  k = 4) + ObservationMethod,
                             family = binomial(link = "logit"),
                             data = dataset)
summary(hopperBIO05_2021_2040)

draw(hopperBIO05_2021_2040, select = "s(BIO05_2021_2040)") +
  geom_hline(yintercept = 0, linetype = "dashed")

grasshopperBIO10_2021_2040 <- gam(grasshopper ~ dev + s(BIO10_2021_2040,  k = 4) + ObservationMethod,
                                  family = binomial(link = "logit"),
                                  data = dataset)
summary(grasshopperBIO10_2021_2040)
draw(grasshopperBIO10_2021_2040, select = "s(BIO10_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")


grasshopperBIO05_2021_2040 <- gam(grasshopper ~ dev + s(BIO05_2021_2040,  k = 4) + ObservationMethod,
                                  family = binomial(link = "logit"),
                                  data = dataset)
summary(grasshopperBIO05_2021_2040)

draw(grasshopperBIO05_2021_2040, select = "s(BIO05_2021_2040)") +
  geom_hline(yintercept = 0, linetype = "dashed")



beetleBIO10_2021_2040 <- gam(beetle ~ dev + s(BIO10_2021_2040,  k = 4) + ObservationMethod,
                             family = binomial(link = "logit"),
                             data = dataset)
summary(beetleBIO10_2021_2040)
draw(beetleBIO10_2021_2040, select = "s(BIO10_2021_2040)")+
  geom_hline(yintercept = 0, linetype = "dashed")


beetleBIO05_2021_2040 <- gam(beetle ~ dev + s(BIO05_2021_2040,  k = 4) + ObservationMethod,
                             family = binomial(link = "logit"),
                             data = dataset)
summary(beetleBIO05_2021_2040)

draw(beetleBIO05_2021_2040, select = "s(BIO05_2021_2040)") +
  geom_hline(yintercept = 0, linetype = "dashed")





p_cat10_2021_2040 <- plot_gam_response(catBIO10_2021_2040, dataset,"BIO10_2021_2040", "Caterpillar")
p_spi10_2021_2040 <- plot_gam_response(spiBIO10_2021_2040, dataset, "BIO10_2021_2040","Spider")
p_ant10_2021_2040 <- plot_gam_response(antBIO10_2021_2040, dataset,"BIO10_2021_2040", "Ant")
p_true10_2021_2040 <- plot_gam_response(truebugBIO10_2021_2040, dataset,"BIO10_2021_2040", "True bug")
p_hopper10_2021_2040 <- plot_gam_response(hopperBIO10_2021_2040, dataset,"BIO10_2021_2040", "Leafhopper")
p_beetle10_2021_2040 <- plot_gam_response(beetleBIO10_2021_2040, dataset,"BIO10_2021_2040", "Beetle")
p_orthopteran10_2021_2040 <- plot_gam_response(grasshopperBIO10_2021_2040, dataset, "BIO10_2021_2040","Orthoptera")

p_cat05_2021_2040 <- plot_gam_response(catBIO05_2021_2040, dataset, "BIO05_2021_2040", "Caterpillar")
p_spi05_2021_2040 <- plot_gam_response(spiBIO05_2021_2040, dataset, "BIO05_2021_2040","Spider")
p_ant05_2021_2040 <- plot_gam_response(antBIO05_2021_2040, dataset, "BIO05_2021_2040","Ant")
p_true05_2021_2040 <- plot_gam_response(truebugBIO05_2021_2040, dataset, "BIO05_2021_2040","True bug")
p_hopper05_2021_2040 <- plot_gam_response(hopperBIO05_2021_2040, dataset, "BIO05_2021_2040","Leafhopper")
p_beetle05_2021_2040 <- plot_gam_response(beetleBIO05_2021_2040, dataset, "BIO05_2021_2040","Beetle")
p_orthopteran05_2021_2040 <- plot_gam_response(grasshopperBIO05_2021_2040, dataset, "BIO05_2021_2040","Orthoptera")

BIO10_2021_2040_panel =
  wrap_plots(
    p_cat10_2021_2040,
    p_orthopteran10_2021_2040,
    p_spi10_2021_2040,
    p_ant10_2021_2040,
    p_beetle10_2021_2040,
    p_hopper10_2021_2040,
    p_true10_2021_2040,
    ncol = 2
  ) +
  plot_annotation(
    title = "Smooth effects of BIO10_2021_2040 on arthropod occurrence")

BIO10_2021_2040_panel


BIO05_2021_2040_panel =
  wrap_plots(
    p_cat05_2021_2040,
    p_orthopteran05_2021_2040,
    p_spi05_2021_2040,
    p_ant05_2021_2040,
    p_beetle05_2021_2040,
    p_hopper05_2021_2040,
    p_true05_2021_2040,
    ncol = 2
  ) +
  plot_annotation(
    title = "Smooth effects of BIO05_2021_2040 on arthropod occurrence")

BIO05_2021_2040_panel







