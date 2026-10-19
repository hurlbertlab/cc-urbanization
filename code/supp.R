# You need to run all script in the main analysis before this one.


library(tidyverse)
library(jsonlite)
library(viridis)
# (1) Read in latest Caterpillars Count! raw dataset from the caterpillars-analysis-public repo----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))  %>% 
  filter(Year <= 2025)

# Read in plant files (up-to-date versions in caterpillars-count-data repo)
officialPlantList = read.csv('data/officialPlantList2024-10-04.csv')
inferredPlantNames = read.csv('data/inferredPlantNames_2024-08-20.csv')
plantOrigin = read.csv('data/plant_origin_status.csv') %>%
  dplyr::select(scientificName, nativeStatus, plantOrigin)

sites = read.csv('data/sites.csv')

ccPlants = fullDataset %>%
  left_join(inferredPlantNames[, c('PlantFK', 'InferredSciName', 'NameConfidence')], by = 'PlantFK') %>%
  left_join(plantOrigin, by = c('sciName' = 'scientificName'))

table(ccPlants$plantOrigin)
sum(is.na(ccPlants$plantOrigin)) # 6931 branch surveys we do not know their plant origin status

# Plants with unknown status and the number of sites they appear
ccPlants %>% 
  filter(is.na(plantOrigin)) %>%
  dplyr::select(Name, PlantSpecies) %>% 
  group_by(PlantSpecies) %>% 
  summarise(nSite = n_distinct(Name)) %>% 
  data.frame()


ccplants_dev = ccPlants %>% 
  left_join(sites %>% dplyr::select(-Name), by = c("Latitude", "Longitude"))%>%
  filter(julianday %in% julianWindow,
         Longitude > -100,
         WetLeaves == 0) 


good_ccplants  = ccplants_dev %>%
  filter(julianday %in% julianWindow,
         Longitude > -100,
         WetLeaves == 0) %>%
  group_by(Latitude, Longitude, ObservationMethod) %>%
  summarize(nSurvs = n_distinct(ID)) %>%
  filter(nSurvs >= minSurveys) 


ccplantDataset = ccplants_dev %>% 
  inner_join(good_ccplants, 
             by = c("Latitude", "Longitude", "ObservationMethod"))

ccdata = ccplantDataset %>%
  filter(
    !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>%
  group_by(Name, Latitude, Longitude, dev, nSurvs, ID, plantOrigin, ObservationMethod) %>%
  summarize(caterpillar = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
            spider = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
            beetle = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
            truebug = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
            hopper = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
            ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0),
            grasshopper = ifelse(sum(Group == "grasshopper", na.rm = TRUE) > 0, 1, 0),
            fly = ifelse(sum(Group == "fly", na.rm = TRUE) > 0, 1, 0),
            daddylonglegs = ifelse(sum(Group == "daddylonglegs", na.rm = TRUE) > 0, 1, 0)) %>% 
  filter(!is.na(dev),
         !is.na(Latitude))
  

scale(ccdata$dev)

cat.Dev.Latitude.origin = glm(caterpillar ~ dev + Latitude + dev*Latitude + ObservationMethod +plantOrigin , 
                       data = ccdata, family = "binomial")
summary(cat.Dev.Latitude.origin)






# datataset filtered for 2018 - 2025

yearWindow = 2018:2025

dataset2018.2025 = inner_join(
fullDataset %>%
  filter(
    julianday %in% julianWindow,
    WetLeaves == 0,
    Longitude > -100,
    Year  %in% c(yearWindow),
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
goodSites,
by = c("Name", "ObservationMethod")) %>% 
  left_join(sites, by = c("Name")) 


cat.Dev.Latitude2018.2025 = glm(caterpillar ~ dev * Latitude + ObservationMethod, 
                                data = dataset2018.2025, family = "binomial")

spi.Dev.Latitude2018.2025 = glm(spider ~ dev * Latitude + ObservationMethod, 
                                data = dataset2018.2025, family = "binomial")

beet.Dev.Latitude2018.2025 = glm(beetle ~ dev * Latitude + ObservationMethod, 
                                 data = dataset2018.2025, family = "binomial")

hop.Dev.Latitude2018.2025 = glm(hopper ~ dev * Latitude + ObservationMethod, 
                                data = dataset2018.2025, family = "binomial")

bug.Dev.Latitude2018.2025 = glm(truebug ~ dev * Latitude + ObservationMethod, 
                                data = dataset2018.2025, family = "binomial")

ant.Dev.Latitude2018.2025 = glm(ant ~ dev * Latitude + ObservationMethod, 
                                data = dataset2018.2025, family = "binomial")

fly.Dev.Latitude2018.2025 = glm(fly ~ dev * Latitude + ObservationMethod, 
                                data = dataset2018.2025, family = "binomial")

grasshopper.Dev.Latitude2018.2025 = glm(grasshopper ~ dev * Latitude + ObservationMethod, 
                                        data = dataset2018.2025, family = "binomial")

daddylonglegs.Dev.Latitude2018.2025 = glm(daddylonglegs ~ dev * Latitude + ObservationMethod, 
                                          data = dataset2018.2025, family = "binomial")



prop.devLatOutput2018.2025 = data.frame(rbind(summary(cat.Dev.Latitude2018.2025)$coefficients, 
                                     summary(spi.Dev.Latitude2018.2025)$coefficients, 
                                     summary(beet.Dev.Latitude2018.2025)$coefficients,
                                     summary(bug.Dev.Latitude2018.2025)$coefficients, 
                                     summary(hop.Dev.Latitude2018.2025)$coefficients,
                                     summary(ant.Dev.Latitude2018.2025)$coefficients,
                                     summary(fly.Dev.Latitude2018.2025)$coefficients,
                                     summary(grasshopper.Dev.Latitude2018.2025)$coefficients,
                                     summary(daddylonglegs.Dev.Latitude2018.2025)$coefficients))
prop.devLatOutput2018.2025$term = rep(c('Intercept', 'dev', 'Latitude', 'dev*Latitude', 'Method'),
                             times = 9)
prop.devLatOutput2018.2025$Group = rep(c('Caterpillar', 'Spider', 'Beetle', 
                                'Leafhopper', 'Truebugs', 'Ant',
                                'Fly', 'Grasshopper', 'Daddylonglegs'), each = 5)

rownames(prop.devLatOutput2018.2025) = NULL


simDev_caterpillar2018.2025 = sim_slopes(cat.Dev.Latitude2018.2025, pred = dev, modx = Latitude,
                                         johnson_neyman = FALSE, digits = 4)

simDev_beetle2018.2025 = sim_slopes(beet.Dev.Latitude2018.2025, pred = dev, modx = Latitude,
                                    johnson_neyman = FALSE, digits = 4)

simDev_spider2018.2025 = sim_slopes(spi.Dev.Latitude2018.2025, pred = dev, modx = Latitude,
                                    johnson_neyman = FALSE, digits = 4)

simDev_hopper2018.2025 = sim_slopes(hop.Dev.Latitude2018.2025, pred = dev, modx = Latitude,
                                    johnson_neyman = FALSE, digits = 4)

simDev_ant2018.2025 = sim_slopes(ant.Dev.Latitude2018.2025, pred = dev, modx = Latitude,
                                 johnson_neyman = FALSE, digits = 4)

simDev_truebug2018.2025 = sim_slopes(bug.Dev.Latitude2018.2025, pred = dev, modx = Latitude,
                                     johnson_neyman = FALSE, digits = 4)

simDev_fly2018.2025 = sim_slopes(fly.Dev.Latitude2018.2025, pred = dev, modx = Latitude,
                                 johnson_neyman = FALSE, digits = 4)

simDev_grasshopper2018.2025 = sim_slopes(grasshopper.Dev.Latitude2018.2025, pred = dev, modx = Latitude,
                                         johnson_neyman = FALSE, digits = 4)

simDev_daddylonglegs2018.2025 = sim_slopes(daddylonglegs.Dev.Latitude2018.2025, pred = dev, modx = Latitude,
                                           johnson_neyman = FALSE, digits = 4)


Dev_sims2018.2025 =
  bind_rows(
    simDev_caterpillar2018.2025$slopes     %>% mutate(Group = "caterpillar"),
    simDev_beetle2018.2025$slopes          %>% mutate(Group = "beetle"),
    simDev_ant2018.2025$slopes             %>% mutate(Group = "ant"),
    simDev_hopper2018.2025$slopes          %>% mutate(Group = "hopper"),
    simDev_truebug2018.2025$slopes         %>% mutate(Group = "truebug"),
    simDev_spider2018.2025$slopes          %>% mutate(Group = "spider"),
    simDev_fly2018.2025$slopes             %>% mutate(Group = "fly"),
    simDev_grasshopper2018.2025$slopes     %>% mutate(Group = "grasshopper"),
    simDev_daddylonglegs2018.2025$slopes   %>% mutate(Group = "daddylonglegs")
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
  )  %>% 
  mutate(
    LatLevel = factor(
      rep(c("Low", "Mid", "High"), times = 9),
      levels = c("Low", "Mid", "High")
    )
  )


  Dev_sims %>% 
    mutate(
      LatLevel = factor(
        rep(c("Low", "Mid", "High"), times = 9),
        levels = c("Low", "Mid", "High")
      )
    ) %>%
  dplyr::select(Group, LatLevel, Est = Est.) %>%
  inner_join(
    Dev_sims2018.2025 %>%
      dplyr::select(Group, LatLevel, Est_2018_2025 = Est.),
    by = c("Group", "LatLevel")) %>% 
ggplot(
       aes(x = Est, y = Est_2018_2025, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~LatLevel) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  labs(
    x = "Slope (Full Dataset)",
    y = "Slope (2018–2025)",
    title = "Comparison of Urbanization Slopes Across Datasets"
  )


  ccdata %>% 
  group_by(Latitude, plantOrigin) %>% 
  summarise(nOrigin = n()) %>% 
    data.frame()

  ccdata %>% 
    group_by(Latitude, dev) %>% 
    summarise(nSurvs = mean(nSurvs),
      alienRatio = sum(plantOrigin == "alien", na.rm = TRUE) /
        sum(plantOrigin == "native", na.rm = TRUE)
    ) %>% 
  ggplot(aes(x= Latitude, y = alienRatio, colour = dev))+
    scale_color_gradientn(
      colours = c("green", "lightyellow", "blue"),
      name = "% Development"
    ) +
    geom_jitter(aes(size = nSurvs),
                width = 0.15,    
                height = 0.005,  
                alpha = 0.8) +
    guides(size = "none")+
    labs(subtitle = "alienRatio = alien/native per site")+
    theme_bw()

  ccdata %>% 
    group_by(Latitude, dev) %>% 
    summarise(
      nSurvs = mean(nSurvs),
      alienProportion = sum(plantOrigin == "alien", na.rm = TRUE) /
        sum(plantOrigin %in% c("alien", "native"), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    ggplot(aes(x = dev, y = alienProportion, colour = Latitude)) +
    scale_color_gradientn(
      colours = c("orange", "green", "blue"),
      name = "Latitude (°)"
    ) +
    geom_jitter(
      aes(size = nSurvs),
      width = 0.5,
      height = 0.005,
      alpha = 0.8
    ) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      se = TRUE,
      colour = "black",
      inherit.aes = TRUE
    ) +
    guides(size = "none") +
    labs(
      x = "% Urban Development",
      y = "Proportion of Alien Plants",
      subtitle = ""
    ) +
    theme_bw(base_size = 20)
  
  
  
  
  
  
  plotData <- ccdata %>% 
    group_by(Latitude, dev) %>% 
    summarise(
      nSurvs = mean(nSurvs),
      alienProportion = sum(plantOrigin == "alien", na.rm = TRUE) /
        sum(plantOrigin %in% c("alien", "native"), na.rm = TRUE),
      .groups = "drop"
    )
  
  # Fit overall linear regression
  model <- lm(alienProportion ~ dev, data = plotData)
  
  # Extract statistics
  r2 <- summary(model)$r.squared
  slope <- coef(model)[["dev"]]
  p <- summary(model)$coefficients["dev", "Pr(>|t|)"]
  
  # Create annotation
  statsLabel <- sprintf(
    "Slope = %.4f\nR² = %.3f\np = %.3g",
    slope, r2, p
  )
  
 alienProp =  ggplot(plotData, aes(x = dev, y = alienProportion, colour = Latitude)) +
    scale_color_gradientn(
      colours = c("#F2B134", "#56B4E9", "#1D4E89"),
      name = "Latitude (°)"
    ) +
    geom_jitter(
      aes(size = nSurvs),
      width = 0.5,
      height = 0.005,
      alpha = 0.8
    ) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      se = TRUE,
      colour = "black"
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = statsLabel,
      hjust = 1.1, vjust = 1.2,
      colour = "black",
      size = 5
    ) +
    guides(size = "none") +
    labs(
      x = "% Urban Development",
      y = "Proportion of Alien Plants",
      subtitle = ""
    ) +
    theme_bw(base_size = 20)
 
 ggsave(
   "images/alienProp.jpeg",
   plot = alienProp,  width = 6.5,
   height = 5.2,  units = "in", dpi = 600)
 
 
  
statusPersite =   ccdata %>% 
    group_by(Latitude, dev) %>% 
    summarise(nSurvs = mean(nSurvs),
              alienRatio = sum(plantOrigin == "alien", na.rm = TRUE) /
                sum(plantOrigin == "native", na.rm = TRUE),
              alienProp = sum(plantOrigin == "alien", na.rm = TRUE)/ nSurvs,
              nativeProp =  sum(plantOrigin == "native", na.rm = TRUE)/ nSurvs) %>% 
    as.data.frame()
  
  
alien.dev = glm(alienProp  ~ dev * Latitude, 
                   data = statusPersite, weights = nSurvs,  family = "binomial")
summary(alien.dev)
