library(tidyverse)
library(jsonlite)

### 1. Read in CC! data files
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))

 

plant = read.csv("data/ccPlants.csv")
plantOrigin = read.csv("data/plantOrigin.csv") %>% filter(!is.na(plantOrigin))



  
fullDataset %>% 
  left_join( # join to a summary of useful plant origin data that is summarised
plant %>%  
  group_by(sciName, plantOrigin) %>% 
  summarise(n = n()) %>% 
  select(-n) %>% 
  filter(!is.na(plantOrigin)), by = c("Species" = "PlantSpecies")) %>% 
  select(PlantSpecies, plantOrigin)
  
  
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>% 
  group_by(Name, ObservationMethod, plantOrigin) %>%
  summarise(surv = n_distinct(ID)) %>% 
  view()
  
  
fullDataset %>% 
  select(Species, PlantSpecies) %>% 
  view()


fullDataset %>% 
  group_by(Name, Code) %>% 
  summarise(nIDs = n_distinct(ID))

  
  summarize(caterpillar = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
            spider = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
            beetle = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
            truebug = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
            hopper = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
            ant = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0)) %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(caterpillar_prop = mean(caterpillar),
            spider_prop = mean(spider),
            beetle_prop = mean(beetle),
            truebug_prop = mean(truebug),
            hopper_prop  = mean(hopper),
            ant_prop = mean(ant),
            Trials = n())


