# Here begins my investigation into how plant composition influences arthropod community composition; abundance and biodiversity


library(geosphere)
library(tidyverse)
require(vegan)



prop_hellinger <- decostand(prop_data[,3:11],
                          method = "hellinger")

propHel_dist <- dist(prop_hellinger, method = "euclidean")
lat_long_dist <- dist(prop_data[,c("Latitude", "Longitude")],
                      method = "euclidean")


coords <- prop_data[, c("Longitude", "Latitude")]
geo_dist <- distm(coords, fun = distHaversine)
geo_dist <- as.dist(geo_dist)


prop_mantel <- mantel(propHel_dist, geo_dist, 
                      method = "pearson", permutations = 999)

###################################################################################



prop_dissimilarity <- as.matrix(vegdist(prop_hellinger, method = "euclidean"))
geo_dist2 <- as.matrix(dist(coords))

prop.dissimilarity_vec <- prop_dissimilarity[lower.tri(prop_dissimilarity)]
geo_distance_vec <- geo_dist2[lower.tri(geo_dist2)]


prop_decay_df <- data.frame(
  Distance = geo_distance_vec,
  Dissimilarity = prop.dissimilarity_vec
)

ggplot(prop_decay_df, aes(x = Distance, y = Dissimilarity)) +
  geom_point(color = "black", size = 2, alpha = 0.7) +   
  geom_smooth(method = "lm", color = "red", se = TRUE) +   
  labs(
    title = "Relationship between geographic distance and site-level arthropod occurence composition dissimilariy, using presence-absence data",
    x = "Geographic Distance",
    y = "Euclidean distance"
  ) +
  theme_minimal()


summary(lm(Dissimilarity ~ Distance, data = prop_decay_df))
# Note that theoretically, the P-value from using lm() on dissimilarity computation, in this case, cannot be valid. This is because pair-wise correlations are not independent. That is what a mantel test corrects for by doing permutations.


#######################################################################################

# plants species that have been sampled at least 200 times.

plant_top_200 = fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>% 
  group_by(PlantSpecies) %>% 
  summarise(Count = n()) %>% 
  filter(Count >= 200, 
         PlantSpecies != "N/A") %>% 
  arrange(desc(Count)) %>% data.frame()



Plant_top_200.site = fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         PlantSpecies %in% plant_top_200$PlantSpecies) %>% 
  group_by(PlantSpecies) %>%
  summarise(nSite = n_distinct(Name),
            nSurv = n(),
            Site.Surv.score = nSite * nSurv) %>%
  arrange(desc(Site.Surv.score)) 
  

minSurveys

fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         PlantSpecies != "N/A") %>% 
  group_by(Name, PlantSpecies) %>% 
  summarise(nSurv = n()) %>% 
  mutate(Count = ifelse(nSurv > 20, 1, 0)) %>% 
  pivot_wider(names_from = PlantSpecies,
              values_from = Count) %>% 
  group_by(Name) %>% 
  summarise(total_plants = sum(across(where(is.numeric) & !matches("nSurv")), 
                               na.rm = TRUE),
            nSurv = sum(nSurv),
           ) %>% 
  arrange(desc(total_plants)) %>% 
  data.frame()

# Each plant is surveyed more than 20 times, and each site must have at lest 5 of these >20 surveyed plants.

site.plantSurv20 = fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         PlantSpecies != "N/A") %>% 
  group_by(Name, PlantSpecies) %>% 
  summarise(nSurv = n(), .groups = "drop") %>% 
  mutate(Count = ifelse(nSurv > 20, 1, 0)) %>% 
  pivot_wider(
    names_from = PlantSpecies,
    values_from = Count,
    values_fill = list(Count = 0)
  ) %>%
  group_by(Name) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), 
            .groups = "drop") %>%
  mutate(total_plants = rowSums(select(., where(is.numeric) & !matches("nSurv")), 
                                na.rm = TRUE)) %>%
  filter(total_plants >= 5) %>% 
  arrange(desc(total_plants)) 



Site20_jac_dist <- vegdist(
  site.plantSurv20 %>% select(-c(nSurv, total_plants, Name)), 
  method = "jaccard", 
  binary = TRUE
)

Site20_jac_sim <- 1 - as.matrix(Site20_jac_dist)


siteJacSim <- as.data.frame(as.matrix(Site20_jac_sim))%>%
  mutate(Total_Jac = rowSums(.)-1)



SiteTopJackSim  <- cbind(site.plantSurv20 %>% select(Name, nSurv, total_plants),
                         Total_Jac = siteJacSim$Total_Jac) %>% 
  as.data.frame() %>% 
  arrange(desc(Total_Jac))

############################### -----Genus level--------------------------------###############


# Plant Genus that have been sampled 200 times or more

plantG_top_200 = fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>% 
  group_by(plantGenus) %>% 
  summarise(Count = n()) %>% 
  filter(Count >= 200, 
         plantGenus != "N/A") %>% 
  arrange(desc(Count)) %>% data.frame()


fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         plantGenus != "N/A") %>% 
  group_by(Name, plantGenus) %>% 
  summarise(nSurv = n()) %>% 
  mutate(Count = ifelse(nSurv > 50, 1, 0)) %>% 
  pivot_wider(names_from = plantGenus,
              values_from = Count) %>% 
  group_by(Name) %>% 
  summarise(total_plants = sum(across(where(is.numeric) & !matches("nSurv")), 
                               na.rm = TRUE),
            nSurv = sum(nSurv),
  ) %>% 
  arrange(desc(total_plants)) %>% 
  data.frame()



site.plantGenusSurv = fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         plantGenus != "N/A") %>% 
  group_by(Name, ObservationMethod, plantGenus) %>% 
  summarise(nSurv = n(), .groups = "drop") %>% 
  mutate(Count = ifelse(nSurv >= 200, 1, 0)) %>% 
  pivot_wider(
    names_from = plantGenus,
    values_from = Count,
    values_fill = list(Count = 0)
  ) %>%
  group_by(Name, ObservationMethod ) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), 
            .groups = "drop") %>%
  mutate(total_plants = rowSums(select(., where(is.numeric) & !matches("nSurv")), 
                                na.rm = TRUE)) %>%
  filter(total_plants >= 1)
 


Goodplant = site.plantGenusSurv %>% 
  pivot_longer(cols = -c("Name", "ObservationMethod", "nSurv", "total_plants"),
               names_to = "PlantGenus",
               values_to = "Count") %>% 
  group_by(PlantGenus) %>% 
  summarise(Count = n_distinct(PlantGenus))


GoodBugs.Plants = fullDataset %>%
  filter(Name %in% site.plantGenusSurv$Name, 
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         plantGenus %in% Goodplant$PlantGenus,
         Group %in% c('caterpillar', 'spider', 'ant', 'leafhopper', 'beetle', 
                      'truebugs', 'fly', 'grasshopper', 
                      'daddylonglegs', 'aphid')) %>%
  group_by(Name, ID, ObservationMethod) %>%
  summarize(
    caterpillar   = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
    spider        = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
    beetle        = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
    truebug       = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
    hopper        = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
    ant           = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0),
    fly           = ifelse(sum(Group == 'fly', na.rm = TRUE) > 0, 1, 0),
    grasshopper   = ifelse(sum(Group == 'grasshopper', na.rm = TRUE) > 0, 1, 0),
    daddylonglegs = ifelse(sum(Group == 'daddylonglegs', na.rm = TRUE) > 0, 1, 0),
    #  aphid         =ifelse(sum(Group == 'aphid', na.rm = TRUE) > 0, 1, 0) high I.D misclass
  ) %>%
  group_by(Name, ObservationMethod) %>%
  summarise(
    Caterpillar   = sum(caterpillar)   / n_distinct(ID),
    Spider        = sum(spider)        / n_distinct(ID),
    Beetle        = sum(beetle)        / n_distinct(ID),
    Truebug       = sum(truebug)       / n_distinct(ID),
    Hopper        = sum(hopper)        / n_distinct(ID),
    Ant           = sum(ant)           / n_distinct(ID),
    fly           = sum(fly)           / n_distinct(ID),
    Grasshopper   = sum(grasshopper)   / n_distinct(ID),
    Daddylongleg = sum(daddylonglegs) / n_distinct(ID),
    #Aphid         = sum(aphid)         / n_distinct(ID),
    surveyNum       = n_distinct(ID))


  
plant.nSurv = fullDataset %>%
  filter(plantGenus %in% Goodplant$PlantGenus,
         Name %in% site.plantGenusSurv$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>% 
  group_by(Name, ObservationMethod, plantGenus) %>% 
  summarise(nSurv = n()) %>% 
  pivot_wider(names_from = plantGenus,
              values_from = nSurv) %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))


# some sites observation 

#  Herbivory EDA ----



fullDataset %>% 
  filter(!HerbivoryScore %in% c(-128, -1)) %>% 
  group_by(Name, HerbivoryScore) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% data.frame()


fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge|NC Botanical Garden"),
    Year == 2024
  ) %>% 
  group_by(Name, julianweek, HerbivoryScore) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  ggplot(aes(y = count, x = julianweek)) +
  geom_point() +
  stat_smooth()+
  facet_grid(Name ~ HerbivoryScore)+
  theme_bw()



fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge|NC Botanical Garden"),
  ) %>% 
  group_by(Name, julianweek, HerbivoryScore) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  ggplot(aes(y = count, x = julianweek)) +
  geom_point() +
  stat_smooth()+
  facet_grid(Name ~ HerbivoryScore)+
  theme_bw()

fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge|NC Botanical Garden|UNC Chapel Hill Campus"),
    Year >= 2021
  ) %>% 
  group_by(Name, Year, julianweek, ID) %>% 
  summarise(Herb  = mean(HerbivoryScore)) %>% 
  group_by(Name, Year, julianweek, Herb) %>% 
  summarise(Herbcount = n (),
            nSurv = n_distinct(ID)) %>% 
  ggplot(aes(x = julianweek, y = Herbcount, color = as.factor(Year), group = Year)) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("#4E79A7", "#F28E2B", "#E15759", "darkgrey", "#59A14F"))+
  facet_grid(Name ~ as.factor(Herb), scales = "free_y") +
  theme_bw() +
  labs(color = "Year")





fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge|NC Botanical Garden|UNC Chapel Hill Campus"),
    Year >= 2021
  ) %>% 
  group_by(Name, Year, julianweek, ID) %>% 
  summarise(Herb  = mean(HerbivoryScore)) %>% 
  group_by(Name, Year, julianweek, Herb) %>% 
  summarise(Herbcount = n (),
            nSurv = n_distinct(ID)) %>% 
  pivot_wider(names_from = Herb,
              values_from = Herbcount)


fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, "Prairie Ridge"),
    Year >= 2024
  ) %>% 
  group_by(Name, Code, julianweek) %>% 
  summarise(nSurv = n_distinct(ID),
            herb = mean(HerbivoryScore)) %>%
  ggplot(aes(y = herb, x = julianweek)) +
  geom_point()+
  facet_wrap( ~ Code)
  
  






fullDataset %>% 
  filter(
    !HerbivoryScore %in% c(-128, -1),
    str_detect(Name, regex("prairie", ignore_case = TRUE )),
    Year == 2017
  ) %>% 
  group_by(PlantSpecies, julianweek) %>% 
  summarise(HerbivoryScore = mean(HerbivoryScore)) %>% 
  ggplot(aes(y = HerbivoryScore, x = julianweek)) +
  geom_point() +
  stat_smooth()+
  facet_wrap(~PlantSpecies)
  
  
  
