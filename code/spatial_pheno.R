minSurveys = 50
Pheno_julianWindow = 140:213
minSiteYear =10



Julian_Arthropod = fullDataset %>%
  filter(julianday %in% Pheno_julianWindow,
         Longitude > -100,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         !Year %in% yearException) %>% 
  group_by(Name, ID, ObservationMethod, Year, julianweek) %>%
  summarize(
    caterpillar   = ifelse(sum(Group == 'caterpillar', na.rm = TRUE) > 0, 1, 0),
    spider        = ifelse(sum(Group == 'spider', na.rm = TRUE) > 0, 1, 0),
    beetle        = ifelse(sum(Group == 'beetle', na.rm = TRUE) > 0, 1, 0),
    truebug       = ifelse(sum(Group == 'truebugs', na.rm = TRUE) > 0, 1, 0),
    hopper        = ifelse(sum(Group == 'leafhopper', na.rm = TRUE) > 0, 1, 0),
    ant           = ifelse(sum(Group == 'ant', na.rm = TRUE) > 0, 1, 0),
    fly           = ifelse(sum(Group == 'fly', na.rm = TRUE) > 0, 1, 0),
    grasshopper   = ifelse(sum(Group == 'grasshopper', na.rm = TRUE) > 0, 1, 0),
    daddylonglegs = ifelse(sum(Group == 'daddylonglegs', na.rm = TRUE) > 0, 1, 0)) %>% 
  group_by(Name, ObservationMethod, Year, julianweek) %>%
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
    surveyNum       = n_distinct(ID)) 

# not so sure: Here it theoretically makes sense that peak period should not differ by what 
# method of data collection was employed. We may not say same for the actual occurence (height)


# good weeks of survey stays. 
goodsurv = Julian_Arthropod %>% 
  filter(surveyNum >=  10)  # each site/year/week should have at least 10 surveys.


Site.julian = goodsurv %>% 
  group_by(Name, ObservationMethod, Year) %>% 
  summarise(nJulianWeek= n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= 5)


JuliandData = inner_join(goodsurv, Site.julian, 
                         by = c("Name", "ObservationMethod", "Year"))

siteyearSpatial =  JuliandData %>% # Good sites with 3 or more years of survey
  group_by(Name, ObservationMethod) %>% 
  summarise(nYear= n_distinct(Year)) %>% 
  filter(nYear >= 1) # Change this if you want to use site that has multiple year survey


JuliSiteDataSpatial = inner_join(JuliandData, siteyearSpatial,
                          by = c("Name", "ObservationMethod"))

JuliSiteDataSpatial %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(nYear = n_distinct(Year)) # all good!


# Return only the maximum value of arthropod group for each year
MaxArthropodSpatial = JuliSiteDataSpatial %>% 
  group_by(Name, ObservationMethod, Year, nJulianWeek, nYear) %>% 
  summarise(across(where(is.numeric)  & !c(julianweek, surveyNum),
                   ~ max(.x, na.rm = TRUE), .names = "{col}")) 



# Caterpillar
JuliSiteData_CaterpillarSpatial= left_join(MaxArthropodSpatial, 
                                           JuliSiteDataSpatial %>% select(Name, ObservationMethod, 
                                                            Year, julianweek, Caterpillar),
                                    by = c("Name", "ObservationMethod", "Year", "Caterpillar")) %>% 
  select(Name, ObservationMethod, Year, julianweek, Caterpillar, nJulianWeek) %>% 
  rename(maxjulianweek = julianweek,
         maxOcc = Caterpillar)




# Spider
JuliSiteData_SpiderSpatial= left_join(MaxArthropodSpatial, 
                               JuliSiteDataSpatial %>% select(Name, ObservationMethod, 
                                                       Year, julianweek, Spider),
                               by = c("Name", "ObservationMethod", "Year", "Spider")) %>% 
  select(Name, ObservationMethod, Year, julianweek, Spider, nJulianWeek) %>% 
  rename(maxjulianweek = julianweek,
         maxOcc = Spider)


# Ant

JuliSiteData_AntSpatial= left_join(MaxArthropodSpatial, 
                            JuliSiteDataSpatial %>% select(Name, ObservationMethod, 
                                                    Year, julianweek, Ant),
                            by = c("Name", "ObservationMethod", "Year", "Ant")) %>% 
  select(Name, ObservationMethod, Year, julianweek, Ant, nJulianWeek) %>% 
  rename(maxjulianweek = julianweek,
         maxOcc = Ant)



spatialPheno = rbind(JuliSiteData_AntSpatial %>% mutate(Group =  "Ant"),
                     JuliSiteData_SpiderSpatial %>% mutate(Group = "Spider"),
                    JuliSiteData_CaterpillarSpatial %>% mutate(Group = "Caterpillar")) %>% 
  as.data.frame() %>% 
  left_join(site.info, by = c("Name", "ObservationMethod"))

spatialPheno %>% 
  filter(if_any(everything(), is.na)) 
        # one location has too much NAs, in three rows, so delete.

spatialPheno = spatialPheno %>% 
  filter(!if_any(everything(), is.na))


GoodSpatialYearObserv  = spatialPheno %>% 
  group_by(Group, Year, ObservationMethod) %>% 
  summarise(nSite = n_distinct(Name)) %>% 
  filter(nSite >= minSiteYear) %>% data.frame()  

# use join function as filter for only the good site

spatialPheno = spatialPheno %>% 
  right_join(GoodSpatialYearObserv, by = c("Year", "ObservationMethod", "Group"))

  spatialPheno %>% 
  filter(Group == "Caterpillar") %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(aes(x = Latitude, y = maxjulianweek, colour = Year, group = Year)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_jitter(width = 0.6, height = 0.3, alpha = 0.7) +
  theme_minimal()
  
  
  


spatialPheno %>% 
  filter(Group == "Ant",
         Year %in% c(2017:max(spatialPheno$Year))) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(aes(x = Latitude, y = maxjulianweek, colour = Year, group = Year)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_jitter(width = 0.6, height = 0.3, alpha = 0.7) +
  theme_minimal()


spatialPheno %>% 
  filter(Group == "Spider",
         Year %in% c(2017: max(spatialPheno$Year))) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(aes(x = Latitude, y = maxjulianweek, colour = Year, group = Year)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(aes(colour = Year))+
  theme_minimal()


# Fit flexible lines

spatialPheno %>% 
  filter(Group == "Caterpillar") %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(aes(x = Latitude, y = maxjulianweek, colour = Year, group = Year)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE) +
  geom_jitter(width = 0.6, height = 0.4, alpha = 0.7, (aes(size = surveyNum))) +
  theme_minimal()


spatialPheno %>% 
  filter(Group == "Caterpillar") %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(aes(x = Latitude, y = maxjulianweek)) +
  facet_wrap(~Year)+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE) +
  geom_jitter(width = 0.6, height = 0.4, alpha = 0.7, (aes(size = surveyNum))) +
  theme_minimal()




# ---------------- Caterpillars
spatialPheno %>% 
  filter(Group == "Caterpillar",
         Latitude >= 33 & Latitude <= 45 ) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(aes(x = Latitude, y = maxjulianweek)) +
  facet_wrap(~Year, nrow = 4, ncol = 2, dir = "v") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, alpha = 0.2,
              colour = "steelblue", linewidth = 0.5) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "darkorange", linewidth = 0.9) +
  geom_jitter(aes(size = surveyNum), 
              width = 0.6, height = 0.4, alpha = 0.5, colour = "darkorange") +
  labs( title = "Caterpillar", subtitle = "Peak Julian Period across site")+
  theme_minimal()


# --------------------- Spiders
spatialPheno %>% 
  filter(Group == "Spider",
         Latitude >= 33 & Latitude <= 45 ) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(aes(x = Latitude, y = maxjulianweek)) +
  facet_wrap(~Year, nrow = 4, ncol = 2, dir = "v") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, alpha = 0.2,
              colour = "steelblue", linewidth = 0.5) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "darkorange", linewidth = 0.9) +
  geom_jitter(aes(size = surveyNum), 
              width = 0.6, height = 0.4, alpha = 0.5, colour = "darkorange") +
  labs( title = "Spider", subtitle = "Peak Julian Period across site")+
  theme_minimal()



# ---------------------- Ants
spatialPheno %>% 
  filter(Group == "Ant",
         Latitude >= 33 & Latitude <= 45 ) %>% 
  mutate(Year = factor(Year)) %>% 
  ggplot(aes(x = Latitude, y = maxjulianweek)) +
  facet_wrap(~Year, nrow = 4, ncol = 2, dir = "v") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, alpha = 0.2,
              colour = "steelblue", linewidth = 0.5) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "darkorange", linewidth = 0.9) +
  geom_jitter(aes(size = surveyNum), 
              width = 0.6, height = 0.4, alpha = 0.5, colour = "darkorange") +
  labs( title = "Ant", subtitle = "Peak Julian Period across site")+
  theme_minimal()


# write.csv(sites, file = "data/urbanDev_forest.csv")
