library(tidyverse)
library(daymetr)
library(geosphere)
require(vegan)
  
  minSurveys = 50
  Pheno_julianWindow = 140:213
  minSurvYears = 3
  yearException = 2025 # because we do not yet have 2025 data for temperature
  TempDayWindow = 110:213   # Adjust
  
  
  # What sites have been surveying arthropods for three or more years?
fullDataset %>%
    filter(julianday %in% Pheno_julianWindow,
           Longitude > -100,
           WetLeaves == 0,
           !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
           !Year %in% yearException) %>%
    group_by(Name, ObservationMethod) %>%
    summarize(nSurvs = n_distinct(ID),
              nYears = n_distinct(Year)) %>%
    filter(nSurvs >= 40,
           nYears >= minSurvYears) %>%
    arrange(desc(nYears)) %>% as.data.frame()


# take the full data and 
# 2 . do necessary filtering,
# 3. Convert abundance data to 1 and 0 (presence and absence)
# summaries data within a given week into proportion of positive data per number of survey in that week
# Also keep information on the number of weekly survey

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


# good weeks of survey stays. 
goodsurv = Julian_Arthropod %>% 
  filter(surveyNum >=  10)  # each site/year/week should have at least 10 surveys.


Site.julian = goodsurv %>% 
  group_by(Name, ObservationMethod, Year) %>% 
  summarise(nJulianWeek= n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= 5)


JuliandData = inner_join(goodsurv, Site.julian, 
                         by = c("Name", "ObservationMethod", "Year"))

siteyear =  JuliandData %>% # Good sites with 3 or more years of survey
  group_by(Name, ObservationMethod) %>% 
  summarise(nYear= n_distinct(Year)) %>% 
  filter(nYear >= 3)


JuliSiteData = inner_join(JuliandData, siteyear,
                          by = c("Name", "ObservationMethod"))

JuliSiteData %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(nYear = n_distinct(Year)) # all good!



# Return only the maximum value of arthropod group for each year
MaxArthropod = JuliSiteData %>% 
  group_by(Name, ObservationMethod, Year, nJulianWeek, nYear) %>% 
  summarise(across(where(is.numeric)  & !c(julianweek, surveyNum),
                   ~ max(.x, na.rm = TRUE), .names = "{col}")) 

 



###################################################################################
# Caterpillar Data

JuliSiteData_Caterpillar= left_join(MaxArthropod, 
                                    JuliSiteData %>% select(Name, ObservationMethod, 
                                                            Year, julianweek, Caterpillar),
                                    by = c("Name", "ObservationMethod", "Year",
                                           "Caterpillar")) %>% 
  select(Name, ObservationMethod, Year, julianweek, Caterpillar) %>% 
  rename(maxjulianweek = julianweek,
         maxOcc = Caterpillar)

# left_join(MaxArthropod, 
#          JuliSiteData %>% select(Name, ObservationMethod, Year, julianweek, Caterpillar),
#          by = c("Name", "ObservationMethod", "Year", "Caterpillar")) %>% 
#  group_by(Name, ObservationMethod, Year) %>% 
#  summarise(count = n()) %>% 
#  arrange(desc(count)) %>% # there are four duplication, which is due to bi-modal peaks within a site/year
#  filter(count > 1)

AbnormCaterpillar= JuliSiteData_Caterpillar %>%  # site/year/occurrence max
  left_join(JuliSiteData_Caterpillar %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(MeanJulWeek = mean(maxjulianweek),  #site mean julWeek
            MeanOccurence = mean(maxOcc))) %>%  # site mean occurrence
  mutate(AnomalJulWeek = maxjulianweek  - MeanJulWeek,
         AnomalOccurence = maxOcc  - MeanOccurence)

 




 



# Spider  data: 

# Return only the maximum value of weekly caterpillar prop. of occurrence for each year and it corresponding Julian week

JuliSiteData_Spider= left_join(MaxArthropod, 
                                    JuliSiteData %>% select(Name, ObservationMethod, 
                                                            Year, julianweek, Spider),
                                    by = c("Name", "ObservationMethod", "Year", "Spider")) %>% 
  select(Name, ObservationMethod, Year, julianweek, Spider) %>% 
  rename(maxjulianweek = julianweek,
         maxOcc = Spider)

AbnormSpider= JuliSiteData_Spider %>%  # site/year/occurrence max
  left_join(JuliSiteData_Spider %>% 
              group_by(Name, ObservationMethod) %>% 
              summarise(MeanJulWeek = mean(maxjulianweek),  #site mean julWeek
                        MeanOccurence = mean(maxOcc))) %>%  # site mean occurrence
  mutate(AnomalJulWeek = maxjulianweek  - MeanJulWeek,
         AnomalOccurence = maxOcc  - MeanOccurence)

# Ant  data: 

# Return only the maximum value of weekly ant prop. of occurrence for each year and it corresponding Julian week

JuliSiteData_Ant= left_join(MaxArthropod, 
                               JuliSiteData %>% select(Name, ObservationMethod, 
                                                       Year, julianweek, Ant),
                               by = c("Name", "ObservationMethod", "Year", "Ant")) %>% 
  select(Name, ObservationMethod, Year, julianweek, Ant) %>% 
  rename(maxjulianweek = julianweek,
         maxOcc = Ant)

AbnormAnt= JuliSiteData_Ant %>%  # site/year/occurrence max
  left_join(JuliSiteData_Ant %>% 
              group_by(Name, ObservationMethod) %>% 
              summarise(MeanJulWeek = mean(maxjulianweek),  #site mean julWeek
                        MeanOccurence = mean(maxOcc))) %>%  # site mean occurrence
  mutate(AnomalJulWeek = maxjulianweek  - MeanJulWeek,
         AnomalOccurence = maxOcc  - MeanOccurence)






# #####################    Temperature data retrieval


 


 

GoodSiteYear = JuliSiteData[, c("Name","Year")] %>% 
  group_by(Name, Year) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))

GoodSiteYearLatLon = left_join(GoodSiteYear, 
                               fullDataset %>% select(Name, Year, Latitude, Longitude), 
                               by = c("Name", "Year")) %>%  # sites is from a different R script
  group_by(Name, Latitude, Longitude) %>% 
  summarise(nYears = n_distinct(Year)) %>% 
  select(-nYears) # really do not need it now

write.csv(GoodSiteYearLatLon %>% 
            rename(
  "site" = "Name",
  "lat" = "Latitude",
  "lon" = "Longitude"),
          file = "Data/PhenoAnomalySites.csv",
  row.names = FALSE)


# 
#    TempSiteData = download_daymet_batch(file_location ="Data/PhenoAnomalySites.csv",
#                                      start = min(GoodSiteYear$Year),
#                                      end = 2024, # this is the most recent available in daymetr
#                                      internal = TRUE)
#  


# The TempSiteData is a list of list. check TempSiteData[[1]]$ 
# so this function would extract the information and add the site information 

TempSiteData_clean <- lapply(TempSiteData, function(x) {
  x$data %>% mutate(site = x$site,
                    Latidue = x$latitude,
                    Longitude = x$longitude)
})

AllTempData <- bind_rows(TempSiteData_clean)




anomalPheno = rbind(AbnormAnt %>% mutate(Group =  "Ant"),
                    AbnormSpider %>% mutate(Group = "Spider"),
                    AbnormCaterpillar %>% mutate(Group = "Caterpillar")) %>% 
  as.data.frame()
 





# join temperature anomaly to phonology anomaly

TempArthropodAnomal = AllTempData %>% 
  filter(yday %in% TempDayWindow) %>% 
  group_by(site, year) %>% 
  summarise(meanTmin = mean(tmin..deg.c.),
            meanTmax = mean(tmax..deg.c.),
            meanPreci = mean(prcp..mm.day.)) %>% 
  left_join(
AllTempData %>% 
  filter(yday %in% TempDayWindow) %>% 
  group_by(site) %>% 
  summarise(AllmeanTmin = mean(tmin..deg.c.),
            AllmeanTmax = mean(tmax..deg.c.),
            AllmeanPreci = mean(prcp..mm.day.)),
by = c("site")) %>% 
  mutate(AnomalTmin = meanTmin - AllmeanTmin,
         AnomalTmax = meanTmax - AllmeanTmax,
         AnomalPreci = meanPreci - AllmeanPreci)%>% 
  right_join(anomalPheno %>%  filter(Year != "2025"), # because daymetr has no 2025 yet.
             by = c("site" = "Name", "year" = "Year")) %>% 
  left_join(GoodSiteYearLatLon, by = c("site" = "Name")) 

 

# Maximum Temperature Anomality

TempArthropodAnomal %>% 
  ggplot(aes(y = AnomalOccurence, x = AnomalTmax, )) +
  geom_point(alpha = 0.6,  aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  labs(
    title = "Anomality in Maximum Daily Temperature",
    x = "Temperature anomaly",
    y = "Occurence anomaly"
  ) +
  facet_wrap(~Group)+
  theme_minimal(base_size = 13)


TempArthropodAnomal %>% 
  ggplot(aes(y = AnomalJulWeek, x = AnomalTmax, )) +
  geom_point(alpha = 0.6,  aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  labs(
    title = "Anomality in Daily Maximum Temperature",
    x = "Temperature Anomaly",
    y = "Peak Julian Week anomaly"
  ) +
  facet_wrap(~Group)+
  theme_minimal(base_size = 13)

# Minimum Temperature Anomality

TempArthropodAnomal %>% 
  ggplot(aes(y = AnomalOccurence, x = AnomalTmin, )) +
  geom_point(alpha = 0.6, aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  labs(
    title = "Anomality in Daily Minimum Temperature",
    x = "Temperature Anomaly",
    y = "Occurence anomaly"
  ) +
  facet_wrap(~Group)+
  theme_minimal(base_size = 13)



TempArthropodAnomal %>% 
  ggplot(aes(y = AnomalJulWeek, x = AnomalTmin, )) +
  geom_point(alpha = 0.6,  aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  labs(
    title = "Anomality in Minimum Daily temperature",
    x = "Temperature Anomaly",
    y = "Peak Julian Week anomaly"
  ) +
  facet_wrap(~Group,  scales = "free_y")+
  theme_minimal(base_size = 13)




catJulTminLat= lm(AnomalOccurence ~ AnomalTmin * scale(Latitude), 
               data = TempArthropodAnomal %>% 
                 filter(Group == "Caterpillar"))
summary(catJulTminLat)

interact_plot(
  catJulTminLat,
  pred = AnomalTmin,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Peak  anomaly",
  x.lab =  "Temperature Anomaly", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 



catJulTmaxLat= lm(AnomalJulWeek~ AnomalTmax, data = TempArthropodAnomal %>% 
                    filter(Group == "Caterpillar"))
summary(catJulTmaxLat)  



# is Phenomenological anomality spatially auto-correlated?



pheno_CV = TempArthropodAnomal %>% select(site, ObservationMethod,
                               Latitude, Longitude, year, maxjulianweek, maxOcc, Group) %>% 
  group_by(site,  Latitude, Longitude, ObservationMethod, Group) %>% 
  summarise(cv_maxjulianweek  = sd(maxjulianweek)/mean(maxjulianweek),
          cv_maxOcc = sd(maxOcc)/mean(maxOcc)) %>% as.data.frame()


# Investigate if the CV for each site is auto-correlated spatially

pheno_LatLonDist = dist(pheno_CV %>% 
                          filter(Group == "Caterpillar") %>% 
                          select(Latitude, Longitude),
                    method = "euclidean")


phenoHel_dist <- dist(pheno_CV%>% 
                               filter(Group == "Caterpillar") %>% 
                               select (cv_maxjulianweek), method = "euclidean")
 
phenoLatLon <- pheno_CV[, c("Longitude", "Latitude")]

phenogeo_dist <- distm(phenoLatLon, fun = distHaversine)
phenogeo_dist <- as.dist(phenogeo_dist)


pheno_mantel <- mantel(phenoHel_dist, phenogeo_dist, 
                      method = "pearson", permutations = 999)






