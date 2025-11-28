# Load library ----
library(tidyverse)
library(daymetr)
library(geosphere)
require(vegan)
library(ggrepel)
library(ggpubr)
library(interactions)
library(broom)

# -- Define data inclusion criteria-----
  minSurveys = 50
  Pheno_julianWindow = 140:213
  minSurvYears = 3 # min number of years of survey at a site.
  min.nSurvWeekYearSite = 10 # minimum number of survey in a week, at a year, at a site.
  yearException = 2025 # because we do not yet have 2025 data for temperature
  TempDayWindow = 90:180   # Adjust
  min.nJulianWeekYearSite = 5
  
  
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
  filter(surveyNum >=  min.nSurvWeekYearSite)  # each site/year/week should have at least 10 surveys.


Site.julian = goodsurv %>% 
  group_by(Name, ObservationMethod, Year) %>% 
  summarise(nJulianWeek= n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= min.nJulianWeekYearSite) # must 5 or more weeks of survey


JuliandData = inner_join(goodsurv, Site.julian, 
                         by = c("Name", "ObservationMethod", "Year"))

siteyear =  JuliandData %>% # Good sites with 3 or more years of survey
  group_by(Name, ObservationMethod) %>% 
  summarise(nYear= n_distinct(Year)) %>% 
  filter(nYear >= minSurvYears)


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

 




## Caterpillar Data----

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

 




 



# Spider  data: ----

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

# Ant  data: ----

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






#  Temperature data retrieval----


 


 

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
 





# Join temperature anomaly to phonology anomaly-----

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
  left_join(GoodSiteYearLatLon, by = c("site" = "Name")) %>% 
  left_join(sites %>% select(Name, dev, forest), by = c( "site" = "Name")) %>% 
  mutate(Group = factor(Group, levels =c("Caterpillar",
                                         "Ant",
                                         "Spider")))

 

# Maximum Temperature Anomality


TempArthropodAnomal %>% 
  ggplot(aes(y = AnomalOccurence, x = AnomalTmax)) +
  geom_point(alpha = 0.6, aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +   
  stat_regline_equation(
    aes(label = paste(..eq.label..,  sep = "~~~")),
    label.x = -1.5,
    label.y = -0.3,
    color = "red",
    size = 3
  ) +
  labs(
    title = "Anomality in Maximum Daily Temperature",
    x = "Temperature anomaly",
    y = "Occurrence anomaly"
  ) +
  facet_wrap(~Group) +
  theme_minimal(base_size = 13)


TempArthropodAnomal %>% 
  ggplot(aes(y = AnomalJulWeek, x = AnomalTmax, )) +
  geom_point(alpha = 0.6,  aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +   
  stat_regline_equation(
    aes(label = paste(..eq.label..,  sep = "~~~")),
    label.x = -1.5,
    label.y = -35,
    color = "red",
    size = 3
  ) +
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
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  stat_regline_equation(
    aes(label = paste(..eq.label..,  sep = "~~~")),
    label.x = -1.5,
    label.y = -0.3,
    color = "red",
    size = 3
  ) +
  labs(
    title = "Anomality in Daily Minimum Temperature",
    x = "Temperature Anomaly",
    y = "Occurence anomaly"
  ) +
  facet_wrap(~Group)+
  theme_minimal(base_size = 13)



TempArthropodAnomal %>% 
  ggplot(aes(y = AnomalJulWeek, x = AnomalTmin)) +
  geom_point(alpha = 0.6, aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  stat_regline_equation(
    aes(label = paste(..eq.label..,  sep = "~~~")),
    label.x = -1.5,
    label.y = -35,
    color = "red",
    size = 3
  ) +
  labs(
    title = "Anomality in Minimum Daily Temperature",
    x = "Temperature Anomaly",
    y = "Peak Julian Week Anomaly"
  ) +
  facet_wrap(~Group) +
  theme_minimal(base_size = 13)





catJulTminLat= lm(AnomalJulWeek ~ AnomalTmin , 
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



### -- testing interaction effects -----

catJulTminDev= lm(AnomalJulWeek ~ AnomalTmin * dev , # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Caterpillar"))

catJulTmaxDev= lm(AnomalJulWeek ~ AnomalTmax * dev , # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Caterpillar"))

catOccTminDev= lm(AnomalOccurence ~ AnomalTmin * dev , # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Caterpillar"))

catOccTmaxDev= lm(AnomalOccurence ~ AnomalTmax * dev , # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Caterpillar"))

antJulTminDev= lm(AnomalJulWeek ~ AnomalTmin * dev , # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Ant"))

antJulTmaxDev= lm(AnomalJulWeek ~ AnomalTmax * dev , # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Ant"))

antOccTminDev= lm(AnomalOccurence ~ AnomalTmin * dev , # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Ant"))

antOccTmaxDev= lm(AnomalOccurence ~ AnomalTmax * dev , # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Ant"))



spiJulTminDev= lm(AnomalJulWeek ~ AnomalTmin * dev , 
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Spider"))

spiJulTmaxDev= lm(AnomalJulWeek ~ AnomalTmax * dev ,  
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Spider"))

spiOccTminDev= lm(AnomalOccurence ~ AnomalTmin * dev ,  # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Spider"))

spiOccTmaxDev= lm(AnomalOccurence ~ AnomalTmax * dev , # no interaction effect
                  data = TempArthropodAnomal %>% 
                    filter(Group == "Spider"))



phenoTempDevOutput = data.frame(rbind(summary(catJulTminDev)$coefficients, 
                                      summary(catJulTmaxDev)$coefficients, 
                                      summary(catOccTminDev)$coefficients,
                                      summary(catOccTmaxDev)$coefficients, 
                                      summary(antJulTminDev)$coefficients, 
                                      summary(antJulTmaxDev)$coefficients, 
                                      summary(antOccTminDev)$coefficients,
                                      summary(antOccTmaxDev)$coefficients, 
                                      summary(spiJulTminDev)$coefficients, 
                                      summary(spiJulTmaxDev)$coefficients, 
                                      summary(spiOccTminDev)$coefficients,
                                      summary(spiOccTmaxDev)$coefficients))
phenoTempDevOutput$term = rep(c('Intercept', 'AnomalTmin', 'dev', 'AnomalTmax:dev',
                            'Intercept', 'AnomalTmax', 'dev', 'AnomalTmax:dev'), times = 6)





spiJulTminDevPlot= interact_plot(
  spiJulTminDev,
  pred = AnomalTmin,
  modx = dev,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Julian week",
  x.lab =  "Min. Temperature", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2)+
  annotation_raster(spiderImage, ymin = 5, ymax = 10, xmin = .2, xmax = .8)



spiJulTmaxDevPlot= interact_plot(
  spiJulTmaxDev,
  pred = AnomalTmax,
  modx = dev,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Julian week",
  x.lab =  "Max. Temperature", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2)+
  annotation_raster(spiderImage, ymin = 5, ymax = 10, xmin = .2, xmax = .8)


antJulTminDevPlot= interact_plot(
  antJulTminDev,
  pred = AnomalTmin,
  modx = dev,
  plot.points = FALSE,
  interval = FALSE,
  y.label =  "Julian week",
  x.lab =  "Min. Temperature", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2)+
  annotation_raster(antImage, ymin = 5, ymax = 12, xmin = .2, xmax = .9)


antOccTminDevPlot= interact_plot(
  antOccTminDev,
  pred = AnomalTmin,
  modx = dev,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "% Occurence",
  x.lab =  "Min. Temperature", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2)+
  annotation_raster(antImage, ymin = .05, ymax = 0.01, xmin = -.9, xmax = -.2)








annotate_figure(
  ggarrange(
    antJulTminDevPlot, antOccTminDevPlot,
    spiJulTminDevPlot, spiJulTmaxDevPlot,
    ncol = 2, nrow = 2,
    common.legend = TRUE, legend = "bottom"
  ),
  top = text_grob(
    "Significant (P < 0.1) effect of interaction with % urban development",
    face = "bold", size = 12
  )
)
# is Phenomenological anomality spatially auto-correlated?



# Investigate if the CV for each site is auto-correlated spatially


ggplot(TempArthropodAnomal, aes(x = "", y = dev)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.1, alpha = 0.4, color = "blue") +  
  labs(
    title = "Distribution of dev",
    y = "dev",
    x = ""
  ) +
  theme_minimal(base_size = 13)





# CENTROIDS  (instead of peak occurrence)----


#
CaterpillarCentroid= JuliSiteData %>%
  group_by(Name, ObservationMethod, Year, nJulianWeek, nYear) %>%
  mutate(across(
    .cols = where(is.numeric) & !matches("^(julianweek|surveyNum)$"),
    .fns = ~ cumsum(replace_na(.x, 0)),
    .names = "{.col}"
  )) %>% 
  select(c(1:4), "Caterpillar") %>% 
  group_by(Name, ObservationMethod, Year, nJulianWeek, nYear) %>%
  mutate(Sum = max(Caterpillar),
         CentroidOcc = 0.5 * Sum) %>%
  filter(Caterpillar >= CentroidOcc) %>%
  slice_min(order_by = julianweek) %>%  #tells R which variable to use to determine the “minimum”
  rename("CentroidJulianWeek" = "julianweek") %>% 
  select(Name, ObservationMethod, nYear, nJulianWeek, Year, 
         CentroidJulianWeek, CentroidOcc, Sum)

AbnormCaterpillarCentroid = CaterpillarCentroid%>% 
  left_join(CaterpillarCentroid %>% 
              group_by(Name, ObservationMethod) %>% 
              summarise(meanCentroidJulWeek = mean(CentroidJulianWeek))) %>% 
  mutate(abnormalCentroid = CentroidJulianWeek  - meanCentroidJulWeek)



#
AntCentroid= JuliSiteData %>%
  group_by(Name, ObservationMethod, Year, nJulianWeek, nYear) %>%
  mutate(across(
    .cols = where(is.numeric) & !matches("^(julianweek|surveyNum)$"),
    .fns = ~ cumsum(replace_na(.x, 0)),
    .names = "{.col}"
  )) %>% 
  select(c(1:4), "Ant") %>% 
  group_by(Name, ObservationMethod, Year, nJulianWeek, nYear) %>%
  mutate(Sum = max(Ant),
         CentroidOcc = 0.5 * Sum) %>%
  filter(Ant >= CentroidOcc) %>%
  slice_min(order_by = julianweek) %>%   
  rename("CentroidJulianWeek" = "julianweek") %>% 
  select(Name, ObservationMethod, nYear, nJulianWeek, Year, 
         CentroidJulianWeek, CentroidOcc, Sum) 

AbnormAntCentroid = AntCentroid%>% 
  left_join(AntCentroid %>% 
              group_by(Name, ObservationMethod) %>% 
              summarise(meanCentroidJulWeek = mean(CentroidJulianWeek))) %>% 
  mutate(abnormalCentroid = CentroidJulianWeek  - meanCentroidJulWeek)


SpiderCentroid= JuliSiteData %>%
  group_by(Name, ObservationMethod, Year, nJulianWeek, nYear) %>%
  mutate(across(
    .cols = where(is.numeric) & !matches("^(julianweek|surveyNum)$"),
    .fns = ~ cumsum(replace_na(.x, 0)),
    .names = "{.col}"
  )) %>% 
  select(c(1:4), "Spider") %>% 
  group_by(Name, ObservationMethod, Year, nJulianWeek, nYear) %>%
  mutate(Sum = max(Spider),
         CentroidOcc = 0.5 * Sum) %>%
  filter(Spider >= CentroidOcc) %>%
  slice_min(order_by = julianweek) %>%   
  rename("CentroidJulianWeek" = "julianweek") %>% 
  select(Name, ObservationMethod, nYear, nJulianWeek, Year, 
         CentroidJulianWeek, CentroidOcc, Sum)

AbnormSpiderCentroid = SpiderCentroid%>% 
  left_join(SpiderCentroid %>% 
              group_by(Name, ObservationMethod) %>% 
              summarise(meanCentroidJulWeek = mean(CentroidJulianWeek))) %>% 
  mutate(abnormalCentroid = CentroidJulianWeek  - meanCentroidJulWeek)



anomalPhenoCentroid = rbind(AbnormAntCentroid %>% mutate(Group =  "Ant"),
                            AbnormSpiderCentroid %>% mutate(Group = "Spider"),
                            AbnormCaterpillarCentroid %>% mutate(Group = "Caterpillar")) %>% 
  as.data.frame()

## -- Combine centroids with Temperature----
TempAnomalPhenoCentroid = AllTempData %>% 
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
  right_join(anomalPhenoCentroid %>%  filter(Year != "2025"), # because daymetr has no 2025 yet.
             by = c("site" = "Name", "year" = "Year")) %>% 
  left_join(GoodSiteYearLatLon, by = c("site" = "Name")) %>% 
  left_join(sites %>% select(Name, dev, forest), by = c( "site" = "Name")) %>% 
  mutate(Group = factor(Group, levels =c("Caterpillar",
                                         "Ant",
                                         "Spider"))) %>% 
  mutate(SiteObserv = paste(site, ObservationMethod, sep = "_"))

# Arthropod images
catImage = readPNG('images/caterpillar.png')
antImage = readPNG('images/ant.png')
beetleImage = readPNG('images/beetle.png')
spiderImage = readPNG('images/spider.png')
hopperImage = readPNG('images/leafhopper.png')
truebugImage = readPNG('images/truebugs.png')

TempAnomalPhenoCentroid %>% 
  ggplot(aes(y = abnormalCentroid, x = AnomalTmax)) +
  geom_point(alpha = 0.6, aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +   
  stat_regline_equation(
    aes(label = paste(..eq.label..,  sep = "~~~")),
    label.x = -1.5,
    label.y = -10.3,
    color = "red",
    size = 3
  ) +
  labs(
    x = "Maximum Daily Temperature anomaly (Degree Celcius)",
    y = "Timing anomaly (Days)"
  ) +
  facet_wrap(~Group) +
  theme_minimal(base_size = 13) 


TempAnomalPhenoCentroid %>% 
  ggplot(aes(y = abnormalCentroid, x = AnomalTmax)) +
  geom_point(alpha = 0.6, aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  labs(
    x = "Maximum Daily Temperature anomaly (Degree Celcius)",
    y = "Timing anomaly (Days)") +
  facet_wrap(~Group) +
  theme_bw(base_size = 13)


TempAnomalPhenoCentroid %>% 
  ggplot(aes(y = abnormalCentroid, x = AnomalTmin)) +
  geom_point(alpha = 0.6, aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  stat_regline_equation(
    aes(label = paste(..eq.label..,  sep = "~~~")),
    label.x = -1.5,
    label.y = -35,
    color = "red",
    size = 3
  ) +
  labs(
    title = "Anomality in Minimum Daily Temperature",
    x = "Temperature Anomaly",
    y = "Centroid Julian week anomaly"
  ) +
  facet_wrap(~Group) +
  theme_minimal(base_size = 13)



TempAnomalPhenoCentroid %>% 
  ggplot(aes(y = abnormalCentroid, x = AnomalTmin)) +
  geom_point(alpha = 0.6, aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  labs(
    x = "Minimum Daily Temperature anomaly (Degree Celcius)",
    y = "Timing anomaly (Days)") +
  facet_wrap(~Group) +
  theme_bw(base_size = 13)




TempAnomalPhenoCentroid %>% 
  ggplot(aes(y = abnormalCentroid, x = AnomalPreci)) +
  geom_point(alpha = 0.6, aes(colour = Latitude)) +  
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              fill = "pink", linewidth = 1) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  stat_regline_equation(
    aes(label = paste(..eq.label..,  sep = "~~~")),
    label.x = -1.5,
    label.y = -35,
    color = "red",
    size = 3
  ) +
  labs(
    title = "Anomality in Precipation",
    x = "Precipiation Anomaly",
    y = "Centroid Julian week anomaly"
  ) +
  facet_wrap(~Group) +
  theme_minimal(base_size = 13)
### 



summary(lm(abnormalCentroid~AnomalTmax, 
           data = TempAnomalPhenoCentroid %>% filter(Group == "Caterpillar")))

summary(lm(abnormalCentroid~AnomalTmax, 
           data = TempAnomalPhenoCentroid %>% filter(Group == "Ant")))

summary(lm(abnormalCentroid~AnomalTmax, 
           data = TempAnomalPhenoCentroid %>% filter(Group == "Spider")))



summary(lm(abnormalCentroid~AnomalTmin, 
           data = TempAnomalPhenoCentroid %>% filter(Group == "Caterpillar")))

summary(lm(abnormalCentroid~AnomalTmin, 
           data = TempAnomalPhenoCentroid %>% filter(Group == "Ant")))

summary(lm(abnormalCentroid~AnomalTmin, 
           data = TempAnomalPhenoCentroid %>% filter(Group == "Spider")))



# What if we fit random slopes-intercept to each site -----




TempAnomalPhenoCentroid %>% 
  filter(nYear >= 4) %>% 
  ggplot(aes(y = abnormalCentroid, 
             x = AnomalTmax, 
             colour = Latitude)) +
  geom_point(alpha = 0.4) +  
  geom_smooth(method = "lm", se = FALSE, aes(group = SiteObserv), alpha = 0.4) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  labs(
    x = "Maximum Daily Temperature anomaly (Degree Celsius)",
    y = "Centroid timing anomaly (Days)",
    title = "Site >= 4 nYears"
  ) +
  facet_wrap(~Group)+
  theme_bw(base_size = 13)


TempAnomalPhenoCentroid %>% 
  filter(nYear >= 4) %>% 
  ggplot(aes(y = abnormalCentroid, 
             x = AnomalTmin, 
             colour = Latitude)) +
  geom_point(alpha = 0.4) +  
  geom_smooth(method = "lm", se = FALSE, aes(group = SiteObserv), alpha = 0.4) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  # geom_text(
  #   data = TempAnomalPhenoCentroid %>%
  #     filter(nYear >= 4) %>%
  #     group_by(SiteObserv, Group) %>%
  #     summarise(
  #       cx = mean(AnomalTmin, na.rm = TRUE),
  #       cy = mean(abnormalCentroid, na.rm = TRUE),
  #       lat = round(mean(Latitude, na.rm = TRUE)) 
  #     ) ,
  #   aes(x = cx, y = cy+1, label = lat), # just to lift the text up a little
  #   size = 5,
  #   inherit.aes = FALSE
  # ) +
  labs(
    x = "Minimum Daily Temperature anomaly (Degree Celsius)",
    y = "Centroid timing anomaly (Days)",
    title = "Site >= 4 nYears"
  ) +
  facet_wrap(~Group)+
  theme_bw(base_size = 13)

## Fitting random effect  for every siteObserv----

library(nlme)
library(lme4)

summary(lme(
  abnormalCentroid ~ AnomalTmin,
  random = ~ AnomalTmin | SiteObserv,
  data = TempAnomalPhenoCentroid %>% 
    filter(Group == "Caterpillar")
))



summary(lm(abnormalCentroid ~ AnomalTmin , 
           data = TempAnomalPhenoCentroid %>% filter(Group == "Caterpillar") ))



summary(lme( # convergence problem
  abnormalCentroid ~ AnomalTmax,
  random = ~ AnomalTmax | SiteObserv,
  data = TempAnomalPhenoCentroid %>% 
    filter(Group == "Caterpillar")
))

summary(lme( # summary estimates change every time it is run due to the autocorrelation accounted for
  abnormalCentroid ~ AnomalTmin,
  random = ~ AnomalTmin | SiteObserv,
  correlation = corExp(form = ~ Latitude_j | SiteObserv),
  data = TempAnomalPhenoCentroid %>% 
    mutate( # Jitter so nlme dosen't have zero distance problem for siteObserv with same coordinate
      Latitude_j = Latitude + runif(n(), -1e-4, 1e-4),
      Longitude_j = Longitude + runif(n(), -1e-4, 1e-4)
    ) %>% 
    filter(Group == "Caterpillar")))


TempAnomalPhenoCentroid %>% 
  mutate(abnormalCentroid.Latitude = abnormalCentroid + meanTmin) %>% 
  filter(nYear >= 4) %>% 
  ggplot(aes(y = abnormalCentroid.Latitude, 
             x = AnomalTmin, 
             colour = Latitude)) +
  #geom_point(alpha = 0.4) +  
  geom_smooth(method = "lm", se = FALSE, aes(group = SiteObserv), alpha = 0.4) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  # geom_text(
  #   data = TempAnomalPhenoCentroid %>%
  #     filter(nYear >= 4) %>%
  #     group_by(SiteObserv, Group) %>%
  #     summarise(
  #       cx = mean(AnomalTmin, na.rm = TRUE),
  #       cy = mean(abnormalCentroid, na.rm = TRUE),
  #       lat = round(mean(Latitude, na.rm = TRUE)) 
  #     ) ,
  #   aes(x = cx, y = cy+1, label = lat), # just to lift the text up a little
  #   size = 5,
#   inherit.aes = FALSE
# ) +
labs(
  x = "Minimum Daily Temperature anomaly (Degree Celsius)",
  y = "Centroid timing anomaly (Days)",
  title = "Site >= 4 nYears"
) +
  facet_wrap(~Group)+
  theme_bw(base_size = 13)


TempAnomalPhenoCentroid %>% 
  filter(Group == "Caterpillar") %>% 
    select(abnormalCentroid, AnomalTmax, site) %>% 
  data.frame()



PhenoTmaxLm = TempAnomalPhenoCentroid %>% 
  group_by(SiteObserv, Group) %>% 
  do({
    model <- lm(abnormalCentroid ~ AnomalTmax, data = .)
    tidy_mod  <- broom::tidy(model)
    glance_mod <- broom::glance(model)
    
    tibble(
      Tmax.intercept = tidy_mod$estimate[ tidy_mod$term == "(Intercept)" ],
      Tmax.slope     = tidy_mod$estimate[ tidy_mod$term == "AnomalTmax" ],
      Tmax.p_value   = tidy_mod$p.value[ tidy_mod$term == "AnomalTmax" ],
      Tmax.r_squared = glance_mod$r.squared
    )
  }) %>% 
  ungroup() 

PhenoTminLm = TempAnomalPhenoCentroid %>% 
  group_by(SiteObserv, Group) %>% 
  do({
    model <- lm(abnormalCentroid ~ AnomalTmin, data = .)
    tidy_mod  <- broom::tidy(model)
    glance_mod <- broom::glance(model)
    
    tibble(
      Tmin.intercept = tidy_mod$estimate[ tidy_mod$term == "(Intercept)" ],
      Tmin.slope     = tidy_mod$estimate[ tidy_mod$term == "AnomalTmin" ],
      Tmin.p_value   = tidy_mod$p.value[ tidy_mod$term == "AnomalTmin" ],
      Tmin.r_squared = glance_mod$r.squared
    )
  }) %>% 
  ungroup() 


 TempPhenoLm = left_join(TempAnomalPhenoCentroid,
      left_join(PhenoTmaxLm, PhenoTminLm, by = c("SiteObserv", "Group")), 
  by = c("SiteObserv", "Group"))
  

 TempPhenoLmSummary =  TempPhenoLm %>% 
   group_by(SiteObserv, Latitude, Group, Tmin.slope, Tmax.slope) %>%
   summarise(n = n()) %>% data.frame()
   

   
   
 TempPhenoLmSummary %>% 
   ggplot(aes(y = abs(Tmin.slope), x = Latitude)) +
   geom_smooth(method = "lm", se = TRUE, color = "red", 
               fill = "pink", linewidth = 1) +   
   stat_regline_equation(
     aes(label = paste(..eq.label..)),
     label.x = 35.5,
     label.y = -0,
     color = "red",
     size = 3
   ) +
   facet_wrap(~Group) +
   theme_minimal(base_size = 13)
   
   TempPhenoLmSummary %>% 
     ggplot(aes(y = abs(Tmax.slope), x = Latitude)) +
     geom_smooth(method = "lm", se = TRUE, color = "red", 
                 fill = "pink", linewidth = 1) +   
     stat_regline_equation(
       aes(label = paste(..eq.label..)),
       label.x = 35.5,
       label.y = -20,
       color = "red",
       size = 3
     ) +
     facet_wrap(~Group) +
     theme_minimal(base_size = 13)
 