 
  
  
  minSurveys = 50
  Pheno_julianWindow = 140:213
  minSurvYears = 3
  
fullDataset %>%
    filter(julianday %in% Pheno_julianWindow,
           Longitude > -100,
           WetLeaves == 0,
           !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK')) %>%
    group_by(Name, ObservationMethod) %>%
    summarize(nSurvs = n_distinct(ID),
              nYears = n_distinct(Year)) %>%
    filter(nSurvs >= 40,
           nYears >= minSurvYears) %>%
    arrange(desc(nYears)) %>% as.data.frame()


Julian_Arthropod = fullDataset %>%
  filter(julianday %in% Pheno_julianWindow,
         Longitude > -100,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         Group %in% c('caterpillar', 'spider', 'ant', 'leafhopper', 'beetle', 
                      'truebugs', 'fly', 'grasshopper', 
                      'daddylonglegs', 'aphid')) %>% 
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



goodsurv = Julian_Arthropod %>% 
  filter(surveyNum >=  10)  # each site/year/week should have at least 10 surveys.


Site.julian = Julian_Arthropod %>% 
  inner_join(goodsurv, by = c("Name", "ObservationMethod", "Year", "julianweek", "surveyNum")) %>% 
  group_by(Name, ObservationMethod, Year) %>% 
  summarise(nJulianWeek= n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= 5)


JuliandData = inner_join(Julian_Arthropod, Site.julian, by = c("Name", "ObservationMethod", "Year"))

siteyear =  JuliandData %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(nYear= n_distinct(Year)) %>% 
  filter(nYear >= 3)


JuliSiteData = inner_join(JuliandData, siteyear,
                          by = c("Name", "ObservationMethod"))




# Return only the maximum value of arthropod group for each year
MaxArthropod = JuliSiteData %>% 
  group_by(Name, ObservationMethod, Year) %>% 
  summarise(across(where(is.numeric)  & !c(julianweek, surveyNum, nJulianWeek, nYear),
                   ~ max(.x, na.rm = TRUE), .names = "{col}")) 


# Caterpillars data: 

# Return only the maximum value of weekly caterpillar prop. of occurrence for each year and it corresponding Julian week


JuliSiteData_uniqueCaterpillar <- JuliSiteData %>%
  group_by(Name, ObservationMethod, Year, Caterpillar) %>% # deals with potential duplicates
  summarise(julianweek = mean(julianweek), .groups = "drop")


MaxJulCaterpillar = MaxArthropod %>% 
  select(Name, ObservationMethod, Year, Caterpillar) %>% 
  inner_join(JuliSiteData_uniqueCaterpillar, 
             by = c("Name", "ObservationMethod", 
                    "Year",
                    "Caterpillar"
                    )) %>% 
  select(Name, ObservationMethod, Year, julianweek, Caterpillar) %>% 
  rename("MaxJulWeek" = "julianweek",
         "MaxCaterpillar" = "Caterpillar"
         )


AbnormCaterpillar =  MaxJulCaterpillar %>% # note that you are joining by the summary statistics.
  left_join(MaxJulCaterpillar %>%          # So, no confusion here!
              group_by(Name, ObservationMethod) %>% 
              summarise(meanJul = mean(MaxJulWeek),
                        meanOccurence = mean(MaxCaterpillar)), 
            by = c("Name", "ObservationMethod")) %>% 
  mutate(AnomalJulWeek = abs(MaxJulWeek - meanJul),
         AnomalOccurence = abs(MaxCaterpillar  - meanOccurence))




# Spider  data: 

# Return only the maximum value of weekly caterpillar prop. of occurrence for each year and it corresponding Julian week


JuliSiteData_uniqueSpider <- JuliSiteData %>%
  group_by(Name, ObservationMethod, Year, Spider) %>%  
  summarise(julianweek = mean(julianweek), .groups = "drop")


MaxJulSpider = MaxArthropod %>% 
  select(Name, ObservationMethod, Year, Spider) %>% 
  inner_join(JuliSiteData_uniqueSpider, 
             by = c("Name", "ObservationMethod", "Year", "Spider")) %>% 
  select(Name, ObservationMethod, Year, julianweek, Spider) %>% 
  rename("MaxJulWeek" = "julianweek",
         "MaxSpider" = "Spider")



AbnormSpider =  MaxJulSpider %>%  
  left_join(MaxJulSpider %>%          
              group_by(Name, ObservationMethod) %>% 
              summarise(meanJul = mean(MaxJulWeek),
                        meanOccurence = mean(MaxSpider)), 
            by = c("Name", "ObservationMethod")) %>% 
  mutate(AnomalJulWeek = abs(MaxJulWeek - meanJul),
         AnomalOccurence = abs(MaxSpider  - meanOccurence))




# Ant  data: 

# Return only the maximum value of weekly ant prop. of occurrence for each year and it corresponding Julian week

MaxJulAnt = MaxArthropod %>% 
  select(Name, ObservationMethod, Year, Ant) %>% 
  left_join(JuliSiteData, 
             by = c("Name", "ObservationMethod", "Year", "Ant")) %>% 
  select(Name, ObservationMethod, Year, julianweek, Ant) %>% 
  rename("MaxJulWeek" = "julianweek",
         "MaxAnt" = "Ant")



AbnormAnt =  MaxJulAnt %>%  
  left_join(MaxJulAnt %>%          
              group_by(Name, ObservationMethod) %>% 
              summarise(meanJul = mean(MaxJulWeek),
                        meanOccurence = mean(MaxAnt)), 
            by = c("Name", "ObservationMethod")) %>% 
  mutate(AnomalJulWeek = abs(MaxJulWeek - meanJul),
         AnomalOccurence = abs(MaxAnt  - meanOccurence))





 

