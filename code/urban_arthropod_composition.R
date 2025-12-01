install.packages("interactions", type = "source")

install.packages("future")
install.packages("interactions", dependencies = TRUE)



fullDataset = read.csv("https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/refs/heads/master/data/fullDataset_2025-07-28.csv")




richness.data = fullDataset %>%
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
            richness = sum(caterpillar, spider, beetle, truebug, hopper, ant))

dataset.rich = left_join(richness.data, sites, by = 'Name')


rich.Dev.Latitude = glm(richness ~ dev*Latitude + ObservationMethod, 
                       data = dataset.rich, family = "poisson")
summary(rich.Dev.Latitude)



# unable to run this successfully because interaction package fails to install.
richDevPlot = interact_plot(rich.Dev.Latitude, pred = dev, modx = Latitude,
                           y.label = "Taxa richness",
                           x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
                           colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2)  


fullDataset %>% 
  filter(ID == 73284) %>% 
  select(Group)


unique(fullDataset$Group)


expertID = read.csv("data/2025-09-13_ExpertIdentification.csv") %>%
  select(-1) %>% 
  rename( "ID" = ArthropodSightingFK)

fullDataset.ID = left_join(fullDataset, expertID, by  = "ID")

fullDataset.ID %>% 
  head(100) %>% 
  View()


fullDataset.ID %>%
  filter(Group == "caterpillar") %>%
  summarise(
    total_caterpillar = n(),
    identified_caterpillar = sum(!is.na(TaxonName)),
    proportion = identified_caterpillar/total_caterpillar)



#############################################################################################

# Summary of expert identifications

require(dplyr)
require(stringr)
require(taxize)

expert = read.csv(list.files()[grepl("ExpertIdentification.csv", list.files())], header = TRUE, stringsAsFactors = FALSE)
arths = read.csv(list.files()[grepl("ArthropodSighting.csv", list.files())], header = TRUE, stringsAsFactors = FALSE)
survey = read.csv(list.files()[grepl("Survey.csv", list.files())], header = TRUE, stringsAsFactors = FALSE)
plant = read.csv(list.files()[grepl("Plant.csv", list.files())], header = TRUE, stringsAsFactors = FALSE)
sites = read.csv(list.files()[grepl("Site.csv", list.files())], header = TRUE, stringsAsFactors = FALSE)

expert_higher_tax = read.csv("classified_expert_identifications.csv", header = T, quote = '\"')

ranks = read.csv("taxon_ranks.csv", header = T)

focalArthGroups = c("ant", "aphid", "bee", "beetle", "caterpillar", 
                    "daddylonglegs", "fly", "grasshopper", "leafhopper",
                    "moths", "spider", "truebugs")

exp = left_join(expert, expert_higher_tax, by = c('TaxonName', 'Rank')) %>%
  left_join(ranks, by = "Rank")


## Basic summaries

totalRecords = nrow(exp)
totalSpecies = length(unique(exp$TaxonName[exp$Rank %in% c("species", "subspecies") ]))

topSpecies = exp %>%
  filter(Rank == "species") %>%
  count(TaxonName, StandardGroup) %>%
  arrange(desc(n))

topSpeciesByNumSites = exp %>%
  left_join(arths[, c("ID", "SurveyFK")], by = c("ArthropodSightingFK" = "ID")) %>%
  left_join(survey[, c("ID", "PlantFK")], by = c("SurveyFK" = "ID")) %>%
  left_join(plant, c("ID", "SiteFK"), by = c("PlantFK" = "ID")) %>%
  filter(Rank == "species") %>%
  group_by(TaxonName, StandardGroup) %>%
  summarize(nSites = length(unique(SiteFK))) %>%
  arrange(desc(nSites))

mostWidespreadByGroup = topSpeciesByNumSites %>%
  filter(StandardGroup %in% focalArthGroups) %>%
  group_by(StandardGroup) %>%
  slice_head(n = 1)

topSitesByNumPhotos = exp %>%
  left_join(arths[, c("ID", "SurveyFK")], by = c("ArthropodSightingFK" = "ID")) %>%
  left_join(survey[, c("ID", "PlantFK")], by = c("SurveyFK" = "ID")) %>%
  left_join(plant[, c("ID", "SiteFK")], by = c("PlantFK" = "ID")) %>%
  left_join(sites[, c("ID", "Name")], by = c("SiteFK" = "ID")) %>%
  count(Name) %>%
  arrange(desc(n))

topSitesByNumSpecies = exp %>%
  left_join(arths[, c("ID", "SurveyFK")], by = c("ArthropodSightingFK" = "ID")) %>%
  left_join(survey[, c("ID", "PlantFK")], by = c("SurveyFK" = "ID")) %>%
  left_join(plant[, c("ID", "SiteFK")], by = c("PlantFK" = "ID")) %>%
  left_join(sites[, c("ID", "Name")], by = c("SiteFK" = "ID")) %>%
  filter(Rank == "species") %>%
  group_by(Name) %>%
  summarize(nSpecies = length(unique(TaxonName))) %>%
  arrange(desc(nSpecies)) %>%
  left_join(topSitesByNumPhotos)




## Breakdown by arthropod groups

groupCount = exp %>%
  mutate(NewGroup = ifelse(StandardGroup %in% focalArthGroups, StandardGroup, "other"),
         NewGroup = ifelse(StandardGroup == "bee", "bee/wasp", NewGroup)) %>%
  group_by(NewGroup) %>%
  summarize(numPhotos = n(),
            numSpecies = length(unique(TaxonName[Rank == "species"]))) %>%
  mutate(pctPhotos = round(100*numPhotos/sum(numPhotos), 2),
         pctSpecies = round(100*numSpecies/sum(numSpecies), 2)) %>%
  arrange(desc(numPhotos))


par(mfrow = c(1,1), mar = c(2, 5, 2, 5))
pie(groupCount$numPhotos, labels = groupCount$NewGroup, col = rainbow(nrow(groupCount)))

#pie(groupCount$numSpecies, labels = groupCount$NewGroup, col = rainbow(nrow(groupCount)))


## Summarizing by taxonomic level of identification

id_rank = exp %>%
  filter(StandardGroup %in% focalArthGroups) %>%
  count(StandardGroup, Rank_group) %>%
  group_by(StandardGroup) %>%
  mutate(tot = sum(n),
         pct = round(100*n/tot,2)) %>%
  arrange(desc(tot), Rank_group)


# Barplots

id_rank_plot = function(dataframe, arthGroup, ...) {
  
  tmp = filter(dataframe, StandardGroup == arthGroup)
  
  allRanks = data.frame(Rank_group = 0:4) %>% left_join(tmp, by = "Rank_group")
  
  allRanks$pct[is.na(allRanks$pct)] = 0
  
  foo = barplot(allRanks$pct ~ allRanks$Rank_group, xaxt = "n", ylab = "Percent", xlab = "", las = 1, ...)
  
  #mtext(c("Unconfirmed", "Order", "Family", "Genus", "Species"), 1, at = foo, cex = 1.2)
  
  return(foo)
}


# All arth groups
par(mfrow = c(3, 4), mar = c(3, 5, 3, 1), cex.lab =1.4, cex.axis = 1.2)
for (a in unique(id_rank$StandardGroup)) {
  
  foo = id_rank_plot(id_rank, a, col = rainbow(5), main = paste(a, ", n =", id_rank$tot[id_rank$StandardGroup == a][1]))
  
  mtext(c("U", "O", "F", "G", "S"), 1, at = foo, cex = 1, line = 0.5)
}


# Arth group subset

par(mfrow = c(2, 3), mar = c(3, 5, 3, 1), cex.lab =1.8, cex.axis = 1.2, cex.main = 1.4)
for (a in c('fly', 'aphid', 'caterpillar', 'beetle', 'grasshopper', 'leafhopper')) {
  
  foo = id_rank_plot(id_rank, a, col = rainbow(5), main = paste(a, ", n =", id_rank$tot[id_rank$StandardGroup == a][1]))
  
  mtext(c("U", "O", "F", "G", "S"), 1, at = foo, cex = 1, line = 0.5)
}
