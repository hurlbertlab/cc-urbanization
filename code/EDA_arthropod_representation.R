# Summary of expert identifications
library(readxl)
library(patchwork)
library(jsonlite)
library(tidyverse)
# (1) Read in latest Caterpillars Count! raw dataset from the caterpillars-analysis-public repo----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file))


# Find a way to get the expert id directly form the github page.
exp.id = read.csv("data/exp.csv") # what generates this data?
# View(exp.id)


fullDataset.ID = left_join(fullDataset %>% 
                             rename("key" = "arthID"), 
                           exp.id %>% rename("key" = "ArthropodSightingFK"), 
                           by = "key")



idClass <- fullDataset.ID %>% 
  filter(
    Name %in% goodSites$Name,
    julianday %in% julianWindow,
    WetLeaves == 0,
    !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
    Group %in% c('caterpillar', 'spider', 'ant', 'leafhopper', 'beetle', 
                 'truebugs', 'fly', 'grasshopper', 'daddylonglegs')
  ) %>%
  filter(Rank %in% c("subfamily", "tribe", "subtribe", "genus", "subgenus", 
                     "species", "subspecies", "complex", "form", "section")) %>%
  select(Group, TaxonName, Latitude, Longitude, Rank) %>%
  mutate(Taxon = word(TaxonName, 1)) %>%
  mutate(
    Rank = case_when(
      Rank %in% c("subgenus", "species", "subspecies", "complex", "form", "section") ~ "genus",
      TRUE ~ Rank
    )
  ) %>%
  left_join(
    sites %>% distinct(Latitude, .keep_all = TRUE),
    by = "Latitude"
  ) %>%
  group_by(Group, Name, Rank, Taxon) %>%
  summarise(
    nOfTaxon = n_distinct(TaxonName, na.rm = TRUE),
    nofCount = n(),
    .groups = "drop"
  )


idClass%>% write.csv(file = "data\\expID1.csv",
                  row.names = FALSE)

total_sites= length(unique(idClass$Name))
spider_sites = length(unique(idClass[idClass$Group=="spider",]$Name))
prop_spider_site = spider_sites/total_sites



# Spider
  spiderID <- idClass %>% 
    filter(Group == "spider") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "spider", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site_p = nSite / n_distinct(idClass[idClass$Group == "spider", ]$Name), # site present with taxon
      prop_site_pa= nSite/n_distinct(idClass$Name),  # site present or absent with taxon
      Taxa_score = Prop_count * prop_site_pa,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
           cum_Taxa_score = cumsum(Taxa_score),
           cum_Taxa_score_rel = cum_Taxa_score / sum(Taxa_score)
    ) %>% 
    as.data.frame()
  
  # Scree plot
    ggplot(spiderID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Spider Taxa",
      x = "",
      y = "Taxa Score"
    )

 top10_spiderID <-  spiderID %>% 
   slice_max(order_by = Taxa_score, n = 20)
 
 spiderID.plot20 <- ggplot(top10_spiderID, aes(x = reorder(Taxon, -Taxa_score), y = Taxa_score)) +
   geom_bar(stat = "identity", alpha = 0.6, aes(colour = Rank, fill = Rank)) +
   theme_minimal(base_size = 12) +
   theme(
     axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
   ) +
   labs(
     title = "Top 20 Spider Taxa by Taxa Score",
     x = "",
     y = "Taxa Score"
   )
 
 spiderID.plot20
 
 
 
  # Ant
  antID = idClass %>% 
    filter(Group == "ant") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "ant", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site_p = nSite / n_distinct(idClass[idClass$Group == "ant", ]$Name), # site present with taxon
      prop_site_pa= nSite/n_distinct(idClass$Name),  # site present or absent with taxon
      Taxa_score = Prop_count * prop_site_pa,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
           cum_Taxa_score = cumsum(Taxa_score),
           cum_Taxa_score_rel = cum_Taxa_score / sum(Taxa_score) 
    ) %>% 
    as.data.frame()
  
  
  top20_antID = antID %>% 
    slice_max(order_by = Taxa_score, n = 20)
  
  # Scree plot
 antID.plot20 =  ggplot(antID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Ants (top 20)",
      x = "",
      y = "Taxa Score"
    )

 antID.plot20
  
 
  
  
  # beetle
  beetleID = idClass %>% 
    filter(Group == "beetle") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "beetle", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site_p = nSite / n_distinct(idClass[idClass$Group == "beetle", ]$Name), # site present with taxon
      prop_site_pa= nSite/n_distinct(idClass$Name),  # site present or absent with taxon
      Taxa_score = Prop_count * prop_site_pa,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
           cum_Taxa_score = cumsum(Taxa_score),
           cum_Taxa_score_rel = cum_Taxa_score / sum(Taxa_score)
    ) %>% 
    as.data.frame()
  
  top20_beetleID = beetleID %>% 
    slice_max(order_by = Taxa_score, n = 20)
  
  
 beetle.plot20= ggplot(top20_beetleID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Beetle (top 20)",
      x = "",
      y = "Taxa Score"
    )
  
 beetle.plot20
  
  # caterpillar
  caterpillarID <- idClass %>% 
    filter(Group == "caterpillar") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "caterpillar", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site_p = nSite / n_distinct(idClass[idClass$Group == "caterpillar", ]$Name), # site present with taxon
      prop_site_pa= nSite/n_distinct(idClass$Name),  # site present or absent with taxon
      Taxa_score = Prop_count * prop_site_pa,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
           cum_Taxa_score = cumsum(Taxa_score),
           cum_Taxa_score_rel = cum_Taxa_score / sum(Taxa_score)
    ) %>% 
    as.data.frame()
  
  top20_caterpillarID = caterpillarID %>% 
    slice_max(order_by = Taxa_score, n = 20)
  
  caterpillar.plot20= ggplot(top20_caterpillarID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Caterpillar (top 20)",
      x = "",
      y = "Taxa Score"
    )
  
  caterpillar.plot20
  
  
  # daddylonglegs
  daddylonglegsID <- idClass %>% 
    filter(Group == "daddylonglegs") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "daddylonglegs", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site_p = nSite / n_distinct(idClass[idClass$Group == "daddylonglegs", ]$Name), # site present with taxon
      prop_site_pa= nSite/n_distinct(idClass$Name),  # site present or absent with taxon
      Taxa_score = Prop_count * prop_site_pa,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
           cum_Taxa_score = cumsum(Taxa_score),
           cum_Taxa_score_rel = cum_Taxa_score / sum(Taxa_score)
    ) %>% 
    as.data.frame()
  
  top20_daddylonglegsID = daddylonglegsID %>% 
    slice_max(order_by = Taxa_score, n = 20)
  
  daddylonglegs.plot20 = ggplot(top20_daddylonglegsID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Daddylonglegs (top 20)",
      x = "",
      y = "Taxa Score"
    )
  
  daddylonglegs.plot20
  
  # truebugs
  truebugsID <- idClass %>% 
    filter(Group == "truebugs") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "truebugs", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site_p = nSite / n_distinct(idClass[idClass$Group == "truebugs", ]$Name), # site present with taxon
      prop_site_pa= nSite/n_distinct(idClass$Name),  # site present or absent with taxon
      Taxa_score = Prop_count * prop_site_pa,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
           cum_Taxa_score = cumsum(Taxa_score),
           cum_Taxa_score_rel = cum_Taxa_score / sum(Taxa_score)
    ) %>% 
    as.data.frame()
  
  top20_truebugsID = truebugsID %>% 
    slice_max(order_by = Taxa_score, n = 20)
  
  truebugs.plot20 = ggplot(top20_truebugsID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Truebugs (top 20)",
      x = "",
      y = "Taxa Score"
    )
  
  truebugs.plot20
  
  
  
  # grasshopper
   grasshopperID = idClass %>% 
    filter(Group == "grasshopper") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "grasshopper", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site_p = nSite / n_distinct(idClass[idClass$Group == "grasshopper", ]$Name), 
      # site present with taxon
      prop_site_pa= nSite/n_distinct(idClass$Name),  # site present or absent with taxon
      Taxa_score = Prop_count * prop_site_pa,
      .groups = "drop"
    ) %>% 
     arrange(desc(Taxa_score)) %>% 
     mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
            cum_Taxa_score = cumsum(Taxa_score),
            cum_Taxa_score_rel = cum_Taxa_score / sum(Taxa_score)
     ) %>% 
     as.data.frame()
  
   top20_grasshopperID = grasshopperID %>% 
     slice_max(order_by = Taxa_score, n = 20)
   
   
   
  grasshopper.plot20 = ggplot(grasshopperID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Grasshopper (top 20)",
      x = "",
      y = "Taxa Score"
    )
  
  grasshopper.plot20  
  
  
  # fly
  flyID = idClass %>% 
    filter(Group == "fly") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "fly", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site_p = nSite / n_distinct(idClass[idClass$Group == "fly", ]$Name), 
      # site present with taxon
      prop_site_pa= nSite/n_distinct(idClass$Name),  # site present or absent with taxon
      Taxa_score = Prop_count * prop_site_pa,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
           cum_Taxa_score = cumsum(Taxa_score),
           cum_Taxa_score_rel = cum_Taxa_score / sum(Taxa_score)
    ) %>% 
    as.data.frame()
  
  top20_flyID = flyID %>% 
    slice_max(order_by = Taxa_score, n = 20)
  
  
fly.plot20 = ggplot(top20_flyID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Fly (top 20)",
      x = "",
      y = "Taxa Score"
    )
fly.plot20
  
  
  # leafhopper
  leafhopperID = idClass %>% 
    filter(Group == "leafhopper") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "leafhopper", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site_p = nSite / n_distinct(idClass[idClass$Group == "leafhopper", ]$Name), 
      # site present with taxon
      prop_site_pa= nSite/n_distinct(idClass$Name),  # site present or absent with taxon
      Taxa_score = Prop_count * prop_site_pa,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)),
           cum_Taxa_score = cumsum(Taxa_score),
           cum_Taxa_score_rel = cum_Taxa_score / sum(Taxa_score)
    ) %>% 
    as.data.frame()
  
  top20_leafhopperID = leafhopperID %>% 
    slice_max(order_by = Taxa_score, n = 20)
  
  leafhopper.plot20 = ggplot(top20_leafhopperID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Leafhopper (top 20)",
      x = "",
      y = "Taxa Score"
    )
  leafhopper.plot20
  
  
  
minScore= 0.9

filter_top_taxa <- function(df, group_name, minScore) {
  df %>%
    { if (nrow(.) < 6) . 
      else filter(., cum_Taxa_score_rel <= minScore) } %>%
    mutate(Group = group_name)
}

arthropodID <- bind_rows(
  filter_top_taxa(antID, "Ant", minScore),
  filter_top_taxa(caterpillarID, "Caterpillar", minScore),
  filter_top_taxa(beetleID, "Beetle", minScore),
  filter_top_taxa(spiderID, "Spider", minScore),
  filter_top_taxa(leafhopperID, "Leafhopper", minScore),
  filter_top_taxa(truebugsID, "Truebugs", minScore),
  filter_top_taxa(flyID, "Fly", minScore),
  filter_top_taxa(grasshopperID, "Grasshopper", minScore),
  filter_top_taxa(daddylonglegsID, "Daddylonglegs", minScore)
) %>%
  mutate(
    Taxa_score_rel.100 = case_when(
      Group == "Ant" ~ Taxa_score / sum(Taxa_score[Group == "Ant"]) * 100,
      Group == "Caterpillar" ~ Taxa_score / sum(Taxa_score[Group == "Caterpillar"]) * 100,
      Group == "Beetle" ~ Taxa_score / sum(Taxa_score[Group == "Beetle"]) * 100,
      Group == "Spider" ~ Taxa_score / sum(Taxa_score[Group == "Spider"]) * 100,
      Group == "Leafhopper" ~ Taxa_score / sum(Taxa_score[Group == "Leafhopper"]) * 100,
      Group == "Truebugs" ~ Taxa_score / sum(Taxa_score[Group == "Truebugs"]) * 100,
      Group == "Fly" ~ Taxa_score / sum(Taxa_score[Group == "Fly"]) * 100,
      Group == "Grasshopper" ~ Taxa_score / sum(Taxa_score[Group == "Grasshopper"]) * 100,
      Group == "Daddylonglegs" ~ Taxa_score / sum(Taxa_score[Group == "Daddylonglegs"]) * 100,
      TRUE ~ NA_real_
    )
  )



arthropodID %>% write.csv(file = "data\\topID.csv",
                          row.names = FALSE)

sum(arthropodID [arthropodID$Group == "Spider",]$Taxa_score_rel.100)
 


feed_top <- read_excel("data/topId feeding guild.xlsx") %>% 
  select(-c(Rank, Group,))

topID_feed <- left_join(arthropodID, feed_top, by = "Taxon") %>% 
  mutate(Herbivore.score = Herbivore *Taxa_score_rel.100,
         Detritivore.score = Detritivore * Taxa_score_rel.100,
         Scavangers.score = Scavangers * Taxa_score_rel.100,
         Predator.score = Predator * Taxa_score_rel.100)

summary.topID_feed  <- topID_feed %>% 
  group_by(Group) %>% 
  summarise(Herbivore.score = sum(Herbivore.score),
            Detritivore.score = sum(Detritivore.score),
            Scavangers.score = sum(Scavangers.score),
            Predator.score = sum(Predator.score)
            ) %>% 
  arrange(desc(Herbivore.score))

topID_feed %>% 
  select(Group, Herbivore, Herbivore.score, Taxa_score_rel.100)

topID_feed %>% select(Group, Herbivore.score)




arthropod_ranks_update = summary.topID_feed %>% 
  select(Group, Herbivore.score) %>% 
  mutate(Rank_update = rownames(summary.topID_feed),
         ) %>% as.data.frame()

arthropod_ranks_update
# prepare data for wilcoxon ranked sum test. 
# I do not trust this part of the analysis.

rank.data <- left_join(arthropod_ranks_update %>% select(Group, Rank_update, Herbivore.score), 
                       species_score.d %>% 
                         select(Group, Rank) %>% 
                         mutate(Group = c("Caterpillar",
                                          "Spider",
                                          "Beetle",
                                          "Truebugs",
                                          "Leafhopper",
                                          "Ant",
                                          "Fly",
                                          "Grasshopper",
                                          "Daddylonglegs"
                                          )),
                       by = "Group") # from the RDA result

wilcox.test(
  as.numeric(rank.data$Rank_update),
  as.numeric(rank.data$Rank),
  exact = FALSE
)




