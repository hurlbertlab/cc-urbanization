# Summary of expert identifications

library(patchwork)

exp.id = read.csv("data/exp.csv")
View(exp.id)


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


idClass%>% write.csv(file = "C:\\Users\\osawe\\Documents\\Git\\cc-urbanization\\data\\expID1.csv",
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
   slice_max(order_by = Taxa_score, n = 10)
 
 spiderID.plot10 <- ggplot(top10_spiderID, aes(x = reorder(Taxon, -Taxa_score), y = Taxa_score)) +
   geom_bar(stat = "identity", alpha = 0.6, aes(colour = Rank, fill = Rank)) +
   theme_minimal(base_size = 12) +
   theme(
     axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
   ) +
   labs(
     title = "Top 10 Spider Taxa by Taxa Score",
     x = "",
     y = "Taxa Score"
   )
 

 print(spiderID.plot+guides(color = "none", fill = "none", title = "none")+
         labs(title = NULL) + 
         theme_minimal()+
         theme(
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.title.y = element_blank())+
         spiderID.plot10 +
         plot_layout(widths = c(1, 4)))
 
 
  # Ant
  antID <- idClass %>% 
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
  
  # Scree plot
  ggplot(antID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Ant Taxa",
      x = "",
      y = "Taxa Score"
    )

  
  ggplot(antID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Ant Taxa",
      x = "",
      y = "Taxa Score"
    )
  
  
  # beetle
  beetleID <- idClass %>% 
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
  
  ggplot(beetleID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Beetle Taxa",
      x = "",
      y = "Taxa Score"
    )
  
  
  
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
  
  ggplot(caterpillarID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Caterpillar Taxa",
      x = "",
      y = "Taxa Score"
    )
  
  
  
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
  
  ggplot(daddylonglegsID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Daddylonglegs Taxa",
      x = "",
      y = "Taxa Score"
    )
  
  
  
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
  
  ggplot(truebugsID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Truebugs Taxa",
      x = "",
      y = "Taxa Score"
    )

  
  
  
  # grasshopper
   grasshopperID <- idClass %>% 
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
  
  ggplot(grasshopperID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Grasshopper Taxa",
      x = "",
      y = "Taxa Score"
    )
  
  
  
  
  # fly
  flyID <- idClass %>% 
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
  
  ggplot(flyID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Fly Taxa",
      x = "",
      y = "Taxa Score"
    )
  
  
  
  # leafhopper
  leafhopperID <- idClass %>% 
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
  
  ggplot(leafhopperID, aes(x = Taxon, y = Taxa_score)) +
    geom_bar(stat = "identity", alpha = 0.6, 
             aes(colour = Rank, fill = Rank), ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)
    ) +
    labs(
      title = "Scree Plot of Leafhopper Taxa",
      x = "",
      y = "Taxa Score"
    )
  
  
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


arthropodID %>% write.csv(file = "C:\\Users\\osawe\\Documents\\Git\\cc-urbanization\\data\\topID.csv",
                          row.names = FALSE)

sum(arthropodID [arthropodID$Group == "Spider",]$Taxa_score_rel.100)
 

library(readxl)
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




