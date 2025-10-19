# Summary of expert identifications


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



spiderID= idClass %>% 
  filter(Group== "spider") %>% 
  group_by(Rank, Taxon) %>% 
  summarise(nofCount = sum(nofCount),
            Prop_count = nofCount/sum(idClass[idClass$Group=="spider",]$nofCount),
            nSite = n_distinct(Name), # calc. proportion of sites with each taxa from total sites
            prop_site = nSite/length(unique(idClass[idClass$Group=="spider",]$Name)),
            Taxa_score = Prop_count * prop_site) %>% 
  arrange(desc(Taxa_score)) %>%
  mutate(Taxon = factor(Taxon, levels = Taxon)) %>% 
  as.data.frame()
  

  ggplot(spiderID, aes(x = Taxon, y = Taxa_score)) +
  geom_bar(stat = "identity",
         #  aes(colour = Rank, fill = Rank), 
           alpha = 0.5) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1)) +
  labs(
    title = "Scree Plot of Spider Taxa",
    x = "",
    y = "Score"
  )


# Spider
  spiderID <- idClass %>% 
    filter(Group == "spider") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "spider", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site = nSite / length(unique(idClass[idClass$Group == "spider", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
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


  
  # Ant
  antID <- idClass %>% 
    filter(Group == "ant") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "ant", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site = nSite / length(unique(idClass[idClass$Group == "ant", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
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


  
  # Ant
  antID <- idClass %>% 
    filter(Group == "ant") %>% 
    group_by(Rank, Taxon) %>% 
    summarise(
      nofCount = sum(nofCount),
      Prop_count = nofCount / sum(idClass[idClass$Group == "ant", ]$nofCount),
      nSite = n_distinct(Name), 
      prop_site = nSite / length(unique(idClass[idClass$Group == "ant", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
    as.data.frame()
  
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
      prop_site = nSite / length(unique(idClass[idClass$Group == "beetle", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
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
      prop_site = nSite / length(unique(idClass[idClass$Group == "caterpillar", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
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
      prop_site = nSite / length(unique(idClass[idClass$Group == "daddylonglegs", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
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
      prop_site = nSite / length(unique(idClass[idClass$Group == "truebugs", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
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
      prop_site = nSite / length(unique(idClass[idClass$Group == "grasshopper", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
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
      prop_site = nSite / length(unique(idClass[idClass$Group == "fly", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
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
      prop_site = nSite / length(unique(idClass[idClass$Group == "leafhopper", ]$Name)),
      Taxa_score = Prop_count * prop_site,
      .groups = "drop"
    ) %>% 
    arrange(desc(Taxa_score)) %>% 
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>% 
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
  
  
  
  
    