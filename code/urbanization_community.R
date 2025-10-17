library(factoextra)
library(FactoMineR)
require(vegan)
require(tidyverse)

proportion <-fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
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

# number of distinct survey IDs for each sites (by observation methods)
#  Where there is at least one observation of arthropods. 

prop_data <- left_join(proportion, sites, by = "Name")





fullDataset %>%
  filter(Name %in% goodSites$Name,
         julianday %in% julianWindow,
         WetLeaves == 0,
         !Name %in% c('Coweeta - BS', 'Coweeta - BB', 'Coweeta - RK'),
         Group %in% c('caterpillar', 'spider', 'ant', 'leafhopper', 'beetle', 
                      'truebugs', 'fly', 'grasshopper', 
                      'daddylonglegs', 'aphid')) %>%
  group_by(ID) %>%
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
    aphid         = ifelse(sum(Group == 'aphid', na.rm = TRUE) > 0, 1, 0),
    other         = ifelse(sum(Group == "other", na.rm = TRUE) > 0, 1, 0),
    Psocodea      = ifelse(sum(Group == "Psocodea", na.rm = TRUE) > 0, 1, 0),
    Trichoptera   = ifelse(sum(Group ==  "Trichoptera", na.rm = TRUE) > 0, 1, 0),
    moths         = ifelse(sum(Group ==  "moths", na.rm = TRUE) > 0, 1, 0),
    bee           = ifelse(sum(Group ==  "bee", na.rm = TRUE) > 0, 1, 0)) %>% str()

  


#   PCA

prop.num <- prop_data[,3:11]
site.info <- prop_data[,c(1,2,12:17)]





### using a RDA
# Read more: https://sites.google.com/site/mb3gustame/constrained-analyses/redundancy-analysis
# We would constrain the arthropod composition by the environmental variable (and maybe covary by latitude)

prop.num <- prop.num %>% as.data.frame()

prop.num_ilogit <- plogis(as.matrix(prop.num)) %>% as.data.frame
prop.num_ilogit


site.z <- scale(prop_data[,c( "Latitude", "dev", "forest")]) %>% 
  as.data.frame()

Observ <- prop_data [, c("ObservationMethod","Longitude")]


part.prop.rda <- rda(prop.num_ilogit, site.z, Observ) # adding observation method as covariate
summary(part.prop.rda)



part.prop.rda <-rda(
  prop.num_ilogit ~ Latitude*dev + forest*Latitude + Condition(ObservationMethod),
  data = cbind(site.z, Observ))


summary(part.prop.rda)

RsquareAdj(part.prop.rda)$adj.r.squared
# Test whether the model is statistically significant


anova.cca(part.prop.rda, step = 999) # good!
anova.cca(part.prop.rda, step = 999, by = "axis")
anova.cca(part.prop.rda, step = 999, by= 'terms')


anova.cca(part.prop.rda, step = 999, by= 'margin') 
# "by = Margin" Tests each term after accounting for all other terms 
#- i.e., its unique effect, also called the partial effect.
# Each term is evaluated as if it were the last term, controlling for everything else.

# Main effects (Latitude, dev, forest) disappear in marginal testing because the interactions   
# account for much of the variation associated with those predictors.

ordiplot(part.prop.rda, scaling = 2, 
         main = "Arthropod Urbanization RDA - Scaling 2")

vif.cca(part.prop.rda) # All vif are less than 3.3, so statistically it makes sense to leave 
# all predictors in the model. Ecologically, it makes sense that some species may be responding
# greater or lesser to urbanization or forestry (i.e, response to both effect are non-equidistant)


# Since we are interested in how species are responding to the urbanization gradient across sites,
# we use the scaling = "species" argument.

# If the question was, 'how does sites differ from each other across urbanization gradient, 
# based on species composition?', we should use scaling = "sites" argument.
scores(part.prop.rda, display = "species", choices = 1:2, scaling = "species")
scores(part.prop.rda, display = "sites", scaling = "species")
scores(part.prop.rda, display = "bp", choices = 1:3, 
       scaling = "species")  # "bp" = biplot arrows



species_score <- scores(part.prop.rda, 
                        display = "species",
                        choices = 1:2) %>% 
  as.data.frame()
site_score <- scores(part.prop.rda, 
                     display = "sites",
                     choices = 1:2) %>% 
  as.data.frame()

site_score.sites  <- cbind(prop_data [,c("Name", "Region", "Longitude", 
                                      "Latitude", "dev", "forest")], site_score)

rda_axis <- scores(part.prop.rda, display = "bp", choices = 1:3) %>% 
  as.data.frame()



# other analysis you may skip
# Begin....

arth_urban_cor <- envfit(prop.num_ilogit, site.z, permutations = 999 ) 
arth_urban_cor

scores(arth_urban_cor, display = "vectors")




correlations <- cor(prop.num_ilogit, site_score)  # correlation with RDA1 and RDA2
head(correlations)


#..... END!







 
 
 
 ggplot() +
   geom_segment(data = rda_axis,
                aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                arrow = arrow(length = unit(0.3, "cm")), 
                color = "darkblue",
                size = 1, alpha = .5) +
   geom_point(data = site_score.sites, 
              aes(x = RDA1, y = RDA2, 
                  color = dev, fill = dev, alpha = 0.5), 
              shape = 21, size = 4) +
  
   geom_text(data = part.prop.rda$CCA$v, 
             aes(x = RDA1, y = RDA2, label = rownames(part.prop.rda$CCA$v)), 
             color = "black") +
   geom_text(data = rda_axis, 
             aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
             color = "darkblue", vjust = -0.5, hjust = 0.1, size =5) +
   scale_fill_gradientn(
     colours = c("green", "yellow", "red"),
     name = "Development"
   ) +
   labs(x = "RDA1", y = "RDA2") +
   guides(color = "none", shape = "none", alpha = "none",
          fill = guide_colorbar(title = "% Development")) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
   geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
   theme_minimal()
 
 
 
 # The sites here should not be interpeted to mean anything, really.
 ggplot() +
   geom_segment(data = rda_axis,
                aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                arrow = arrow(length = unit(0.3, "cm")), 
                color = "darkblue",
                size = 1, alpha = .5) +
   geom_text(data = site_score.sites, 
             aes(x = RDA1, y = RDA2, label = Region, colour = dev), 
            # color = "blue", 
             size = 3) + 
   geom_text(data = part.prop.rda$CCA$v, 
             aes(x = RDA1, y = RDA2, label = rownames(part.prop.rda$CCA$v)), 
             color = "black") + 
   geom_text(data = rda_axis, 
             aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
             color = "darkblue", vjust = -0.5, hjust = 0.1, size =5) +
   scale_color_gradientn(
     colours = c("green", "yellow", "red")) +
   theme(
     text = element_text(family = "Times New Roman", size = 20)
   ) + 
   labs(x = "RDA1", y = "RDA2") +
   theme_minimal() +
   guides(
     color = guide_legend(title = "Urban cover"), 
     fill = "none"  
   )+  
   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
   geom_vline(xintercept = 0, linetype = "dashed", color = "black")
 


 ggplot() +
   geom_point(data = site_score.sites, 
              aes(x = RDA1, y = RDA2, 
                  color = Latitude,   
                  fill = Latitude,
                  alpha = 0.5), 
              shape = 21,  # important: to use both color and fill
              size = 4) + 
   geom_text(data = part.prop.rda$CCA$v, 
             aes(x = RDA1, y = RDA2, label = rownames(part.prop.rda$CCA$v)), 
             color = "black") + 
   scale_fill_gradientn(
     colours = c("green", "yellow", "red"),
   ) +
   theme(
     text = element_text(family = "Times New Roman", size = 20)
   ) + 
   labs(x = "RDA1", y = "RDA2") +
   theme_minimal() +
   guides(
     color = "none",   
     shape = "none",
     alpha = "none",
     fill = guide_colorbar(title = "Latitude")  # keep a nice fill legend
   )+  
   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
   geom_vline(xintercept = 0, linetype = "dashed", color = "black")
 

 
 
prop_num.site.score <- cbind(site_score.sites, prop.num )


ggplot() +
  geom_segment(data = rda_axis,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black",
               size = 1, alpha = .5) +
  geom_point(data = prop_num.site.score, 
             aes(x = RDA1, y = RDA2, 
                 color = Hopper,   
                 fill = Hopper,
                 alpha = 0.8), 
             shape = 21,  # important: to use both color and fill
             size = 4) + 
  geom_text(data = rda_axis, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("yellow", "darkgreen"),
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + 
  labs(x = "RDA1", y = "RDA2", title = "Hopper") +
  theme_minimal() +
  guides(
    color = "none",   
    shape = "none",
    alpha = "none",
    fill = guide_colorbar(title = "% occurence")  # keep a nice fill legend
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  annotation_raster(hopperImage, ymin = -1.2, ymax = -0.6, xmin = -1.1, xmax = -0.4)







ggplot() +
  geom_segment(data = rda_axis,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black",
               size = 1, alpha = .5) +
  geom_point(data = prop_num.site.score, 
             aes(x = RDA1, y = RDA2, 
                 color = Caterpillar,   
                 fill = Caterpillar,
                 alpha = 0.8), 
             shape = 21,  # important: to use both color and fill
             size = 4) + 
  geom_text(data = rda_axis, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("yellow", "darkgreen"),
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + 
  labs(x = "RDA1", y = "RDA2", title = "Caterpillar") +
  theme_minimal() +
  guides(
    color = "none",   
    shape = "none",
    alpha = "none",
    fill = guide_colorbar(title = "% occurence")  # keep a nice fill legend
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  annotation_raster(catImage, ymin = -1.2, ymax = -0.7, xmin = .2, xmax =1.1)





ggplot() +
  geom_segment(data = rda_axis,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black",
               size = 1, alpha = .5) +
  geom_point(data = prop_num.site.score, 
             aes(x = RDA1, y = RDA2, 
                 color = Ant,   
                 fill = Ant,
                 alpha = 0.8), 
             shape = 21,  # important: to use both color and fill
             size = 4) + 
  geom_text(data = rda_axis, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("yellow", "darkgreen"),
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + 
  labs(x = "RDA1", y = "RDA2", title = "Ant") +
  theme_minimal() +
  guides(
    color = "none",   
    shape = "none",
    alpha = "none",
    fill = guide_colorbar(title = "% occurence")  # keep a nice fill legend
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  annotation_raster(antImage, ymin = 1.1, ymax = 1.5, xmin = -0.7, xmax =-0.1)




ggplot() +
  geom_segment(data = rda_axis,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black",
               size = 1, alpha = .5) +
  geom_point(data = prop_num.site.score, 
             aes(x = RDA1, y = RDA2, 
                 color = Beetle,   
                 fill = Beetle,
                 alpha = 0.8), 
             shape = 21,  # important: to use both color and fill
             size = 4) + 
  geom_text(data = rda_axis, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("yellow", "darkgreen"),
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + 
  labs(x = "RDA1", y = "RDA2", title = "Beetle") +
  theme_minimal() +
  guides(
    color = "none",   
    shape = "none",
    alpha = "none",
    fill = guide_colorbar(title = "% occurence")  # keep a nice fill legend
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  annotation_raster(beetleImage, ymin = -0.8, ymax = -0.4, xmin = -0.5, xmax = -0.1)



ggplot() +
  geom_segment(data = rda_axis,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black",
               size = 1, alpha = .5) +
  geom_point(data = prop_num.site.score, 
             aes(x = RDA1, y = RDA2, 
                 color = Spider,   
                 fill = Spider,
                 alpha = 0.8), 
             shape = 21,  # important: to use both color and fill
             size = 4) + 
  geom_text(data = rda_axis, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("yellow", "darkgreen"),
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + 
  labs(x = "RDA1", y = "RDA2", title = "Spider") +
  theme_minimal() +
  guides(
    color = "none",   
    shape = "none",
    alpha = "none",
    fill = guide_colorbar(title = "% occurence")  
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  annotation_raster(spiderImage, ymin = -0.7, ymax = 0.0, xmin = -1.4, xmax = -0.8)




ggplot() +
  geom_segment(data = rda_axis,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black",
               size = 1, alpha = .5) +
  geom_point(data = prop_num.site.score, 
             aes(x = RDA1, y = RDA2, 
                 color = Truebug,   
                 fill = Truebug,
                 alpha = 0.8), 
             shape = 21,  
             size = 4) + 
  geom_text(data = rda_axis, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("yellow", "darkgreen"),
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + 
  labs(x = "RDA1", y = "RDA2", title = "Truebug") +
  theme_minimal() +
  guides(
    color = "none",   
    shape = "none",
    alpha = "none",
    fill = guide_colorbar(title = "% occurence")  
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  annotation_raster(truebugImage, ymin = 0.4, ymax = 0.9, xmin = 0, xmax = 0.4)




ggplot() +
  geom_segment(data = rda_axis,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black",
               size = 1, alpha = .5) +
  geom_point(data = cbind(prop_data[,c("surveyNum")], prop_num.site.score), 
             aes(x = RDA1, y = RDA2, 
                 color = Hopper,   
                 fill = Hopper,
                 alpha = 0.8,
                 size = log10(surveyNum+10)), 
             shape = 21,  # important: to use both color and fill
             ) + 
  geom_text(data = rda_axis, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("yellow", "darkgreen"),
  ) +
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + 
  labs(x = "RDA1", y = "RDA2", title = "Hopper") +
  theme_minimal() +
  guides(
    color = "none",   
    shape = "none",
    alpha = "none",
    size = "none",
    fill = guide_colorbar(title = "% occurence")  # keep a nice fill legend
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  annotation_raster(hopperImage, ymin = -1.2, ymax = -0.6, xmin = -1.1, xmax = -0.4)




# How does site differ based on arthropod composition response to urbanization?
# what arthropod groups are most associated with sites based on urbanization level and latitude

scores(part.prop.rda, display = "species", choices = 1:2, scaling = "sites")
scores(part.prop.rda, display = "sites", scaling = "sites")
scores(part.prop.rda, display = "bp", choices = 1:3, 
       scaling = "sites")  # "bp" = biplot arrows


# site = scale 1; species = scale 2; no scaling preference = scale 3


species_score3 <- scores(part.prop.rda, display = "species", 
                         choices = 1:2, scaling = 3) %>% 
  as.data.frame()
site_score3 <- scores(part.prop.rda, 
                     display = "sites",
                     choices = 1:2, scaling = 3) %>% 
  as.data.frame()

site_score.sites3  <- cbind(prop_data [,c("Name", "Region", "Longitude", 
                                         "Latitude", "dev", "forest")], site_score3)

rda_axis3 <- scores(part.prop.rda, display = "bp", choices = 1:3, scaling = 3) %>% 
  as.data.frame()

ggplot() +
  geom_segment(data = rda_axis3,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "darkblue",
               size = 1, alpha = .5) +
  geom_text(data = site_score.sites3, 
            aes(x = RDA1, y = RDA2, label = Region, colour = dev), 
            # color = "blue", 
            size = 3) + 
  geom_text(data = species_score3, 
            aes(x = RDA1, y = RDA2, label = rownames(species_score3)), 
            color = "black") + 
  geom_text(data = rda_axis3, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis3)), 
            color = "darkblue", vjust = -0.5, hjust = 0.1, size =5) +
  scale_color_gradientn(
    colours = c("green", "yellow", "red")) +
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + 
  labs(x = "RDA1", y = "RDA2", title = "Scaling = 3") +
  theme_minimal() +
  guides(
    color = guide_legend(title = "Urban cover"), 
    fill = "none"  
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")




####################################################################################



species_score1 <- scores(part.prop.rda, display = "species", 
                         choices = 1:2, scaling = 1) %>% 
  as.data.frame()
site_score1 <- scores(part.prop.rda, 
                      display = "sites",
                      choices = 1:2, scaling = 1) %>% 
  as.data.frame()

site_score.sites1  <- cbind(prop_data [,c("Name", "Region", "Longitude", 
                                          "Latitude", "dev", "forest")], site_score1)

rda_axis1 <- scores(part.prop.rda, display = "bp", choices = 1:3, scaling = 1) %>% 
  as.data.frame()

ggplot() +
  geom_segment(data = rda_axis1,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "darkblue",
               size = 1, alpha = .5) +
  geom_text(data = site_score.sites1, 
            aes(x = RDA1, y = RDA2, label = Region, colour = dev), 
            # color = "blue", 
            size = 3) + 
  geom_text(data = species_score1, 
            aes(x = RDA1, y = RDA2, label = rownames(species_score1)), 
            color = "black") + 
  geom_text(data = rda_axis1, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis1)), 
            color = "darkblue", vjust = -0.5, hjust = 0.1, size =5) +
  scale_color_gradientn(
    colours = c("green", "yellow", "red")) +
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + 
  labs(x = "RDA1", y = "RDA2", title = "Scaling = 1") +
  theme_minimal() +
  guides(
    color = guide_legend(title = "Urban cover"), 
    fill = "none"  
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")




# Kmeans clustering of species and sites


rda_sp.matrix <-species_score %>% 
  remove_rownames()

rda_sp.matrix <- as.data.frame(rda_sp.matrix)
rda_sp.matrix <- dplyr::select_if(rda_sp.matrix, is.numeric)

fviz_nbclust(rda_sp.matrix, kmeans, method = "wss", k.max =9) # decide optimum clusters?

km.sp. <- kmeans(rda_sp.matrix, 4, nstart = 25)
print(km.res)

sp.cluster <- cbind(species_score, cluster = km.sp.$cluster)
sp.cluster %>% arrange(cluster)

sp.cluster.unscaled = cbind(part.prop.rda$CCA$v,  cluster = km.sp.$cluster)

# The unscaled looks better for vizual
# Actually, I am not sure. it seems to do some thing to the interpretation
# Compare both


# unscaled
ggplot() +
  geom_segment(data = rda_axis,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "darkblue",
               size = 1, alpha = .5)+
  geom_point(data = sp.cluster.unscaled, aes(x = RDA1, y = RDA2, 
                                   color = factor(cluster), 
                                   fill = factor(cluster)), 
             size = 10) + 
  
  geom_text(data = sp.cluster.unscaled, 
            aes(x = RDA1, y = RDA2, label = rownames(sp.cluster.unscaled)))+ 
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + labs(x = "RDA1", y = "RDA2")+
  theme_minimal()+
  guides(
    color = "none",   
    shape = "none", 
    fill = "none"   
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_text(data = rda_axis, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) 


# scaled
ggplot() +
  geom_segment(data = rda_axis,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "darkblue",
               size = 1, alpha = .5)+
  geom_point(data = sp.cluster, aes(x = RDA1, y = RDA2, 
                                             color = factor(cluster), 
                                             fill = factor(cluster)), 
             size = 10) + 
  
  geom_text(data = sp.cluster, 
            aes(x = RDA1, y = RDA2, label = rownames(sp.cluster)))+ 
  theme(
    text = element_text(family = "Times New Roman", size = 20)
  ) + labs(x = "RDA1", y = "RDA2")+
  theme_minimal()+
  guides(
    color = "none",   
    shape = "none", 
    fill = "none"   
  )+  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_text(data = rda_axis, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis)), 
            color = "black", vjust = -0.5, hjust = 0.1, size =5) 



set.seed(123)
km_res <- kmeans(rda_sp.matrix, centers = 4, nstart = 25)

# Get distances of each species to its assigned cluster center
dist_to_center <- sqrt(rowSums((rda_sp.matrix - km_res$centers[km_res$cluster, ])^2))

# Combine results into a dataframe
species_strength <- data.frame(
  species = rownames(rda_sp.matrix),
  cluster = km_res$cluster,
  dist_to_center = dist_to_center
)

# Rank species within each cluster by their distance (low = strong, high = weak)
species_strength <- species_strength %>%
  group_by(cluster) %>%
  arrange(dist_to_center, .by_group = TRUE)
#Does local urbanization influence arthropod community composition independently of where the 
# sites are located and how they were sampled?



part.prop.rda2 <- rda(
  prop.num ~ dev  +forest + Condition(ObservationMethod + Latitude),
  data = cbind(site.z, Observ))

summary(part.prop.rda2)

RsquareAdj(part.prop.rda2)$adj.r.squared
# Test whether the model is statistically significant
anova.cca(part.prop.rda2, step = 999) # good!
anova.cca(part.prop.rda2, step = 999, by = "axis")
anova.cca(part.prop.rda2, step = 999, by= 'terms')
anova.cca(part.prop.rda2, step = 999, by= 'margin')
ordiplot(part.prop.rda2, scaling = 2, 
         main = "Arthropod Urbanization RDA - Scaling 2")





part.prop.rda.df <- rda(
  prop.num ~ dev * forest + Condition(ObservationMethod + Latitude),
  data = cbind(site.z, Observ))

summary(part.prop.rda.df)
RsquareAdj(part.prop.rda.df)$adj.r.squared

anova.cca(part.prop.rda.df, step = 999) # good!
anova.cca(part.prop.rda.df, step = 999, by = "axis")
anova.cca(part.prop.rda.df, step = 999, by= 'terms')
anova.cca(part.prop.rda.df, step = 999, by= 'margin')
ordiplot(part.prop.rda.df, scaling = 2, 
         main = "Arthropod Urbanization RDA - Scaling 2")
vif.cca(part.prop.rda.df)



part.prop.rda.df$CCA$v



species_score.df <- scores(part.prop.rda.df, 
                           display = "species",
                           choices = 1:2) %>% 
  as.data.frame()
site_score.df <- scores(part.prop.rda.df, 
                        display = "sites",
                        choices = 1:2) %>% 
  as.data.frame()

site_score.sites.df  <- cbind(prop_data [,c("Name", "Region", "Longitude", 
                                            "Latitude", "dev", "forest")], site_score.df)

rda_axis.df <- scores(part.prop.rda.df, display = "bp", choices = 1:2) %>% 
  as.data.frame()




ggplot() +
  geom_segment(data = rda_axis.df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "darkblue",
               size = 1, alpha = .5) +
  geom_point(data = site_score.sites.df, 
             aes(x = RDA1, y = RDA2, 
                 color = dev, fill = dev, alpha = 0.5), 
             shape = 21, size = 4) +
  
  geom_text(data = part.prop.rda.df$CCA$v, 
            aes(x = RDA1, y = RDA2, label = rownames(part.prop.rda.df$CCA$v)), 
            color = "black") +
  geom_text(data = rda_axis.df, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis.df)), 
            color = "darkblue", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("green", "yellow", "red"),
    name = "Development"
  ) +
  labs(x = "RDA1", y = "RDA2") +
  guides(color = "none", shape = "none", alpha = "none",
         fill = guide_colorbar(title = "% Development")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_minimal()

