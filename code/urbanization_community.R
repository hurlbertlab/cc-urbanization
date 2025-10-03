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
    aphid         = ifelse(sum(Group == 'aphid', na.rm = TRUE) > 0, 1, 0)) %>%
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
    Aphid         = sum(aphid)         / n_distinct(ID),
    surveyNum       = n_distinct(ID)) 

# number of distinct survey IDs for each sites (by observation methods)
#  Where there is at least one observation of arthropods. 

prop_data <- left_join(proportion, sites, by = "Name")






#   PCA

prop.num <- prop_data[,3:12]
site.info <- prop_data[,c(1,2,13:18)]



# rownames(omo_physic_PCA) <- omo_physic_PCA$Location


prop.pca <- PCA(prop.num, graph = FALSE)

eig.val_prop <- get_eigenvalue(prop.pca)
eig.val_prop

fviz_eig(prop.pca, addlabels = TRUE, ylim = c(0, 50)) # scree plot

fviz_pca_var(prop.pca, col.var = "black", repel = TRUE)

fviz_pca_var(prop.pca, col.var = "cos2",
             gradient.cols = c("red", "#E7B800", "darkgreen"),
             repel = TRUE # Avoid text overlapping
)

fviz_pca_ind(prop.pca, col.ind = "cos2", pointsize = "cos2",
             pointshape = 21, fill = "#E7B800", 
             repel = TRUE) # Avoid text overlapping (slow if many points)

fviz_pca_var(prop.pca, col.var = "contrib",
             gradient.cols = c("red", "#E7B800", "darkgreen"),
             repel = TRUE # Avoid text overlapping
)
fviz_contrib(prop.pca, choice = "var", axes = 1)
fviz_contrib(prop.pca, choice = "var", axes = 2)
fviz_contrib(prop.pca, choice = "var", axes = 3)



coordinate_ind <- as.data.frame(omo.pca$ind$coord)
coordinate_ind12<- coordinate_ind[,1:2]   

coordinate_ind123 <- as.data.frame(cbind((omo[,1:3]),
                                         coordinate_ind12))
head(coordinate_ind123)

quality_of_rep_ind <- as.data.frame(omo.pca$ind$cos2)
quality_of_rep_ind12<- quality_of_rep_ind[,1:2]   
quality_of_rep_ind12<- quality_of_rep_ind12 %>% 
  rename(CO2.Dim.1 = Dim.1, 
         CO2.Dim.2 = Dim.2)
ind_viz_data <- cbind(coordinate_ind123, quality_of_rep_ind12)
head(ind_viz_data)


co2_var_dataframe12 <-  as.data.frame(omo.pca$var$coord[,1:2])

co2_var<- cbind(row.names(co2_var_dataframe12),co2_var_dataframe12) 
row.names(co2_var)  <- NULL     
head(co2_var) 

co2_var <- co2_var %>% 
  rename(parameter = "row.names(co2_var_dataframe12)",
         PCA1 = Dim.1, 
         PCA2 = Dim.2)
head(co2_var) 
head(co2_ind_viz_data)
co2_var<- as.data.frame(co2_var)
prop.pca$var$cos2
prop.pca$var$coord
prop.pca$ind$contrib
prop.pca$var$contrib
length(prop.pca$Location)


### using a RDA

# We would constrain the arthropod composition by the environmental variable (and maybe covary by latitude)

prop.num <- prop.num %>% as.data.frame()

prop.num_ilogit <- plogis(as.matrix(prop.num)) %>% as.data.frame
prop.num_ilogit <- 1 / (1 + exp(-prop.num))


site.z <- scale(prop_data[,c( "Latitude", "dev", "forest")]) %>% 
  as.data.frame()

Observ <- prop_data [, c("ObservationMethod","Longitude")]


prop.rda <- rda(prop.num_ilogit ~ Latitude * dev + forest, data = site.z)
summary(prop.rda)

anova.cca(prop.rda, step = 1000) # model is significant when compared to a NULL
anova.cca(prop.rda, step = 1000, by = 'axis') # only RDA1 and RDA2 is significant 
anova.cca(prop.rda, step = 1000, by = 'terms')

part.prop.rda <- rda(prop.num_ilogit, site.z, Observ) # adding observation method as covariate
summary(part.prop.rda)



part.prop.rda <- rda(
  prop.num ~ Latitude * dev + forest + Condition(ObservationMethod + Longitude),
  data = cbind(site.z, Observ))

summary(part.prop.rda)

RsquareAdj(part.prop.rda)$adj.r.squared
# Test whether the model is statistically significant
anova.cca(part.prop.rda, step = 999) # good!
anova.cca(part.prop.rda, step = 999, by = "axis")
anova.cca(part.prop.rda, step = 999, by= 'terms')
anova.cca(part.prop.rda, step = 999, by= 'margin')
ordiplot(part.prop.rda, scaling = 2, 
         main = "Arthropod Urbanization RDA - Scaling 2")

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





arth_urban_cor <- envfit(prop.num_ilogit, site.z, permutations = 999 ) 
arth_urban_cor

scores(arth_urban_cor, display = "vectors")



site_scores <- scores(mod, display = "sites")
correlations <- cor(prop.num_ilogit, site_score)  # correlation with RDA1 and RDA2
head(correlations)


# make an pRDA for each species and calculate its adjusted R squared.
species_var <- apply(prop.num, 2, function(sp) {
  mod_sp <- rda(sp ~ Latitude * dev + forest + Condition(ObservationMethod + Longitude),
                data = cbind(site.z, Observ))
  RsquareAdj(mod_sp)$adj.r.squared
})
species_var <- sort(species_var, decreasing = TRUE)










 
 
 
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
