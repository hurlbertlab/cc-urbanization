


part.prop.rda.LD <-rda(
  prop.num_ilogit ~ Latitude*dev +  Condition(ObservationMethod),
  data = cbind(site.z, Observ))

summary(part.prop.rda.LD)
RsquareAdj(part.prop.rda.LD)$adj.r.squared
# Test whether the model is statistically significant


anova.cca(part.prop.rda.LD, step = 999) # good!
anova.cca(part.prop.rda.LD, step = 999, by = "axis")
anova.cca(part.prop.rda.LD, step = 999, by= 'terms')

anova.cca(part.prop.rda.LD, step = 999, by= 'margin')

vif.cca(part.prop.rda.LD) # Below 1.1


vpart.prop.LD <- varpart(prop.num_ilogit, site.z %>%select(-forest),  Observ)
vpart.prop.LD


species_score.ld <- scores(part.prop.rda.LD, 
                        display = "species",
                        choices = 1:2) %>% 
  as.data.frame()
site_score.ld <- scores(part.prop.rda.LD, 
                     display = "sites",
                     choices = 1:2) %>% 
  as.data.frame()

site_score.sites.ld  <- cbind(prop_data [,c("Name", "Region", "Longitude", 
                                         "Latitude", "dev", "forest")], site_score.ld)

rda_axis.ld <- scores(part.prop.rda.LD, display = "bp", choices = 1:3) %>% 
  as.data.frame()




ggplot() +
  geom_segment(data = rda_axis.ld,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "darkblue",
               size = 1, alpha = .5) +
  geom_point(data = site_score.sites.ld, 
             aes(x = RDA1, y = RDA2, 
                 color = dev, fill = dev, alpha = 0.5), 
             shape = 21, size = 4) +
  
  geom_text(data = part.prop.rda.LD$CCA$v, 
            aes(x = RDA1, y = RDA2, label = rownames(part.prop.rda.LD$CCA$v)), 
            color = "black") +
  geom_text(data = rda_axis.ld, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis.ld)), 
            color = "darkblue", vjust = -0.5, hjust = 0.1, size =5) +
  scale_fill_gradientn(
    colours = c("green", "lightyellow", "skyblue"),
    name = "Development"
  ) +
  labs(x = "RDA1", y = "RDA2") +
  guides(color = "none", shape = "none", alpha = "none",
         fill = guide_colorbar(title = "% Development")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_minimal()

##################################################################

Group <- c("Caterpillar", "Hopper", "Grasshopper",
           "Beetle", "fly", "Truebug",
           "Ant", "Spider", "Daddylongleg")

FeedingGuild <- c("Herbivore", "Herbivore", "Herbivore",
                  "Mixed",  "Mixed",
                  "Mixed", "Predator",
                  "Predator", "Predator")

# Assign ranks from herbivore → predator
Rank <- 1:9  # 1 = most herbivorous, 9 = most predatory

arthropod_ranks <- data.frame(Group, FeedingGuild, Rank)


###########################################################

# if latitude was partialed out, what would be the effect of urbanization?


part.prop.rda.D <- rda(
  prop.num_ilogit ~ dev + Condition(Latitude + ObservationMethod),
  data = cbind(site.z, Observ)
)


summary(part.prop.rda.D)
RsquareAdj(part.prop.rda.D)$adj.r.squared
# Test whether the model is statistically significant


anova.cca(part.prop.rda.D, step = 999) # good!
anova.cca(part.prop.rda.D, step = 999, by = "axis")
anova.cca(part.prop.rda.D, step = 999, by= 'terms')
anova.cca(part.prop.rda.D, step = 999, by= 'margin')
vif.cca(part.prop.rda.D) # Below 1.1


species_score.d <- scores(part.prop.rda.D, 
                           display = "species",
                           choices = 1) %>% 
  as.data.frame()




species_score.d$Group <- rownames(species_score.d)

species_score.d= remove_rownames(species_score.d) %>% 
  mutate(Rank = rownames(species_score.d)) %>% as.data.frame()


################################################################

###########################################################

left_join(species_score.d, arthropod_ranks, by = "Group") %>% 
  ggplot(aes(x = reorder(Group, RDA1), y = RDA1)) +
  geom_point(size = 4, aes( fill = FeedingGuild, colour = FeedingGuild)) +
  geom_segment(aes(xend = Group, y = 0, yend = RDA1),
               linewidth = 1) +
  coord_flip() +
  scale_fill_manual(
    values = c("red", "steelblue", "green"),
    name = "Feeding guild"
  )+
  scale_color_manual(
    values = c("red", "steelblue", "green"),
    name = "Feeding guild"
  )+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(x = "Arthropod group", y = "RDA1 (urbanization response)") +
  theme_minimal(base_size = 14)



left_join(species_score.d, arthropod_ranks, by = "Group") %>% 
  ggplot(aes(x = reorder(Group, RDA1), y = RDA1)) +
  geom_point(size = 4, aes( fill = Rank, colour = RDA1)) +
  geom_segment(aes(xend = Group, y = 0, yend = RDA1),
               linewidth = 1) +
  coord_flip() +
  scale_fill_gradientn(
    colours = c("blue", "grey", "blue"),
    name = "Response"
  )+
  scale_color_gradientn(
    colours = c("blue", "grey", "blue"),
    name = "Response"
  )+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(x = "Arthropod group", y = "RDA1 (urbanization response)") +
  theme_minimal(base_size = 14)


#################################################################

rda_sp.matrix <-species_score.ld %>%  # Lat * dev model
  remove_rownames()

rda_sp.matrix <- as.data.frame(rda_sp.matrix)
rda_sp.matrix <- dplyr::select_if(rda_sp.matrix, is.numeric)

fviz_nbclust(rda_sp.matrix, kmeans, method = "wss", k.max =8) # decide optimum clusters?

km.sp. <- kmeans(rda_sp.matrix, 4, nstart = 25)
print(km.res)

sp.cluster <- cbind(species_score.ld, cluster = km.sp.$cluster)
sp.cluster %>% arrange(cluster)

sp.cluster.unscaled = cbind(part.prop.rda.LD$CCA$v,  cluster = km.sp.$cluster)

# The unscaled looks better for visual
# Actually, I am not sure. it seems to do some thing to the interpretation
# Compare both




 

# unscaled
ggplot() +
  geom_segment(data = rda_axis.ld,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "darkblue",
               size = 1, alpha = .5)+
  geom_point(data = sp.cluster.unscaled, aes(x = RDA1, y = RDA2, 
                                            # color = factor(cluster), 
                                            # fill = factor(cluster)
                                            ), 
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
  geom_text(data = rda_axis.ld, 
            aes(x = RDA1, y = RDA2, label = rownames(rda_axis.ld)), 
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



##########################################################################################################
testtt = sim_slopes(cat.Dev.Latitude, pred = dev, modx = Latitude,modx.values = seq(34, 44, by = 0.2),
                                  johnson_neyman = FALSE, digits = 4)$slopes %>% as.data.frame() 


ggplot(testtt, aes(x = Est., 
                      y = "Value of Latitude" , 
                      xmin = X2.5., 
                      xmax = X97.5.)) +
  geom_point(size = 2, aes(colour = Est.)) +
  geom_errorbarh(width = 0.02, size = 1, aes(colour = Est.)) +   # horizontal error bars
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2)+
  scale_color_gradientn(
    colours = c("red", "#424242", "green"),
  ) +
  labs(
    x = "Estimated effect of urban cover",
    y = "Latitude",
   # title = "Effect of Urban Cover on Arthropod Groups by Latitude",
   # subtitle = "Points = estimates. Estimate bars (=95% CI) touching the dashed line are not statistically significant."
  ) +
  guides(color = "none")+
  scale_y_continuous(limits = c(34, 45))+
  scale_x_continuous(limits = c(min(Dev_sims$X2.5.), max(Dev_sims$X97.5.)))+
  theme_bw()



testsimDev_5caterpillar <- sim_slopes(cat.Dev.Latitude, pred = dev, 
                                      modx = Latitude,modx.values = seq(34, 44, by = 1),
                                  johnson_neyman = FALSE, digits = 4)

testsimDev_5beetle <-sim_slopes(beet.Dev.Latitude, pred = dev, 
                                modx = Latitude,modx.values = seq(34, 44, by = 1),
                            johnson_neyman = FALSE, digits = 4)

testsimDev_5spider <-sim_slopes(spi.Dev.Latitude, pred = dev, 
                                modx = Latitude,modx.values = seq(34, 44, by = 1),
                            johnson_neyman = FALSE, digits = 4)

testsimDev_5hopper <-sim_slopes(hop.Dev.Latitude, pred = dev, 
                                modx = Latitude,modx.values = seq(34, 44, by = 1),
                            johnson_neyman = FALSE, digits = 4)

testsimDev_5ant <-sim_slopes(ant.Dev.Latitude, pred = dev, 
                             modx = Latitude,modx.values = seq(34, 44, by = 1),
                         johnson_neyman = FALSE, digits = 4)

testsimDev_5truebug <-sim_slopes(ant.Dev.Latitude, pred = dev, 
                                 modx = Latitude,modx.values = seq(34, 44, by = 1),
                             johnson_neyman = FALSE, digits = 4)


testDev_sims5 = data.frame(rbind(testsimDev_5caterpillar$slopes %>% mutate(Group = "caterpillar"), 
                                 testsimDev_5hopper$slopes %>% mutate(Group = "hopper"), 
                                 testsimDev_5beetle$slopes %>% mutate(Group = "beetle"), 
                                 testsimDev_5truebug$slopes %>% mutate(Group = "truebug"),
                                 testsimDev_5ant$slopes %>% mutate(Group = "ant"), 
                                 testsimDev_5spider$slopes %>% mutate(Group = "spider"))) %>% 
  mutate(Group = factor(Group, levels = c("caterpillar", "beetle", "ant",
                                          "hopper", "truebug", "spider")))



ggplot(testDev_sims5, aes(x = Est., 
                      y = Value.of.Latitude, 
                      xmin = X2.5., 
                      xmax = X97.5.)) +
  
  geom_errorbarh(width = 0, size =1, aes(colour = Est.), alpha = 0.5) +   # horizontal error bars
  geom_point(size = 1, color = "black", alpha = 0.3,) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2)+
  scale_color_gradientn(
    colours = c("red", "darkblue", "green"),
  ) +
  facet_wrap(~ Group, scales = "free_x") +  
  labs(
    x = "Estimated effect of urban cover",
    y = "Latitude",
    title = "Effect of Urban Cover on Arthropod Groups by Latitude",
    subtitle = "Points = estimates. Estimate bars (=95% CI) touching the dashed line are not statistically significant."
  ) +
  guides(color = "none")+
  scale_y_continuous(limits = c(33, 45))+
  scale_x_continuous(limits = c(min(testDev_sims5$X2.5.), max(testDev_sims5$X97.5.)))+
  theme_bw()





testsimFor_5caterpillar <- sim_slopes(cat.For.Latitude, pred = forest, 
                                      modx = Latitude,modx.values = seq(34, 44, by = 1),
                                      johnson_neyman = FALSE, digits = 4)

testsimFor_5beetle <-sim_slopes(beet.For.Latitude, pred = forest, 
                                modx = Latitude,modx.values = seq(34, 44, by = 1),
                                johnson_neyman = FALSE, digits = 4)

testsimFor_5spider <-sim_slopes(spi.For.Latitude, pred = forest, 
                                modx = Latitude,modx.values = seq(34, 44, by = 1),
                                johnson_neyman = FALSE, digits = 4)

testsimFor_5hopper <-sim_slopes(hop.For.Latitude, pred = forest, 
                                modx = Latitude,modx.values = seq(34, 44, by = 1),
                                johnson_neyman = FALSE, digits = 4)

testsimFor_5ant <-sim_slopes(ant.For.Latitude, pred = forest, 
                             modx = Latitude,modx.values = seq(34, 44, by = 1),
                             johnson_neyman = FALSE, digits = 4)

testsimFor_5truebug <-sim_slopes(ant.For.Latitude, pred = forest, 
                                 modx = Latitude,modx.values = seq(34, 44, by = 1),
                                 johnson_neyman = FALSE, digits = 4)


testFor_sims5 = data.frame(rbind(testsimFor_5caterpillar$slopes %>% mutate(Group = "caterpillar"), 
                                 testsimFor_5hopper$slopes %>% mutate(Group = "hopper"), 
                                 testsimFor_5beetle$slopes %>% mutate(Group = "beetle"), 
                                 testsimFor_5truebug$slopes %>% mutate(Group = "truebug"),
                                 testsimFor_5ant$slopes %>% mutate(Group = "ant"), 
                                 testsimFor_5spider$slopes %>% mutate(Group = "spider"))) %>% 
  mutate(Group = factor(Group, levels = c("caterpillar", "beetle", "ant",
                                          "hopper", "truebug", "spider")))


ggplot(testFor_sims5, aes(x = Est., 
                          y = Value.of.Latitude, 
                          xmin = X2.5., 
                          xmax = X97.5.)) +
  
  geom_errorbarh(width = 0, size =1, aes(colour = Est.), alpha = 0.5) +   # horizontal error bars
  geom_point(size = 1, color = "black", alpha = 0.3,) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2)+
  scale_color_gradientn(
    colours = c("red", "darkblue", "green"),
  ) +
  facet_wrap(~ Group, scales = "free_x") +  
  labs(
    x = "Estimated effect of Forest cover",
    y = "Latitude",
    title = "Effect of Urban Cover on Arthropod Groups by Latitude",
    subtitle = "Points = estimates. Estimate bars (=95% CI) touching the dashed line are not statistically significant."
  ) +
  guides(color = "none")+
  scale_y_continuous(limits = c(34, 45))+
 scale_x_continuous(limits = c(min(testFor_sims5$X2.5.), max(testFor_sims5$X97.5.)))+
  theme_bw()




###########################################################################################

########### ----------------- Combined-------------------------------------------------


testsimDev_6caterpillar <- sim_slopes(cat.Dev.Latitude, pred = dev, 
                                      modx = Latitude,modx.values = seq(34.2, 44.2, by = 1),
                                      johnson_neyman = FALSE, digits = 4)

testsimDev_6beetle <-sim_slopes(beet.Dev.Latitude, pred = dev, 
                                modx = Latitude,modx.values = seq(34.2, 44.2, by = 1),
                                johnson_neyman = FALSE, digits = 4)

testsimDev_6spider <-sim_slopes(spi.Dev.Latitude, pred = dev, 
                                modx = Latitude,modx.values = seq(34.2, 44.2, by = 1),
                                johnson_neyman = FALSE, digits = 4)

testsimDev_6hopper <-sim_slopes(hop.Dev.Latitude, pred = dev, 
                                modx = Latitude,modx.values = seq(34.2, 44.2, by = 1),
                                johnson_neyman = FALSE, digits = 4)

testsimDev_6ant <-sim_slopes(ant.Dev.Latitude, pred = dev, 
                             modx = Latitude,modx.values = seq(34.2, 44.2, by = 1),
                             johnson_neyman = FALSE, digits = 4)

testsimDev_6truebug <-sim_slopes(ant.Dev.Latitude, pred = dev, 
                                 modx = Latitude,modx.values = seq(34.2, 44.2, by = 1),
                                 johnson_neyman = FALSE, digits = 4)



testDev_sims6 = data.frame(rbind(testsimDev_6caterpillar$slopes %>% mutate(Group = "caterpillar"), 
                                 testsimDev_6hopper$slopes %>% mutate(Group = "hopper"), 
                                 testsimDev_6beetle$slopes %>% mutate(Group = "beetle"), 
                                 testsimDev_6truebug$slopes %>% mutate(Group = "truebug"),
                                 testsimDev_6ant$slopes %>% mutate(Group = "ant"), 
                                 testsimDev_6spider$slopes %>% mutate(Group = "spider"))) %>% 
  mutate(Group = factor(Group, levels = c("caterpillar", "beetle", "ant",
                                          "hopper", "truebug", "spider")))


testDev_sims7 = testDev_sims6 %>% 
  mutate(Est. = 0-Est.,
         X2.5. = 0- X2.5.,
         X97.5. = 0 - X97.5.) %>% 
  mutate(Cover = "Urban cover")


ggplot(testDev_sims7, aes(x = Est., 
                          y = Value.of.Latitude, 
                          xmin = X2.5., 
                          xmax = X97.5.)) +
  geom_errorbarh(width = 0, size =1, color = "blue", alpha = 0.5) +   # horizontal error bars
  geom_point(size = 1, color = "blue", alpha = 0.3,) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2)+
  facet_wrap(~ Group, 
             scales = "free_x"
             ) +  
  labs(
    x = "Estimated effect of Forest cover",
    y = "Latitude",
    title = "Effect of Urban Cover on Arthropod Groups by Latitude",
    subtitle = "Points = estimates. Estimate bars (=95% CI) touching the dashed line are not statistically significant."
  ) +
  guides(color = "none")+
  scale_y_continuous(limits = c(34, 45))+
  #scale_x_continuous(limits = c(max(testDev_sims7$X2.5.), min(testDev_sims7$X97.5.)))+
  theme_bw()




Dev.For = rbind(testDev_sims6%>% mutate(Cover = "Urban cover"),
                testFor_sims5%>% 
                  mutate(Est. = 0-Est.,
                         X2.5. = 0- X2.5.,
                         X97.5. = 0 - X97.5.) %>% mutate(Cover = "Forest cover"))



ggplot(Dev.For, aes(x = Est., 
                          y = Value.of.Latitude, 
                          xmin = X2.5., 
                          xmax = X97.5.)) +
  geom_errorbarh(width = 0, size =1.1, aes(colour = Cover), alpha = 0.5) +   # horizontal error bars
  geom_point(size = 1, aes(colour = Cover), alpha = 0.3,) +
  geom_vline(xintercept = 0,  color = "black", linewidth = 0.2)+
  coord_flip()+
  facet_wrap(~ Group, scales = "free_x") +  
  scale_colour_manual(values = c("Forest cover" = "green", 
                                 "Urban cover"= "blue"))+
  labs(
    x = "Estimated effect",
    y = "Latitude"
  ) +
  guides(color = "none")+
  scale_y_continuous(limits = c(34, 45))+
  scale_x_continuous(limits = c(min(Dev.For$X2.5. -0.005), max(Dev.For$X97.5.)))+
  theme_bw()





ggplot(Dev.For, aes(x = Value.of.Latitude, y = Est.)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = Cover),
              alpha = 0.05, colour = NA) +
  geom_line(aes(colour = Cover), size = 1) +
  scale_colour_manual(values = c("Forest cover" = "green", "Urban cover" = "blue")) +
  scale_fill_manual(values = c("Forest cover" = "green", "Urban cover" = "blue")) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(x = "Latitude", y = "Estimated effect") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
  theme_bw()


# add background colour to differentiate response direction
ggplot(Dev.For, aes(x = Value.of.Latitude, y = Est.)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),
            fill = "mistyrose", alpha = 0.03) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "honeydew", alpha = 0.05) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = Cover),
              alpha = 0.05, colour = NA) +
  geom_line(aes(colour = Cover), size = 1) +
  scale_colour_manual(values = c("Forest cover" = "green", "Urban cover" = "blue")) +
  scale_fill_manual(values = c("Forest cover" = "green", "Urban cover" = "blue")) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(x = "Latitude", y = "Estimated effect of landcover change") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
  theme_bw()


##############3


Dev.For.2 <- rbind(testFor_sims5 %>% mutate(Cover = "Forest cover"),
                   testDev_sims5 %>% mutate(Cover = "Urban cover")) %>% as.data.frame()


ggplot(Dev.For.2, aes(x = Value.of.Latitude, y = Est.)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),
            fill = "#E6F0FF", alpha = 0.03) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "#E6F4E6", alpha = 0.05) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = Cover),
              alpha = 0.05, colour = NA) +
  geom_line(aes(colour = Cover), size = 1) +
  scale_colour_manual(values = c("Forest cover" = "green", "Urban cover" = "blue")) +
  scale_fill_manual(values = c("Forest cover" = "green", "Urban cover" = "blue")) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(x = "Latitude", y = "Estimated effect of landcover change") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.2) +
  theme_bw()









ggplot(testDev_sims5, aes(x = Est., 
                          y = Value.of.Latitude, 
                          xmin = X2.5., 
                          xmax = X97.5.)) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "mistyrose", alpha = 0.025) +
  
  # Your main layers
  #geom_errorbarh(width = 0, size =1, color = "blue", alpha = 0.5) +
  geom_line(size = 1, color = "blue", alpha = 0.7) +
  geom_ribbon(aes(xmin = X2.5., xmax = X97.5.),
              alpha = 0.05, colour = NA) +
  geom_vline(xintercept = 0,  color = "black", linewidth = 0.3) +
  
  facet_wrap(~ Group, scales = "free_y") +
  labs(
    x = "Estimated effect",
    y = "Latitude",
    title = "Effect of Urbanization on Arthropod Groups",
    subtitle = "Red = + Response; White = - Response"
  ) +
  guides(color = "none") +
 # scale_y_continuous(limits = c(34, 45)) +
  theme_bw()+
  coord_flip()

 