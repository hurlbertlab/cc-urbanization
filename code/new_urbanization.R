library(ggrepel)

# Developed cover models
cat.Dev = glm(caterpillar ~ dev +  ObservationMethod, 
                       data = dataset, family = "binomial")

spi.Dev = glm(spider ~ dev + ObservationMethod, 
                       data = dataset, family = "binomial")

beet.Dev = glm(beetle ~ dev + ObservationMethod, 
                        data = dataset, family = "binomial")

hop.Dev = glm(hopper ~ dev + ObservationMethod, 
                       data = dataset, family = "binomial")

bug.Dev = glm(truebug ~ dev + ObservationMethod, 
                       data = dataset, family = "binomial")

ant.Dev = glm(ant ~ dev + ObservationMethod, 
                       data = dataset, family = "binomial")


# Forest cover models
cat.For = glm(caterpillar ~ forest+ ObservationMethod, 
                       data = dataset, family = "binomial")

spi.For = glm(spider ~ forest + ObservationMethod, 
                       data = dataset, family = "binomial")

beet.For = glm(beetle ~ forest + ObservationMethod, 
                        data = dataset, family = "binomial")

hop.For = glm(hopper ~ forest + ObservationMethod, 
                       data = dataset, family = "binomial")

bug.For = glm(truebug ~ forest + ObservationMethod, 
                       data = dataset, family = "binomial")

ant.For = glm(ant ~ forest + ObservationMethod, 
                       data = dataset, family = "binomial")



# GLM output

devOnlyOutput = data.frame(rbind(summary(cat.Dev)$coefficients, 
                             summary(spi.Dev)$coefficients, 
                             summary(beet.Dev)$coefficients,
                             summary(hop.Dev)$coefficients, 
                             summary(bug.Dev)$coefficients,
                             summary(ant.Dev)$coefficients))
devOnlyOutput$term = rep(c('Intercept', 'dev', 'Method'), times = 6)
devOnlyOutput$Group = rep(c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'), each = 3)
rownames(devOnlyOutput) = NULL

forOnlyOutput = data.frame(rbind(summary(cat.For)$coefficients, 
                             summary(spi.For)$coefficients, 
                             summary(beet.For)$coefficients,
                             summary(hop.For)$coefficients, 
                             summary(bug.For)$coefficients,
                             summary(ant.For)$coefficients))
forOnlyOutput$term = rep(c('Intercept', 'forest', 'Method'), times = 6)
forOnlyOutput$Group = rep(c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'), each = 3)
rownames(forOnlyOutput) = NULL




#########################################################################################################


prop_fullDataset<- fullDataset %>%
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
            grasshopper = ifelse(sum(Group == "grasshopper", na.rm = TRUE) > 0, 1, 0),
            fly = ifelse(sum(Group == "fly", na.rm = TRUE) > 0, 1, 0),
            daddylonglegs = ifelse(sum(Group == "daddylonglegs", na.rm = TRUE) > 0, 1, 0)) %>% 
  group_by(Name, ObservationMethod) %>% 
  summarise(caterpillar_prop = mean(caterpillar),
            spider_prop = mean(spider),
            beetle_prop = mean(beetle),
            truebug_prop = mean(truebug),
            hopper_prop  = mean(hopper),
            ant_prop = mean(ant),
            grasshopper_prop = mean(grasshopper),
            fly_prop = mean(fly),
            daddylonglegs_prop = mean(daddylonglegs),
            Trials = n())  

prop_dataset = left_join(prop_fullDataset, sites, by = 'Name')



prop.cat.Dev = glm(caterpillar_prop ~ dev + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")

prop.spi.Dev = glm(spider_prop ~ dev + ObservationMethod,
                            data = prop_dataset, weights = Trials,  family = "binomial")

prop.bet.Dev= glm(beetle_prop ~ dev + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")

prop.bug.Dev = glm(truebug_prop ~ dev + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")


prop.hop.Dev = glm(hopper_prop ~ dev + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")


prop.ant.Dev = glm(ant_prop ~ dev + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")

prop.fly.Dev = glm(fly_prop ~ dev + ObservationMethod, 
                   data = prop_dataset, weights = Trials,  family = "binomial")

prop.grasshopper.Dev = glm(grasshopper_prop ~ dev + ObservationMethod, 
                   data = prop_dataset, weights = Trials,  family = "binomial")

prop.daddylonglegs.Dev = glm(daddylonglegs_prop ~ dev + ObservationMethod, 
                             data = prop_dataset, weights = Trials,  family = "binomial")






prop.devOnlyOutput = data.frame(rbind(summary(prop.cat.Dev)$coefficients, 
                                 summary(prop.spi.Dev)$coefficients, 
                                 summary(prop.bet.Dev)$coefficients,
                                 summary(prop.bug.Dev)$coefficients, 
                                 summary(prop.hop.Dev)$coefficients,
                                 summary(prop.ant.Dev)$coefficients,
                                 summary(prop.fly.Dev)$coefficients,
                                 summary(prop.grasshopper.Dev)$coefficients,
                                 summary(prop.daddylonglegs.Dev)$coefficients))
prop.devOnlyOutput$term = rep(c('Intercept', 'dev', 'Method'), times = 9)
prop.devOnlyOutput$Group = rep(c('Caterpillar', 'Spider', 'Beetle', 
                                 'Leafhopper', 'Truebugs', 'Ant',
                                 'Fly', 'Grasshopper', 'Daddylonglegs'), each = 3)

rownames(prop.devOnlyOutput) = NULL




# forest cover

prop.cat.For = glm(caterpillar_prop ~ forest + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")


prop.spi.For = glm(spider_prop ~ forest + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")
 

prop.bet.For = glm(beetle_prop ~ forest + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")


prop.bug.For = glm(truebug_prop ~ forest + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")

prop.hop.For = glm(hopper_prop ~ forest + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")


prop.ant.For = glm(ant_prop ~ forest + ObservationMethod, 
                            data = prop_dataset, weights = Trials,  family = "binomial")
 
prop.fly.For = glm(fly_prop ~ forest + ObservationMethod, 
                   data = prop_dataset, weights = Trials,  family = "binomial")

prop.grasshopper.For = glm(grasshopper_prop ~ forest + ObservationMethod, 
                           data = prop_dataset, weights = Trials,  family = "binomial")

prop.daddylonglegs.For = glm(daddylonglegs_prop ~ forest + ObservationMethod, 
                             data = prop_dataset, weights = Trials,  family = "binomial")

prop.forOnlyOutput = data.frame(rbind(summary(prop.cat.For)$coefficients, 
                                 summary(prop.spi.For)$coefficients, 
                                 summary(prop.bet.For)$coefficients,
                                 summary(prop.hop.For)$coefficients, 
                                 summary(prop.bug.For)$coefficients,
                                 summary(prop.ant.For)$coefficients,
                                 summary(prop.grasshopper.For)$coefficients,
                                 summary(prop.daddylonglegs.For)$coefficients,
                                 summary(prop.daddylonglegs.For)$coefficients))
prop.forOnlyOutput$term = rep(c('Intercept', 'forest', 'Method'), times = 9)
prop.forOnlyOutput$Group = rep(c('Caterpillar', 'Spider', 'Beetle', 
                                 'Leafhopper', 'Truebugs', 'Ant',
                                 'Fly', 'Grasshopper', 'Daddylonglegs'), each = 3)

rownames(prop.forOnlyOutput) = NULL




prop_dataset %>% 
  pivot_longer(cols = -c("Name", "ObservationMethod", "Region", 
                         "Longitude", "Latitude", "dev", "forest", "Trials"),
               names_to = "Group",
               values_to = "Occurence") %>% 
ggplot(aes(x = dev, y = Occurence)) +
  geom_smooth(
    method = "glm",
    method.args = list(family = binomial),
    aes(weight = Trials),
    se = TRUE,
    color = "blue"
  )+
  facet_wrap(~Group,  scales = "free_y")+
  theme_bw()







# Question 2

#########################################################################


arthropod_ranks_update
prop.devOnlyOutput

ArthRank.devEffect = left_join(arthropod_ranks_update, prop.devOnlyOutput, by = c("Group"))

ArthRank.devEffect %>% 
  filter(term == "dev") %>% 
  ggplot(aes(x = Herbivore.score, y = Estimate, colour = Herbivore.score)) +
  geom_point() +
  geom_text(aes(label = Group), vjust = -0.5, size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 1)+
  theme_bw()




ArthRank.devEffect %>% 
  filter(term == "dev") %>% 
  ggplot(aes(x = Herbivore.score, y = Estimate, colour = Herbivore.score)) +
  geom_point(size = 3) +
  geom_text(aes(label = Group), vjust = -0.5, size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 1) +
  stat_regline_equation(
    aes(label = paste(..rr.label.., sep = "~~~~")),
    formula = y ~ x,
    label.x = 65,
    label.y = 0.005,
    size = 4,
    color = "black"
  ) +
  labs(
    y = "Est. effect of % Urban Development",
    x = "Herbivore score",
    colour = "Herbivore score"
  ) +
  scale_colour_gradientn(
    colours = c("blue", "yellow", "green")
  ) +
  coord_cartesian(clip = "off") +
  theme_bw()



estDevFor.Abs = left_join(
prop.devOnlyOutput %>% filter(term == "dev") %>% 
  select(Group, Estimate) %>% rename("Est_Dev" = "Estimate") %>% 
  mutate(Est_Dev = abs(Est_Dev)),
prop.forOnlyOutput %>% filter(term == "forest")%>% 
  select(Group, Estimate) %>% rename("Est_For" = "Estimate")%>% 
  mutate(Est_For = abs(Est_For)),
by = c("Group"))



ggplot(estDevFor.Abs, aes(y = Est_Dev, x = Est_For)) +
  geom_point(size = 3, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_text_repel(aes(label = Group), size = 3.5, max.overlaps = Inf) +  
  coord_equal() +
  labs(
    y = "Absolute Estimated effect of Urban Development",
    x = "Absolute Estimated effect of Forest Cover"
  ) +
  theme_bw()



prop.cat.DevLat = glm(caterpillar_prop ~ dev * Latitude  + ObservationMethod, 
                   data = prop_dataset, weights = Trials,  family = "binomial")

prop.spi.DevLat = glm(spider_prop ~ dev * Latitude + ObservationMethod,
                   data = prop_dataset, weights = Trials,  family = "binomial")

prop.bet.DevLat= glm(beetle_prop ~ dev * Latitude + ObservationMethod, 
                  data = prop_dataset, weights = Trials,  family = "binomial")

prop.bug.DevLat = glm(truebug_prop ~ dev * Latitude + ObservationMethod, 
                   data = prop_dataset, weights = Trials,  family = "binomial")


prop.hop.DevLat = glm(hopper_prop ~ dev * Latitude + ObservationMethod, 
                   data = prop_dataset, weights = Trials,  family = "binomial")


prop.ant.DevLat = glm(ant_prop ~ dev * Latitude + ObservationMethod, 
                   data = prop_dataset, weights = Trials,  family = "binomial")

prop.fly.DevLat = glm(fly_prop ~ dev * Latitude + ObservationMethod, 
                   data = prop_dataset, weights = Trials,  family = "binomial")

prop.grasshopper.DevLat = glm(grasshopper_prop ~ dev * Latitude + ObservationMethod, 
                           data = prop_dataset, weights = Trials,  family = "binomial")

prop.daddylonglegs.DevLat = glm(daddylonglegs_prop ~ dev * Latitude + ObservationMethod, 
                             data = prop_dataset, weights = Trials,  family = "binomial")






prop.devLatOutput = data.frame(rbind(summary(prop.cat.DevLat)$coefficients, 
                                      summary(prop.spi.DevLat)$coefficients, 
                                      summary(prop.bet.DevLat)$coefficients,
                                      summary(prop.bug.DevLat)$coefficients, 
                                      summary(prop.hop.DevLat)$coefficients,
                                      summary(prop.ant.DevLat)$coefficients,
                                      summary(prop.fly.DevLat)$coefficients,
                                      summary(prop.grasshopper.DevLat)$coefficients,
                                      summary(prop.daddylonglegs.DevLat)$coefficients))
prop.devLatOutput$term = rep(c('Intercept', 'dev', 'Latitude', 'dev*Latitude', 'Method'),
                             times = 9)
prop.devLatOutput$Group = rep(c('Caterpillar', 'Spider', 'Beetle', 
                                 'Leafhopper', 'Truebugs', 'Ant',
                                 'Fly', 'Grasshopper', 'Daddylonglegs'), each = 5)

rownames(prop.devLatOutput) = NULL


