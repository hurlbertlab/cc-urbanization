library(ggrepel)
library(lme4)
library(pscl)
library(tidyverse)


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

prop_dataset = left_join(prop_fullDataset, sites, by = 'Name') %>% 
  filter(Trials >= 50)



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
  geom_point(aes(colour = ObservationMethod, size = Trials, alpha = 0.2))+
  geom_smooth(
    method = "glm",
    method.args = list(family = binomial),
    aes(weight = Trials, colour = ObservationMethod ),
    se = TRUE
  )+
  facet_wrap(~Group,  scales = "free_y")+
  theme_bw()




# long format data

prop_dataset_LONG = prop_dataset %>% 
  pivot_longer(cols = -c("Name", "ObservationMethod", "Region", 
                         "Longitude", "Latitude", "dev", "forest", "Trials"),
               names_to = "Group",
               values_to = "Occurence")  

# Caterpillar predictions
prop.cat.Dev.pred = predict(prop.cat.Dev, 
                            newdata = prop_dataset_LONG %>% 
                              filter(Group == "caterpillar_prop"), 
                            type = "response") %>% as.data.frame()%>% 
  setNames("Predicted") %>% 
  cbind(prop_dataset_LONG %>% 
          filter(Group == "caterpillar_prop"))

# Spiders
prop.spi.Dev.pred = predict(prop.spi.Dev, 
                            newdata = prop_dataset_LONG %>% 
                              filter(Group == "spider_prop"), 
                            type = "response") %>% as.data.frame()%>% 
  setNames("Predicted") %>% 
  cbind(prop_dataset_LONG %>% 
          filter(Group == "spider_prop"))

# Beetles
prop.bet.Dev.pred = predict(prop.bet.Dev, 
                            newdata = prop_dataset_LONG %>% 
                              filter(Group == "beetle_prop"), 
                            type = "response") %>% as.data.frame()%>% 
  setNames("Predicted") %>% 
  cbind(prop_dataset_LONG %>% 
          filter(Group == "beetle_prop"))

# True bugs
prop.bug.Dev.pred = predict(prop.bug.Dev, 
                            newdata = prop_dataset_LONG %>% 
                              filter(Group == "truebug_prop"), 
                            type = "response") %>% as.data.frame()%>% 
  setNames("Predicted") %>% 
  cbind(prop_dataset_LONG %>% 
          filter(Group == "truebug_prop"))


# Leaf hopper
prop.hop.Dev.pred = predict(prop.hop.Dev, 
                            newdata = prop_dataset_LONG %>% 
                              filter(Group == "hopper_prop"), 
                            type = "response") %>% as.data.frame()%>% 
  setNames("Predicted") %>% 
  cbind(prop_dataset_LONG %>% 
          filter(Group == "hopper_prop"))


# ant
prop.ant.Dev.pred = predict(prop.ant.Dev, 
                            newdata = prop_dataset_LONG %>% 
                              filter(Group == "ant_prop"), 
                            type = "response") %>% as.data.frame()%>% 
  setNames("Predicted") %>% 
  cbind(prop_dataset_LONG %>% 
          filter(Group == "ant_prop"))


# Grasshopper
prop.grasshopper.Dev.pred = predict(prop.grasshopper.Dev, 
                            newdata = prop_dataset_LONG %>% 
                              filter(Group == "grasshopper_prop"), 
                            type = "response") %>% as.data.frame()%>% 
  setNames("Predicted") %>% 
  cbind(prop_dataset_LONG %>% 
          filter(Group == "grasshopper_prop"))


# fly 
prop.fly.Dev.pred = predict(prop.fly.Dev, 
                                    newdata = prop_dataset_LONG %>% 
                                      filter(Group == "fly_prop"), 
                                    type = "response") %>% as.data.frame()%>% 
  setNames("Predicted") %>% 
  cbind(prop_dataset_LONG %>% 
          filter(Group == "fly_prop"))

# daddylonglegs 
prop.daddylonglegs.Dev.pred = predict(prop.daddylonglegs.Dev, 
                            newdata = prop_dataset_LONG %>% 
                              filter(Group == "daddylonglegs_prop"), 
                            type = "response") %>% as.data.frame()%>% 
  setNames("Predicted") %>% 
  cbind(prop_dataset_LONG %>% 
          filter(Group == "daddylonglegs_prop"))

prop_dataset.prediction = rbind(prop.cat.Dev.pred,
                                prop.spi.Dev.pred,
                                prop.bet.Dev.pred,
                                prop.bet.Dev.pred,
                                prop.bug.Dev.pred,
                                prop.hop.Dev.pred,
                                prop.ant.Dev.pred,
                                prop.grasshopper.Dev.pred,
                                prop.fly.Dev.pred,
                                prop.daddylonglegs.Dev.pred) %>% 
                                data.frame()


prop_dataset.prediction %>%
  ggplot() +
  geom_point( aes(Occurence, x = dev) ) +
  geom_smooth(aes(y = Predicted, x = dev))+
  facet_wrap(~Group,  scales = "free_y")+
  theme_bw()+
  labs(subtitle = "Here, slope is made from fitting the predicted occurence (% occ ~ dev + Observ), not the observed occurence")

prop_dataset.prediction %>%
  ggplot() +
  #geom_point( aes(Occurence, x = dev) ) +
  geom_smooth(aes(y = Predicted, x = dev))+
  facet_wrap(~Group,  scales = "free_y")+
  theme_bw()+
  labs(subtitle = "Here, slope is made from fitting the predicted occurence (% occ ~ dev + Observ), not the observed occurence")

prop_dataset.prediction %>%
  ggplot() +
  geom_point( aes(Occurence, x = Latitude) ) +
  geom_smooth(aes(y = Predicted, x = Latitude), color = "steelblue")+
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

###########################################################################################

# - trying interaction effect in multilevel model with interaction in fixed effects

prop_datasetScaled = prop_dataset %>%
  mutate(
    dev_sc = scale(dev),
    forest_sc = scale(forest),
    Latitude_sc = scale(Latitude),
    ObservationMethod = factor(ObservationMethod))




prop.cat.DevLat = glm(caterpillar_prop ~ dev * Latitude  + ObservationMethod, 
                      data = prop_dataset, weights = Trials,  family = "binomial")

catDevPlotM = interact_plot(
  prop.cat.DevLat,
  pred = dev,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Prop. of surveys with caterpillar",
  x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 


prop.bug.DevLat.mixed <- glmer(
  truebug_prop ~ scale(dev) * scale(Latitude) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)

bugDevPlotM = interact_plot(
  prop.bug.DevLat.mixed,
  pred = dev,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Prop. of surveys with truebug",
  x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 



prop.spi.DevLat.mixed <- glmer(
  spider_prop ~ scale(dev) * scale(Latitude) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial
)


spiDevPlotM = interact_plot(
  prop.spi.DevLat.mixed,
  pred = dev,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Prop. of surveys with spider",
  x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 


prop.bet.DevLat= glm(beetle_prop ~ dev * Latitude + ObservationMethod, 
                     data = prop_dataset, weights = Trials,  family = "binomial")

prop.bet.DevLat.mixed <- glmer(
  beetle_prop ~ scale(dev) * scale(Latitude) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial
)


betDevPlotM = interact_plot(
  prop.bet.DevLat.mixed,
  pred = dev,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Prop. of surveys with beetle",
  x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 

 


prop.grasshopper.DevLat.mixed <- glmer(
  grasshopper_prop ~ scale(dev) * scale(Latitude) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial
)

graDevPlotM = interact_plot(
  prop.grasshopper.DevLat.mixed,
  pred = dev,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Prop. of surveys with grasshopper",
  x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 



prop.ant.DevLat.mixed <- glmer(
  ant_prop ~ scale(dev) * scale(Latitude) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)

antDevPlotM = interact_plot(
  prop.ant.DevLat.mixed,
  pred = dev,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Prop. of surveys with ant",
  x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 



prop.fly.DevLat.mixed = glmer(
  fly_prop ~ scale(dev) * scale(Latitude) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial
)

flyDevPlotM = interact_plot(
  prop.fly.DevLat.mixed,
  pred = dev,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Prop. of surveys with fly",
  x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 



prop.hop.DevLat = glm(hopper_prop ~ dev * Latitude + ObservationMethod, 
                              data = prop_dataset, weights = Trials,  family = "binomial")

hopDevPlotM = interact_plot(
  prop.hop.DevLat,
  pred = dev,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Prop. of surveys with hopper",
  x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 



prop.daddylonglegs.DevLat.mixed = glmer(
  daddylonglegs_prop ~ scale(dev) * scale(Latitude) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial
)

daddylonglegsDevPlotM = interact_plot(
  prop.daddylonglegs.DevLat.mixed,
  pred = dev,
  modx = Latitude,
  plot.points = FALSE,
  interval = FALSE,
  y.label = "Prop. of surveys with daddylonglegs",
  x.lab = "% developed cover", cex.lab = 2, vary.lty = FALSE,
  colors = c('darkblue', 'blue', 'powderblue'), line.thickness = 2) 


ggarrange(catDevPlotM, spiDevPlotM, betDevPlotM, bugDevPlotM, hopDevPlotM, antDevPlotM, 
          flyDevPlotM, graDevPlotM, daddylonglegsDevPlotM,
          ncol=3, nrow=3, common.legend = TRUE, legend="bottom")
##########################################################################################



prop.cat.Dev = glm(caterpillar_prop ~ dev + ObservationMethod, 
                      data = prop_dataset, weights = Trials,  family = "binomial")


caterpillar_prop.mixed <- glmer(
  caterpillar_prop ~ scale(dev) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)

summary(caterpillar_prop.mixed)



prop.bug.Dev.mixed <- glmer(
  truebug_prop ~ scale(dev) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)

prop.spi.Dev.mixed <- glmer(
  spider_prop ~ scale(dev) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)


prop.hop.Dev.mixed <- glmer(
  hopper_prop ~ scale(dev) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)

prop.fly.Dev.mixed = glmer(
  fly_prop ~ scale(dev) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)


prop.daddylonglegs.Dev.mixed = glmer(
  daddylonglegs_prop ~ scale(dev) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)


prop.ant.Dev.mixed = glmer(
  ant_prop ~ scale(dev) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)

prop.grasshopper.Dev.mixed <- glmer(
  grasshopper_prop ~ scale(dev) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial
)


prop.bet.Dev.mixed <- glmer(
  beetle_prop ~ scale(dev) + (1 | ObservationMethod ),
  data = prop_dataset,
  weights = Trials,
  family = binomial)


##########################################################################


# ------  compare the r-squred of development cover and forest cover models




# Collect pseudo R² for development models
dev_r2 <- data.frame(
  Group = c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'),
  R2 = c(
    pscl::pR2(cat.Dev)['McFadden'],
    pscl::pR2(spi.Dev)['McFadden'],
    pscl::pR2(beet.Dev)['McFadden'],
    pscl::pR2(hop.Dev)['McFadden'],
    pscl::pR2(bug.Dev)['McFadden'],
    pscl::pR2(ant.Dev)['McFadden']
  ),
  ModelType = 'Development'
)

# Collect pseudo R² for forest models
for_r2 <- data.frame(
  Group = c('caterpillar', 'spider', 'beetle', 'leafhopper', 'truebugs', 'ant'),
  R2 = c(
    pscl::pR2(cat.For)['McFadden'],
    pscl::pR2(spi.For)['McFadden'],
    pscl::pR2(beet.For)['McFadden'],
    pscl::pR2(hop.For)['McFadden'],
    pscl::pR2(bug.For)['McFadden'],
    pscl::pR2(ant.For)['McFadden']
  ),
  ModelType = 'Forest'
)

# Combine both
r2_all <- rbind(dev_r2, for_r2)
r2_all = r2_all %>% 
  pivot_wider(names_from = ModelType, 
              values_from = R2)


ggplot(r2_all, aes(y = Development, x = Forest)) +
  geom_point(size = 3, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_text_repel(aes(label = Group), size = 3.5, max.overlaps = Inf) +  
  coord_equal() +
  labs(
    y = "R2 of Urban Development",
    x = "R2 of Forest Cover"
  ) +
  theme_bw()

