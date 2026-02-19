
# load 'dataset' object first

# Load libraries ----
library(tidyverse)
library(tidybayes)
library(viridis)
library(raster)
library(ggpubr)
library(coda)

# Load the rjags model fits saved as rds files ----

caterpillarFit = readRDS("caterpillarFit.rds")
spiderFit      = readRDS("spiderFit.rds")
beetleFit      = readRDS("beetleFit.rds")
truebugFit     = readRDS("truebugFit.rds")
hopperFit      = readRDS("hopperFit.rds")
antFit         = readRDS("antFit.rds")
grasshopperFit = readRDS("grasshopperFit.rds")
flyFit         = readRDS("flyFit.rds")
daddylonglegsFit = readRDS("daddylonglegsFit.rds")

# Arthropod silhouette images ----



### Arthropod images ----
catImage = readPNG('images/caterpillar.png')
antImage = readPNG('images/ant.png')
beetleImage = readPNG('images/beetle.png')
spiderImage = readPNG('images/spider.png')
hopperImage = readPNG('images/leafhopper.png')
truebugImage = readPNG('images/truebugs.png')

# https://www.phylopic.org/images/7174527d-3060-4c78-ab59-c3ccf76075de/opilio-saxatilis
daddylonglegImage = image_read('images/daddylongleg.png') %>% 
  image_convert(type = "truecolor") %>% 
  as.raster()

# https://www.phylopic.org/images/2bec28d1-3d13-4512-8b4d-00cd0fcef6e4/clogmia-albipunctata
flyImage = image_read('images/fly.png') %>% 
  image_convert(type = "truecolor") %>% 
  as.raster()

# https://www.phylopic.org/images/0de35750-d9ba-472a-b4f2-3c630777fcc3/acrididae
grasshopperImage = image_read('images/grasshopper.png') %>% 
  image_convert(type = "truecolor") %>% 
  as.raster()


# simulate hypothesis plot------------------------------------------------------------------
set.seed(123)
simulate_lines <- function(beta_dev, beta_int, label) {
  
  n <- 500
  
  dev <- rnorm(n, 0, 0.4)
  Latitude <- rnorm(n, 0, 0.5)
  ObservationMethod <- sample(c("Visual","Beating sheet"), n, replace = TRUE)
  
  method_num <- ifelse(ObservationMethod == "Visual", 0, 1)
  
  beta0  <- 0.5
  beta_lat <- 0.5
  beta_method <- 0.05
  
  y <- beta0 +
    beta_dev*dev +
    beta_lat*Latitude +
    beta_int*dev*Latitude +
    beta_method*method_num +
    rnorm(n, 0, 0.3)
  
  sim_dataset <- data.frame(y, dev, Latitude, ObservationMethod)
  
  mod <- lm(y ~ dev * Latitude + ObservationMethod, data = sim_dataset)
  
  dev_seq <- seq(min(dev), max(dev), length = 100)
  
  lat_vals <- c(-0.7, 0, 0.7)
  lat_labels <- c("Low","Mid","High")
  
  # predict response, hold latitude constant, fix method to visual
  preds <- lapply(1:3, function(i){
    predict(mod,
            newdata = data.frame(dev = dev_seq,
                                 Latitude = lat_vals[i],
                                 ObservationMethod = "Visual"))
  })
  
  data.frame(
    dev = rep(dev_seq, 3),
    response = unlist(preds),
    Latitude = factor(rep(lat_labels, each = length(dev_seq))),
    Scenario = label
  )
}


plot_data <- rbind(
  simulate_lines( 1.3, 1.25, "Positive effect, UHI"),
  simulate_lines( 1.3, 0,    "Positive effect, no UHI"),
  simulate_lines(-1.3, 1.25, "Negative effect, UHI"),
  simulate_lines(-1.3, 0,    "Negative effect, no UHI")
) %>% mutate(Latitude = factor(Latitude, levels= c("High", "Mid", "Low")))


hypothesisPlot = ggplot(plot_data, aes(dev, response, color = Latitude)) +
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  facet_wrap(~Scenario, ncol = 2) +
  theme_bw(base_size = 13) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +labs(y = "Proportion of surey with arthropod",
          x = "% Urban development")

hypothesisPlot

panel_labels <- data.frame(
  Scenario = c("Positive effect, UHI",
               "Positive effect, no UHI",
               "Negative effect, UHI",
               "Negative effect, no UHI"),
  label = c("d", "c", "b", "a"),
  x = -Inf,
  y = Inf
)

hypothesisPlot +
  geom_text(data = panel_labels,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = -0.5,
            vjust = 1.5,
            size = 5,
            fontface = "bold")

################################################################################

dataset$dev_c = as.numeric(scale(dataset$dev))
dataset$lat_c  = as.numeric(scale(dataset$Latitude))
dataset$dev_lat = dataset$dev_c * dataset$lat_c
dataset$method = as.numeric(dataset$ObservationMethod == "Visual")


# Three latitude levels: -1 SD, mean (0), +1 SD
lat_levels = c(-1, 0, 1)

pred_grid = expand.grid(
  dev_c = seq(min(dataset$dev_c), max(dataset$dev_c), length.out = 50),
  lat_c = lat_levels,
  method = 0  # fix ObservationMethod, as Visual. Change to 1 if you prefer beating sheet.
)

################################################################################
# Caterpillar
################################################################################


# Extract posterior draws
caterpillarPost_draws = as.mcmc(do.call(rbind, caterpillarFit)) %>%
  as.data.frame()

# For each posterior sample, calculate predicted probability
caterpillarPred_grid_long <- pred_grid %>%
  crossing(post_draw = 1:nrow(caterpillarPost_draws)) %>%
  mutate(
    beta_0      = caterpillarPost_draws$beta_0[post_draw],
    beta_dev    = caterpillarPost_draws$beta_dev[post_draw],
    beta_lat    = caterpillarPost_draws$beta_lat[post_draw],
    beta_int    = caterpillarPost_draws$beta_int[post_draw],
    beta_method = caterpillarPost_draws$beta_method[post_draw],
    logit_p = beta_0 + beta_dev*dev_c + beta_lat*lat_c + beta_int*(dev_c*lat_c) + beta_method*method,
    p       = plogis(logit_p)
  )

caterpillarPred_summary = caterpillarPred_grid_long %>%
  group_by(dev_c, lat_c) %>%
  summarise(
    p_median = median(p),
    p_lower  = quantile(p, 0.025),
    p_upper  = quantile(p, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(Latitude = case_when(
    lat_c == -1 ~ "Low",
    lat_c == 0  ~ "Mid",
    lat_c == 1  ~ "High"
  )) %>%
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High"))) %>% 
  mutate(dev = dev_c * sd(dataset$dev) + mean(dataset$dev)) # converting back to raw values.





CaterpillarPlot = ggplot(caterpillarPred_summary, aes(x = dev, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "% Urban Development", y = "Prop. of surveys with spider presence") +
  annotation_raster(catImage, ymin = .085, ymax = .1, xmin = 60, xmax = 100) +
  theme_minimal()

CaterpillarPlot

################################################################################
# spider
################################################################################


# Extract posterior draws
spiderPost_draws = as.mcmc(do.call(rbind, spiderFit)) %>%
  as.data.frame()

# For each posterior sample, calculate predicted probability
spiderPred_grid_long <- pred_grid %>%
  crossing(post_draw = 1:nrow(spiderPost_draws)) %>%
  mutate(
    beta_0      = spiderPost_draws$beta_0[post_draw],
    beta_dev    = spiderPost_draws$beta_dev[post_draw],
    beta_lat    = spiderPost_draws$beta_lat[post_draw],
    beta_int    = spiderPost_draws$beta_int[post_draw],
    beta_method = spiderPost_draws$beta_method[post_draw],
    logit_p = beta_0 + beta_dev*dev_c + beta_lat*lat_c + beta_int*(dev_c*lat_c) + beta_method*method,
    p       = plogis(logit_p)
  )

spiderPred_summary = spiderPred_grid_long %>%
  group_by(dev_c, lat_c) %>%
  summarise(
    p_median = median(p),
    p_lower  = quantile(p, 0.025),
    p_upper  = quantile(p, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(Latitude = case_when(
    lat_c == -1 ~ "Low",
    lat_c == 0  ~ "Mid",
    lat_c == 1  ~ "High"
  )) %>%
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High"))) %>% 
  mutate(dev = dev_c * sd(dataset$dev) + mean(dataset$dev))


spiderPlot = ggplot(spiderPred_summary, aes(x = dev, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "% Urban Development", y = "prop. of surveys with spider presence") +
  annotation_raster(spiderImage, ymin = .29, ymax = .35, xmin = 60, xmax = 100) +
  theme_minimal()


################################################################################
# beetle
################################################################################


# Extract posterior draws
beetlePost_draws = as.mcmc(do.call(rbind, beetleFit)) %>%
  as.data.frame()

# For each posterior sample, calculate predicted probability
beetlePred_grid_long <- pred_grid %>%
  crossing(post_draw = 1:nrow(beetlePost_draws)) %>%
  mutate(
    beta_0      = beetlePost_draws$beta_0[post_draw],
    beta_dev    = beetlePost_draws$beta_dev[post_draw],
    beta_lat    = beetlePost_draws$beta_lat[post_draw],
    beta_int    = beetlePost_draws$beta_int[post_draw],
    beta_method = beetlePost_draws$beta_method[post_draw],
    logit_p = beta_0 + beta_dev*dev_c + beta_lat*lat_c + beta_int*(dev_c*lat_c) + beta_method*method,
    p       = plogis(logit_p)
  )

beetlePred_summary = beetlePred_grid_long %>%
  group_by(dev_c, lat_c) %>%
  summarise(
    p_median = median(p),
    p_lower  = quantile(p, 0.025),
    p_upper  = quantile(p, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(Latitude = case_when(
    lat_c == -1 ~ "Low",
    lat_c == 0  ~ "Mid",
    lat_c == 1  ~ "High"
  )) %>%
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High"))) %>% 
  mutate(dev = dev_c * sd(dataset$dev) + mean(dataset$dev))


beetlePlot = ggplot(beetlePred_summary, aes(x = dev, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "% Urban Development", y = "Prop. of surveys with beetle presence") +
  annotation_raster(beetleImage, ymin = .28, ymax = .30, xmin = 60, xmax = 100)  +
  theme_minimal()




################################################################################
# truebug
################################################################################


# Extract posterior draws
truebugPost_draws = as.mcmc(do.call(rbind, truebugFit)) %>%
  as.data.frame()

# For each posterior sample, calculate predicted probability
truebugPred_grid_long <- pred_grid %>%
  crossing(post_draw = 1:nrow(truebugPost_draws)) %>%
  mutate(
    beta_0      = truebugPost_draws$beta_0[post_draw],
    beta_dev    = truebugPost_draws$beta_dev[post_draw],
    beta_lat    = truebugPost_draws$beta_lat[post_draw],
    beta_int    = truebugPost_draws$beta_int[post_draw],
    beta_method = truebugPost_draws$beta_method[post_draw],
    logit_p = beta_0 + beta_dev*dev_c + beta_lat*lat_c + beta_int*(dev_c*lat_c) + beta_method*method,
    p       = plogis(logit_p)
  )

truebugPred_summary = truebugPred_grid_long %>%
  group_by(dev_c, lat_c) %>%
  summarise(
    p_median = median(p),
    p_lower  = quantile(p, 0.025),
    p_upper  = quantile(p, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(Latitude = case_when(
    lat_c == -1 ~ "Low",
    lat_c == 0  ~ "Mid",
    lat_c == 1  ~ "High"
  )) %>%
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High"))) %>% 
  mutate(dev = dev_c * sd(dataset$dev) + mean(dataset$dev))


truebugPlot = ggplot(truebugPred_summary, aes(x = dev, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "% Urban Development", y = "Prop. of surveys with truebug presence") +
  annotation_raster(truebugImage, ymin = .09, ymax = .11, xmin = 0, xmax = 40)  +
  theme_minimal()


################################################################################
# hopper
################################################################################

# Extract posterior draws
hopperPost_draws = as.mcmc(do.call(rbind, hopperFit)) %>%
  as.data.frame()

# For each posterior sample, calculate predicted probability
hopperPred_grid_long <- pred_grid %>%
  crossing(post_draw = 1:nrow(hopperPost_draws)) %>%
  mutate(
    beta_0      = hopperPost_draws$beta_0[post_draw],
    beta_dev    = hopperPost_draws$beta_dev[post_draw],
    beta_lat    = hopperPost_draws$beta_lat[post_draw],
    beta_int    = hopperPost_draws$beta_int[post_draw],
    beta_method = hopperPost_draws$beta_method[post_draw],
    logit_p = beta_0 + beta_dev*dev_c + beta_lat*lat_c + beta_int*(dev_c*lat_c) + beta_method*method,
    p       = plogis(logit_p)
  )

hopperPred_summary = hopperPred_grid_long %>%
  group_by(dev_c, lat_c) %>%
  summarise(
    p_median = median(p),
    p_lower  = quantile(p, 0.025),
    p_upper  = quantile(p, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(Latitude = case_when(
    lat_c == -1 ~ "Low",
    lat_c == 0  ~ "Mid",
    lat_c == 1  ~ "High"
  )) %>%
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High"))) %>% 
  mutate(dev = dev_c * sd(dataset$dev) + mean(dataset$dev))


hopperPlot = ggplot(hopperPred_summary, aes(x = dev, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "% Urban Development", y = "Prop. of surveys with truehopper presence") +
  annotation_raster(hopperImage, ymin = .07, ymax = .095, xmin = 60, xmax = 100)  +
  theme_minimal()


################################################################################
# ant
################################################################################


# Extract posterior draws
antPost_draws = as.mcmc(do.call(rbind, antFit)) %>%
  as.data.frame()

# For each posterior sample, calculate predicted probability
antPred_grid_long <- pred_grid %>%
  crossing(post_draw = 1:nrow(antPost_draws)) %>%
  mutate(
    beta_0      = antPost_draws$beta_0[post_draw],
    beta_dev    = antPost_draws$beta_dev[post_draw],
    beta_lat    = antPost_draws$beta_lat[post_draw],
    beta_int    = antPost_draws$beta_int[post_draw],
    beta_method = antPost_draws$beta_method[post_draw],
    logit_p = beta_0 + beta_dev*dev_c + beta_lat*lat_c + beta_int*(dev_c*lat_c) + beta_method*method,
    p       = plogis(logit_p)
  )

antPred_summary = antPred_grid_long %>%
  group_by(dev_c, lat_c) %>%
  summarise(
    p_median = median(p),
    p_lower  = quantile(p, 0.025),
    p_upper  = quantile(p, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(Latitude = case_when(
    lat_c == -1 ~ "Low",
    lat_c == 0  ~ "Mid",
    lat_c == 1  ~ "High"
  )) %>%
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High"))) %>% 
  mutate(dev = dev_c * sd(dataset$dev) + mean(dataset$dev))


antPlot = ggplot(antPred_summary, aes(x = dev, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "% Urban Development", y = "Prop. of surveys with ant presence") +
  annotation_raster(antImage, ymin = .14, ymax = .16, xmin = 60, xmax = 99)  +
  theme_minimal()



################################################################################
# grasshopper
################################################################################


# Extract posterior draws
grasshopperPost_draws = as.mcmc(do.call(rbind, grasshopperFit)) %>%
  as.data.frame()

# For each posterior sample, calculate predicted probability
grasshopperPred_grid_long <- pred_grid %>%
  crossing(post_draw = 1:nrow(grasshopperPost_draws)) %>%
  mutate(
    beta_0      = grasshopperPost_draws$beta_0[post_draw],
    beta_dev    = grasshopperPost_draws$beta_dev[post_draw],
    beta_lat    = grasshopperPost_draws$beta_lat[post_draw],
    beta_int    = grasshopperPost_draws$beta_int[post_draw],
    beta_method = grasshopperPost_draws$beta_method[post_draw],
    logit_p = beta_0 + beta_dev*dev_c + beta_lat*lat_c + beta_int*(dev_c*lat_c) + beta_method*method,
    p       = plogis(logit_p)
  )

grasshopperPred_summary = grasshopperPred_grid_long %>%
  group_by(dev_c, lat_c) %>%
  summarise(
    p_median = median(p),
    p_lower  = quantile(p, 0.025),
    p_upper  = quantile(p, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(Latitude = case_when(
    lat_c == -1 ~ "Low",
    lat_c == 0  ~ "Mid",
    lat_c == 1  ~ "High"
  )) %>%
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High"))) %>% 
  mutate(dev = dev_c * sd(dataset$dev) + mean(dataset$dev))


grasshopperPlot = ggplot(grasshopperPred_summary, aes(x = dev, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "% Urban Development", y = "Prop. of surveys with grasshopper presence") +
  annotation_raster(grasshopperImage, ymin = .12, ymax = .16, xmin = 15, xmax = 60)  +
  theme_minimal()



################################################################################
# fly
################################################################################

# Extract posterior draws
flyPost_draws = as.mcmc(do.call(rbind, flyFit)) %>%
  as.data.frame()

# For each posterior sample, calculate predicted probability
flyPred_grid_long <- pred_grid %>%
  crossing(post_draw = 1:nrow(flyPost_draws)) %>%
  mutate(
    beta_0      = flyPost_draws$beta_0[post_draw],
    beta_dev    = flyPost_draws$beta_dev[post_draw],
    beta_lat    = flyPost_draws$beta_lat[post_draw],
    beta_int    = flyPost_draws$beta_int[post_draw],
    beta_method = flyPost_draws$beta_method[post_draw],
    logit_p = beta_0 + beta_dev*dev_c + beta_lat*lat_c + beta_int*(dev_c*lat_c) + beta_method*method,
    p       = plogis(logit_p)
  )

flyPred_summary = flyPred_grid_long %>%
  group_by(dev_c, lat_c) %>%
  summarise(
    p_median = median(p),
    p_lower  = quantile(p, 0.025),
    p_upper  = quantile(p, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(Latitude = case_when(
    lat_c == -1 ~ "Low",
    lat_c == 0  ~ "Mid",
    lat_c == 1  ~ "High"
  )) %>%
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High"))) %>% 
  mutate(dev = dev_c * sd(dataset$dev) + mean(dataset$dev))


flyPlot = ggplot(flyPred_summary, aes(x = dev, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "% Urban Development", y = "Prop. of surveys with fly presence") +
  annotation_raster(flyImage, ymin = .14, ymax = .16, xmin = 20, xmax = 60)  +
  theme_minimal()


################################################################################
# daddylonglegs
################################################################################


# Extract posterior draws
daddylonglegsPost_draws = as.mcmc(do.call(rbind, daddylonglegsFit)) %>%
  as.data.frame()

# For each posterior sample, calculate predicted probability
daddylonglegsPred_grid_long <- pred_grid %>%
  crossing(post_draw = 1:nrow(daddylonglegsPost_draws)) %>%
  mutate(
    beta_0      = daddylonglegsPost_draws$beta_0[post_draw],
    beta_dev    = daddylonglegsPost_draws$beta_dev[post_draw],
    beta_lat    = daddylonglegsPost_draws$beta_lat[post_draw],
    beta_int    = daddylonglegsPost_draws$beta_int[post_draw],
    beta_method = daddylonglegsPost_draws$beta_method[post_draw],
    logit_p = beta_0 + beta_dev*dev_c + beta_lat*lat_c + beta_int*(dev_c*lat_c) + beta_method*method,
    p       = plogis(logit_p)
  )

daddylonglegsPred_summary = daddylonglegsPred_grid_long %>%
  group_by(dev_c, lat_c) %>%
  summarise(
    p_median = median(p),
    p_lower  = quantile(p, 0.025),
    p_upper  = quantile(p, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(Latitude = case_when(
    lat_c == -1 ~ "Low",
    lat_c == 0  ~ "Mid",
    lat_c == 1  ~ "High"
  )) %>%
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High"))) %>% 
  mutate(dev = dev_c * sd(dataset$dev) + mean(dataset$dev))


daddylonglegsPlot = ggplot(daddylonglegsPred_summary, aes(x = dev, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "% Urban Development", y = "Prop. of surveys with harvestman presence") +
  annotation_raster(daddylonglegImage, ymin = .04, ymax = .06, xmin = 60, xmax = 90)  +
  theme_minimal()



ggarrange(CaterpillarPlot +
            annotate("text", x = I(0.05), y = I(0.95), label = "a", size = 8,
                     fontface = "bold"), 
          spiderPlot +
            annotate("text", x = I(0.05), y = I(0.95), label = "b", size = 8,
                     fontface = "bold"), 
          beetlePlot +
            annotate("text", x = I(0.05), y = I(0.95), label = "c", size = 8,
                     fontface = "bold"), 
          truebugPlot +
            annotate("text", x = I(0.05), y = I(0.95), label = "d", size = 8,
                     fontface = "bold"), 
          hopperPlot +
            annotate("text", x = I(0.05), y = I(0.95), label = "e", size = 8,
                     fontface = "bold"), 
          antPlot +
            annotate("text", x = I(0.05), y = I(0.95), label = "f", size = 8,
                     fontface = "bold"),
           grasshopperPlot +
            annotate("text", x = I(0.05), y = I(0.95), label = "g", size = 8,
                     fontface = "bold"), 
          ncol= 2, nrow= 4, common.legend = TRUE, legend = "right")








##################################################################################

 









































