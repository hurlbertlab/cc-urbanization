library(rjags)
library(coda)
library(tidyverse)
library(tidybayes)
library(viridis)
library(raster)
library(ggpubr)

# Load the rjags model fits

caterpillarFit = readRDS("caterpillarFit.rds")
spiderFit      = readRDS("spiderFit.rds")
beetleFit      = readRDS("beetleFit.rds")
truebugFit     = readRDS("truebugFit.rds")
hopperFit      = readRDS("hopperFit.rds")
antFit         = readRDS("antFit.rds")
grasshopperFit = readRDS("grasshopperFit.rds")
flyFit         = readRDS("flyFit.rds")
daddylonglegsFit = readRDS("daddylonglegsFit.rds")


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
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High")))


CaterpillarPlot = ggplot(caterpillarPred_summary, aes(x = dev_c, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "Scaled Urban Development", y = "Predicted Probability of Caterpillar Presence") +
  annotation_raster(catImage, ymin = .085, ymax = .1, xmin = 0, xmax = 1.5) +
  theme_minimal()


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
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High")))


spiderPlot = ggplot(spiderPred_summary, aes(x = dev_c, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "Scaled Urban Development", y = "Predicted Probability of spider Presence") +
  annotation_raster(spiderImage, ymin = .29, ymax = .35, xmin = 0, xmax = 1.5) +
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
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High")))


beetlePlot = ggplot(beetlePred_summary, aes(x = dev_c, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "Scaled Urban Development", y = "Predicted Probability of beetle Presence") +
  annotation_raster(beetleImage, ymin = .29, ymax = .31, xmin = 0, xmax = 1.5)  +
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
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High")))


truebugPlot = ggplot(truebugPred_summary, aes(x = dev_c, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "Scaled Urban Development", y = "Predicted Probability of truebug Presence") +
  annotation_raster(truebugImage, ymin = .1, ymax = .12, xmin = -0.8, xmax = 0.6)  +
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
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High")))


hopperPlot = ggplot(hopperPred_summary, aes(x = dev_c, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "Scaled Urban Development", y = "Predicted Probability of True hopper Presence") +
  annotation_raster(hopperImage, ymin = .07, ymax = .095, xmin = 0, xmax = 1.3)  +
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
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High")))


antPlot = ggplot(antPred_summary, aes(x = dev_c, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "Scaled Urban Development", y = "Predicted Probability of ant Presence") +
  annotation_raster(antImage, ymin = .14, ymax = .16, xmin = 0.2, xmax = 1.3)  +
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
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High")))


grasshopperPlot = ggplot(grasshopperPred_summary, aes(x = dev_c, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "Scaled Urban Development", y = "Predicted Probability of grasshopper Presence") +
  annotation_raster(grasshopperImage, ymin = .12, ymax = .16, xmin = -0.8, xmax = 0.8)  +
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
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High")))


flyPlot = ggplot(flyPred_summary, aes(x = dev_c, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "Scaled Urban Development", y = "Predicted Probability of fly Presence") +
  annotation_raster(flyImage, ymin = .14, ymax = .17, xmin = -0.8, xmax = 0.6)  +
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
  mutate(Latitude = factor(Latitude, levels = c("Low", "Mid", "High")))


daddylonglegsPlot = ggplot(daddylonglegsPred_summary, aes(x = dev_c, y = p_median, color = factor(Latitude))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper, fill = Latitude), alpha = 0.2, color = NA) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(x = "Scaled Urban Development", y = "Predicted Probability of daddylonglegs Presence") +
  annotation_raster(daddylonglegImage, ymin = .04, ymax = .07, xmin = -0.9, xmax = 0.7)  +
  theme_minimal()



ggarrange(CaterpillarPlot, spiderPlot, beetlePlot, truebugPlot, hopperPlot, antPlot,
          flyPlot, grasshopperPlot, daddylonglegsPlot,
          ncol=3, nrow=3, common.legend = TRUE, legend="bottom")








##################################################################################

mean_dev = mean(dataset$dev)
sd_dev  = sd(dataset$dev)

mean_lat = mean(dataset$Latitude)
sd_lat = sd(dataset$Latitude)
post = as.matrix(caterpillarFit)
beta_dev_orig = post[, "beta_dev"] / sd_dev
beta_lat_orig = post[, "beta_lat"] / sd_lat
beta_int_orig = post[, "beta_int"] / (sd_dev * sd_lat)

beta_0_orig =
  post[, "beta_0"] -
  post[, "beta_dev"] * mean_dev / sd_dev -
  post[, "beta_lat"] * mean_lat / sd_lat +
  post[, "beta_int"] * mean_dev * mean_lat / (sd_dev * sd_lat)

beta_method_orig = post[, "beta_method"]
median(beta_dev_orig)
quantile(beta_dev_orig, c(0.025, 0.975))


OR_dev <- exp(beta_dev_orig)
median(OR_dev)
quantile(OR_dev, c(0.025, 0.975))











































