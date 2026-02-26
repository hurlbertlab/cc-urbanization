
# load 'dataset' object first

# Load libraries ----
library(tidyverse)
library(tidybayes)
library(viridis)
library(raster)
library(ggpubr)
library(coda)
library(gt) # to create tables
library(posterior)
library(ggdist)

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


get_mcmc_diagnostics = function(fit, model_name) {
  
  rhat <- gelman.diag(fit, autoburnin = FALSE)$psrf[,1]
  
  ess <- effectiveSize(fit)

  
  
  summary_stats <- summary(fit)$statistics
  quantiles <- summary(fit)$quantiles
  
  tibble(
    model = model_name,
    Parameter = names(rhat),
    Mean = summary_stats[, "Mean"],
    SE = summary_stats[, "SD"],      # posterior SD
    CI_lower = quantiles[, "2.5%"],
    CI_upper = quantiles[, "97.5%"],
    Rhat = as.numeric(rhat),
    ESS = as.numeric(ess)
  )
}

fits = list(
  Caterpillar = caterpillarFit,
  Spider = spiderFit,
  Beetle = beetleFit,
  TrueBug = truebugFit,
  Hopper = hopperFit,
  Ant = antFit,
  Grasshopper = grasshopperFit,
  Fly = flyFit,
  DaddyLonglegs = daddylonglegsFit
)

diagnostics_table = bind_rows(
  lapply(names(fits), function(name) {
    get_mcmc_diagnostics(fits[[name]], name)
  })
) %>% as.data.frame()

diagnostics_table






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

caterpillarSlopes = caterpillarPost_draws %>%
  as.data.frame() %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    mean_low  = mean(slope_low),
    se_low    = sd(slope_low),
    ci_low_l  = quantile(slope_low, 0.025),
    ci_low_u  = quantile(slope_low, 0.975),
    
    mean_mid  = mean(slope_mid),
    se_mid    = sd(slope_mid),
    ci_mid_l  = quantile(slope_mid, 0.025),
    ci_mid_u  = quantile(slope_mid, 0.975),
    
    mean_high = mean(slope_high),
    se_high   = sd(slope_high),
    ci_high_l = quantile(slope_high, 0.025),
    ci_high_u = quantile(slope_high, 0.975),
    
    # Bayesian tail-area p-values
    p_low  = 2 * min(mean(slope_low  > 0), mean(slope_low  < 0)),
    p_mid  = 2 * min(mean(slope_mid  > 0), mean(slope_mid  < 0)),
    p_high = 2 * min(mean(slope_high > 0), mean(slope_high < 0))
  )


caterpillarSlopes_long = tibble(
  Latitude = factor(c("low", "mid", "high"),
                    levels = c("low","mid","high")),
  Estimate = c(caterpillarSlopes$mean_low,
               caterpillarSlopes$mean_mid,
               caterpillarSlopes$mean_high),
  SE       = c(caterpillarSlopes$se_low,
               caterpillarSlopes$se_mid,
               caterpillarSlopes$se_high),
  CI_lower = c(caterpillarSlopes$ci_low_l,
               caterpillarSlopes$ci_mid_l,
               caterpillarSlopes$ci_high_l),
  CI_upper = c(caterpillarSlopes$ci_low_u,
               caterpillarSlopes$ci_mid_u,
               caterpillarSlopes$ci_high_u),
  p_value  = c(caterpillarSlopes$p_low,
               caterpillarSlopes$p_mid,
               caterpillarSlopes$p_high)
)


caterpillarSlopes_long

caterpillar_slope_draws = caterpillarPost_draws %>%
  as.data.frame() %>%
  mutate(
    Low  = beta_dev + beta_int * -1,
    Mid  = beta_dev + beta_int * 0,
    High = beta_dev + beta_int * 1
  ) %>%
  dplyr::select(Low, Mid, High) %>%   
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Latitude",
    values_to = "Slope"
  ) %>% 
  mutate(Latitude = factor(Latitude, 
                                  levels= c("Low","Mid","High")))


ggplot(caterpillar_slope_draws, aes(x = Slope, y = Latitude, fill = Latitude)) +
  stat_halfeye(.width = 0.95) +
  scale_color_viridis_d(name = "Latitude") +
  scale_fill_viridis_d(name = "Latitude") +
  labs(
    x = "Urbanization effect",
    y = "",
    title = "Caterpillars: Posterior slopes by latitude"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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

spiderSlopes = spiderPost_draws %>%
  as.data.frame() %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    mean_low  = mean(slope_low),
    se_low    = sd(slope_low),
    ci_low_l  = quantile(slope_low, 0.025),
    ci_low_u  = quantile(slope_low, 0.975),
    
    mean_mid  = mean(slope_mid),
    se_mid    = sd(slope_mid),
    ci_mid_l  = quantile(slope_mid, 0.025),
    ci_mid_u  = quantile(slope_mid, 0.975),
    
    mean_high = mean(slope_high),
    se_high   = sd(slope_high),
    ci_high_l = quantile(slope_high, 0.025),
    ci_high_u = quantile(slope_high, 0.975),
    
    # Bayesian p-values (two-tailed)
    p_low  = 2 * min(mean(slope_low  > 0), mean(slope_low  < 0)),
    p_mid  = 2 * min(mean(slope_mid  > 0), mean(slope_mid  < 0)),
    p_high = 2 * min(mean(slope_high > 0), mean(slope_high < 0)))

spiderSlopes_long = tibble(
  Latitude = factor(c("low", "mid", "high"),
                    levels = c("low","mid","high")),
  Estimate = c(spiderSlopes$mean_low,
               spiderSlopes$mean_mid,
               spiderSlopes$mean_high),
  SE       = c(spiderSlopes$se_low,
               spiderSlopes$se_mid,
               spiderSlopes$se_high),
  CI_lower = c(spiderSlopes$ci_low_l,
               spiderSlopes$ci_mid_l,
               spiderSlopes$ci_high_l),
  CI_upper = c(spiderSlopes$ci_low_u,
               spiderSlopes$ci_mid_u,
               spiderSlopes$ci_high_u),
  p_value  = c(spiderSlopes$p_low,
               spiderSlopes$p_mid,
               spiderSlopes$p_high))

spiderSlopes_long


spider_slope_draws = spiderPost_draws %>%
  as.data.frame() %>%
  mutate(
    Low  = beta_dev + beta_int * -1,
    Mid  = beta_dev + beta_int * 0,
    High = beta_dev + beta_int * 1
  ) %>%
  dplyr::select(Low, Mid, High) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Latitude",
    values_to = "Slope"
  ) %>%
  mutate(Latitude = factor(Latitude,
                           levels = c("Low","Mid","High")))


ggplot(spider_slope_draws,
       aes(x = Slope, y = Latitude, fill = Latitude)) +
  stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(name = "Latitude") +
  labs(
    x = "Urbanization effect",
    y = "",
    title = "Spiders: Posterior slopes by latitude"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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



beetlePost_draws %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    P_low  = mean(slope_low  > 0),
    P_mid  = mean(slope_mid  > 0),
    P_high = mean(slope_high > 0)
  )




beetleSlopes = beetlePost_draws %>%
  as.data.frame() %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    mean_low  = mean(slope_low),
    se_low    = sd(slope_low),
    ci_low_l  = quantile(slope_low, 0.025),
    ci_low_u  = quantile(slope_low, 0.975),
    
    mean_mid  = mean(slope_mid),
    se_mid    = sd(slope_mid),
    ci_mid_l  = quantile(slope_mid, 0.025),
    ci_mid_u  = quantile(slope_mid, 0.975),
    
    mean_high = mean(slope_high),
    se_high   = sd(slope_high),
    ci_high_l = quantile(slope_high, 0.025),
    ci_high_u = quantile(slope_high, 0.975),
    
    # Bayesian p-values (two-tailed)
    p_low  = 2 * min(mean(slope_low  > 0), mean(slope_low  < 0)),
    p_mid  = 2 * min(mean(slope_mid  > 0), mean(slope_mid  < 0)),
    p_high = 2 * min(mean(slope_high > 0), mean(slope_high < 0)))

beetleSlopes_long = tibble(
  Latitude = factor(c("low", "mid", "high"),
                    levels = c("low","mid","high")),
  Estimate = c(beetleSlopes$mean_low,
               beetleSlopes$mean_mid,
               beetleSlopes$mean_high),
  SE       = c(beetleSlopes$se_low,
               beetleSlopes$se_mid,
               beetleSlopes$se_high),
  CI_lower = c(beetleSlopes$ci_low_l,
               beetleSlopes$ci_mid_l,
               beetleSlopes$ci_high_l),
  CI_upper = c(beetleSlopes$ci_low_u,
               beetleSlopes$ci_mid_u,
               beetleSlopes$ci_high_u),
  p_value  = c(beetleSlopes$p_low,
               beetleSlopes$p_mid,
               beetleSlopes$p_high))

beetleSlopes_long


beetle_slope_draws = beetlePost_draws %>%
  as.data.frame() %>%
  mutate(
    Low  = beta_dev + beta_int * -1,
    Mid  = beta_dev + beta_int * 0,
    High = beta_dev + beta_int * 1
  ) %>%
  dplyr::select(Low, Mid, High) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Latitude",
    values_to = "Slope"
  ) %>%
  mutate(Latitude = factor(Latitude,
                           levels = c("Low","Mid","High")))


ggplot(beetle_slope_draws,
       aes(x = Slope, y = Latitude, fill = Latitude)) +
  stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(name = "Latitude") +
  labs(
    x = "Urbanization effect",
    y = "",
    title = "beetles: Posterior slopes by latitude"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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


truebugSlopes = truebugPost_draws %>%
  as.data.frame() %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    mean_low  = mean(slope_low),
    se_low    = sd(slope_low),
    ci_low_l  = quantile(slope_low, 0.025),
    ci_low_u  = quantile(slope_low, 0.975),
    
    mean_mid  = mean(slope_mid),
    se_mid    = sd(slope_mid),
    ci_mid_l  = quantile(slope_mid, 0.025),
    ci_mid_u  = quantile(slope_mid, 0.975),
    
    mean_high = mean(slope_high),
    se_high   = sd(slope_high),
    ci_high_l = quantile(slope_high, 0.025),
    ci_high_u = quantile(slope_high, 0.975),
    
    # Bayesian p-values (two-tailed)
    p_low  = 2 * min(mean(slope_low  > 0), mean(slope_low  < 0)),
    p_mid  = 2 * min(mean(slope_mid  > 0), mean(slope_mid  < 0)),
    p_high = 2 * min(mean(slope_high > 0), mean(slope_high < 0)))

truebugSlopes_long = tibble(
  Latitude = factor(c("low", "mid", "high"),
                    levels = c("low","mid","high")),
  Estimate = c(truebugSlopes$mean_low,
               truebugSlopes$mean_mid,
               truebugSlopes$mean_high),
  SE       = c(truebugSlopes$se_low,
               truebugSlopes$se_mid,
               truebugSlopes$se_high),
  CI_lower = c(truebugSlopes$ci_low_l,
               truebugSlopes$ci_mid_l,
               truebugSlopes$ci_high_l),
  CI_upper = c(truebugSlopes$ci_low_u,
               truebugSlopes$ci_mid_u,
               truebugSlopes$ci_high_u),
  p_value  = c(truebugSlopes$p_low,
               truebugSlopes$p_mid,
               truebugSlopes$p_high))

truebugSlopes_long


truebug_slope_draws = truebugPost_draws %>%
  as.data.frame() %>%
  mutate(
    Low  = beta_dev + beta_int * -1,
    Mid  = beta_dev + beta_int * 0,
    High = beta_dev + beta_int * 1
  ) %>%
  dplyr::select(Low, Mid, High) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Latitude",
    values_to = "Slope"
  ) %>%
  mutate(Latitude = factor(Latitude,
                           levels = c("Low","Mid","High")))


ggplot(truebug_slope_draws,
       aes(x = Slope, y = Latitude, fill = Latitude)) +
  stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(name = "Latitude") +
  labs(
    x = "Urbanization effect",
    y = "",
    title = "truebugs: Posterior slopes by latitude"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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


hopperSlopes = hopperPost_draws %>%
  as.data.frame() %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    mean_low  = mean(slope_low),
    se_low    = sd(slope_low),
    ci_low_l  = quantile(slope_low, 0.025),
    ci_low_u  = quantile(slope_low, 0.975),
    
    mean_mid  = mean(slope_mid),
    se_mid    = sd(slope_mid),
    ci_mid_l  = quantile(slope_mid, 0.025),
    ci_mid_u  = quantile(slope_mid, 0.975),
    
    mean_high = mean(slope_high),
    se_high   = sd(slope_high),
    ci_high_l = quantile(slope_high, 0.025),
    ci_high_u = quantile(slope_high, 0.975),
    
    # Bayesian p-values (two-tailed)
    p_low  = 2 * min(mean(slope_low  > 0), mean(slope_low  < 0)),
    p_mid  = 2 * min(mean(slope_mid  > 0), mean(slope_mid  < 0)),
    p_high = 2 * min(mean(slope_high > 0), mean(slope_high < 0)))

hopperSlopes_long = tibble(
  Latitude = factor(c("low", "mid", "high"),
                    levels = c("low","mid","high")),
  Estimate = c(hopperSlopes$mean_low,
               hopperSlopes$mean_mid,
               hopperSlopes$mean_high),
  SE       = c(hopperSlopes$se_low,
               hopperSlopes$se_mid,
               hopperSlopes$se_high),
  CI_lower = c(hopperSlopes$ci_low_l,
               hopperSlopes$ci_mid_l,
               hopperSlopes$ci_high_l),
  CI_upper = c(hopperSlopes$ci_low_u,
               hopperSlopes$ci_mid_u,
               hopperSlopes$ci_high_u),
  p_value  = c(hopperSlopes$p_low,
               hopperSlopes$p_mid,
               hopperSlopes$p_high))

hopperSlopes_long


hopper_slope_draws = hopperPost_draws %>%
  as.data.frame() %>%
  mutate(
    Low  = beta_dev + beta_int * -1,
    Mid  = beta_dev + beta_int * 0,
    High = beta_dev + beta_int * 1
  ) %>%
  dplyr::select(Low, Mid, High) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Latitude",
    values_to = "Slope"
  ) %>%
  mutate(Latitude = factor(Latitude,
                           levels = c("Low","Mid","High")))


ggplot(hopper_slope_draws,
       aes(x = Slope, y = Latitude, fill = Latitude)) +
  stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(name = "Latitude") +
  labs(
    x = "Urbanization effect",
    y = "",
    title = "hoppers: Posterior slopes by latitude"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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



antSlopes = antPost_draws %>%
  as.data.frame() %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    mean_low  = mean(slope_low),
    se_low    = sd(slope_low),
    ci_low_l  = quantile(slope_low, 0.025),
    ci_low_u  = quantile(slope_low, 0.975),
    
    mean_mid  = mean(slope_mid),
    se_mid    = sd(slope_mid),
    ci_mid_l  = quantile(slope_mid, 0.025),
    ci_mid_u  = quantile(slope_mid, 0.975),
    
    mean_high = mean(slope_high),
    se_high   = sd(slope_high),
    ci_high_l = quantile(slope_high, 0.025),
    ci_high_u = quantile(slope_high, 0.975),
    
    # Bayesian p-values (two-tailed)
    p_low  = 2 * min(mean(slope_low  > 0), mean(slope_low  < 0)),
    p_mid  = 2 * min(mean(slope_mid  > 0), mean(slope_mid  < 0)),
    p_high = 2 * min(mean(slope_high > 0), mean(slope_high < 0)))

antSlopes_long = tibble(
  Latitude = factor(c("low", "mid", "high"),
                    levels = c("low","mid","high")),
  Estimate = c(antSlopes$mean_low,
               antSlopes$mean_mid,
               antSlopes$mean_high),
  SE       = c(antSlopes$se_low,
               antSlopes$se_mid,
               antSlopes$se_high),
  CI_lower = c(antSlopes$ci_low_l,
               antSlopes$ci_mid_l,
               antSlopes$ci_high_l),
  CI_upper = c(antSlopes$ci_low_u,
               antSlopes$ci_mid_u,
               antSlopes$ci_high_u),
  p_value  = c(antSlopes$p_low,
               antSlopes$p_mid,
               antSlopes$p_high))

antSlopes_long


ant_slope_draws = antPost_draws %>%
  as.data.frame() %>%
  mutate(
    Low  = beta_dev + beta_int * -1,
    Mid  = beta_dev + beta_int * 0,
    High = beta_dev + beta_int * 1
  ) %>%
  dplyr::select(Low, Mid, High) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Latitude",
    values_to = "Slope"
  ) %>%
  mutate(Latitude = factor(Latitude,
                           levels = c("Low","Mid","High")))


ggplot(ant_slope_draws,
       aes(x = Slope, y = Latitude, fill = Latitude)) +
  stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(name = "Latitude") +
  labs(
    x = "Urbanization effect",
    y = "",
    title = "ants: Posterior slopes by latitude"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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



grasshopperSlopes = grasshopperPost_draws %>%
  as.data.frame() %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    mean_low  = mean(slope_low),
    se_low    = sd(slope_low),
    ci_low_l  = quantile(slope_low, 0.025),
    ci_low_u  = quantile(slope_low, 0.975),
    
    mean_mid  = mean(slope_mid),
    se_mid    = sd(slope_mid),
    ci_mid_l  = quantile(slope_mid, 0.025),
    ci_mid_u  = quantile(slope_mid, 0.975),
    
    mean_high = mean(slope_high),
    se_high   = sd(slope_high),
    ci_high_l = quantile(slope_high, 0.025),
    ci_high_u = quantile(slope_high, 0.975),
    
    # Bayesian p-values (two-tailed)
    p_low  = 2 * min(mean(slope_low  > 0), mean(slope_low  < 0)),
    p_mid  = 2 * min(mean(slope_mid  > 0), mean(slope_mid  < 0)),
    p_high = 2 * min(mean(slope_high > 0), mean(slope_high < 0)))

grasshopperSlopes_long = tibble(
  Latitude = factor(c("low", "mid", "high"),
                    levels = c("low","mid","high")),
  Estimate = c(grasshopperSlopes$mean_low,
               grasshopperSlopes$mean_mid,
               grasshopperSlopes$mean_high),
  SE       = c(grasshopperSlopes$se_low,
               grasshopperSlopes$se_mid,
               grasshopperSlopes$se_high),
  CI_lower = c(grasshopperSlopes$ci_low_l,
               grasshopperSlopes$ci_mid_l,
               grasshopperSlopes$ci_high_l),
  CI_upper = c(grasshopperSlopes$ci_low_u,
               grasshopperSlopes$ci_mid_u,
               grasshopperSlopes$ci_high_u),
  p_value  = c(grasshopperSlopes$p_low,
               grasshopperSlopes$p_mid,
               grasshopperSlopes$p_high))

grasshopperSlopes_long


grasshopper_slope_draws = grasshopperPost_draws %>%
  as.data.frame() %>%
  mutate(
    Low  = beta_dev + beta_int * -1,
    Mid  = beta_dev + beta_int * 0,
    High = beta_dev + beta_int * 1
  ) %>%
  dplyr::select(Low, Mid, High) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Latitude",
    values_to = "Slope"
  ) %>%
  mutate(Latitude = factor(Latitude,
                           levels = c("Low","Mid","High")))


ggplot(grasshopper_slope_draws,
       aes(x = Slope, y = Latitude, fill = Latitude)) +
  stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(name = "Latitude") +
  labs(
    x = "Urbanization effect",
    y = "",
    title = "grasshoppers: Posterior slopes by latitude"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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


flySlopes = flyPost_draws %>%
  as.data.frame() %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    mean_low  = mean(slope_low),
    se_low    = sd(slope_low),
    ci_low_l  = quantile(slope_low, 0.025),
    ci_low_u  = quantile(slope_low, 0.975),
    
    mean_mid  = mean(slope_mid),
    se_mid    = sd(slope_mid),
    ci_mid_l  = quantile(slope_mid, 0.025),
    ci_mid_u  = quantile(slope_mid, 0.975),
    
    mean_high = mean(slope_high),
    se_high   = sd(slope_high),
    ci_high_l = quantile(slope_high, 0.025),
    ci_high_u = quantile(slope_high, 0.975),
    
    # Bayesian p-values (two-tailed)
    p_low  = 2 * min(mean(slope_low  > 0), mean(slope_low  < 0)),
    p_mid  = 2 * min(mean(slope_mid  > 0), mean(slope_mid  < 0)),
    p_high = 2 * min(mean(slope_high > 0), mean(slope_high < 0)))

flySlopes_long = tibble(
  Latitude = factor(c("low", "mid", "high"),
                    levels = c("low","mid","high")),
  Estimate = c(flySlopes$mean_low,
               flySlopes$mean_mid,
               flySlopes$mean_high),
  SE       = c(flySlopes$se_low,
               flySlopes$se_mid,
               flySlopes$se_high),
  CI_lower = c(flySlopes$ci_low_l,
               flySlopes$ci_mid_l,
               flySlopes$ci_high_l),
  CI_upper = c(flySlopes$ci_low_u,
               flySlopes$ci_mid_u,
               flySlopes$ci_high_u),
  p_value  = c(flySlopes$p_low,
               flySlopes$p_mid,
               flySlopes$p_high))

flySlopes_long


fly_slope_draws = flyPost_draws %>%
  as.data.frame() %>%
  mutate(
    Low  = beta_dev + beta_int * -1,
    Mid  = beta_dev + beta_int * 0,
    High = beta_dev + beta_int * 1
  ) %>%
  dplyr::select(Low, Mid, High) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Latitude",
    values_to = "Slope"
  ) %>%
  mutate(Latitude = factor(Latitude,
                           levels = c("Low","Mid","High")))


ggplot(fly_slope_draws,
       aes(x = Slope, y = Latitude, fill = Latitude)) +
  stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(name = "Latitude") +
  labs(
    x = "Urbanization effect",
    y = "",
    title = "flys: Posterior slopes by latitude"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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



daddylonglegsSlopes = daddylonglegsPost_draws %>%
  as.data.frame() %>%
  mutate(
    slope_low  = beta_dev + beta_int * -1,
    slope_mid  = beta_dev + beta_int * 0,
    slope_high = beta_dev + beta_int * 1
  ) %>%
  summarise(
    mean_low  = mean(slope_low),
    se_low    = sd(slope_low),
    ci_low_l  = quantile(slope_low, 0.025),
    ci_low_u  = quantile(slope_low, 0.975),
    
    mean_mid  = mean(slope_mid),
    se_mid    = sd(slope_mid),
    ci_mid_l  = quantile(slope_mid, 0.025),
    ci_mid_u  = quantile(slope_mid, 0.975),
    
    mean_high = mean(slope_high),
    se_high   = sd(slope_high),
    ci_high_l = quantile(slope_high, 0.025),
    ci_high_u = quantile(slope_high, 0.975),
    
    # Bayesian p-values (two-tailed)
    p_low  = 2 * min(mean(slope_low  > 0), mean(slope_low  < 0)),
    p_mid  = 2 * min(mean(slope_mid  > 0), mean(slope_mid  < 0)),
    p_high = 2 * min(mean(slope_high > 0), mean(slope_high < 0)))

daddylonglegsSlopes_long = tibble(
  Latitude = factor(c("low", "mid", "high"),
                    levels = c("low","mid","high")),
  Estimate = c(daddylonglegsSlopes$mean_low,
               daddylonglegsSlopes$mean_mid,
               daddylonglegsSlopes$mean_high),
  SE       = c(daddylonglegsSlopes$se_low,
               daddylonglegsSlopes$se_mid,
               daddylonglegsSlopes$se_high),
  CI_lower = c(daddylonglegsSlopes$ci_low_l,
               daddylonglegsSlopes$ci_mid_l,
               daddylonglegsSlopes$ci_high_l),
  CI_upper = c(daddylonglegsSlopes$ci_low_u,
               daddylonglegsSlopes$ci_mid_u,
               daddylonglegsSlopes$ci_high_u),
  p_value  = c(daddylonglegsSlopes$p_low,
               daddylonglegsSlopes$p_mid,
               daddylonglegsSlopes$p_high))

daddylonglegsSlopes_long


daddylonglegs_slope_draws = daddylonglegsPost_draws %>%
  as.data.frame() %>%
  mutate(
    Low  = beta_dev + beta_int * -1,
    Mid  = beta_dev + beta_int * 0,
    High = beta_dev + beta_int * 1
  ) %>%
  dplyr::select(Low, Mid, High) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Latitude",
    values_to = "Slope"
  ) %>%
  mutate(Latitude = factor(Latitude,
                           levels = c("Low","Mid","High")))


ggplot(daddylonglegs_slope_draws,
       aes(x = Slope, y = Latitude, fill = Latitude)) +
  stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(name = "Latitude") +
  labs(
    x = "Urbanization effect",
    y = "",
    title = "daddylonglegss: Posterior slopes by latitude"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
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




posteriorSlopes = bind_rows(
  Caterpillar = caterpillarSlopes_long,
  Beetle      = beetleSlopes_long,
  Spider      = spiderSlopes_long,
  TrueBug     = truebugSlopes_long,
  Hopper      = hopperSlopes_long,
  Ant         = antSlopes_long,
  Grasshopper = grasshopperSlopes_long,
  Daddylonglegs = daddylonglegsSlopes_long,
  Fly        = flySlopes_long,
  .id = "Arthropod"
)

posteriorSlopes %>% 
  filter(!Arthropod %in% c("Fly", "Daddylonglegs")) %>% 
  rename(
         "2.5% CI" =  "CI_lower",
         "97.5% CI" =  "CI_upper",
  ) %>% 
  gt() %>% 
  tab_header(
    title = "Table 1. Effect of urbanization on arthropod occurence at low, mid, and high latitudes"
  ) %>% 
  fmt_number(
    columns = where(is.numeric),
    decimals = 3
  ) %>% 
  gtsave("images/tableS2.png")

ggplot(posteriorSlopes, aes(Estimate, Arthropod, color = Latitude)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper),
                 position = position_dodge(width = 0.5)) +
  scale_color_viridis_d(name = "Latitude") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal()



diagnostics_table %>% 
  filter(!model %in% c("Fly", "DaddyLonglegs")) %>% 
  rename("Arthropod" = "model",
         "2.5% CI" =  "CI_lower",
         "97.5% CI" =  "CI_upper",
  ) %>% 
  gt() %>% 
  tab_header(
    title = "Table 1. Model summaries"
  ) %>% 
  fmt_number(
    columns = where(is.numeric),
    decimals = 3
  ) %>% 
  gtsave("images/table S1.png")

##################################################################################

 









































