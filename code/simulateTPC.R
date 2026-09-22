
rm(list = ls())

library(tidyverse)
temp <- seq(0, 60, length.out = 300)
# Broad thermal optimum with left skew
broad <- ifelse(
  temp <= 41,
  exp(-((temp - 41)/19)^2),  
  exp(-((temp - 41)/5)^2))

broad <- (broad / max(broad)) * 0.7 # scale the height

# Slightly left-skewed Gaussian-like curve
normal <- ifelse(
  temp <= 35,
  exp(-((temp - 35)/15)^2),
  exp(-((temp - 35)/5)^2))

normal <- normal / max(normal)

# Strongly left-skewed curve
left <- dlnorm(55 - temp, meanlog = 2.35, sdlog = 0.65)
left <- left / max(left)

df <- bind_rows(
  tibble(temp, performance = broad, curve = "Broad optimum"),
  tibble(temp, performance = normal, curve = "Gaussian-like"),
  tibble(temp, performance = left, curve = "Left-skewed"))

df <- df %>%
  mutate(curve = factor(curve,
                        levels = c("Gaussian-like",  "Broad optimum", "Left-skewed"))) %>%
  arrange(curve)

ggplot(df, aes(temp, performance, fill = curve, colour = curve)) +
  geom_area(alpha = 0.5, position = "identity") +
  geom_line(linewidth = 0.7, color = "NA") +
  scale_fill_manual(values = c(
    "Broad optimum" = "#009E73",     
    "Left-skewed" =  "#D55E00",       
    "Gaussian-like" = "#0072B2" )) +
  scale_colour_manual(values = c(
    "Broad optimum" = "#009E73",     
    "Left-skewed" =  "#D55E00",       
    "Gaussian-like" = "#0072B2" )) +
  # theme_void()
  theme_classic()

ggplot(df, aes(temp, performance, fill = curve, colour = curve)) +
  annotate("rect",
           xmin = -Inf, xmax = (min(df$temp) + max(df$temp)) / 1/1.65,
           ymin = -Inf, ymax = Inf,
           fill = "#0072B2", alpha = 0.12) +
  annotate("rect",
           xmin = (min(df$temp) + max(df$temp)) / 1/1.65, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "orange", alpha = 0.12) +
  geom_area(alpha = 0.7, position = "identity") + # alpha = 0.5
  geom_line(linewidth = 0.7) +
  #facet_wrap(~ curve, nrow = 1) +
  scale_fill_manual(values = c(
    "Broad optimum" = "#009E73",
    "Left-skewed" = "#D55E00",
    "Gaussian-like" = "#0072B2")) +
  scale_colour_manual(values = c(
    "Broad optimum" = "#009E73",
    "Left-skewed" = "#D55E00",
    "Gaussian-like" = "#0072B2")) +
  theme_void() +
  theme(legend.position = "none")


ggplot(df %>% filter(curve == "Gaussian-like"), aes(temp, performance, fill = curve, colour = curve)) +
  geom_area(alpha = 0.5, position = "identity") +
  geom_line(linewidth = 0.7, color = "white") +
  scale_fill_manual(values = c(
    "Gaussian-like" = "#0072B2"  )) +
  scale_colour_manual(values = c(
    "Gaussian-like" = "#0072B2" )) +
  theme_classic(base_size = 14) +
  labs(x = "Temperature (°C)", y = "Thermal performance")


####################################################################################

# Extract Gaussian-like curve
gaussian_df <- tibble(
  temp = temp,
  performance = normal
)

# Create modified segments
cool_segment <- gaussian_df %>%
  filter(temp >= 23, temp <= 32) %>%
  mutate(
    performance = performance - 0.1,
    curve = "Cool segment (-0.1)"
  )

warm_segment <- gaussian_df %>%
  filter(temp >= 31, temp <= 42) %>%
  mutate(
    performance = performance + 0.1,
    curve = "Warm segment (+0.1)"
  )

# Original curve
original <- gaussian_df %>%
  mutate(curve = "Gaussian-like")

# Combine
three_lines <- bind_rows(
  original,
  cool_segment,
  warm_segment)



ggplot() +
  geom_area(
    data = original,
    aes(temp, performance),
    fill = "#0072B2",
    alpha = 0.35) +
  # Original Gaussian-like outline
  geom_line(
    data = original,
    aes(temp, performance),
    colour = "#0072B2",
    linewidth = 1.2) +
  # Cool segment shifted downward
  geom_line(
    data = cool_segment,
    aes(temp, performance),
    colour = "#0072B2",
    linewidth = 2,
    linetype = 2) +
  # Warm segment shifted upward
  geom_line(
    data = warm_segment,
    aes(temp, performance),
    colour = "#D55E00",
    linewidth = 2,
    linetype = 2) +
  theme_classic(base_size = 14) +
  labs(
    x = "Temperature (°C)",
    y = "Thermal performance")

ggplot() +
  geom_area(
    data = original,
    aes(temp, performance),
    fill = "#0072B2",
    alpha = 0.35) +
  geom_line(
    data = original,
    aes(temp, performance),
    colour = "#0072B2",
    linewidth = 1.2) +
  geom_line(
    data = cool_segment,
    aes(temp, performance),
    colour = "#0072B2",
    linewidth = 2,
    linetype = 2) +
  geom_line(
    data = warm_segment,
    aes(temp, performance),
    colour = "#D55E00",
    linewidth = 2,
    linetype = 2) +
  theme_void()


#################################################################################


library(ggplot2)
library(dplyr)

df <- tibble(
  temp = seq(5, 35, length.out = 500)
) %>%
  mutate(
    performance = 1 / (1 + exp(-(temp - 25) / 4)))

df <- df %>%
  mutate(
    dydx = c(NA, diff(performance) / diff(temp))
  )

ggplot(df, aes(temp, performance)) +
  geom_line(linewidth = 1) +
  theme_classic()

