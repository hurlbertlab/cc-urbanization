
rm(list = ls())

library(tidyverse)


temp <- seq(0, 40, length.out = 500)

# Broad thermal optimum with left skew
broad <- ifelse(
  temp <= 27,
  exp(-((temp - 27)/13)^4),  
  exp(-((temp - 27)/8)^4)
)

broad <- broad / max(broad)
broad <- broad * 0.4


# Slightly left-skewed Gaussian-like curve
normal <- ifelse(
  temp <= 25,
  exp(-((temp - 25)/9)^2),
  exp(-((temp - 25)/6)^2)
)

normal <- normal / max(normal)


# Strongly left-skewed curve
left <- dlnorm(40 - temp, meanlog = 2.35, sdlog = 0.35)
left <- left / max(left)


df <- bind_rows(
  tibble(temp, performance = broad, curve = "Broad optimum"),
  tibble(temp, performance = normal, curve = "Gaussian-like"),
  tibble(temp, performance = left, curve = "Left-skewed"))


df <- df %>%
  mutate(curve = factor(curve,
                        levels = c("Gaussian-like", "Left-skewed", "Broad optimum"))) %>%
  arrange(curve)

ggplot(df, aes(temp, performance, fill = curve, colour = curve)) +
  geom_area(alpha = 0.5, position = "identity") +
  geom_line(linewidth = 0.7, color = "NA") +
  scale_fill_manual(values = c(
    "Broad optimum" = "#0072B2",    # vermillion/orange-red
    "Left-skewed" =  "#D55E00",      # green
    "Gaussian-like" = "#009E73"    
  )) +
  scale_colour_manual(values = c(
    "Broad optimum" = "#0072B2",    # vermillion/orange-red
    "Left-skewed" =  "#D55E00",      # green
    "Gaussian-like" = "#009E73" 
  )) +
  theme_void()



ggplot(df %>% filter(curve == "Gaussian-like"), aes(temp, performance, fill = curve, colour = curve)) +
  geom_area(alpha = 0.5, position = "identity") +
  geom_line(linewidth = 0.7, color = "white") +
  scale_fill_manual(values = c(
    "Gaussian-like" = "#009E73"    
  )) +
  scale_colour_manual(values = c(
    "Gaussian-like" = "#009E73" 
  )) +
  theme_classic(base_size = 14) +
  labs(
    x = "Temperature (°C)",
    y = "Thermal performance")
####################################################################################

# Extract Gaussian-like curve
gaussian_df <- tibble(
  temp = temp,
  performance = normal
)

# Create modified segments
cool_segment <- gaussian_df %>%
  filter(temp >= 17, temp <= 26) %>%
  mutate(
    performance = performance - 0.1,
    curve = "Cool segment (-0.1)"
  )

warm_segment <- gaussian_df %>%
  filter(temp >= 22, temp <= 32) %>%
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
    fill = "#009E73",
    alpha = 0.35
  ) +
  # Original Gaussian-like outline
  geom_line(
    data = original,
    aes(temp, performance),
    colour = "#009E73",
    linewidth = 1.2
  ) +
  # Cool segment shifted downward
  geom_line(
    data = cool_segment,
    aes(temp, performance),
    colour = "#0072B2",
    linewidth = 2,
    linetype = 2
  ) +
  # Warm segment shifted upward
  geom_line(
    data = warm_segment,
    aes(temp, performance),
    colour = "#D55E00",
    linewidth = 2,
    linetype = 2
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Temperature (°C)",
    y = "Thermal performance")

ggplot() +
  geom_area(
    data = original,
    aes(temp, performance),
    fill = "#009E73",
    alpha = 0.35
  ) +
  geom_line(
    data = original,
    aes(temp, performance),
    colour = "#009E73",
    linewidth = 1.2
  ) +
  geom_line(
    data = cool_segment,
    aes(temp, performance),
    colour = "#0072B2",
    linewidth = 2,
    linetype = 2
  ) +
  geom_line(
    data = warm_segment,
    aes(temp, performance),
    colour = "#D55E00",
    linewidth = 2,
    linetype = 2
  ) +
  theme_void()
