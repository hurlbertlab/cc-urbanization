# extract land cover for year 2020 in usa and canada.
# urban scale = 1000 (1km), but you may change it

usethis::use_git_ignore("largeFile/")

library(rnaturalearth)
library(terra)
library(sf)
library(tidyverse)
library(viridis)

# (1) Read in latest Caterpillars Count! raw dataset from the caterpillars-analysis-public repo----
options(timeout = 300)  

api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)

dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]

# pick the latest one
latest_file <- dataset_file[1]

github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"

fullDataset <- read.csv(paste0(github_raw, latest_file)) %>% filter(Year <= 2025)


nlcd = rast("largeFile/Annual_NLCD_LndCov_2020_CU_C1V0.tif")
clcd = rast("largeFile/Canada-landcover-2020-classification.tif")





# -----------------------------
# 1. CREATE SITES
# -----------------------------
sites <- fullDataset %>%
  distinct(Name, Region, Longitude, Latitude) %>%
  mutate(ID = row_number())


# -----------------------------
# NORTH AMERICA MASK
# -----------------------------
northAmerica <- rnaturalearth::ne_countries(
  country = c("United States of America", "Canada"),
  returnclass = "sf")

na_union <- st_union(northAmerica)
northAmerica_buff <- sf::st_buffer(northAmerica, 0.02)


sites_sf <- sites %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  sf::st_join(
    northAmerica_buff["admin"],
    join = sf::st_intersects,
    left = TRUE
  )

sites_sf %>%
  filter(is.na(admin)) %>%
  st_drop_geometry()
# sites_sf <- st_join(
#   sites_sf,
#   northAmerica["admin"],
#   join = st_within)


# -----------------------------
# FILTER SITES
# -----------------------------

sites_us <- sites_sf %>%
  filter(admin == "United States of America")

sites_can <- sites_sf %>%
  filter(admin == "Canada")

# -----------------------------
# STATES + PROVINCES
# -----------------------------
states <- rnaturalearth::ne_states(
  country = "United States of America",
  returnclass = "sf")

provinces <- rnaturalearth::ne_states(
  country = "Canada",
  returnclass = "sf")


ggplot() +
  geom_sf(data = northAmerica, fill = "gray95", color = "white") +
  geom_sf(data = states, fill = NA, color = "gray60", linewidth = 0.25) +
  geom_sf(data = provinces, fill = NA, color = "gray80", linewidth = 0.4) +
  geom_sf(data = sites_sf, color = "blue", size = 2) +
  theme_minimal()


# -----------------------------
# 3. NLCD PIPELINE
# -----------------------------
sites_nlcd <- st_transform(sites_us, crs(nlcd))
buffers_nlcd <- st_buffer(sites_nlcd, dist = 500)

nlcd_vals <- terra::extract(nlcd, vect(buffers_nlcd))

# rename raster column safely
names(nlcd_vals)[2] <- "landcover"

# REMOVE NA 
nlcd_vals <- nlcd_vals %>%
  filter(!is.na(landcover))

nlcd_vals$ID <- buffers_nlcd$ID[nlcd_vals$ID]

# -----------------------------
# 4. CANADA PIPELINE
# -----------------------------
sites_can <- st_transform(sites_can, crs(clcd))
buffers_can <- st_buffer(sites_can, dist = 500)

can_vals <- terra::extract(clcd, vect(buffers_can))

names(can_vals)[2] <- "landcover"

can_vals <- can_vals %>%
  filter(!is.na(landcover))

can_vals$ID <- buffers_can$ID[can_vals$ID]
# -----------------------------
# 5. URBAN CLASS DEFINITIONS
# -----------------------------
urban_classes_nlcd <- c(22, 23, 24)
urban_class_can <- 17


# -----------------------------
# 6. NLCD SUMMARY
# https://www.usgs.gov/media/images/annual-nlcd-land-cover-change-legend
# -----------------------------
npix_nlcd <- nlcd_vals %>%
  group_by(ID) %>%
  summarise(numpix = n(), .groups = "drop")

# For each site (ID), count how many pixels belong to each
# landcover category (e.g., forest, urban, agriculture).

land_nlcd <- nlcd_vals %>%
  group_by(ID, landcover) %>%
  summarise(count = n(), .groups = "drop") %>% 
  right_join(npix_nlcd, by = "ID") %>%
  mutate(percent = 100 * count / numpix)

# Filter only urban landcover classes (low, medium, high intensity)
urban_nlcd <- land_nlcd %>% 
  filter(landcover %in% urban_classes_nlcd) %>%
  group_by(ID) %>%
  summarise(urban_percent = sum(percent), .groups = "drop")

 
urban_nlcd <- sites %>% right_join(npix_nlcd, by = "ID") %>% 
  dplyr::select(ID) %>%
  left_join(urban_nlcd, by = "ID") %>%
  mutate(urban_percent = ifelse(is.na(urban_percent), 0, urban_percent),
         source = "NLCD")


# -----------------------------
# 7. CANADA SUMMARY
# -----------------------------
npix_can <- can_vals %>%
  group_by(ID) %>%
  summarise(numpix = n(), .groups = "drop")

land_can <- can_vals %>%
  group_by(ID, landcover) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(npix_can, by = "ID") %>%
  mutate(percent = 100 * count / numpix)

urban_can_summary <- land_can %>%
  filter(landcover == urban_class_can) %>%
  group_by(ID) %>%
  summarise(urban_percent = sum(percent), .groups = "drop")

urban_can <- sites %>%
  right_join(npix_can, by = "ID") %>%
  dplyr::select(ID) %>%
  left_join(urban_can_summary, by = "ID") %>%
  mutate(
    urban_percent = coalesce(urban_percent, 0),
    source = "Canada")

# -----------------------------
# 8. COMBINE RESULTS
# -----------------------------
urban_all <- bind_rows(urban_nlcd, urban_can)

sites_final <- sites %>%
  right_join(urban_all, by = "ID")


# -----------------------------
# 9. FINAL CHECKS
# -----------------------------
sum(is.na(sites_final$urban_percent))

sites_final %>%
  filter(is.na(urban_percent))

  write.csv(x = sites_final, file = "data/urbanization_USA_Canada_500m.csv")

#########################################################

states <- rnaturalearth::ne_states(
  country = "United States of America",
  returnclass = "sf")

provinces <- rnaturalearth::ne_states(
  country = "Canada",
  returnclass = "sf")

crop_bbox <- st_bbox(c(
  xmin = -100,   # west cutoff
  xmax = -65,    # east Atlantic
  ymin = 25,     # southern extent 
  ymax = 50),      # latitude cutoff
 crs = st_crs(northAmerica))

crop_poly = st_as_sfc(crop_bbox)

# -----------------------------
#  CROP EVERYTHING

sites_urban_sf <- sites_sf %>%
  right_join(urban_all, by = "ID") # adding urbanization info
states_crop <- st_intersection(states, crop_poly) # Clip US state polygons to the defined crop polygon
provinces_crop <- st_intersection(provinces, crop_poly)
sites_crop <- st_intersection(sites_urban_sf, crop_poly) # Clip site locations to the study area
na_crop <- st_intersection(northAmerica, crop_poly) # Clip the N.America country polygons to the same extent

# -----------------------------
ggplot() +
  geom_sf(data = na_crop, fill = "gray95", color = "white") +
  geom_sf(data = states_crop, fill = NA, color = "gray70", linewidth = 0.3) +
  geom_sf(data = provinces_crop, fill = NA, color = "gray50", linewidth = 0.3) +
  geom_sf(data = sites_crop, color = "blue", size = 3) +
  coord_sf(
    xlim = c(-100, -65),
    ylim = c(25, 50),
    expand = FALSE) +
  theme_minimal() +
  labs(title = "")


ggplot() +
  geom_sf(data = na_crop, fill = "gray95", color = "white") +
  geom_sf(data = states_crop, fill = NA, color = "gray70", linewidth = 0.3) +
  geom_sf(data = provinces_crop, fill = NA, color = "gray50", linewidth = 0.3) +
  geom_sf(
    data = sites_crop,
    aes(color = urban_percent),
    size = 3) +
  scale_color_viridis_c(option = "viridis", na.value = "transparent", direction = 1) +
  coord_sf(
    xlim = c(-100, -65),
    ylim = c(25, 50),
    expand = FALSE) +
  theme_minimal() +
  labs(
    color = "Urban (%)",
    title = "")

sum(!is.na(sites_crop$urban_percent))
table(is.na(sites_crop$urban_percent))

lost_sites <- sites_urban_sf %>%
  anti_join(
    st_drop_geometry(sites_crop),
    by = "ID"
  )

lost_sites %>%
  dplyr::select(Name, Region)




