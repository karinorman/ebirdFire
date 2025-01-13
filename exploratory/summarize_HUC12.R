library(dplyr)
library(terra)
library(ggplot2)
library(tidyterra)

boundary <- vect(here::here("data/study_boundary.shp"))
cbi <- rast(here::here("data/cbi.tif"))

mtbs <- vect(here::here("raw_data/mtbs_perimeter_data/mtbs_perims_DD.shp")) %>%
  project("epsg:4326") %>%
  crop(boundary)


cbi_mask <- cbi %>%
  mask(mtbs %>% select(Map_ID))

ggplot() +
  geom_spatraster(data = cbi_mask) +
  geom_spatvector(data = mtbs)

HUC12 <- vect(here::here("raw_data/WBD_HUC12_CONUS_pulled10262020/WBD_HUC12_CONUS_pulled10262020.shp"))

# Get UT watersheds for fun
ut_huc12 <- HUC12 %>% filter(stringr::str_detect(states, "UT"))
ut_huc12_rast <- ut_huc12 %>% rasterize(hotspot_rast)

hotspot_rast <- rast(here::here("data/hotspot_rasts.tif"))

huc_raster <- HUC12 %>% crop(boundary) %>%
  rasterize(hotspot_rast, field = "huc12")

huc_hotspots <- c(huc_raster, hotspot_rast) %>%
  mutate(across(-huc12, ~ifelse(.x == 1, 1, 0))) %>%
  group_by(huc12) %>%
  summarize(across(everything(), mean))

huc_raster %>%
  group_b

ut_hotspot <- hotspot_rast %>%
  select(breeding_richness) %>%
  mutate(breeding_richness = ifelse(breeding_richness == 1, 1, 0)) %>%
  crop(ut_huc12) %>%
  extract(ut_huc12, fun = "mean", na.rm = TRUE)
