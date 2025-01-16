library(dplyr)
library(terra)
library(ggplot2)
library(tidyterra)

# read in data
boundary <- vect(here::here("data/study_boundary.shp"))
cbi <- rast(here::here("data/cbi.tif"))

HUC12 <- vect(here::here("raw_data/WBD_HUC12_CONUS_pulled10262020/WBD_HUC12_CONUS_pulled10262020.shp")) %>%
  project("epsg:4326") %>%
  crop(cbi)

metric_rast <- rast(here::here("data/metric_rast.tiff"))
fd_rast <- rast(here::here("data/fd_rast.tiff"))
biodiv_rast <- c(metric_rast, fd_rast)

ecoregion_shp <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  crop(., boundary)

#huc12_ecoregion <- intersect(HUC12 %>% select(huc12), ecoregion_shp %>% select(ECO_NAME, ECO_NUM))

# mtbs <- vect(here::here("raw_data/mtbs_perimeter_data/mtbs_perims_DD.shp")) %>%
#   project("epsg:4326") %>%
#   crop(boundary)
#
# cbi_mask <- cbi %>%
#   mask(mtbs %>% select(Map_ID))
#
# ggplot() +
#   geom_spatraster(data = cbi_mask) +
#   geom_spatvector(data = mtbs)

# get HUC12 level means of biodiversity metrics
biodiv_zonal <- zonal(biodiv_rast, HUC12, fun = "mean", na.rm = TRUE) %>%
  cbind(huc12 = HUC12$huc12, .)

# rejoin with polygon labels and geometry
biodiv_zonal_vec <- HUC12 %>% select(huc12) %>%
  left_join(biodiv_zonal)

# get spatvector of hotspots for each metric
hotspot_zonal_vec <- biodiv_zonal_vec %>%
  mutate(across(-c(huc12),
                ~ifelse(.x > quantile(.x, probs = 0.95, na.rm = TRUE), 1, NA)))

# get frequency of CBI values within polygons to figure out what kind of hotspot they are
huc12_cbi_modal_vect <- zonal(cbi, HUC12, fun = "modal", na.rm = TRUE, as.polygons = TRUE)

huc_cbi_perc <- extract(cbi, HUC12) %>%
  select(ID, predict.high.severity.fire.final) %>%
  count(ID, predict.high.severity.fire.final) %>%
  pivot_wider(names_from = predict.high.severity.fire.final, values_from = n) %>%
  #select(-c(`0`, `NA`)) %>%
  rename(unforested = `0`, low_sev = `1`, high_sev = `2`, no_data = `NA`) %>%
  rowwise() %>%
  mutate(across(-c(ID), ~ .x/sum(low_sev,high_sev, unforested, no_data, na.rm = TRUE))) %>%
  cbind(HUC12$huc12) %>%
  rename(huc12 = `HUC12$huc12`) %>%
  # we can add back in the model results from zonal to double check it's ok to cbind the huc12 ID's
  right_join(values(huc12_cbi_modal_vect) %>% select(huc12, predict.high.severity.fire.final)) %>%
  # and remove them when we're done
  select(-c(ID, predict.high.severity.fire.final)) %>%
  filter(no_data != 1)

biodiv_cbi_huc12 <- biodiv_zonal_vec %>% filter(huc12 %in% huc_cbi_perc$huc12) %>% left_join(huc_cbi_perc)

test_plot <- ggplot() +
  geom_spatvector(data = biodiv_zonal_vec, aes(fill = breeding_richness), color = "transparent")

ggsave(here::here("figures/test_map.jpeg"), test_plot)

hotspot_plot <- ggplot() +
  geom_spatvector(data = hotspot_zonal_vec, aes(fill = breeding_richness), color = "transparent")

ggsave(here::here("figures/huc12_hotspot.jpeg"), hotspot_plot)
