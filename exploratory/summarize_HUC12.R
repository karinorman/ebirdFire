library(dplyr)
library(terra)
library(ggplot2)
library(tidyterra)

boundary <- vect(here::here("data/study_boundary.shp"))
cbi <- rast(here::here("data/cbi.tif"))

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

HUC12 <- vect(here::here("raw_data/WBD_HUC12_CONUS_pulled10262020/WBD_HUC12_CONUS_pulled10262020.shp")) %>%
  project("epsg:4326")

# Get UT watersheds for fun
# ut_huc12 <- HUC12 %>% filter(stringr::str_detect(states, "UT"))
# ut_huc12_rast <- ut_huc12 %>% rasterize(hotspot_rast)

metric_rast <- rast(here::here("data/metric_rast.tiff"))
test_zonal <- zonal(metric_rast, ut_huc12, fun = "mean", na.rm = TRUE, as.polygons = TRUE)

# get frequency of values within polygons
huc12_cbi_modal_vect <- zonal(cbi, ut_huc12, fun = "modal", na.rm = TRUE, as.polygons = TRUE)

huc_cbi_perc <- extract(cbi, ut_huc12) %>%
  select(ID, predict.high.severity.fire.final) %>%
  count(ID, predict.high.severity.fire.final) %>%
  pivot_wider(names_from = predict.high.severity.fire.final, values_from = n) %>%
  #select(-c(`0`, `NA`)) %>%
  rename(unforested = `0`, low_sev = `1`, high_sev = `2`, no_data = `NA`) %>%
  rowwise() %>%
  mutate(across(-c(ID), ~ .x/sum(low_sev,high_sev, unforested, no_data, na.rm = TRUE))) %>%
  cbind(ut_huc12$huc12) %>%
  rename(huc12 = `ut_huc12$huc12`) %>%
  # we can add back in the model results from zonal to double check it's ok to cbind the huc12 ID's
  right_join(values(huc12_cbi_modal_vect) %>% select(huc12, predict.high.severity.fire.final)) %>%
  # and remove them when we're done
  select(-c(ID, predict.high.severity.fire.final))
