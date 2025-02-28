library(dplyr)
library(terra)
library(ggplot2)
library(tidyterra)

# read in data
boundary <- vect(here::here("data/study_boundary.shp"))
cbi <- rast(here::here("data/cbi.tif"))

HUC12_uncrop <- vect(here::here("raw_data/WBD_HUC12_CONUS_pulled10262020/WBD_HUC12_CONUS_pulled10262020.shp")) %>%
  project("epsg:4326") %>%
  crop(boundary)
# #
# HUC12_forest <- vect(here::here("raw_data/WBD_HUC12_CONUS_pulled10262020/WBD_HUC12_CONUS_pulled10262020.shp")) %>%
#   project("epsg:4326") %>%
#   crop(cbi_forest)

metric_rast <- rast(here::here("data/metric_rast.tiff"))
fd_rast <- rast(here::here("data/fd_rast.tiff"))
biodiv_rast <- c(metric_rast, fd_rast)

ecoregion_shp <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  crop(., boundary)

ecoregion_rast <- ecoregion_shp %>%
  select("ECO_NAME") %>%
  terra::rasterize(., cbi, field = "ECO_NAME")

HUC12 <- vect(here::here("raw_data/WBD_HUC12_CONUS_pulled10262020/WBD_HUC12_CONUS_pulled10262020.shp")) %>%
  project("epsg:4326") %>%
  crop(ecoregion_shp)

# get ecoregion labels for each HUC12
# ecoregion_labels <- extract(ecoregion_rast, HUC12) %>%
#   # get number of pixels of each type intersecting each zone
#   count(ID, ECO_NAME) %>%
#   group_by(ID) %>%
#   # get rid of any zones that have NA's
#   filter(!any(is.na(ECO_NAME))) %>%
#   group_by(ID) %>%
#   # get label associated with the modal intersection
#   slice(which.max(n)) %>%
#   ungroup() %>%.
#   left_join(tibble(huc12 = HUC12$huc12) %>% mutate(ID = row_number())) %>%
#   select(huc12, ECO_NAME)

ecoregion_labels <- zonal(ecoregion_rast, HUC12, fun = "modal", na.rm = TRUE) %>%
  cbind(huc12 = HUC12$huc12, .) %>%
  filter(!is.na(ECO_NAME))

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

# get frequency of CBI values within polygons to figure out what kind of hotspot they are
# want to exclude ecoregions that are majority non-forested (either cbi == 0 or NA)
huc12_cbi_modal_vect <- zonal(cbi, HUC12, fun = "modal", na.rm = TRUE, as.polygons = TRUE)

# get the break down of CBI types within each huc12 to classify what kind of hotspot it is
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
  # right_join(values(huc12_cbi_modal_vect) %>% select(huc12, predict.high.severity.fire.final)) %>%
  # # and remove them when we're done
  # select(-c(ID, predict.high.severity.fire.final)) %>%
  filter(no_data != 1 | is.na(no_data))

#biodiv_cbi_huc12 <- biodiv_zonal_vec %>% filter(huc12 %in% huc_cbi_perc$huc12) %>% left_join(huc_cbi_perc)

# get HUC12 level means of biodiversity metrics
biodiv_zonal <- zonal(biodiv_rast, HUC12_uncrop, fun = "mean", na.rm = TRUE) %>%
  cbind(huc12 = HUC12_uncrop$huc12, .) %>%
  left_join(ecoregion_labels) %>%
  left_join(values(huc12_cbi_modal_vect) %>% select(huc12, fire = predict.high.severity.fire.final)) %>%
  left_join(huc_cbi_perc) %>%
  mutate(max_severity = case_when(
    low_sev == pmax(low_sev, high_sev, na.rm = TRUE) ~ "low_sev",
    high_sev == pmax(low_sev, high_sev, na.rm = TRUE) ~ "high_sev"
  )) %>%
  rowwise() %>%
  mutate(forest_total = sum(low_sev, high_sev, na.rm = TRUE))

# rejoin with polygon labels and geometry
biodiv_zonal_vec <- HUC12_uncrop %>% select(huc12) %>%
  right_join(biodiv_zonal)

# get hotspots after filtering out non-forest pixels
noforest_zonal <- biodiv_rast %>%
  crop(cbi, mask = TRUE) %>%
  c(., cbi) %>%
  filter(predict.high.severity.fire.final %in% c(1,2)) %>%
  select(-predict.high.severity.fire.final) %>%
  zonal(., HUC12_uncrop, fun = "mean", na.rm = TRUE) %>%
  cbind(huc12 = HUC12_uncrop$huc12, .) %>%
  left_join(ecoregion_labels) %>%
  left_join(values(huc12_cbi_modal_vect) %>% select(huc12, fire = predict.high.severity.fire.final)) %>%
  left_join(huc_cbi_perc) %>%
  mutate(max_severity = case_when(
    low_sev == pmax(low_sev, high_sev, na.rm = TRUE) ~ "low_sev",
    high_sev == pmax(low_sev, high_sev, na.rm = TRUE) ~ "high_sev"
  )) %>%
  rowwise() %>%
  mutate(forest_total = sum(low_sev, high_sev, na.rm = TRUE)) %>%
  filter(!is.na(max_severity), !is.na(ECO_NAME))

noforest_ecoregion_hotspot <- noforest_zonal %>%
  #filter(forest_total > 0.20) %>%
  group_by(ECO_NAME) %>%
  mutate(across(-c(huc12, fire, low_sev, high_sev, no_data, unforested, max_severity, forest_total),
                ~ifelse(.x > quantile(.x, probs = 0.90, na.rm = TRUE), 1, NA)),
         hotspot_type = case_when(
           ## THIS IS THE SENSITIVITY FOR BEING CONSIDERED "MIXED"
           abs(low_sev - high_sev) < 0.1 ~ "Mixed",
           max_severity == "high_sev" ~ "Area of Concern",
           max_severity == "low_sev" ~ "Refugia"
         )) %>%
  ungroup()

ecoregion_hotspot_zonal_vec <- HUC12 %>% select(huc12) %>%
  right_join(noforest_ecoregion_hotspot)

# # get spatvector of hotspots for each metric
# hotspot_zonal_vec <- biodiv_zonal_vec %>%
#   # get hotspots only in forested areas
#   filter(!is.na(ECO_NAME), !is.na(max_severity), forest_total > 0.3) %>%
#   mutate(across(-c(huc12, ECO_NAME, fire, low_sev, high_sev, no_data, unforested, max_severity, forest_total),
#                 ~ifelse(.x > quantile(.x, probs = 0.95, na.rm = TRUE), 1, NA)),
#          hotspot_type = case_when(
#            abs(low_sev - high_sev) < 0.1 ~ "Mixed",
#            fire == 2 ~ "Area of Concern",
#            fire == 1 ~ "Refugia"
#          ))
#
# ecoregion_hotspot_zonal_vec <- HUC12 %>% select(huc12) %>%
#   right_join(biodiv_zonal %>%
#   ## A WATERSHED CAN BE A HOTSPOT IF AT LEAST 20% IS FORESTED ##
#   filter(!is.na(ECO_NAME), !is.na(max_severity), forest_total > 0.20) %>%
#   group_by(ECO_NAME) %>%
#   mutate(across(-c(huc12, fire, low_sev, high_sev, no_data, unforested, max_severity, forest_total),
#                 ~ifelse(.x > quantile(.x, probs = 0.95, na.rm = TRUE), 1, NA)),
#          hotspot_type = case_when(
#            ## THIS IS THE SENSITIVITY FOR BEING CONSIDERED "MIXED"
#            abs(low_sev - high_sev) < 0.1 ~ "Mixed",
#            max_severity == "high_sev" ~ "Area of Concern",
#            max_severity == "low_sev" ~ "Refugia"
#          )) %>%
#   ungroup())

writeVector(biodiv_zonal_vec, here::here("data/biodiv_zonal.shp"))
writeVector(ecoregion_hotspot_zonal_vec, here::here("data/ecoregion_hotspots_huc12.shp"))

# how many hotspots fall in different fire severity categories for each metric
hotspot_cat <- ecoregion_hotspot_zonal_vec %>%
  as.data.frame() %>%
  select(huc12, breeding_richness, ecoregion_breeding_lcbd,
         FRic_breeding, hotspot_type) %>%
  pivot_longer(-c(huc12, hotspot_type), names_to = "metric", values_to = "type") %>%
  filter(!is.na(type)) %>%
  left_join(as.data.frame(HUC12) %>% select(huc12, areasqkm)) %>%
  select(-type) %>%
  group_by(hotspot_type, metric) %>%
  summarize(total_area = sum(areasqkm)) %>%
  group_by(metric) %>%
  mutate(perc = total_area/sum(total_area))

# plot an exsummarise()# plot an example watershed for debugging
watershed = "180201530301"
ggplot() +
  geom_spatraster(data = cbi %>% mutate(predict.high.severity.fire.final = as.factor(predict.high.severity.fire.final)) %>% crop(HUC12_uncrop %>% filter(huc12 == watershed))) +
  scale_fill_manual(values = list(`1` = "#D3B073", `2` =  "#BC4749"), na.value = "transparent", na.translate = FALSE,
                    labels = c("Low Severity", "High Severity")) +
  geom_spatvector(data = HUC12 %>% filter(huc12 == watershed), fill = "transparent") +
  theme_void()

