library(dplyr)
library(terra)
library(tidyterra)
library(purrr)
library(ggplot2)
library(cowplot)

#boundary to clip shapefiles
boundary <- vect(here::here("data/study_boundary.shp"))

cbi <- rast(here::here("data/cbi.tif")) %>%
  crop(boundary, mask = TRUE)

# forest_poly <- cbi %>% filter(predict.high.severity.fire.final %in% c(1,2)) %>%
#   as.polygons(extent = FALSE, na.rm = TRUE, crs = "epsg:4326")

# get ecoregions as shapefile
ecoregion_shp <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  crop(., boundary)

ecoregion_rast <- ecoregion_shp %>%
  select("ECO_NAME") %>%
  terra::rasterize(., cbi, field = "ECO_NAME")# %>%
#resample(cbi, method = "near")

metric_rast <- rast(here::here("data/metric_rast.tiff")) %>%
  crop(ecoregion_shp, mask = TRUE) %>%
  c(., rast(here::here("data/fd_rast.tiff"))%>%
      crop(ecoregion_shp, mask = TRUE)) %>%
  c(., ecoregion_rast)

ecoregion_hotspots <- metric_rast %>%
  terra::as.data.frame(., xy = TRUE) %>%
  group_by(ECO_NAME) %>%
  mutate(across(-c(x,y),
                ~ifelse(.x > quantile(.x, probs = 0.95, na.rm = TRUE), 1, NA))) %>%
  ungroup() %>%
  select(-ECO_NAME) %>%
  rast(.,  type="xyz", crs = "epsg:4326")

forest_hotspots <- c(cbi, metric_rast) %>%
  filter(predict.high.severity.fire.final %in% c(1,2)) %>%
  terra::as.data.frame(., xy = TRUE) %>%
  group_by(ECO_NAME) %>%
  mutate(across(-c(x,y),
                ~ifelse(.x > quantile(.x, probs = 0.95, na.rm = TRUE), 1, NA))) %>%
  ungroup() %>%
  select(-ECO_NAME) %>%
  rast(.,  type="xyz", crs = "epsg:4326") %>%
  rename_with(~paste0("forest_", .x))

forest_quantile_cutoffs <- c(cbi, metric_rast) %>%
  filter(predict.high.severity.fire.final %in% c(1,2)) %>%
  terra::as.data.frame(., xy = FALSE) %>%
  select(-predict.high.severity.fire.final) %>%
  group_by(ECO_NAME) %>%
  summarize(across(everything(), ~quantile(.x, probs = 0.95, na.rm = TRUE))) %>%
  ungroup()

hotspot_rasts <- c(ecoregion_hotspots, forest_hotspots) %>%
  select(-forest_predict.high.severity.fire.final)

writeRaster(hotspot_rasts, here::here("data/hotspot_rasts.tif"))

# # Let's get polygons for low and high severity fire
# low_sev <- cbi %>% filter(predict.high.severity.fire.final == 1) %>%
#   as.polygons()
#
# high_sev <- cbi %>% filter(predict.high.severity.fire.final == 2) %>%
#   as.polygons()
#
# # get overlap between high severity and hotspots
# hotspot_highsev <- lapply(1:length(names(hotspots_poly)), function(x) intersect(hotspots_poly[[1]], high_sev))



# # get metric layers in the right order
# metrics <- c(metric_rast, fd_rast) %>%
#   select(-ends_with("_quantile_cutoff")) %>%
#   mutate(across(c(breeding_richness, nonbreeding_richness, starts_with("FRic"), starts_with("FEve"), starts_with("FDiv")), .names = "ecoregion_{col}")) %>%
#   #mutate(ecoregion_breeding_richness = breeding_richness, ecoregion_nonbreeding_richness = nonbreeding_richness) %>%
#   # select(breeding_lcbd, nonbreeding_lcbd, breeding_richness, nonbreeding_richness,
#   #        ecoregion_breeding_lcbd, ecoregion_nonbreeding_lcbd) %>%
#   select(names(hotspot_poly))
#
# # list of limits for the color ramp based on the plotted metric
# limits <- map(names(metrics), ~if (stringr::str_detect(.x, "lcbd")) c(0,1) else c(0,160))
#
# # this should be true
# length(names(metrics)) == length(hotspot_poly)
#
# # get a US boundary base map
# US_boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
#   vect() %>%
#   project("epsg:4326") %>%
#   filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
#                      "Montana", "Arizona", "Utah", "Wyoming", "Texas", "Colorado", "New Mexico",
#                      "Kansas", "Oklahoma")) %>%
#   crop(ext(c(-130, -98, 18, 50)))
#
# US_boundary <- US_boundary_states %>% aggregate()
#
# #this sort of works for visualization but the scales are way off
# hotspot_maps <- map(1:length(names(metrics)), ~ ggplot() +
#                       geom_spatvector(data = US_boundary_states, color = "black", fill = "transparent", alpha = 0.5) +
#                       geom_spatraster(data = metrics[[.x]], maxcell = 2500000) +
#                       scale_fill_continuous(type = "viridis", na.value = "transparent", limits = limits[[.x]]) +
#                       geom_spatvector(data = hotspot_poly[[.x]], fill = "#FDE725FF", color = "#FDE725FF", linewidth = 0.2) +
#                       #geom_spatvector(data = highsev_int_poly[[.x]], fill = "black", color = "black") +
#                       theme_void() +
#                       theme(legend.title=element_blank(), plot.title = element_blank()) +
#                       ggtitle(label = names(metrics[[.x]])))
#
# cbi_hotspot_poly <-  map(hotspot_poly, ~crop(.x, cbi))
# cbi_boundary <- as.polygons(cbi) %>% aggregate() %>% fillHoles()
#
# hotspot_overlap_maps <- map(1:length(names(metrics)), ~ ggplot() +
#                       geom_spatvector(data = US_boundary_states, color = "black", fill = "transparent", alpha = 0.5) +
#                       #geom_spatraster(data = cbi, maxcell = 2500000) +
#                       #geom_spatvector(data = cbi_boundary, color = "black",fill = "transparent",) +
#                       #scale_fill_continuous(type = "viridis", na.value = "transparent", limits = limits[[.x]]) +
#                       #geom_spatvector(data = cbi_hotspot_poly[[.x]], fill = "#FDE725FF", color = "#FDE725FF", linewidth = 0.2) +
#                       geom_spatvector(data = highsev_int_poly[[.x]], fill = "red", color = "red") +
#                       theme_void() +
#                       theme(legend.title=element_blank(), plot.title = element_blank()) +
#                       ggtitle(label = names(metrics[[.x]])))
#
# hotspot_plot <- plot_grid(plotlist = hotspot_maps[5:8], nrow = 2)
# hotspot_overlap_plot <- plot_grid(plotlist = hotspot_overlap_maps[5:8], nrow = 2)
# ggsave(here::here("figures/hotspot_maps_big_legend.jpg"), hotspot_plot, width = 10, height = 10, dpi = "retina")
# ggsave(here::here("figures/hotspot_maps_small_legend.jpg"), hotspot_plot, width = 30, height = 30, dpi = "retina")
#

area <- hotspot_rasts  %>%
  expanse(usenames = TRUE) %>%
  rename(metric = layer)

lowsev_int <- hotspot_rasts %>%
  c(., cbi) %>%
  filter(predict.high.severity.fire.final == 1) %>%
  expanse(usenames = TRUE) %>%
  rename(metric = layer, lowsev_area = area)

highsev_int <-  hotspot_rasts %>%
  c(., cbi) %>%
  filter(predict.high.severity.fire.final == 2) %>%
  expanse(usenames = TRUE) %>%
  rename(metric = layer, highsev_area = area)

# why don't the forest hotspot percentages sum to zero? Should be any other area in the CBI...
overlap_df <- full_join(area, lowsev_int) %>%
  full_join(highsev_int) %>%
  mutate(percent_lowsev = (lowsev_area/area)*100, percent_highsev = (highsev_area/area)*100) %>%
  filter(metric != "predict.high.severity.fire.final")

usethis::use_data(overlap_df)


##################################################################
#### test hotspot distributions against landscape expectation ####
##################################################################

### Perform binomial hypothesis test for each tail, where high severity is defined as a "success"
### Alternative hypothesis: there is significantly more high severity (success) than expected by the ecoregion distribution &
### there is significantly less than expected by the ecoregion distribution

# get ratio of high to low severity for each ecoregion to parameterize "true" distribution
ecoregion_names <- unique(ecoregion_shp$ECO_NAME)

# ecoregion_values <- map_dfr(ecoregion_names, function(x) {
#   browser()
#   filter_ecoregion <- ecoregion_shp %>% filter(ECO_NAME == x)
#   ecoregion_cbi <- crop(cbi, filter_ecoregion, mask = TRUE)
#
#   freq(ecoregion_cbi) %>%
#     select(-layer) %>%
#     pivot_wider(names_from = value, values_from = count) %>%
#     mutate(ECO_NAME = x)
#
# }) %>% mutate(p = `2`/(`1` + `2`))

ecoregion_values <- map_dfr(ecoregion_names, function(x) {

  #browser()
  filter_ecoregion <- ecoregion_shp %>% filter(ECO_NAME == x)
  ecoregion_cbi <- extract(cbi, filter_ecoregion)

  ecoregion_cbi %>%
    count(predict.high.severity.fire.final) %>%
    pivot_wider(names_from = predict.high.severity.fire.final, values_from = n) %>%
    mutate(ECO_NAME = x)


}) %>% mutate(p = `2`/(`1` + `2`))

# get number of success (cells == 2) and number of trials (cells %in% c(1,2)) for each
# ecoregion for hotspots of each metric
hotspot_names <- names(hotspot_rasts)

hotspot_values <- map_dfr(hotspot_names, function(x){

  map_dfr(ecoregion_names, hotspot = hotspot_rasts[[x]], function(name, hotspot) {
    #browser()

    metric_name <- names(hotspot)

    filter_ecoregion <- ecoregion_shp %>% filter(ECO_NAME == name)

    ecoregion_hotspot <- crop(hotspot, filter_ecoregion, mask = TRUE)

    ecoregion_cbi <- crop(cbi, filter_ecoregion, mask = TRUE) %>%
                            c(., ecoregion_hotspot) %>%
      filter(!!as.symbol(metric_name) == 1)

    ecoregion_cbi %>%
      select(predict.high.severity.fire.final) %>%
      freq() %>%
      pivot_wider(names_from = value, values_from = count) %>%
      mutate(ECO_NAME = name) %>%
      select(-layer)

  }) %>% mutate(hotspot_type = x)
})

# dataframe of values for binomial test
binomial_df <- hotspot_values %>%
  rowwise() %>%
  mutate(n = sum(`1`,`2`, na.rm = TRUE)) %>%
  # x is the argument for "successes"
  select(hotspot_type, ECO_NAME, x = `2`, n) %>%
  left_join(ecoregion_values %>% select(ECO_NAME,p)) %>%
  mutate(x = ifelse(is.na(x), 0, x))

# perform a test for upper and lower tail
binomial_test <- binomial_df %>%
  group_by(ECO_NAME, hotspot_type) %>%
  mutate(greater.p.value = binom.test(x = x, n = n, p = p, alternative = "greater") %>% .$p.value,
         lower.p.value = binom.test(x = x, n = n, p = p, alternative = "less") %>% .$p.value) %>%
  ungroup() %>%
  tidyr::separate(hotspot_type, into = "type", sep = "_", remove = FALSE) %>%
  mutate(type = ifelse(type == "forest", "forest", "ecoregion"),
         metric = stringr::str_remove(hotspot_type, "forest_")) %>%
  select(-hotspot_type)

# ecoregions that are significant for upper tail test
greater_sig_ecoregions <- binomial_test %>% filter(greater.p.value < 0.05) %>%
  filter(!metric %in% c("breeding_lcbd", "nonbreeding_lcbd")) %>%
  # let's just look at forest hotspots for now
  filter(type == "forest")

# ecoregions that are significant for lower tail test
lower_sig_ecoregions <- binomial_test %>% filter(lower.p.value < 0.05) %>%
  filter(!metric %in% c("breeding_lcbd", "nonbreeding_lcbd")) %>%
  # let's just look at forest hotspots for now
  filter(type == "forest")

sig_ecoregions_df <- greater_sig_ecoregions %>%
  mutate(test_tail = "upper") %>%
  bind_rows(lower_sig_ecoregions %>% mutate(test_tail = "lower"))

# Let's map some stuff!

# get a US boundary base map
US_boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming", "Texas", "Colorado", "New Mexico")) %>%
  crop(ext(c(-130, -104, 18, 50)))

US_boundary <- US_boundary_states %>% aggregate()

map(c("ecoregion_breeding_lcbd", "breeding_richness", "FRic_breeding"), function(metric_name) {
  #browser()

  greater_sig_metrics <- sig_ecoregions_df %>%
  filter(test_tail == "upper", metric == metric_name) %>%
  pull(ECO_NAME)

  lower_sig_metrics <- sig_ecoregions_df %>%
    filter(test_tail == "lower", metric == metric_name) %>%
    pull(ECO_NAME)

indicator_shp <- ecoregion_shp %>%
  mutate(significant = case_when(
    ECO_NAME %in% greater_sig_metrics ~ "more higher severity",
    ECO_NAME %in% lower_sig_metrics ~ "less high severity",
    .default = "nonsignificant"
  )) %>%
  crop(., boundary)

plot <- ggplot() +
  geom_spatvector(mapping = aes(fill = as.factor(significant)), data = indicator_shp, color = "black", linewidth = 0.4) +
  geom_spatvector(data = US_boundary, fill = NA) +
  scale_fill_discrete(type = c("#117733", "#e1ad01", "white"), name = "") +
  theme_map() +
  ggtitle(metric_name)

ggsave(here::here(paste0("figures/", metric_name, "_binomial_test.png")), plot)
})
