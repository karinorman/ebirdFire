library(dplyr)
library(terra)
library(tidyterra)
library(purrr)
library(ggplot2)
library(cowplot)

cbi <- rast(here::here("data/cbi.tif"))

forest_poly <- cbi %>% filter(predict.high.severity.fire.final %in% c(1,2)) %>%
  as.polygons() #%>%
  # aggregate() %>%
  # fillHoles()

# get ecoregions as raster
ecoregions <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  select("ECO_NAME") %>%
  terra::rasterize(., cbi, field = "ECO_NAME")# %>%
  #resample(cbi, method = "near")

metric_rast <- rast(here::here("data/metric_rast.tiff")) %>%
  crop(ecoregions, mask = TRUE) %>%
  c(., ecoregions, rast(here::here("data/fd_rast.tiff"))%>%
      crop(ecoregions, mask = TRUE)) # %>%
  #crop(cbi, mask = TRUE)

ecoregion_hotspots <- metric_rast %>%
  terra::as.data.frame(., xy = TRUE) %>%
  group_by(ECO_NAME) %>%
  mutate(across(-c(x,y),
                ~ifelse(.x > quantile(.x, probs = 0.95, na.rm = TRUE), 1, NA))) %>%
  ungroup() %>%
  select(-ECO_NAME) %>%
  rast(.,  type="xyz", crs = "epsg:4326")

forest_hotspots <- metric_rast %>%
  crop(forest_poly, mask = TRUE) %>%
  terra::as.data.frame(., xy = TRUE) %>%
  group_by(ECO_NAME) %>%
  mutate(across(-c(x,y),
                ~ifelse(.x > quantile(.x, probs = 0.95, na.rm = TRUE), 1, NA))) %>%
  ungroup() %>%
  select(-ECO_NAME) %>%
  rast(.,  type="xyz", crs = "epsg:4326") %>%
  rename_with(~paste0("forest_", .x))

hotspots_poly <- c(lapply(1:length(names(ecoregion_hotspots)), function(x) as.polygons(ecoregion_hotspots[[x]])),
                      lapply(1:length(names(forest_hotspots)), function(x) as.polygons(forest_hotspots[[x]])))

names(hotspots_poly) <- c(names(ecoregion_hotspots), names(forest_hotspots))

# Let's get polygons for low and high severity fire
low_sev <- cbi %>% filter(predict.high.severity.fire.final == 1) %>%
  as.polygons()

high_sev <- cbi %>% filter(predict.high.severity.fire.final == 2) %>%
  as.polygons()

# get overlap between high severity and hotspots
hotspot_highsev <- lapply(1:length(names(hotspots_poly)), function(x) intersect(hotspots_poly[[1]], high_sev))



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

area <-  map(hotspots_poly, ~expanse(.x)) %>%
  #map(., expanse) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "hotspot_area")

lowsev_int <- map(hotspots_poly, ~ terra::intersect(.x, low_sev) %>% expanse()) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "lowsev_intersect")

highsev_int <- map(hotspots_poly, ~ intersect(.x, high_sev) %>% expanse()) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "highsev_intersect")

# why don't the forest hotspot percentages sum to zero? Should be any other area in the CBI...
overlap_df <- full_join(area, lowsev_int) %>%
  full_join(highsev_int) %>%
  mutate(percent_lowsev = (lowsev_intersect/hotspot_area)*100, percent_highsev = (highsev_intersect/hotspot_area)*100)

usethis::use_data(overlap_df)
