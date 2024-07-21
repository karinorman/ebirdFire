library(dplyr)
library(terra)
library(tidyterra)
library(purrr)
library(ggplot2)

#boundary <- vect(here::here("data/study_boundary.shp"))

cbi <- rast(here::here("raw_data/predict.high.severity.fire.draft.tif"))
#cbi_tfw <- rast(here::here("raw_data/predict.high.severity.fire.draft.tfw"))

richness_rast <- rast(here::here("data/richness_rast.tif"))
lcbd_rast <- rast(here::here("data/lcbd_rast.tiff"))

# combine and crop raster
metric_rast <- c(richness_rast, lcbd_rast %>% select(breeding_lcbd, nonbreeding_lcbd, ecoregion_breeding_lcbd, ecoregion_nonbreeding_lcbd)) %>% crop(cbi, mask = TRUE)

# # let's visualize it
# ggplot() +
#   geom_spatraster(data = metric_rast) +
#   scale_fill_continuous(type = "viridis", na.value = "transparent") +
#   facet_wrap(~lyr)

# get raster of hotspots for non-ecoregion layers
hotspot_rast <- metric_rast %>%
  select(breeding_lcbd, nonbreeding_lcbd) %>%
  mutate(across(everything(),
                              ~ifelse(.x > quantile(.x, probs = 0.95, na.rm = TRUE), "a", NA)))

# do it for ecoregions which already have quantile calculated for each group
ecoregion_hotspot_rast <- lcbd_rast %>%
  crop(cbi, mask = TRUE) %>%
  select(-c(breeding_lcbd, nonbreeding_lcbd)) %>%
  mutate(ecoregion_breeding_lcbd = ifelse(ecoregion_breeding_lcbd > breeding_quantile_cutoff, "a", NA),
         ecoregion_nonbreeding_lcbd = ifelse(ecoregion_nonbreeding_lcbd > nonbreeding_quantile_cutoff, "a", NA)) %>%
  select(-c(breeding_quantile_cutoff, nonbreeding_quantile_cutoff))

#hotspot_rast <- c(hotspot_rast, ecoregion_hotspot_rast)

## danger! this depends on layer order! ##
hotspot_poly <- c(map(1:4, ~as.polygons(hotspot_rast[[.x]])),
                  map(1:2, ~as.polygons(ecoregion_hotspot_rast[[.x]])))

# this should be true
length(names(metric_rast)) == length(hotspot_poly)

#this sort of works for visualization but the scales are way off
hotspot_maps <- map(1:length(names(metric_rast)), ~ ggplot() +
      geom_spatraster(data = metric_rast[[.x]]) +
      scale_fill_continuous(type = "viridis", na.value = "transparent") +
      geom_spatvector(data = hotspot_poly[[.x]], color = "black", fill = "#FDE725FF") +
      theme_void() +
  ggtitle(label = names(metric_rast[[.x]])))

hotspot_plot <- plot_grid(plotlist = hotspot_maps, nrow = 3)
ggsave(here::here("figures/hotspot_maps.jpg"), hotspot_plot, width = 10, height = 10)

# Let's get polygons for low and high severity fire
low_sev <- cbi %>% filter(predict.high.severity.fire.draft == 1) %>%
  as.polygons()

high_sev <- cbi %>% filter(predict.high.severity.fire.draft == 2) %>%
  as.polygons()

# Now let's see the overlap with different hotspots
layer_names <- c(names(hotspot_rast), names(ecoregion_hotspot_rast))

area <- map(hotspot_poly, expanse) %>%
  set_names(layer_names) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "hotspot_area")

lowsev_int <- map(hotspot_poly, ~ intersect(.x, low_sev) %>% expanse()) %>%
  set_names(layer_names) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "lowsev_intersect")

highsev_int <- map(hotspot_poly, ~ intersect(.x, high_sev) %>% expanse()) %>%
  set_names(layer_names) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "highsev_intersect")

overlap_df <- full_join(area, lowsev_int) %>%
  full_join(highsev_int) %>%
  mutate(percent_lowsev = (lowsev_intersect/hotspot_area)*100, percent_highsev = (highsev_intersect/hotspot_area)*100)
