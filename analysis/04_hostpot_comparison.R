library(dplyr)
library(terra)
library(tidyterra)
library(purrr)
library(ggplot2)
library(cowplot)

#boundary <- vect(here::here("data/study_boundary.shp"))

cbi <- rast(here::here("raw_data/predict.high.severity.fire.draft.tif"))
#cbi_tfw <- rast(here::here("raw_data/predict.high.severity.fire.draft.tfw"))

metric_rast <- rast(here::here("data/metric_rast.tiff"))
fd_rast <- rast(here::here("data/fd_rast.tiff"))

# get raster of hotspots for non-ecoregion layers
hotspot_rast <- c(metric_rast, fd_rast) %>%
  select(breeding_lcbd, nonbreeding_lcbd, breeding_richness, nonbreeding_richness,
         FRic_breeding, FEve_breeding, FDiv_breeding,
         FRic_nonbreeding, FEve_nonbreeding, FDiv_nonbreeding) %>%
  mutate(across(everything(),
                              ~ifelse(.x > quantile(.x, probs = 0.95, na.rm = TRUE), "a", NA)))

# do it for ecoregions which already have quantile calculated for each group
ecoregion_hotspot_rast <- metric_rast %>%
  select(-c(breeding_lcbd, nonbreeding_lcbd)) %>%
  mutate(ecoregion_breeding_lcbd = ifelse(ecoregion_breeding_lcbd > lcbd_breeding_quantile_cutoff, "a", NA),
         ecoregion_nonbreeding_lcbd = ifelse(ecoregion_nonbreeding_lcbd > lcbd_nonbreeding_quantile_cutoff, "a", NA),
         ecoregion_breeding_richness = ifelse(breeding_richness > breeding_quantile_cutoff, "a", NA),
         ecoregion_nonbreeding_richness = ifelse(nonbreeding_richness > nonbreeding_quantile_cutoff, "a", NA)) %>%
  select(-ends_with("cutoff"), -c(breeding_richness, nonbreeding_richness))

# Ecoregion quantiles for fd metrics
ecoregion_fd_hotspot_rast <- fd_rast %>%
  mutate(across(-starts_with("quantile"), ~ifelse(.x > get(paste0("quantile_", cur_column())), "a", NA))) %>%
  select(-starts_with("quantile_")) %>%
  rename_with(~paste0("ecoregion_", .x))

hotspot_poly <- c(map(1:length(names(hotspot_rast)), ~as.polygons(hotspot_rast[[.x]])),
                  map(1:length(names(ecoregion_hotspot_rast)), ~as.polygons(ecoregion_hotspot_rast[[.x]])),
                  map(1:length(names(ecoregion_fd_hotspot_rast)), ~as.polygons(ecoregion_fd_hotspot_rast[[.x]])))

names(hotspot_poly) <- c(names(hotspot_rast), names(ecoregion_hotspot_rast), names(ecoregion_fd_hotspot_rast))

# Let's get polygons for low and high severity fire
low_sev <- cbi %>% filter(predict.high.severity.fire.draft == 1) %>%
  as.polygons()

high_sev <- cbi %>% filter(predict.high.severity.fire.draft == 2) %>%
  as.polygons()

# get overlap between high severity and hotspots
highsev_int_poly <- map(hotspot_poly, ~ intersect(.x, high_sev))

# get metric layers in the right order
metrics <- c(metric_rast, fd_rast) %>%
  select(-ends_with("_quantile_cutoff")) %>%
  mutate(across(c(breeding_richness, nonbreeding_richness, starts_with("FRic"), starts_with("FEve"), starts_with("FDiv")), .names = "ecoregion_{col}")) %>%
  #mutate(ecoregion_breeding_richness = breeding_richness, ecoregion_nonbreeding_richness = nonbreeding_richness) %>%
  # select(breeding_lcbd, nonbreeding_lcbd, breeding_richness, nonbreeding_richness,
  #        ecoregion_breeding_lcbd, ecoregion_nonbreeding_lcbd) %>%
  select(names(hotspot_poly))

# list of limits for the color ramp based on the plotted metric
limits <- map(names(metrics), ~if (stringr::str_detect(.x, "lcbd")) c(0,1) else c(0,160))

# this should be true
length(names(metrics)) == length(hotspot_poly)

# get a US boundary base map
US_boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming", "Texas", "Colorado", "New Mexico",
                     "North Dakota", "South Dakota", "Kansas", "Nebraska", "Oklahoma")) %>%
  crop(ext(c(-130, -98, 18, 50)))

US_boundary <- US_boundary_states %>% aggregate()

#this sort of works for visualization but the scales are way off
hotspot_maps <- map(1:length(names(metrics)), ~ ggplot() +
                      geom_spatvector(data = US_boundary_states, color = "black", fill = "transparent", alpha = 0.5) +
                      geom_spatraster(data = metrics[[.x]], maxcell = 2500000) +
                      scale_fill_continuous(type = "viridis", na.value = "transparent", limits = limits[[.x]]) +
                      geom_spatvector(data = hotspot_poly[[.x]], fill = "#FDE725FF", color = "#FDE725FF", linewidth = 0.2) +
                      #geom_spatvector(data = highsev_int_poly[[.x]], fill = "black", color = "black") +
                      theme_void() +
                      theme(legend.title=element_blank(), plot.title = element_blank()) +
                      ggtitle(label = names(metrics[[.x]])))

cbi_hotspot_poly <-  map(hotspot_poly, ~crop(.x, cbi))
cbi_boundary <- as.polygons(cbi) %>% aggregate() %>% fillHoles()

hotspot_overlap_maps <- map(1:length(names(metrics)), ~ ggplot() +
                      geom_spatvector(data = US_boundary_states, color = "black", fill = "transparent", alpha = 0.5) +
                      #geom_spatraster(data = cbi, maxcell = 2500000) +
                      #geom_spatvector(data = cbi_boundary, color = "black",fill = "transparent",) +
                      #scale_fill_continuous(type = "viridis", na.value = "transparent", limits = limits[[.x]]) +
                      #geom_spatvector(data = cbi_hotspot_poly[[.x]], fill = "#FDE725FF", color = "#FDE725FF", linewidth = 0.2) +
                      geom_spatvector(data = highsev_int_poly[[.x]], fill = "red", color = "red") +
                      theme_void() +
                      theme(legend.title=element_blank(), plot.title = element_blank()) +
                      ggtitle(label = names(metrics[[.x]])))

hotspot_plot <- plot_grid(plotlist = hotspot_maps[5:8], nrow = 2)
hotspot_overlap_plot <- plot_grid(plotlist = hotspot_overlap_maps[5:8], nrow = 2)
ggsave(here::here("figures/hotspot_maps_big_legend.jpg"), hotspot_plot, width = 10, height = 10, dpi = "retina")
ggsave(here::here("figures/hotspot_maps_small_legend.jpg"), hotspot_plot, width = 30, height = 30, dpi = "retina")


# we want the area relative to that in the CBI extent
# this is NOT RIGHT for the non-ecoregion hotspots, since those quantiles were calculated relative
# to the whole study extent, not just CBI, so only use the ecoregion overlap!
area <-  map(hotspot_poly, ~crop(.x, cbi)) %>%
  map(., expanse) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "hotspot_area")

lowsev_int <- map(hotspot_poly, ~ terra::intersect(.x, low_sev) %>% expanse()) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "lowsev_intersect")

highsev_int <- map(hotspot_poly, ~ intersect(.x, high_sev) %>% expanse()) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "highsev_intersect")

overlap_df <- full_join(area, lowsev_int) %>%
  full_join(highsev_int) %>%
  mutate(percent_lowsev = (lowsev_intersect/hotspot_area)*100, percent_highsev = (highsev_intersect/hotspot_area)*100)

usethis::use_data(overlap_df)
