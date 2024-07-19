library(dplyr)
library(terra)
library(tidyterra)
library(ggnewscale)
library(viridis)
library(scales)
library(cowplot)

# testing things out
#
# breeding_97 <-quantile(values(test$breeding_lcbd), probs = 0.975, na.rm = TRUE)
#
# test <- lcbd_rast %>%
#   mutate(breeding_hotspot = ifelse(breeding_lcbd > breeding_97, "a", NA)) %>%
#   select(breeding_lcbd, breeding_hotspot)
#
# hotspot_poly <- as.polygons(test$breeding_hotspot)
#
# ggplot() +
#   geom_spatraster(data = test$breeding_lcbd) +
#   #scale_fill_viridis_c(limits = c(0.2, 1), oob = scales::squish, na.value = "transparent") +
#   scale_fill_continuous(type = "viridis", na.value = "transparent") +
#   geom_spatvector(data = hotspot_poly, color = "black", fill = "#FDE725FF") +
#   theme_void()
#   #scale_fill_discrete("red")
#   #scale_fill_manual(values = "red", na.value = "transparent")


# maps of metrics for Western US
richness_rast <- rast(here::here("data/richness_rast.tif"))
lcbd_rast <- rast(here::here("data/lcbd_rast.tiff"))

# let's visualize it
richness_plot <- ggplot() +
  geom_spatraster(data = richness_rast) +
  scale_fill_continuous(type = "viridis", na.value = "transparent", name = "species richness") +
  #geom_spatvector(data = boundary, fill = "transparent") +
  facet_wrap(~lyr, labeller = as_labeller(c(`breeding_richness` = "breeding",
                                            `nonbreeding_richness` = "nonbreeding"))) +
  theme_void()

lcbd_plot <- ggplot() +
  geom_spatraster(data = lcbd_rast %>% select(-c("breeding_quantile_cutoff", "nonbreeding_quantile_cutoff"))) +
  scale_fill_continuous(type = "viridis", na.value = "transparent", name = "lcbd") +
  #geom_spatvector(data = boundary, fill = "transparent") +
  facet_wrap(~lyr, labeller = as_labeller(c(`breeding_lcbd` = "",
                                            `nonbreeding_lcbd` = ""))) +
  theme_void()

wus_metric_plot <- plot_grid(richness_plot, lcbd_plot, nrow = 2)

ggsave(here::here("figures/wus_metrics_map.jpg"), wus_metric_plot, width = 10, height = 10)
