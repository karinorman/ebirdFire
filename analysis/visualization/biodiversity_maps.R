library(dplyr)
library(ggplot2)
library(terra)
library(tidyterra)
library(ggnewscale)
library(viridis)
library(scales)
library(cowplot)
library(purrr)
library(ggridges)

# maps of metrics for Western US
# richness_rast <- rast(here::here("data/richness_rast.tif"))
# lcbd_rast <- rast(here::here("data/lcbd_rast.tiff"))


boundary <- terra::vect(here::here("data/study_boundary.shp"))

metric_rast <- rast(here::here("data/metric_rast.tiff")) %>%
  crop(boundary, mask = TRUE)

fd_rast <- rast(here::here("data/fd_rast.tiff")) %>%
  crop(boundary, mask = TRUE)

biodiv_rast <- c(metric_rast, fd_rast)

cbi <- rast(here::here("data/cbi.tif")) %>%
  resample(., biodiv_rast, "near")

ecoregions <- terra::vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326")

# let's visualize it
# richness_plot <- ggplot() +
#   geom_spatraster(data = metric_rast %>% select(breeding_richness, nonbreeding_richness)) +
#   scale_fill_continuous(type = "viridis", na.value = "transparent", name = "species richness") +
#   #geom_spatvector(data = boundary, fill = "transparent") +
#   facet_wrap(~lyr, labeller = as_labeller(c(`breeding_richness` = "breeding",
#                                             `nonbreeding_richness` = "nonbreeding"))) +
#   theme_void()
#
# lcbd_plot <- ggplot() +
#   geom_spatraster(data = metric_rast %>% select(breeding_lcbd, nonbreeding_lcbd)) +
#   scale_fill_continuous(type = "viridis", na.value = "transparent", name = "lcbd") +
#   #geom_spatvector(data = boundary, fill = "transparent") +
#   facet_wrap(~lyr, labeller = as_labeller(c(`breeding_lcbd` = "",
#                                             `nonbreeding_lcbd` = ""))) +
#   theme_void()
#
# wus_metric_plot <- plot_grid(richness_plot, lcbd_plot, nrow = 2)
#
# ggsave(here::here("figures/wus_metrics_map.jpg"), wus_metric_plot, width = 10, height = 10)

#################################
######### metric map ############
#################################

plot_metric_map <- function(data, metric_col){

  col = sym(metric_col)
  ggplot() +
  geom_spatraster(data = data %>% select(!!col),
                  maxcell = 1000e+05) +
  scale_fill_continuous(type = "viridis", na.value = "transparent") +
  geom_spatvector(data = boundary, fill = "transparent", colour = "black") +
  theme_void() +
  theme(legend.title=element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = NA, color = NA),
        legend.box.background = element_rect(fill = NA, color = NA),
        legend.key = element_rect(fill = "transparent"),
        legend.box = element_blank()
        )
}

breeding_richness_plot <- plot_metric_map(metric_rast, "breeding_richness")
ggsave(here::here("figures/wus_bredding_richness_map.png"), breeding_richness_plot, bg = "transparent")

breeding_lcbd_plot <- plot_metric_map(metric_rast, "breeding_lcbd")
ggsave(here::here("figures/wus_breeding_lcbd_map.png"), breeding_lcbd_plot, bg = "transparent")

ecoreg_breeding_lcbd_plot <- plot_metric_map(metric_rast, "ecoregion_breeding_lcbd")
ggsave(here::here("figures/wus_ecoreg_breedgin_lcbd_map.png"), ecoreg_breeding_lcbd_plot, bg = "transparent")

breeding_fric_plot <- plot_metric_map(fd_rast, "FRic_breeding")
ggsave(here::here("figures/wus_breeding_fric_map.png"), breeding_fric_plot, bg = "transparent")


#################################
########## CBI map ##############
#################################

# Plot of just CBI
cbi_plot <- ggplot() +
  geom_spatvector(data = boundary, color = "black", fill = "white") +
  #geom_spatvector(data = ecoregions %>% crop(boundary) %>% terra::aggregate(), color = "black", fill = "white", linewidth = 0.4) +
  geom_spatraster(data = cbi %>%
                    filter(predict.high.severity.fire.final %in% c(1,2)) %>%
                    mutate(predict.high.severity.fire.final = as.factor(predict.high.severity.fire.final))) +
  scale_fill_manual(values = list(`1` = "#D3B073", `2` =  "#BC4749"), na.value = "transparent", na.translate = FALSE,
                    labels = c("Low Severity", "High Severity")) +
  theme_void() +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = NA, color = NA),
        legend.box.background = element_rect(fill = NA, color = NA),
        legend.key = element_rect(fill = "transparent"),
        legend.box = element_blank()
  )

ggsave(here::here("figures/cbi.png"), cbi_plot, width = 10, height = 10, bg = "transparent")

#################################
####### hotspot maps ###########
#################################

hotspot_rast <- rast(here::here("data/hotspot_rasts.tif"))

plot_cbi <- cbi %>%
  mutate(predict.high.severity.fire.final =
           as.factor(ifelse(predict.high.severity.fire.final == 0, NA, predict.high.severity.fire.final)))

plot_refugia <- function(metric_col, boundary_shp, mask = FALSE){

  plot <- ggplot() +
    geom_spatvector(data = boundary_shp, color = "black", fill = "white") +
    theme_void() +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "transparent",
                                          colour = NA_character_), # necessary to avoid drawing panel outline
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          plot.background = element_rect(fill = "transparent",
                                         colour = NA_character_), # necessary to avoid drawing plot outline
          legend.background = element_rect(fill = NA, color = NA),
          legend.box.background = element_rect(fill = NA, color = NA),
          legend.key = element_rect(fill = "transparent"),
          legend.box = element_blank()
    )

  if (mask == TRUE){
    plot <- plot +
      geom_spatvector(data = cbi %>% crop(boundary_shp, mask = TRUE) %>% as.polygons() %>% aggregate(), color = "#DDDDDDBF", alpha = 0.5) +
      geom_spatraster(data = plot_cbi %>% crop(hotspot_rast[[metric_col]], mask = TRUE) %>% crop(boundary_shp, mask = TRUE)) +
      scale_fill_manual(values = list(`1` = "#82A6B1", `2` =  "#BC4749"), na.value = "#DDDDDDBF", na.translate = FALSE,
                        labels = c("Refugia", "Areas of Concern"))
  } else{
    plot <- plot +
      geom_spatvector(data = as.polygons(cbi) %>% aggregate(), color = "#DDDDDDBF", alpha = 0.5) +
      geom_spatraster(data = plot_cbi %>% crop(hotspot_rast[[metric_col]], mask = TRUE)) +
      scale_fill_manual(values = list(`1` = "#82A6B1", `2` =  "#BC4749"), na.value = "#DDDDDDBF", na.translate = FALSE,
                        labels = c("Refugia", "Areas of Concern"))
  }
}

breeding_lcbd_refugia <- plot_refugia("forest_ecoregion_breeding_lcbd", boundary_shp = boundary)
ggsave(here::here("figures/refugia_eco_breeding_lcbd.png"), breeding_lcbd_refugia)

breeding_richness_refugia <- plot_refugia("forest_breeding_richness", boundary_shp = boundary)
ggsave(here::here("figures/refugia_breeding_richness.png"), breeding_richness_refugia)

breeding_fric_refugia <- plot_refugia("forest_FRic_breeding", boundary_shp = boundary)
ggsave(here::here("figures/refugia_breeding_fric.png"), breeding_fric_refugia)

### Zoom in on Arizona
arizona <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name == "Arizona")

montana <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name == "Montana")

arizona_breeding_lcbd_refugia <- plot_refugia("forest_ecoregion_breeding_lcbd", boundary_shp = arizona, mask = TRUE)
ggsave(here::here("figures/refugia_pullouts/arizona_refugia_eco_breeding_lcbd.png"), arizona_breeding_lcbd_refugia)

arizona_breeding_richness_refugia <- plot_refugia("forest_breeding_richness", boundary_shp = arizona, mask = TRUE)
ggsave(here::here("figures/refugia_pullouts/arizona_refugia_breeding_richness.png"), arizona_breeding_richness_refugia)

arizona_breeding_fric_refugia <- plot_refugia("forest_FRic_breeding", boundary_shp = arizona, mask = TRUE)
ggsave(here::here("figures/refugia_pullouts/arizona_refugia_breeding_fric.png"), arizona_breeding_fric_refugia)

montana_richness_zoom <- plot_refugia("forest_breeding_richness", boundary_shp = montana, mask = TRUE)
ggsave(here::here("figures/refugia_pullouts/montana_refugia_breeding_richness.png"), montana_richness_zoom, bg = "#DDDDDDBF")

#################################
### metric + hotspot maps #######
#################################

plot_metric_hotspot_map <- function(data, metric_col){

  col = sym(metric_col)

  poly_name = paste0("forest_", metric_col)

  ggplot() +
    geom_spatraster(data = data %>% select(!!col) %>% crop(., ecoregions, mask = TRUE),
                    maxcell = 1000e+05) +
    scale_fill_continuous(type = "viridis", na.value = "transparent") +
    geom_spatvector(data = boundary, fill = "transparent", colour = "black") +
    #geom_spatvector(data = ecoregions, fill = "transparent", colour = "white") +
    geom_spatvector(data = as.polygons(hotspot_rast[[poly_name]]), fill = "#FDE725FF", color = "#FDE725FF") +
    theme_void() +
    theme(legend.title=element_blank(),
          panel.background = element_rect(fill = "transparent",
                                          colour = NA_character_), # necessary to avoid drawing panel outline
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          plot.background = element_rect(fill = "transparent",
                                         colour = NA_character_), # necessary to avoid drawing plot outline
          legend.background = element_rect(fill = NA, color = NA),
          legend.box.background = element_rect(fill = NA, color = NA),
          legend.key = element_rect(fill = "transparent"),
          legend.box = element_blank()
    )
}

breeding_richness_hotspot <- plot_metric_hotspot_map(metric_rast, "breeding_richness")
ggsave(here::here("figures/breeding_richness_hotspot.png"), breeding_richness_hotspot)

breeding_lcbd_hotspot <- plot_metric_hotspot_map(metric_rast, "ecoregion_breeding_lcbd")
ggsave(here::here("figures/breeding_lcbd_hotspot.png"), breeding_lcbd_hotspot)

breeding_fric_hotspot <- plot_metric_hotspot_map(fd_rast, "FRic_breeding")
ggsave(here::here("figures/breeding_fric_hotspot.png"), breeding_fric_hotspot)
