library(dplyr)
library(terra)
library(tidyterra)
library(ggnewscale)
library(viridis)
library(scales)
library(cowplot)
library(purrr)
library(ggplot2)
library(patchwork)

# read in CBI raster
cbi <- rast(here::here("data/cbi.tif"))
cbi_forest <- cbi %>%
  filter(predict.high.severity.fire.final %in% c(1,2)) %>%
  as.polygons() %>%
  aggregate()

# get basemap for plotting
boundary <- vect(here::here("data/study_boundary.shp"))
boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming", "Colorado", "New Mexico", "Texas")) %>%
  crop(boundary)

# shapefiles of biodiversity metrics and hotspots
biodiv_zonal_vec <- vect(here::here("data/biodiv_zonal.shp"))
names(biodiv_zonal_vec) <- c("huc12", "breeding_lcbd", "nonbreeding_lcbd", "ecoregion_breeding_lcbd", "ecoregion_nonbreeding_lcbd", "breeding_richness",
                             "nonbreeding_richness", "FRic_breeding", "FEve_breeding", "FDiv_breeding", "FRic_nonbreeding", "FEve_nonbreeding", "FDiv_nonbreeding",
                             "ECO_NAME", "fire", "ID", "low_sev", "high_sev", "unforested", "no_data", "max_severity", "forest_total")

ecoregion_hotspot_zonal_vec <- vect(here::here("data/ecoregion_hotspots_huc12.shp"))
names(ecoregion_hotspot_zonal_vec) <- c("huc12", "breeding_lcbd", "nonbreeding_lcbd", "ecoregion_breeding_lcbd", "ecoregion_nonbreeding_lcbd", "breeding_richness",
                                        "nonbreeding_richness", "FRic_breeding", "FEve_breeding", "FDiv_breeding", "FRic_nonbreeding", "FEve_nonbreeding", "FDiv_nonbreeding",
                                        "ECO_NAME", "fire", "ID", "low_sev", "high_sev", "unforested", "no_data", "max_severity", "forest_total", "hotspot_type")
###############################
### Biodiversity Metric Map ###
###############################

# co_bound <- rnaturalearth::ne_states(iso_a2 = "US") %>%
#   vect() %>%
#   project("epsg:4326") %>%
#   filter(name == "Utah")
#
# test_rast <- biodiv_zonal_vec %>% crop(co_bound)

# ggplot() +
#   #geom_spatvector(data = boundary, color = "black", fill = "white") +
#   geom_spatvector(data = biodiv_zonal_vec, aes(fill = breeding_richness, color = breeding_richness)) +
#   theme_void() +
#   viridis::scale_fill_viridis() +
#   viridis::scale_color_viridis()

plot_metric_map <- function(data, metric_col, legend_title){

  col = sym(metric_col)

  ggplot() +
    geom_spatvector(data = data, aes(fill = !!col, color = !!col)) +
    theme_void() +
    viridis::scale_fill_viridis(name = legend_title) +
    viridis::scale_color_viridis(guide = "none") +
    theme(panel.background = element_rect(fill = "transparent",
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


metric_plot_df <- tibble(metric_col = c("breeding_richness", "ecoregion_breeding_lcbd", "FRic_breeding"), legend_title = c("Richness", "LCBD", "Functional \nRichness"))

metric_plot_list <- pmap(metric_plot_df, plot_metric_map, data = biodiv_zonal_vec)

metric_join_plot <- wrap_plots(metric_plot_list, nrow = 1) + plot_annotation(tag_levels = "A")

ggsave(here::here("figures/huc12_metric_join.jpeg"), metric_join_plot)

#####################
### Hotspot Types ###
#####################

plot_hotspot_types <- function(data, metric_col){

  col = sym(metric_col)

  ggplot() +
    geom_spatvector(data = boundary_states, color = "black", fill = "white") +
    geom_spatvector(data = cbi_forest, fill = "#DDDDDDBF", color = "transparent", alpha = 0.4) +
    geom_spatvector(data = data %>%
                      filter(!!sym(metric_col) == 1),
                    aes(fill = hotspot_type, color = hotspot_type)) +
    theme_void() +
    scale_fill_manual(values = list(`Refugia` = "#82A6B1", `Area of Concern` =  "#BC4749", `Mixed` = "#D5A021"), na.value = "transparent") +
    scale_color_manual(values = list(`Refugia` = "#82A6B1", `Area of Concern` =  "#BC4749", `Mixed` = "#D5A021"), na.value = "transparent") +
    theme(legend.title=element_blank(), panel.background = element_rect(fill = "transparent",
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


hotspot_plot_list <- map(c("breeding_richness", "ecoregion_breeding_lcbd", "FRic_breeding"), plot_hotspot_types, data = ecoregion_hotspot_zonal_vec)

hotspot_type_figure <- patchwork::wrap_plots(hotspot_plot_list, nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(here::here("figures/huc12_hotspot_join.jpeg"), hotspot_type_figure)

