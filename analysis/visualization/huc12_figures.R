library(dplyr)
library(terra)
library(tidyterra)
library(ggnewscale)
library(viridis)
library(scales)
library(cowplot)
library(purrr)

###############################
### Biodiversity Metric Map ###
###############################

co_bound <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name == "Utah")

test_rast <- biodiv_zonal_vec %>% crop(co_bound)

ggplot() +
  #geom_spatvector(data = boundary, color = "black", fill = "white") +
  geom_spatvector(data = biodiv_zonal_vec, aes(fill = breeding_richness, color = breeding_richness)) +
  theme_void() +
  viridis::scale_fill_viridis() +
  viridis::scale_color_viridis()

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

breeding_rich_plot <- plot_metric_map(data = test_rast, metric_col = "breeding_richness", legend_title = "Richness")

metric_plot_df <- tibble(metric_col = c("breeding_richness", "breeding_lcbd", "FRic_breeding"), legend_title = c("Richness", "LCBD", "Functional Richness"))

metric_plot_list <- pmap(metric_plot_df, plot_metric_map, data = test_rast)

metric_join_plot <- plot_grid(plot_list = metric_plot_list, nrow = 1, labels = "AUTO", align = 'vh',
                              hjust = -1, axis = "l", scale = 0.9)

save_plot(here::here("figures/huc12_metric_join.jpeg"), metric_join_plot, nrow = 1)

#####################
### Hotspot Types ###
#####################

ecoregion_hotspot_plot <- ggplot() +
  geom_spatvector(data = boundary_states, color = "black", fill = "white") +
  geom_spatvector(data = cbi_forest, fill = "#DDDDDDBF", color = "transparent", alpha = 0.4) +
  geom_spatvector(data = ecoregion_hotspot_zonal_vec %>%
                    filter(breeding_richness == 1),
                  aes(fill = hotspot_type, color = hotspot_type)) +
  theme_void() +
  scale_fill_manual(values = list(`Refugia` = "#82A6B1", `Area of Concern` =  "#BC4749", `Mixed` = "#D5A021"), na.value = "transparent") +
  scale_color_manual(values = list(`Refugia` = "#82A6B1", `Area of Concern` =  "#BC4749", `Mixed` = "#D5A021"), na.value = "transparent") +
  theme(legend.title=element_blank())

ggsave(here::here("figures/huc12_ecoregion_hotspot.jpeg"), ecoregion_hotspot_plot)

# ggplot() +
#   geom_spatvector(data = cbi %>% crop(boundary, mask = TRUE) %>% as.polygons() %>% aggregate(), color = "#DDDDDDBF", alpha = 0.5) +
#   geom_spatvector(data = huc12_cbi_modal_vect %>% filter(predict.high.severity.fire.final %in% c(1,2)), color = "black") +
#   theme_void()



