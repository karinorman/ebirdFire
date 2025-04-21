library(dplyr)
library(terra)
library(tidyterra)

boundary <- terra::vect(here::here("data/study_boundary.shp")) %>% project("epsg:5070")
#cbi <- rast(here::here("data/cbi.tif"))

regime_group <- rast(here::here("raw_data/US_105_FRG/Tif/us_105frg.tif")) %>%
  crop(., boundary)

## This exhuasts vector memory
# severity_types <- regime_group %>%
#   mutate(severity = case_when(
#     Value_1 %in% c(1,3) ~ "low",
#     Value_1 %in% c(2,4) ~ "high",
#     .default = NA
#   ))

ecoregion_hotspot_zonal_vec <- vect(here::here("data/ecoregion_hotspots_huc12.shp"))
names(ecoregion_hotspot_zonal_vec) <- c("huc12", "breeding_lcbd", "nonbreeding_lcbd", "ecoregion_breeding_lcbd", "ecoregion_nonbreeding_lcbd", "breeding_richness",
                                        "nonbreeding_richness", "FRic_breeding", "FEve_breeding", "FDiv_breeding", "FRic_nonbreeding", "FEve_nonbreeding", "FDiv_nonbreeding",
                                        "ECO_NAME", "fire", "ID", "low_sev", "high_sev", "unforested", "no_data", "max_severity", "forest_total", "hotspot_type")

high_sev_huc12 <- ecoregion_hotspot_zonal_vec %>%
  filter(max_severity == "high_sev") %>%
  select(huc12, max_severity)

high_sev_regime <- extract(regime_group, high_sev_huc12 %>% project("epsg:5070"), fun = table) %>%
  bind_cols(huc12 = high_sev_huc12$huc12) %>%
  select(-ID)

# save out because it takes a long time to run
usethis::use_data(high_sev_regime)

hist_regime <- high_sev_regime %>%
  mutate(total_pixel = rowSums(across(-huc12))) %>%
  rowwise() %>%
  mutate(low_sev = sum(`1`, `3`),
         high_sev = sum(`2`, `4`)) %>%
  select(huc12, total_pixel, low_sev, high_sev) %>%
  mutate(across(-c(huc12, total_pixel), ~.x/total_pixel)) %>%
  left_join(as.data.frame(high_sev_huc12)) %>%
  mutate(historical_type = case_when(
    abs(low_sev - high_sev) < 0.1 ~ "mixed",
    low_sev > high_sev ~ "low_sev",
    high_sev > low_sev ~ "high_sev"
  )) %>% left_join(ecoregion_hotspot_zonal_vec %>%
                     as.data.frame() %>%
                     select(huc12, ecoregion_breeding_lcbd, breeding_richness, FRic_breeding)) %>%
  filter(if_any(c(ecoregion_breeding_lcbd, breeding_richness, FRic_breeding), ~!is.na(.x)))

hotspot_dist <- hist_regime %>%
  select(historical_type, ecoregion_breeding_lcbd, breeding_richness, FRic_breeding) %>%
  group_by(historical_type) %>%
  summarize(across(c(ecoregion_breeding_lcbd, breeding_richness, FRic_breeding), ~sum(!is.na(.x)))) %>%
  pivot_longer(-historical_type, names_to = "metric", values_to = "count") %>%
  group_by(metric) %>%
  mutate(percent = count/sum(count))

load(here::here("data/hotspot_cat.rda"))

pal_sev <- list(`Refugia` = "#82A6B1",  `Mixed` = "#D5A021", `Area of Concern` =  "#BC4749")
pal_hist <- list(`low_sev` = "#82A6B1",  `mixed` = "#D5A021", `high_sev` =  "#BC4749")

sev_bar <- ggplot(hotspot_cat %>% mutate(hotspot_type = factor(hotspot_type, levels = c("Refugia", "Mixed", "Area of Concern")),
                                         metric = factor(metric, levels = c("FRic_breeding",  "ecoregion_breeding_lcbd", "breeding_richness"))),
                  aes(fill=hotspot_type, y=percent, x=metric)) +
  geom_bar(position="dodge", stat="identity", color = "white") +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=10),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 10)) +
  scale_fill_manual(values = pal_sev,
                    labels = list(`Refugia` = "low severity",  `Mixed` = "mixed", `Area of Concern` =  "high severity"), guide = "none") +
  ylab("percent of hotspots") +
  xlab("") +
  scale_x_discrete(labels = list("breeding_richness" = "Richness", "ecoregion_breeding_lcbd" = "LCBD", "FRic_breeding" =  "Functional Richness")) +
  coord_flip() +
  ggtitle("Predicted")

hist_sev_bar <- ggplot(hotspot_dist %>% mutate(historical_type = factor(historical_type, levels = c("low_sev", "mixed", "high_sev")),
                                               metric = factor(metric, levels = c("FRic_breeding",  "ecoregion_breeding_lcbd", "breeding_richness"))),
                       aes(fill=historical_type, y=percent, x=metric)) +
  geom_bar(position="dodge", stat="identity", width = 0.5, color = "white") +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 10)) +
  scale_fill_manual(values = pal_hist, labels = list("mixed" = "mixed", "high_sev" = "high severity", "low_sev" = "low severity")) +
  ylab("percent of hotspots") +
  xlab("") +
  scale_x_discrete(labels = list("breeding_richness" = "Richness", "ecoregion_breeding_lcbd" = "LCBD", "FRic_breeding" =  "Functional Richness")) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) +
  ggtitle("Historical")


hotspot_dist_join <- cowplot::plot_grid(sev_bar, hist_sev_bar, nrow = 1)

cowplot::save_plot(here::here("figures/hotspot_dist.jpeg"), hotspot_dist_join, nrow = 1, base_asp = 1.9)

# hist_regime_vect <- ecoregion_hotspot_zonal_vec %>%
#   select(huc12) %>%
#   left_join(hist_regime)

####################################
### Map of historical fire regime ##
####################################

library(ggplot2)

hist_hotspots <- ecoregion_hotspot_zonal_vec %>%
  left_join(hist_regime %>% select(huc12, historical_type)) %>%
  filter(!is.na(historical_type)) %>%
  mutate(historical_plot = case_when(
    historical_type == "low_sev" ~ "low severity",
    historical_type == "high_sev" ~ "high severity",
    historical_type == "mixed" ~ "mixed severity"
  ))

# get basemap stuff
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

plot_hist_hotspot_types <- function(data, metric_col, legend = TRUE){

  col = sym(metric_col)

  plot <- ggplot() +
    geom_spatvector(data = boundary_states, color = "black", fill = "white") +
    geom_spatvector(data = cbi_forest, fill = "#DDDDDDBF", color = "transparent", alpha = 0.4) +
    geom_spatvector(data = data %>%
                      filter(!!sym(metric_col) == 1, hotspot_type == "Area of Concern"),
                    aes(fill = historical_plot, color = historical_plot)) +
    theme_void() +
    scale_fill_manual(values = list(`low severity` = "#82A6B1", `high severity` =  "#BC4749", `mixed severity` = "#D5A021"), na.value = "transparent") +
    scale_color_manual(values = list(`low severity` = "#82A6B1", `high severity` =  "#BC4749", `mixed severity` = "#D5A021"), na.value = "transparent") +
    theme(legend.title=element_blank(), panel.background = element_rect(fill = "transparent",
                                                                        colour = NA_character_), # necessary to avoid drawing panel outline
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          plot.background = element_rect(fill = "transparent",
                                         colour = NA_character_), # necessary to avoid drawing plot outline
          legend.background = element_rect(fill = NA, color = NA),
          legend.box.background = element_rect(fill = NA, color = NA),
          legend.key = element_rect(fill = "transparent"),
          legend.box = element_blank(),
          legend.text = element_text(size = 14),
          legend.key.size = unit(3,"line")
          )

  if (legend == "FALSE"){
    plot <- plot + theme(legend.position = "none")
  }

  return(plot)
}


hotspot_hist_list <- map(c("breeding_richness", "ecoregion_breeding_lcbd", "FRic_breeding"), plot_hist_hotspot_types, data = hist_hotspots, legend = TRUE)

# cowplot attempt
# legend_plot <- plot_hist_hotspot_types(data = hist_hotspots, "breeding_richness")
# legend <- cowplot::get_legend(legend_plot + theme(legend.box.margin = margin(0, 10, 10, 0)))
# hotspot_hist_figure <- cowplot::plot_grid(hotspot_hist_list[[1]], hotspot_hist_list[[2]], hotspot_hist_list[[3]], nrow = 1)
#
# hotspot_hist_figure <- cowplot::plot_grid(hotspot_hist_figure, legend, rel_widths = c(3, .4))
# cowplot::save_plot(here::here("figures/huc12_hotspot_hist_join.jpeg"), hotspot_hist_figure, nrow = 1)

hotspot_hist_figure <- hotspot_hist_list[[1]] + hotspot_hist_list[[2]] + hotspot_hist_list[[3]] + plot_layout(nrow = 1, guides = "collect")

ggsave(here::here("figures/huc12_hotspot_hist_join.jpeg"), hotspot_hist_figure, width = 420, height = 250, unit = "mm", dpi = 1000)

