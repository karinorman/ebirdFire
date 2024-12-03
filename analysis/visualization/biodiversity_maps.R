library(dplyr)
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

metric_rast <- rast(here::here("data/metric_rast.tiff"))
fd_rast <- rast(here::here("data/fd_rast.tiff"))

cbi <- rast(here::here("data/cbi.tif"))

boundary <- terra::vect(here::here("data/study_boundary.shp"))

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
####### bivariate map ###########
#################################

# we need quantiles for diversity metrics
rich_quant_b <- quantile(values(metric_rast$breeding_richness),c(0.33,0.66,1), na.rm = TRUE)
rich_quant_nb <- quantile(values(metric_rast$nonbreeding_richness),c(0.33,0.66,1), na.rm = TRUE)

richness_quant <- metric_rast %>% mutate(breeding_quantile = ifelse(breeding_richness<rich_quant_b[1],1,ifelse(breeding_richness<rich_quant_b[2],2,3)),
                         nonbreeding_quantile = ifelse(nonbreeding_richness<rich_quant_nb[1],1,ifelse(nonbreeding_richness<rich_quant_nb[2],2,3)))

### Let's try this excluding the non-forest cells
# we need one raster with a "group" label for each cell
quant_rast <- c(richness_quant %>% crop(cbi, mask = TRUE),
                cbi %>% mutate(fire_quant = ifelse(predict.high.severity.fire.final == 0, NA, ifelse(predict.high.severity.fire.final == 1, 2, 3)))) %>%
  #biodiversity metric is x, fire is y
  mutate(plot_group = paste0(as.character(breeding_quantile), ",", as.character(fire_quant)),
         plot_group = ifelse(plot_group %in% c("NA,NA", "1,NA", "2,NA", "3,NA"), NA, plot_group)) %>%
  select(plot_group)

#color_assign <- setNames(pals::stevens.greenblue(n = 9)[-c(4,5,6)], c("1,2", "2,2", "3,2", "1,3", "2,3", "3,3"))
color_assign <- setNames(c("#fff2c6", "#fed755", "#e1ad01", "#f2dadb", "#d79192", "#BC4749"), c("1,2", "2,2", "3,2", "1,3", "2,3", "3,3"))

color_assign_df <- tibble::enframe(color_assign, name = "plot_group", value = "fill") %>%
  tidyr::separate_wider_delim(plot_group, names = c("metric", "fire"), delim = ",") %>%
  mutate(metric = as.integer(metric),
         fire = as.integer(fire))

legend <- ggplot() +
  geom_tile(
    data = color_assign_df,
    mapping = aes(
      x = metric,
      y = fire,
      fill = fill)
  ) +
  scale_fill_identity() +
  theme_classic() +
  labs(x="Species Richness \U2192",y="Fire Severity \U2192") +
  # make font small enough
  theme(
    axis.title = element_text(size = 10),
    axis.title.y = element_text(angle=90),
    axis.line=element_blank(),
    axis.ticks=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
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
  ) +
  # quadratic tiles
  coord_fixed()


bivar <- ggplot() +
  geom_spatvector(data = boundary, color = "black", fill = "white", alpha = 0.5) +
  #geom_spatvector(data = ecoregions %>% filter(ECO_NAME %in% plot_ecoregions) %>% terra::aggregate(), color = "black", fill = "white", linewidth = 0.4) +
  #geom_spatvector(data = ecoregions %>% terra::aggregate(), color = "black", fill = "transparent", linetype = "dotted", linewidth = .4) +
  geom_spatraster(data = quant_rast, maxcell = 2500000) +
  scale_fill_discrete(type = color_assign, na.value = "transparent", name = "species richness", guide = "none") +
  theme_void() +
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


richness_cbi_bivar <- ggdraw() +
  draw_plot(bivar, 0, 0, .95, .95) +
  draw_plot(legend, 0.75, 0.3, 0.25, 0.25)

ggsave(here::here("figures/richness_cbi_bivar.png"), richness_cbi_bivar, width = 10, height = 10, bg = "transparent")

#################################
########## CBI map ##############
#################################

# Plot of just CBI
cbi_plot <- ggplot() +
  geom_spatvector(data = US_boundary, color = "black", fill = "white") +
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

plot_refugia <- function(metric_col){
ggplot() +
  geom_spatvector(data = US_boundary, color = "black", fill = "white") +
  geom_spatvector(data = as.polygons(cbi) %>% aggregate(), color = "#DDDDDDBF", alpha = 0.5) +
  geom_spatraster(data = plot_cbi %>% crop(hotspot_rast[[metric_col]], mask = TRUE)) +
  scale_fill_manual(values = list(`1` = "#82A6B1", `2` =  "#BC4749"), na.value = "#DDDDDDBF", na.translate = FALSE,
                    labels = c("Refugia", "Areas of Concern")) +
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
}

breeding_lcbd_refugia <- plot_refugia("forest_ecoregion_breeding_lcbd")
ggsave(here::here("figures/refugia_eco_breeding_lcbd.png"), breeding_lcbd_refugia)

breeding_richness_refugia <- plot_refugia("forest_breeding_richness")
ggsave(here::here("figures/refugia_breeding_richness.png"), breeding_richness_refugia)

breeding_fric_refugia <- plot_refugia("forest_FRic_breeding")
ggsave(here::here("figures/refugia_breeding_fric.png"), breeding_fric_refugia)

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
