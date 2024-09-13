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

cbi <- rast(here::here("raw_data/predict.high.severity.fire.draft.tif"))

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

#################################
######### metric map ############
#################################

plot_metric_map <- function(data, metric_col){

  col = sym(metric_col)
breeding_richness_plot <- ggplot() +
  geom_spatraster(data = data %>% select(!!col),
                  maxcell = 1000e+05) +
  scale_fill_continuous(type = "viridis", na.value = "transparent", name = "species richness") +
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

ggsave(here::here("figures/wus_metrics_map.jpg"), wus_metric_plot, width = 10, height = 10)

#################################
####### bivariate map ###########
#################################

# # we need quantiles for diversity metrics
# rich_quant_b <- quantile(values(metric_rast$breeding_richness),c(0.33,0.66,1), na.rm = TRUE)
# rich_quant_nb <- quantile(values(metric_rast$nonbreeding_richness),c(0.33,0.66,1), na.rm = TRUE)
#
# richness_quant <- metric_rast %>% mutate(breeding_quantile = ifelse(breeding_richness<rich_quant_b[1],1,ifelse(breeding_richness<rich_quant_b[2],2,3)),
#                          nonbreeding_quantile = ifelse(nonbreeding_richness<rich_quant_nb[1],1,ifelse(nonbreeding_richness<rich_quant_nb[2],2,3)))
#
# # we need one raster with a "group" label for each cell
# quant_rast <- c(richness_quant %>% crop(cbi, mask = TRUE),
#                 cbi %>% mutate(fire_quant = ifelse(predict.high.severity.fire.draft == 0, 1, ifelse(predict.high.severity.fire.draft == 1, 2, 3)))) %>%
#   #biodiversity metric is x, fire is y
#   mutate(plot_group = paste0(as.character(breeding_quantile), ",", as.character(fire_quant)),
#          plot_group = ifelse(plot_group == "NA,NA", NA, plot_group)) %>%
#   select(plot_group)
#
# # create a get quantile function
#
# # bivariate color scale
# bivariate_color_scale <- tibble(
#   "3 - 3" = "#3F2949", # high inequality, high income
#   "2 - 3" = "#435786",
#   "1 - 3" = "#4885C1", # low inequality, high income
#   "3 - 2" = "#77324C",
#   "2 - 2" = "#806A8A", # medium inequality, medium income
#   "1 - 2" = "#89A1C8",
#   "3 - 1" = "#AE3A4E", # high inequality, low income
#   "2 - 1" = "#BC7C8F",
#   "1 - 1" = "#CABED0" # low inequality, low income
# ) %>%
#   pivot_longer(everything(), names_to = "group", values_to = "fill")
#
#
# legend<-
#   ggplot(d, aes(x,y,fill=atan(y/x),alpha=x+y))+
#   geom_tile()+
#   #geom_text(alpha=1)+
#   scale_fill_viridis()+
#   theme_void()+
#   theme(legend.position="none",
#         panel.background=element_blank(),
#         plot.margin=margin(t=10,b=10,l=10),
#         axis.title.y = element_text(angle=90))+
#   labs(x="Species Richness \U2192",y="Fire Severity \U2192")+
#   #labs(x = "Species Richness ⟶️",
#        #y = "Fire Severity ⟶️") +
#   theme(axis.title=element_text(color="black"))#+
#   #Draw some arrows:
#   # geom_segment(aes(x=1, xend = 3 , y=0, yend = 0), size=1.5,
#   #              arrow = arrow(length = unit(0.6,"cm"))) +
#   # geom_segment(aes(x=0, xend = 0 , y=1, yend = 3), size=1.5,
#   #              arrow = arrow(length = unit(0.6,"cm")))
# legend
#
# pal_viridis <- ggplot_build(g.legend)$data[[1]] %>%
#   select(x, y, fill, alpha, label) %>%
#   mutate(plot_group = paste0(x, ",", y),
#          new_fill = scales::alpha(fill, alpha))
#
#
# # plotting_rast <- quant_rast %>%
# #   as.data.frame(xy = TRUE) %>%
# #   left_join(pal_viridis) %>%
# #   select(x, y, everything()) %>%
# #   rast(., type = "xyz", crs = "epsg:4326")
#
# color_assign <- setNames(as.character(pal_viridis$new_fill), pal_viridis$plot_group)
#
# bivar <- ggplot() +
#   geom_spatraster(data = quant_rast) +
#   scale_fill_discrete(type = color_assign, na.value = "transparent", name = "species richness", guide = "none") +
#   theme_void()
#
# ggdraw() +
#   draw_plot(bivar, 0, 0, 1, 1) +
#   draw_plot(legend, 0.1, 0.3, 0.2, 0.2)



### Let's try this excluding the non-forest cells
# we need one raster with a "group" label for each cell
quant_rast <- c(richness_quant %>% crop(cbi, mask = TRUE),
                cbi %>% mutate(fire_quant = ifelse(predict.high.severity.fire.draft == 0, NA, ifelse(predict.high.severity.fire.draft == 1, 2, 3)))) %>%
  #biodiversity metric is x, fire is y
  mutate(plot_group = paste0(as.character(breeding_quantile), ",", as.character(fire_quant)),
         plot_group = ifelse(plot_group %in% c("NA,NA", "1,NA", "2,NA", "3,NA"), NA, plot_group)) %>%
  select(plot_group)

color_assign <- setNames(pals::stevens.greenblue(n = 9)[-c(4,5,6)], c("1,2", "2,2", "3,2", "1,3", "2,3", "3,3"))

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
    panel.background = element_rect(color = "white", fill = "white", linewidth = 10)
  ) +
  # quadratic tiles
  coord_fixed()

#background map
US_boundary <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming", "Texas", "Colorado", "New Mexico")) %>%
  crop(ext(c(-130, -103.5, 18, 50)))

plot_ecoregions <- c("Canadian Rocky Mountains", "Middle Rockies - Blue Mountains", "Utah-Wyoming Rocky Mountains", "Southern Rocky Mountains",
                     "Utah High Plateaus", "Colorado Plateau", "Apache Highlands", "Arizona-New Mexico Mountains")

bivar <- ggplot() +
  geom_spatvector(data = US_boundary, color = "black", fill = "white", alpha = 0.5) +
  geom_spatvector(data = ecoregions %>% filter(ECO_NAME %in% plot_ecoregions) %>% terra::aggregate(), color = "black", fill = "white", linewidth = 0.4) +
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
  draw_plot(bivar, 0, 0, 1, 1) +
  draw_plot(legend, 0.25, 0.4, 0.2, 0.2)

ggsave(here::here("figures/richness_cbi_bivar.png"), richness_cbi_bivar, width = 10, height = 10, bg = "transparent")

#################################
########## CBI map ##############
#################################

# Plot of just CBI
cbi_plot <- ggplot() +
  geom_spatvector(data = US_boundary, color = "black", fill = "white") +
  geom_spatraster(data = cbi %>%
                    filter(predict.high.severity.fire.draft %in% c(1,2)) %>%
                    mutate(predict.high.severity.fire.draft = as.factor(predict.high.severity.fire.draft))) +
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

hotspot_poly <- vect(here::here("data/hotspots.shp"))

plot_cbi <- cbi %>%
  mutate(predict.high.severity.fire.draft =
           as.factor(ifelse(predict.high.severity.fire.draft == 0, NA, predict.high.severity.fire.draft)))

plot_refugia <- function(metric_col){
ggplot() +
  geom_spatvector(data = US_boundary, color = "black", fill = "white") +
  geom_spatvector(data = as.polygons(cbi) %>% aggregate(), color = "#DDDDDDBF", alpha = 0.5) +
  geom_spatraster(data = plot_cbi %>% crop(hotspot_poly[[metric_col]], mask = TRUE)) +
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

breeding_lcbd_refugia <- plot_refugia("ecoregion_breeding_lcbd")
ggsave(here::here("figures/refugia_eco_breeding_lcbd.png", breeding_lcbd_refugia))

breeding_richness_refugia <- plot_refugia("breeding_richness")
ggsave(here::here("figures/refugia_breeding_richness.png"), breeding_richness_refugia)

breeding_fric_refugia <- plot_refugia("FRic_breeding")
ggsave(here::here("figures/refugia_breeding_fric.png"), breeding_fric_refugia)


