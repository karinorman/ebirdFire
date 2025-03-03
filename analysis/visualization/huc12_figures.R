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
    geom_spatvector(data = data, aes(fill = !!col, color = !!col), na.rm = TRUE) +
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
          legend.box = element_blank(),
          legend.key.size = unit(.5, 'cm')
    ) #+
    # theme(legend.key.size=unit(1,'mm'),
    #       legend.text=element_text(size=6),
    #       legend.title=element_text(size=7))
}


metric_plot_df <- tibble(metric_col = c("breeding_richness", "ecoregion_breeding_lcbd", "FRic_breeding"), legend_title = c("Richness", "LCBD", "Functional \nRichness"))

metric_plot_list <- pmap(metric_plot_df, plot_metric_map, data = biodiv_zonal_vec)

#metric_join_plot <- wrap_plots(metric_plot_list, nrow = 1) + plot_annotation(tag_levels = "A")
metric_join_plot <- plot_grid(metric_plot_list[[1]], metric_plot_list[[2]], metric_plot_list[[3]], nrow = 1)

#ggsave(here::here("figures/huc12_metric_join.jpeg"), metric_join_plot)
save_plot(here::here("figures/huc12_metric_join.jpeg"), metric_join_plot, nrow = 1, dpi = 800, base_width = 12)

# ecoregion_lcbd <- plot_metric_map(data = biodiv_zonal_vec, metric_col = "ecoregion_breeding_lcbd", legend_title = "LCBD")
# ggsave(here::here("figures/huc12_ecoreg_lcbd.jpeg"), ecoregion_lcbd)

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
          legend.box = element_blank(),
          legend.text = element_text(size = 14),
          legend.key.size = unit(3,"line")
    )
}


hotspot_plot_list <- map(c("breeding_richness", "ecoregion_breeding_lcbd", "FRic_breeding"), plot_hotspot_types, data = ecoregion_hotspot_zonal_vec)

#hotspot_type_figure <- plot_grid(hotspot_plot_list[[1]], hotspot_plot_list[[2]], hotspot_plot_list[[3]], nrow =1)

hotspot_type_figure <- patchwork::wrap_plots(hotspot_plot_list, nrow = 1, guides = "collect")# +
  #plot_annotation(tag_levels = "A")

ggsave(here::here("figures/huc12_hotspot_join.jpeg"), hotspot_type_figure, width = 180, height = 100, unit = "mm", dpi = 1000)


##### Congruence map between richness, lcbd, fric
cong_hotspot <- ecoregion_hotspot_zonal_vec %>%
  mutate(cong_hotspot = ifelse(breeding_richness ==1 & ecoregion_breeding_lcbd == 1 & FRic_breeding == 1, 1, NA))

cong_hotspot_map <- plot_hotspot_types(cong_hotspot, "cong_hotspot")
ggsave(here::here("figures/huc12_cong_hotspot.jpeg"), cong_hotspot_map, width = 180, height = 100, unit = "mm", dpi = 1000)

#################################
####### bivariate map ###########
#################################

metric_rast <- rast(here::here("data/metric_rast.tiff")) %>%
  crop(boundary, mask = TRUE)

fd_rast <- rast(here::here("data/fd_rast.tiff")) %>%
  crop(boundary, mask = TRUE)

biodiv_rast <- c(metric_rast, fd_rast)

cbi <- rast(here::here("data/cbi.tif")) %>%
  resample(., biodiv_rast, "near")

# exclude unforested areas first
biodiv_rast <- crop(biodiv_rast, cbi, mask = TRUE) %>%
  c(., cbi) %>%
  filter(predict.high.severity.fire.final %in% c(1,2))

# we need quantiles for diversity metrics
rich_quant <- quantile(values(biodiv_rast$breeding_richness),c(0.33,0.66,1), na.rm = TRUE)
lcbd_quant <- quantile(values(biodiv_rast$ecoregion_breeding_lcbd),c(0.33,0.66,1), na.rm = TRUE)
fric_quant <- quantile(values(biodiv_rast$FRic_breeding),c(0.33,0.66,1), na.rm = TRUE)

biodiv_quant <- biodiv_rast %>% mutate(breeding_quantile = ifelse(breeding_richness<rich_quant[1],1,ifelse(breeding_richness<rich_quant[2],2,3)),
                                       lcbd_quantile = ifelse(ecoregion_breeding_lcbd<lcbd_quant[1],1,ifelse(ecoregion_breeding_lcbd<lcbd_quant[2],2,3)),
                                       fric_quantile = ifelse(FRic_breeding<fric_quant[1],1,ifelse(FRic_breeding<fric_quant[2],2,3)))

biodiv_quant <- biodiv_rast %>%
  mutate(across(-predict.high.severity.fire.final, ~factor(findInterval(.x, c(-Inf, quantile(.x, probs=c(0.33, .67), na.rm = TRUE), Inf)),
                                                           labels = c(1,2,3))))

### Let's try this excluding the non-forest cells
# we need one raster with a "group" label for each cell
forest_biodiv_rast <- c(biodiv_quant %>% crop(resample(cbi, biodiv_quant, method = "near"), mask = TRUE),
                        cbi %>% mutate(fire_quant = ifelse(predict.high.severity.fire.final == 0, NA, ifelse(predict.high.severity.fire.final == 1, 2, 3)))) #%>%
#biodiversity metric is x, fire is y
# mutate(plot_group = paste0(as.character(breeding_quantile), ",", as.character(fire_quant)),
#        plot_group = ifelse(plot_group %in% c("NA,NA", "1,NA", "2,NA", "3,NA"), NA, plot_group)) %>%
# select(plot_group)

#color_assign <- setNames(pals::stevens.greenblue(n = 9)[-c(4,5,6)], c("1,2", "2,2", "3,2", "1,3", "2,3", "3,3"))
color_assign <- setNames(c("#fff2c6", "#fed755", "#e1ad01", "#f2dadb", "#d79192", "#BC4749"), c("1,1", "2,1", "3,1", "1,2", "2,2", "3,2"))

color_assign_df <- tibble::enframe(color_assign, name = "plot_group", value = "fill") %>%
  tidyr::separate_wider_delim(plot_group, names = c("metric", "fire"), delim = ",") %>%
  mutate(metric = as.integer(metric),
         fire = as.integer(fire))

# make legend
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
  labs(x="Biodiversity Metric \U2192",y="Fire Severity \U2192") +
  # make font small enough
  theme(
    axis.title = element_text(size = 14),
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

bivar_plot_list <- map(c("breeding_richness", "ecoregion_breeding_lcbd", "FRic_breeding"), function(x){

  plot_data <- biodiv_quant %>%
    mutate(plot_group = paste0(as.character(!!sym(x)), ",", as.character(predict.high.severity.fire.final)),
           plot_group = ifelse(plot_group %in% c("NA,NA", "1,NA", "2,NA", "3,NA", "NA,1", "NA,2"), NA, plot_group)) %>%
    select(plot_group)

  bivar <- ggplot() +
    geom_spatvector(data = boundary, color = "black", fill = "white", alpha = 0.5) +
    #geom_spatvector(data = ecoregions %>% filter(ECO_NAME %in% plot_ecoregions) %>% terra::aggregate(), color = "black", fill = "white", linewidth = 0.4) +
    #geom_spatvector(data = ecoregions %>% terra::aggregate(), color = "black", fill = "transparent", linetype = "dotted", linewidth = .4) +
    geom_spatraster(data = plot_data, maxcell = 2500000) +
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
})

bivar_legend <- ggdraw() +
  draw_plot(bivar_plot_list[[3]], 0, 0, .95, .95) +
  draw_plot(legend, 0.75, 0.3, 0.25, 0.25)

bivar_plot_join <- patchwork::wrap_plots(bivar_plot_list) + legend +
  plot_layout(nrow = 1, ncol = 4, width = c(1,1, 1,.4))

#bivar_plot_list[[1]] + bivar_plot_list[[2]] + bivar_legend + plot_layout(nrow = 1, width = c(1, 1, 1))

ggsave(here::here("figures/metric_bivar.png"), bivar_plot_join, bg = "transparent")

## Calculate percentage in each landscape
# quant <- as.data.frame(biodiv_quant, cells = TRUE)
#
# quant %>% select(cell, breeding_richness, ecoregion_breeding_lcbd,
#                  FRic_breeding, predict.high.severity.fire.final) %>%
#   pivot_longer(-c(cell, predict.high.severity.fire.final), names_to = "metric", values_to = "type") %>%
#   group_by(predict.high.severity.fire.final, metric, type) %>%
#   count() %>%
#   group_by(metric) %>%
#   mutate(perc = n/sum(n)) %>% View()

#######################################
##### master bivar/hostpot figure #####
#######################################

# patchwork::wrap_plots(hotspot_plot_list, nrow = 1, guides = "collect") / (patchwork::wrap_plots(bivar_plot_list) + legend +
#   plot_layout(nrow = 1, ncol = 4, width = c(1.2,1.2, 1.2,.4)))
#
# ((patchwork::wrap_plots(hotspot_plot_list, nrow = 1, guides = "collect") / patchwork::wrap_plots(bivar_plot_list, nrow = 1)) | cong_hotspot_map) #+
#   #plot_layout(widths = c(2,1))


bivar_hotspot <-  (bivar_plot_list[[1]] + bivar_plot_list[[2]] + bivar_plot_list[[3]] +
     (plot_spacer() + legend + plot_spacer() + plot_layout(widths = c(1,4,3))) +
     plot_layout(nrow = 1, width = c(1,1,1,1), heights = c(1,1,1,.5))) /
  (hotspot_plot_list[[1]] + hotspot_plot_list[[2]] + hotspot_plot_list[[3]] + guide_area() + plot_layout(guides = "collect", nrow = 1))

ggsave(here::here("figures/bivar_hotspot_join.jpeg"), bivar_hotspot, dpi = 800, width = 20)

#bivar_hotspot + cong_hotspot_map + plot_layout(widths = c(2.5, 1), heights = c(1, .8))

#(plot_spacer() + cong_hotspot_map + plot_spacer() + plot_layout(heights = c(1,4,1)))
