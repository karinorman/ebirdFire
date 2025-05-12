library(dplyr)
library(terra)
library(tidyterra)
library(purrr)
library(ggplot2)
library(cowplot)

#boundary to clip shapefiles
boundary <- vect(here::here("data/study_boundary.shp"))

cbi <- rast(here::here("data/cbi.tif")) %>%
  crop(boundary, mask = TRUE)

# forest_poly <- cbi %>% filter(predict.high.severity.fire.final %in% c(1,2)) %>%
#   as.polygons(extent = FALSE, na.rm = TRUE, crs = "epsg:4326")

# get ecoregions as shapefile
ecoregion_shp <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  crop(., boundary)

ecoregion_rast <- ecoregion_shp %>%
  select("ECO_NAME") %>%
  terra::rasterize(., cbi, field = "ECO_NAME")# %>%
#resample(cbi, method = "near")

metric_rast <- rast(here::here("data/metric_rast.tiff")) %>%
  crop(ecoregion_shp, mask = TRUE) %>%
  c(., rast(here::here("data/fd_rast.tiff"))%>%
      crop(ecoregion_shp, mask = TRUE)) %>%
  c(., ecoregion_rast)

ecoregion_hotspots <- metric_rast %>%
  terra::as.data.frame(., xy = TRUE) %>%
  group_by(ECO_NAME) %>%
  mutate(across(-c(x,y),
                ~ifelse(.x > quantile(.x, probs = 0.90, na.rm = TRUE), 1, NA))) %>%
  ungroup() %>%
  select(-ECO_NAME) %>%
  rast(.,  type="xyz", crs = "epsg:4326")

forest_hotspots <- c(cbi, metric_rast) %>%
  filter(predict.high.severity.fire.final %in% c(1,2)) %>%
  terra::as.data.frame(., xy = TRUE) %>%
  group_by(ECO_NAME) %>%
  mutate(across(-c(x,y),
                ~ifelse(.x > quantile(.x, probs = 0.90, na.rm = TRUE), 1, NA))) %>%
  ungroup() %>%
  select(-ECO_NAME) %>%
  rast(.,  type="xyz", crs = "epsg:4326") %>%
  rename_with(~paste0("forest_", .x))

forest_quantile_cutoffs <- c(cbi, metric_rast) %>%
  filter(predict.high.severity.fire.final %in% c(1,2)) %>%
  terra::as.data.frame(., xy = FALSE) %>%
  select(-predict.high.severity.fire.final) %>%
  group_by(ECO_NAME) %>%
  summarize(across(everything(), ~quantile(.x, probs = 0.90, na.rm = TRUE))) %>%
  ungroup()

hotspot_rasts <- c(ecoregion_hotspots, forest_hotspots) %>%
  select(-forest_predict.high.severity.fire.final)

writeRaster(hotspot_rasts, here::here("data/hotspot_rasts.tif"))

# get distribution of hotspots across severity types
area <- hotspot_rasts  %>%
  expanse(usenames = TRUE) %>%
  rename(metric = layer)

lowsev_int <- hotspot_rasts %>%
  c(., cbi) %>%
  filter(predict.high.severity.fire.final == 1) %>%
  expanse(usenames = TRUE) %>%
  rename(metric = layer, lowsev_area = area)

highsev_int <-  hotspot_rasts %>%
  c(., cbi) %>%
  filter(predict.high.severity.fire.final == 2) %>%
  expanse(usenames = TRUE) %>%
  rename(metric = layer, highsev_area = area)

# why don't the forest hotspot percentages sum to zero? Should be any other area in the CBI...
overlap_df <- full_join(area, lowsev_int) %>%
  full_join(highsev_int) %>%
  mutate(percent_lowsev = (lowsev_area/area)*100, percent_highsev = (highsev_area/area)*100) %>%
  filter(metric != "predict.high.severity.fire.final")

usethis::use_data(overlap_df)

##################################################################
#### test hotspot distributions against landscape expectation ####
##################################################################

### Perform binomial hypothesis test for each tail, where high severity is defined as a "success"
### Alternative hypothesis: there is significantly more high severity (success) than expected by the ecoregion distribution &
### there is significantly less than expected by the ecoregion distribution

# get ratio of high to low severity for each ecoregion to parameterize "true" distribution
ecoregion_names <- unique(ecoregion_shp$ECO_NAME)

# ecoregion_values <- map_dfr(ecoregion_names, function(x) {
#   browser()
#   filter_ecoregion <- ecoregion_shp %>% filter(ECO_NAME == x)
#   ecoregion_cbi <- crop(cbi, filter_ecoregion, mask = TRUE)
#
#   freq(ecoregion_cbi) %>%
#     select(-layer) %>%
#     pivot_wider(names_from = value, values_from = count) %>%
#     mutate(ECO_NAME = x)
#
# }) %>% mutate(p = `2`/(`1` + `2`))

ecoregion_values <- map_dfr(ecoregion_names, function(x) {

  #browser()
  filter_ecoregion <- ecoregion_shp %>% filter(ECO_NAME == x)
  ecoregion_cbi <- extract(cbi, filter_ecoregion)

  ecoregion_cbi %>%
    count(predict.high.severity.fire.final) %>%
    pivot_wider(names_from = predict.high.severity.fire.final, values_from = n) %>%
    mutate(ECO_NAME = x)


}) %>% mutate(p = `2`/(`1` + `2`))

# get number of success (cells == 2) and number of trials (cells %in% c(1,2)) for each
# ecoregion for hotspots of each metric
hotspot_names <- names(hotspot_rasts)

hotspot_values <- map_dfr(hotspot_names, function(x){

  map_dfr(ecoregion_names, hotspot = hotspot_rasts[[x]], function(name, hotspot) {
    #browser()

    metric_name <- names(hotspot)

    filter_ecoregion <- ecoregion_shp %>% filter(ECO_NAME == name)

    ecoregion_hotspot <- crop(hotspot, filter_ecoregion, mask = TRUE)

    ecoregion_cbi <- crop(cbi, filter_ecoregion, mask = TRUE) %>%
                            c(., ecoregion_hotspot) %>%
      filter(!!as.symbol(metric_name) == 1)

    ecoregion_cbi %>%
      select(predict.high.severity.fire.final) %>%
      freq() %>%
      pivot_wider(names_from = value, values_from = count) %>%
      mutate(ECO_NAME = name) %>%
      select(-layer)

  }) %>% mutate(hotspot_type = x)
})

# dataframe of values for binomial test
binomial_df <- hotspot_values %>%
  rowwise() %>%
  mutate(n = sum(`1`,`2`, na.rm = TRUE)) %>%
  # x is the argument for "successes"
  select(hotspot_type, ECO_NAME, x = `2`, n) %>%
  left_join(ecoregion_values %>% select(ECO_NAME,p)) %>%
  mutate(x = ifelse(is.na(x), 0, x))

# perform a test for upper and lower tail
binomial_test <- binomial_df %>%
  group_by(ECO_NAME, hotspot_type) %>%
  mutate(greater.p.value = binom.test(x = x, n = n, p = p, alternative = "greater") %>% .$p.value,
         lower.p.value = binom.test(x = x, n = n, p = p, alternative = "less") %>% .$p.value) %>%
  ungroup() %>%
  tidyr::separate(hotspot_type, into = "type", sep = "_", remove = FALSE) %>%
  mutate(type = ifelse(type == "forest", "forest", "ecoregion"),
         metric = stringr::str_remove(hotspot_type, "forest_")) %>%
  select(-hotspot_type)

# ecoregions that are significant for upper tail test
greater_sig_ecoregions <- binomial_test %>% filter(greater.p.value < 0.05) %>%
  filter(!metric %in% c("breeding_lcbd", "nonbreeding_lcbd")) %>%
  # let's just look at forest hotspots for now
  filter(type == "forest")

# ecoregions that are significant for lower tail test
lower_sig_ecoregions <- binomial_test %>% filter(lower.p.value < 0.05) %>%
  filter(!metric %in% c("breeding_lcbd", "nonbreeding_lcbd")) %>%
  # let's just look at forest hotspots for now
  filter(type == "forest")

sig_ecoregions_df <- greater_sig_ecoregions %>%
  mutate(test_tail = "upper") %>%
  bind_rows(lower_sig_ecoregions %>% mutate(test_tail = "lower"))

usethis::use_data(sig_ecoregions_df)


########################
####### Figures ########
########################

# get the tile figure showing test results

sig_ecoregions_plot_df <- sig_ecoregions_df %>%
  group_by(ECO_NAME) %>%
  mutate(ECO_NUM = cur_group_id()) %>%
  ungroup() %>%
  mutate(plot_eco_name = paste(ECO_NUM, ECO_NAME, sep = " "))

# horizontal
# sig_ecoregions_fig <- sig_ecoregions_plot_df %>%
#   filter(metric %in% c("breeding_richness", "FRic_breeding", "ecoregion_breeding_lcbd")) %>%
#   select(plot_eco_name, ECO_NUM, test_tail, metric) %>%
#   ggplot(aes(x = reorder(plot_eco_name, ECO_NUM), y = metric, fill = test_tail)) +
#   geom_tile(alpha = 0.8, color = "white") +
#   theme_minimal() +
#   scale_fill_manual(values = list("lower" = "#5296A5", "upper" = "#BC4749"), labels = list("lower" = "lower tail", "upper" = "upper tail"),
#                     na.value = "lightgrey") +
#   scale_x_discrete(expand = c(0,0),
#                    position = "top") +
#   scale_y_discrete(labels = list("breeding_richness" = "Species Richness", "ecoregion_breeding_lcbd" = "LCBD",
#                                  "FRic_breeding" = "Functional Richness"),
#                    expand = c(0,0)) +
#   theme(panel.grid = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 0),
#         panel.border=element_rect(fill = NA, colour=alpha('lightgrey', .5),size=1)) +
#   coord_fixed()

# vertical
sig_ecoregions_fig <- sig_ecoregions_plot_df %>%
  filter(metric %in% c("breeding_richness", "FRic_breeding", "ecoregion_breeding_lcbd")) %>%
  select(plot_eco_name, ECO_NUM, test_tail, metric) %>%
  ggplot(aes(y = reorder(plot_eco_name, ECO_NUM, decreasing = TRUE), x = metric, fill = test_tail)) +
  geom_tile(alpha = 0.8, color = "white") +
  theme_minimal() +
  scale_fill_manual(values = list("lower" = "#5296A5", "upper" = "#BC4749"), labels = list("lower" = "low severity", "upper" = "high severity"),
                    na.value = "lightgrey") +
  scale_y_discrete(expand = c(0,0),
                   position = "left") +
  scale_x_discrete(labels = list("breeding_richness" = "Species Richness", "ecoregion_breeding_lcbd" = "Uniqueness",
                                 "FRic_breeding" = "Functional Richness"), position = "top",
                   expand = c(0,0)) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        panel.border=element_rect(fill = NA, colour=alpha('lightgrey', .5),size=1)) +
  coord_fixed()


#ggsave(here::here("figures/binomial_heatmap.jpeg"), sig_ecoregions_fig)

# get map of ecoregions

# pal <- colorspace::terrain_hcl(n = 20)
#pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))\

pal <- colorRampPalette(c("#d7e8d2", "#3F6634"))

# get a US boundary base map
US_boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming", "Texas", "Colorado", "New Mexico")) %>%
  crop(ext(c(-130, -104, 18, 50)))

US_boundary <- US_boundary_states %>% aggregate()

# get centroids of ecoregions for labels
eco_cent <- centroids(ecoregion_shp, inside = TRUE) %>%
  select(ECO_NAME) %>%
  left_join(sig_ecoregions_plot_df %>% select(ECO_NAME, ECO_NUM) %>% as.data.frame()) %>%
  sf::st_as_sf() %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  as.data.frame() %>%
  mutate(lon = case_when(
    ECO_NUM == 1 ~ -109.7,
    ECO_NUM == 9 ~ -121,
    ECO_NUM == 14 ~ -118.7,
    ECO_NUM == 8 ~ -118,
    ECO_NUM == 13 ~ -106.5,
    ECO_NUM == 18 ~ -110.3,
    ECO_NUM == 6 ~ -115,
    ECO_NUM == 15 ~ -119.6,
    .default = lon
  ),
  lat = case_when(
    ECO_NUM == 1 ~ 32,
    ECO_NUM == 9 ~ 42,
    ECO_NUM == 14 ~ 48.6,
    ECO_NUM == 8 ~ 43,
    ECO_NUM == 13 ~ 47,
    ECO_NUM == 18 ~ 44.5,
    ECO_NUM == 6 ~ 48,
    ECO_NUM == 15 ~ 37.85,
    .default = lat
  ))

set.seed(45)

eco_plt <- ggplot() +
  geom_spatvector(data = US_boundary_states %>% crop(boundary), fill = NA, linetype = "dashed") +
  geom_spatvector(data = boundary, fill = NA, color = "black") +
  geom_spatvector(data = ecoregion_shp, aes(fill = ECO_NAME), alpha = 0.75) +
  geom_point(data = eco_cent, aes(x = lon, y = lat), size = 8, shape = 21, color = "black", fill = "white") +
  geom_text(data = eco_cent %>% select(ECO_NUM, lon, lat) %>% distinct(), aes(label = ECO_NUM, x = lon, y = lat)) +
  scale_fill_manual(values = sample(pal(19))) +
  ggthemes::theme_map() +
  theme(legend.position = "none")

sig_tail_join <- sig_ecoregions_fig + eco_plt + plot_layout(nrow = 1, widths = c(.3, 1.2))
ggsave(here::here("figures/sig_tail_ecoregion_join.jpeg"), sig_tail_join, height = 200, width = 270, units = "mm", dpi = 800)


# map(c("ecoregion_breeding_lcbd", "breeding_richness", "FRic_breeding"), function(metric_name) {
#   #browser()
#
#   greater_sig_metrics <- sig_ecoregions_df %>%
#   filter(test_tail == "upper", metric == metric_name) %>%
#   pull(ECO_NAME)
#
#   lower_sig_metrics <- sig_ecoregions_df %>%
#     filter(test_tail == "lower", metric == metric_name) %>%
#     pull(ECO_NAME)
#
# indicator_shp <- ecoregion_shp %>%
#   mutate(significant = case_when(
#     ECO_NAME %in% greater_sig_metrics ~ "more higher severity",
#     ECO_NAME %in% lower_sig_metrics ~ "less high severity",
#     .default = "nonsignificant"
#   )) %>%
#   crop(., boundary)
#
# plot <- ggplot() +
#   geom_spatvector(mapping = aes(fill = as.factor(significant)), data = indicator_shp, color = "black", linewidth = 0.4) +
#   geom_spatvector(data = US_boundary, fill = NA) +
#   scale_fill_discrete(type = c("#117733", "#e1ad01", "white"), name = "") +
#   theme_map() +
#   ggtitle(metric_name)
#
# ggsave(here::here(paste0("figures/", metric_name, "_binomial_test.png")), plot)
# })
