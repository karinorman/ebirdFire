library(dplyr)
library(terra)
library(tidyterra)
library(ggnewscale)
library(viridis)
library(scales)
library(cowplot)


# maps of metrics for Western US
richness_rast <- rast(here::here("data/richness_rast.tif"))
lcbd_rast <- rast(here::here("data/lcbd_rast.tiff"))

cbi <- rast(here::here("raw_data/predict.high.severity.fire.draft.tif"))

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

### Let's try a bivariate map ###

# we need quantiles for diversity metrics
rich_quant_b <- quantile(values(richness_rast$breeding_richness),c(0.33,0.66,1), na.rm = TRUE)
rich_quant_nb <- quantile(values(richness_rast$nonbreeding_richness),c(0.33,0.66,1), na.rm = TRUE)

richness_quant <- richness_rast %>% mutate(breeding_quantile = ifelse(breeding_richness<rich_quant_b[1],1,ifelse(breeding_richness<rich_quant_b[2],2,3)),
                         nonbreeding_quantile = ifelse(nonbreeding_richness<rich_quant_nb[1],1,ifelse(nonbreeding_richness<rich_quant_nb[2],2,3)))

# we need one raster with a "group" label for each cell
quant_rast <- c(richness_quant %>% crop(cbi, mask = TRUE),
                cbi %>% mutate(fire_quant = ifelse(predict.high.severity.fire.draft == 0, 1, ifelse(predict.high.severity.fire.draft == 1, 2, 3)))) %>%
  #biodiversity metric is x, fire is y
  mutate(plot_group = paste0(as.character(breeding_quantile), ",", as.character(fire_quant)),
         plot_group = ifelse(plot_group == "NA,NA", NA, plot_group)) %>%
  select(plot_group)

# create a get quantile function

# bivariate color scale
bivariate_color_scale <- tibble(
  "3 - 3" = "#3F2949", # high inequality, high income
  "2 - 3" = "#435786",
  "1 - 3" = "#4885C1", # low inequality, high income
  "3 - 2" = "#77324C",
  "2 - 2" = "#806A8A", # medium inequality, medium income
  "1 - 2" = "#89A1C8",
  "3 - 1" = "#AE3A4E", # high inequality, low income
  "2 - 1" = "#BC7C8F",
  "1 - 1" = "#CABED0" # low inequality, low income
) %>%
  pivot_longer(everything(), names_to = "group", values_to = "fill")


legend<-
  ggplot(d, aes(x,y,fill=atan(y/x),alpha=x+y))+
  geom_tile()+
  #geom_text(alpha=1)+
  scale_fill_viridis()+
  theme_void()+
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.margin=margin(t=10,b=10,l=10),
        axis.title.y = element_text(angle=90))+
  labs(x="Species Richness \U2192",y="Fire Severity \U2192")+
  #labs(x = "Species Richness ⟶️",
       #y = "Fire Severity ⟶️") +
  theme(axis.title=element_text(color="black"))#+
  #Draw some arrows:
  # geom_segment(aes(x=1, xend = 3 , y=0, yend = 0), size=1.5,
  #              arrow = arrow(length = unit(0.6,"cm"))) +
  # geom_segment(aes(x=0, xend = 0 , y=1, yend = 3), size=1.5,
  #              arrow = arrow(length = unit(0.6,"cm")))
legend

pal_viridis <- ggplot_build(g.legend)$data[[1]] %>%
  select(x, y, fill, alpha, label) %>%
  mutate(plot_group = paste0(x, ",", y),
         new_fill = scales::alpha(fill, alpha))


# plotting_rast <- quant_rast %>%
#   as.data.frame(xy = TRUE) %>%
#   left_join(pal_viridis) %>%
#   select(x, y, everything()) %>%
#   rast(., type = "xyz", crs = "epsg:4326")

color_assign <- setNames(as.character(pal_viridis$new_fill), pal_viridis$plot_group)

bivar <- ggplot() +
  geom_spatraster(data = quant_rast) +
  scale_fill_discrete(type = color_assign, na.value = "transparent", name = "species richness", guide = "none") +
  theme_void()

ggdraw() +
  draw_plot(bivar, 0, 0, 1, 1) +
  draw_plot(legend, 0.1, 0.3, 0.2, 0.2)
