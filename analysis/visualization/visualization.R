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
richness_rast <- rast(here::here("data/richness_rast.tif"))
lcbd_rast <- rast(here::here("data/lcbd_rast.tiff"))

cbi <- rast(here::here("raw_data/predict.high.severity.fire.draft.tif"))

boundary <- terra::vect(here::here("data/study_boundary.shp"))

ecoregions <- terra::vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326")

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
  geom_spatvector(data = US_boundary, color = "black", fill = "transparent", alpha = 0.5) +
  geom_spatvector(data = ecoregions %>% filter(ECO_NAME %in% plot_ecoregions) %>% terra::aggregate(), color = "black", fill = "white", linewidth = 0.4) +
  #geom_spatvector(data = ecoregions %>% terra::aggregate(), color = "black", fill = "transparent", linetype = "dotted", linewidth = .4) +
  geom_spatraster(data = quant_rast, maxcell = 2500000) +
  scale_fill_discrete(type = color_assign, na.value = "transparent", name = "species richness", guide = "none") +
  theme_void()

richness_cbi_bivar <- ggdraw() +
  draw_plot(bivar, 0, 0, 1, 1) +
  draw_plot(legend, 0.25, 0.4, 0.2, 0.2)

ggsave(here::here("figures/richness_cbi_bivar.jpg"), richness_cbi_bivar, width = 10, height = 10)


####################################
######## Density Plots #############
####################################

metric_rast <- rast(here::here("data/metric_rast.tiff"))
fd_rast <- rast(here::here("data/fd_rast.tiff"))

#### Quadrants of Concern ####
ecoregions <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  rasterize(., metric_rast, field = "ECO_NAME")

metric_stack <- c(metric_rast, fd_rast, ecoregions) %>%
  crop(., cbi) %>%
  c(., cbi) %>%
  select(cbi = predict.high.severity.fire.draft, breeding_lcbd, nonbreeding_lcbd, breeding_richness, nonbreeding_richness,
         FRic_breeding, FEve_breeding, FDiv_breeding,
         FRic_nonbreeding, FEve_nonbreeding, FDiv_nonbreeding, ecoregion = ECO_NAME)

metric_stack_df <- as.data.frame(metric_stack, xy = TRUE)
#
# density_plt <- metric_stack_df %>%
#   filter(cbi %in% c(1,2)) %>%
#   ggplot(aes(x = breeding_lcbd, fill = as.factor(cbi))) +
#   geom_density(alpha = 0.25) +
#   theme_classic()
#
# build <- ggplot2::ggplot_build(density_plt)
#
# df_breaks <- build$data[[1]] %>%
#   group_by(fill) %>%
#   mutate(mean_val = mean(x),
#          severity = ifelse(fill == "#00BFC4", "low severity", "high severity"),
#          status = ifelse(x < mean_val, "low lcbd", "high lcbd")) %>%
#   ungroup() %>%
#   mutate(plot_group = paste(severity, status))
#
# pal <- list(`low severity low lcbd` = "#c9d5c5", `low severity high lcbd` = "#77976e",
#             `high severity low lcbd` = "#fee48e", `high severity high lcbd` = "#e1ad01")
# df_breaks %>%
#   ggplot() +
#   geom_area(
#     aes(x = x, y = y, fill = plot_group)
#   ) +
#   scale_fill_manual(values = pal) +
#   guides(fill=guide_legend(title=element_blank())) +
#   theme_classic() +
#   xlab("Breeding LCBD") +
#   theme(axis.line.y=element_blank(),
#         axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#         axis.title.y=element_blank())


plot_metric_density <- function(metric_col){
  col = sym(metric_col)

  means <- metric_stack_df %>%
    filter(cbi %in% c(1,2)) %>%
    group_by(cbi) %>%
    summarize(mean = mean(!!col))

  hold_plt <- metric_stack_df %>%
    filter(cbi %in% c(1,2)) %>%
    ggplot(aes(x = !!col, fill = as.factor(cbi))) +
    geom_density(alpha = 0.25)

  build <- ggplot2::ggplot_build(hold_plt)

  df_breaks <- build$data[[1]] %>%
    group_by(fill) %>%
    mutate(cbi = ifelse(fill == "#00BFC4", 1, 2),
           severity = ifelse(fill == "#00BFC4", "low severity", "high severity")) %>%
    ungroup() %>%
    left_join(means) %>%
    mutate(status = ifelse(x < mean, "low biodiv value", "high biodiv value")) %>%
    mutate(plot_group = paste(severity, status))

  pal <- list(`low severity low biodiv value` = "#fee48e", `low severity high biodiv value` = "#e1ad01",
              `high severity low biodiv value` = "#f4afae", `high severity high biodiv value` = "#bd1b19")
  pal_line <- list(`1` = "#e1ad01", `2` =  "#bd1b19" )

  df_breaks %>%
    ggplot() +
    geom_area(
      aes(x = x, y = density, fill = plot_group), position = "identity", alpha = 0.60, color = "black"
    ) +
    scale_fill_manual(values = pal) +
    guides(fill=guide_legend(title=element_blank())) +
    theme_classic() +
    xlab(metric_col) +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          axis.title.y=element_blank()) +
    geom_vline(aes(xintercept = mean, color = as.factor(cbi))) +
    scale_color_manual(values = pal_line, guide = "none")
}

plot_metric_density("breeding_lcbd")

metric_names <- colnames(metric_stack_df)[!colnames(metric_stack_df) %in% c("x", "y", "cbi")]

map(metric_names,plot_metric_density)

plot_ridglines <- function(metric_col){
  col = sym(metric_col)

  ## ridgline plot approach
  metric_stack_df %>%
    filter(cbi == 2) %>%
    ggplot(aes(x = !!col, y = ecoregion, fill = stat(quantile))) +
    stat_density_ridges(quantile_lines = TRUE,
                        calc_ecdf = TRUE,
                        geom = "density_ridges_gradient",
                        quantiles = c(0.95)) +
    scale_fill_manual(name = "Prob.", values = c("#f4afae", "#bd1b19"),
                      labels = c("low biodiversity", "area of concern")) +
    theme_classic()


  metric_stack_df %>%
    filter(cbi == 1) %>%
    ggplot(aes(x = !!col, y = ecoregion, fill = stat(quantile))) +
    stat_density_ridges(quantile_lines = TRUE,
                        calc_ecdf = TRUE,
                        geom = "density_ridges_gradient",
                        quantiles = c(0.95)) +
    scale_fill_manual(name = "Prob.", values = c("#fee48e", "#e1ad01"),
                      labels = c("low biodiversity", "refugia")) +
    theme_classic()
}

map(metric_names[metric_names != "ecoregion"], plot_ridglines)


# can we manually change the fill when we have a grouped ridgeline plot
# doesn't work because it's not expecting the aesthetic to change
# grouped_ridgeline <- metric_stack_df %>%
#   filter(cbi %in% c(1,2)) %>%
#   ggplot(aes(x = breeding_richness, y = ecoregion, fill = as.factor(cbi))) +
#   geom_density_ridges(quantile_lines = TRUE, alpha = 0.75,
#                       calc_ecdf = TRUE,
#                       #geom = "density_ridges_gradient",
#                       quantiles = c(0.95)) +
#   theme_classic() +
#   ylab(element_blank()) +
#   guides(fill = guide_legend(title = "Fire Severity"))
#
# plot_parts <- ggplot_build(grouped_ridgeline)
#
# plt_data <- plot_parts$data[[1]] %>%
#   mutate(fill = case_when(
#   fill == "#F8766D" & quantile == 1 ~ "#fee48e",
#   fill == "#F8766D" & quantile == 2 ~ "#e1ad01",
#   fill == "#00BFC4" & quantile == 1 ~ "#f4afae",
#   fill == "#00BFC4" & quantile == 2 ~ "#bd1b19"
# ))
#
# plot_parts$data[[1]] <- plt_data
#
# ggplot_gtable(plot_parts)


# let's just got overlapping ridgeline plots, without the shaded upper region but with
# the quantile line
plot_grouped_ridgelines <- function(metric_col){
  col = sym(metric_col)

  grouped_ridgeline <- metric_stack_df %>%
    filter(cbi %in% c(1,2)) %>%
    ggplot(aes(x = !!col, y = ecoregion, fill = factor(cbi, levels = c("2", "1")))) +
    geom_density_ridges(quantile_lines = TRUE, alpha = 0.75,
                        calc_ecdf = TRUE,
                        #geom = "density_ridges_gradient",
                        quantiles = c(0.95)) +
    theme_classic() +
    ylab(element_blank()) +
    guides(fill = guide_legend(title = "Fire Severity")) +
    scale_fill_manual(values = list(`1` = "#e1ad01", `2` = "#bd1b19"),
                      labels = c("low", "high"))
}

map(metric_names[metric_names != "ecoregion"], plot_grouped_ridgelines)


## Plots with custom palette
pal <- list(`Canadian Rocky Mountains 2` = "#CC6677", `Canadian Rocky Mountains 1` = "#ebc2c9",
            `Middle Rockies - Blue Mountains 1` = "#c6beef", `Middle Rockies - Blue Mountains 2` = "#332288",
            `Utah-Wyoming Rocky Mountains 1` = "#3ee375", `Utah-Wyoming Rocky Mountains 2` = "#117733",
            `Southern Rocky Mountains 1` = "#d2ecf9", `Southern Rocky Mountains 2` = "#88CCEE",
            `Utah High Plateaus 1` = "#db6fa5", `Utah High Plateaus 2` = "#882255",
            `Colorado Plateau 1` = "#b7e2db", `Colorado Plateau 2` = "#44AA99",
            `Arizona-New Mexico Mountains 1` = "#e3e3aa", `Arizona-New Mexico Mountains 2` = "#999933",
            `Apache Highlands 1` = "#e9c9e4", `Apache Highlands 2` = "#AA4499")

levels_order <- c("Canadian Rocky Mountains 2","Canadian Rocky Mountains 1",
                  "Middle Rockies - Blue Mountains 2","Middle Rockies - Blue Mountains 1",
                  "Utah-Wyoming Rocky Mountains 2", "Utah-Wyoming Rocky Mountains 1",
                  "Southern Rocky Mountains 2", "Southern Rocky Mountains 1",
                  "Utah High Plateaus 2", "Utah High Plateaus 1",
                  "Colorado Plateau 2", "Colorado Plateau 1",
                  "Arizona-New Mexico Mountains 2" ,"Arizona-New Mexico Mountains 1",
                  "Apache Highlands 2", "Apache Highlands 1")

# Example with colors by ecoregion with shaded tails
xlab_df <- data.frame(metric_col = metric_names[metric_names != "ecoregion"],
                      xlab = c("Breeding LCBD", "Nonbreeding LCBD", "Breeding Richness", "Nonbreeding Richness",
                        "Breeding Functional Richness", "Breeding Functional Evenness", "Breeding Functional Divergence",
                        "Nonbreeding Functional Richness", "Nonbreeding Functional Evenness", "Nonbreeding Functional Divergence"))

pmap(xlab_df, function(metric_col, xlab){
  col = sym(metric_col)

  metric_stack_df %>%
    filter(cbi %in% c(1, 2)) %>%
    mutate(fill_var = paste(ecoregion, cbi)) %>%
    filter(cbi %in% c(1,2)) %>%
    ggplot(aes(x = !!col, y = ecoregion, fill = factor(fill_var, levels = levels_order))) +
    geom_density_ridges(quantile_lines = TRUE, alpha = 0.75,
                        calc_ecdf = TRUE,
                        #geom = "density_ridges_gradient",
                        quantiles = c(0.95)) +
    theme_classic() +
    ylab(element_blank()) +
    scale_fill_manual(values = pal, guide = "none") +
    # new_scale_fill() +
    # stat_density_ridges(aes(fill = stat(quantile)), quantile_lines = TRUE,
    #                     calc_ecdf = TRUE,
    #                     geom = "density_ridges_gradient",
    #                     quantiles = c(0.95)) +
    # scale_fill_manual(name = "Prob.",
    #                   values = c("#BD1B1900", "#DDDDDDBF"),
    #                   labels = c("low biodiversity", "area of concern"),
    #                   guide = "none") +
    xlab(xlab)
})

### Plot Map
ecoregions <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326")

inc_ecoregion <- metric_stack_df %>%
  filter(cbi %in% c(1, 2)) %>%
  pull(ecoregion) %>% unique()

ggplot() +
  geom_spatvector(data = US_boundary, color = "black", fill = "transparent", alpha = 0.5) +
  geom_spatvector(data = ecoregions %>% filter(ECO_NAME %in% inc_ecoregion),
                  aes(fill = factor(ECO_NAME)),
                  color = "black", linewidth = 0.4) +
  scale_fill_manual(values = list(`Canadian Rocky Mountains` = "#CC6677", `Middle Rockies - Blue Mountains` = "#332288",
                                  `Utah-Wyoming Rocky Mountains` = "#117733", `Southern Rocky Mountains` = "#88CCEE",
                                  `Utah High Plateaus` = "#882255", `Colorado Plateau` = "#44AA99", `Arizona-New Mexico Mountains` = "#999933",
                                  `Apache Highlands` = "#AA4499"), drop = F, name = NULL) +
  theme_void()
