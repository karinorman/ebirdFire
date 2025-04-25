library(dplyr)
library(terra)
library(tidyterra)
library(ggnewscale)
library(viridis)
library(scales)
library(purrr)
library(ggridges)
library(ggplot2)


####################################
######## Density Plots #############
####################################

metric_rast <- rast(here::here("data/metric_rast.tiff"))
fd_rast <- rast(here::here("data/fd_rast.tiff"))

cbi <- rast(here::here("data/cbi.tif"))

#### Quadrants of Concern ####
ecoregion_rast <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  rasterize(., metric_rast, field = "ECO_NAME")

metric_stack <- c(metric_rast, fd_rast, ecoregion_rast) %>%
  crop(., cbi) %>%
  c(., cbi) %>%
  select(cbi = predict.high.severity.fire.final, breeding_lcbd, nonbreeding_lcbd, breeding_richness, nonbreeding_richness,
         FRic_breeding, FEve_breeding, FDiv_breeding,
         FRic_nonbreeding, FEve_nonbreeding, FDiv_nonbreeding, ecoregion = ECO_NAME)

metric_stack_df <- as.data.frame(metric_stack, xy = TRUE)

### Plots a single aggreagate density plot, with tails shaded with darker color
# plot_metric_density <- function(metric_col){
#   col = sym(metric_col)
#
#   means <- metric_stack_df %>%
#     filter(cbi %in% c(1,2)) %>%
#     group_by(cbi) %>%
#     summarize(mean = mean(!!col))
#
#   hold_plt <- metric_stack_df %>%
#     filter(cbi %in% c(1,2)) %>%
#     ggplot(aes(x = !!col, fill = as.factor(cbi))) +
#     geom_density(alpha = 0.25)
#
#   build <- ggplot2::ggplot_build(hold_plt)
#
#   df_breaks <- build$data[[1]] %>%
#     group_by(fill) %>%
#     mutate(cbi = ifelse(fill == "#00BFC4", 1, 2),
#            severity = ifelse(fill == "#00BFC4", "low severity", "high severity")) %>%
#     ungroup() %>%
#     left_join(means) %>%
#     mutate(status = ifelse(x < mean, "low biodiv value", "high biodiv value")) %>%
#     mutate(plot_group = paste(severity, status))
#
#   pal <- list(`low severity low biodiv value` = "#fee48e", `low severity high biodiv value` = "#e1ad01",
#               `high severity low biodiv value` = "#f4afae", `high severity high biodiv value` = "#bd1b19")
#   pal_line <- list(`1` = "#e1ad01", `2` =  "#bd1b19" )
#
#   df_breaks %>%
#     ggplot() +
#     geom_area(
#       aes(x = x, y = density, fill = plot_group), position = "identity", alpha = 0.60, color = "black"
#     ) +
#     scale_fill_manual(values = pal) +
#     guides(fill=guide_legend(title=element_blank())) +
#     theme_classic() +
#     xlab(metric_col) +
#     theme(axis.line.y=element_blank(),
#           axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#           axis.title.y=element_blank()) +
#     geom_vline(aes(xintercept = mean, color = as.factor(cbi))) +
#     scale_color_manual(values = pal_line, guide = "none")
# }
#
 metric_names <- colnames(metric_stack_df)[!colnames(metric_stack_df) %in% c("x", "y", "cbi", "ecoregion")]
#
# map(metric_names,plot_metric_density)


### Plot two seperate plots for high and low severity, density plots by ecoregion
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

map(metric_names, plot_ridglines)

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

map(metric_names, plot_grouped_ridgelines)

#### Custom palette for each ecoregion ####

#get palette values
pal <- Polychrome::kelly.colors(22)
pal <- within(as.list(pal), rm(white, black))

ecoregion_names <- na.omit(unique(metric_stack_df$ecoregion))

# create dataframe with ecoregion names, palette assignment, and programatically get lighter value of palette
pal_df <- tibble(ecoregion_names) %>%
  # alphabetical order
  arrange(ecoregion_names) %>%
  # assign palette in order
  mutate(pal = unlist(pal[row_number()])) %>%
  rowwise() %>%
  # get lighter palette value
  mutate(pal_light = colorspace::lighten(pal, 0.6, space = "HLS"),
         ecoregion_dark = paste(ecoregion_names, "2"),
         ecoregion_light = paste(ecoregion_names, "1"))

# dataframe to named list that can be passed to the manual fill in ggplot
pal_list <- c(tibble::deframe(pal_df %>% select(ecoregion_dark, pal)),
              tibble::deframe(pal_df %>% select(ecoregion_light, pal_light)))

# want to maintain order or fill assignment, get list to enforce order
levels_order <- pal_df %>% select(ecoregion_light, ecoregion_dark) %>% t %>% as.list()

# dataframe with metrics and their x label
xlab_df <- data.frame(metric_col = metric_names[metric_names != "ecoregion"],
                      xlab = c("Breeding LCBD", "Nonbreeding LCBD", "Breeding Richness", "Nonbreeding Richness",
                               "Breeding Functional Richness", "Breeding Functional Evenness", "Breeding Functional Divergence",
                               "Nonbreeding Functional Richness", "Nonbreeding Functional Evenness", "Nonbreeding Functional Divergence"))

density_plot_list <- pmap(xlab_df %>%
                            # modify arguments slightly for final figure, add which plot get's y axis
                            filter(metric_col %in% c("breeding_richness", "breeding_lcbd", "FRic_breeding")) %>%
                            arrange(match(metric_col, c("breeding_richness", "breeding_lcbd", "FRic_breeding"))) %>%
                            bind_cols(y_axis = c(TRUE, FALSE, FALSE)) %>%
                            mutate(xlab = c("Species Richness", "LCBD", "Functional Richness")),
                          function(metric_col, xlab, y_axis){
  col = sym(metric_col)

  plt <- metric_stack_df %>%
    #filter(!ecoregion %in% c("West Cascades", "Sierra Nevada", "Okanagan", "Klamath Mountains", "Great Basin",
                             # "East Cascades - Modoc Plateau", "Columbia Plateau",
                             # "California South Coast", "California Central Coast",
                             # "Northern Great Plains Steppe", "California North Coast")) %>%
    filter(cbi %in% c(1, 2)) %>%
    mutate(fill_var = paste(ecoregion, cbi)) %>%
    filter(cbi %in% c(1,2)) %>%
    ggplot(aes(x = !!col, y = forcats::fct_rev(factor(ecoregion, levels = pal_df$ecoregion_names)), fill = factor(fill_var, levels = levels_order))) +
    geom_density_ridges(#quantile_lines = TRUE,
      alpha = 0.75,
                        #calc_ecdf = TRUE,
                        #geom = "density_ridges_gradient",
                        #quantiles = c(0.95)
      ) +
    theme_classic() +
    ylab(element_blank()) +
    scale_fill_manual(values = pal_list, guide = "none") +
    ## code below adds a color to the shaded tails
    # new_scale_fill() +
    # stat_density_ridges(aes(fill = stat(quantile)), quantile_lines = TRUE,
    #                     calc_ecdf = TRUE,
    #                     geom = "density_ridges_gradient",
    #                     quantiles = c(0.95)) +
    # scale_fill_manual(name = "Prob.",
    #                   values = c("#BD1B1900", "#DDDDDDBF"),
    #                   labels = c("low biodiversity", "area of concern"),
    #                   guide = "none") +
    xlab(xlab) +
    theme(text = element_text(size = 10))

  if (y_axis == FALSE){
  plt <- plt +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  }

  return(plt)

  #ggsave(here::here("figures", paste0(metric_col, "density.png")))
})

density_plt_join <- patchwork::wrap_plots(density_plot_list, nrow = 1, axes = "collect_y")

# Get Cascades as main text example

cascades_plots <- pmap(xlab_df %>%
                         # modify arguments slightly for final figure, add which plot get's y axis
                         filter(metric_col %in% c("breeding_richness", "breeding_lcbd", "FRic_breeding")) %>%
                         arrange(match(metric_col, c("breeding_richness", "breeding_lcbd", "FRic_breeding"))) %>%
                         bind_cols(y_axis = c(TRUE, FALSE, FALSE)) %>%
                         mutate(xlab = c("Species Richness", "Uniqueness", "Functional Richness")),
                       function(metric_col, xlab, y_axis){
                         col = sym(metric_col)

                         plt <- metric_stack_df %>%
                           filter(cbi %in% c(1, 2),
                                  ecoregion %in% c("West Cascades", "East Cascades - Modoc Plateau")) %>%
                           mutate(fill_var = paste(ecoregion, cbi)) %>%
                           filter(cbi %in% c(1,2)) %>%
                           ggplot(aes(x = !!col, y = forcats::fct_rev(factor(ecoregion, levels = pal_df$ecoregion_names)),
                                      fill = factor(cbi, levels = c("2", "1"))#,
                                      #height = ..ndensity..
                                      )) +
                           geom_density_ridges(alpha = 0.75) +
                           theme_classic() +
                           ylab(element_blank()) +
                           scale_fill_manual(values = list(`2` = "#bd1b19", `1` = "#e1ad01"),
                                             labels = c("high severity", "low severity")) +
                           xlab(xlab) +
                           theme(text = element_text(size = 15),
                                 legend.title = element_blank()) +
                           scale_y_discrete(expand = c(0,0),
                                            labels = list("West Cascades", "East Cascades"))
                           #ylim(c(0.5, 4))

                         if (y_axis == FALSE){
                           plt <- plt +
                             theme(axis.line.y = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())
                         }

                         return(plt)

                         #ggsave(here::here("figures", paste0(metric_col, "density.png")))
                       })

cascades_plt_join <- patchwork::wrap_plots(cascades_plots, nrow = 1, axes = "collect_y", guides = "collect")

ggsave(here::here("figures/cascades_density.jpeg"), cascades_plt_join, width = 15, height = 5.5, units = "in", dpi = 800)

### Ecoregion Map with color assignments ###
US_boundary <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming", "Texas", "Colorado", "New Mexico")) %>%
  crop(ext(c(-130, -103.5, 18, 50)))

ecoregion_vect <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  crop(US_boundary)

ecoregion_map <- ggplot() +
  geom_spatvector(data = US_boundary, color = "black", fill = "transparent", alpha = 0.5) +
  geom_spatvector(data = ecoregion_vect,
                  aes(fill = factor(ECO_NAME)),
                  color = "black", linewidth = 0.4) +
  scale_fill_manual(values = tibble::deframe(pal_df %>% select(ecoregion_names, pal)), drop = F, name = NULL, guide = "none") +
  theme_void()

#ggsave(here::here("figures/ecoregion_map.png"), ecoregion_map)

# Add density plots and map, save
density_map <- density_plt_join + ecoregion_map + plot_layout(heights = c(2, 1.6))
ggsave(here::here("figures/density_map.jpeg"), density_map, width = 15, height = 8)
