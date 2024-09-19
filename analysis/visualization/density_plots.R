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

cbi <- rast(here::here("raw_data/predict.high.severity.fire.draft.tif"))

#### Quadrants of Concern ####
ecoregion_rast <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  rasterize(., metric_rast, field = "ECO_NAME")

metric_stack <- c(metric_rast, fd_rast, ecoregion_rast) %>%
  crop(., cbi) %>%
  c(., cbi) %>%
  select(cbi = predict.high.severity.fire.draft, breeding_lcbd, nonbreeding_lcbd, breeding_richness, nonbreeding_richness,
         FRic_breeding, FEve_breeding, FDiv_breeding,
         FRic_nonbreeding, FEve_nonbreeding, FDiv_nonbreeding, ecoregion = ECO_NAME)

metric_stack_df <- as.data.frame(metric_stack, xy = TRUE)

### Plots a single aggreagate density plot, with tails shaded with darker color
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

metric_names <- colnames(metric_stack_df)[!colnames(metric_stack_df) %in% c("x", "y", "cbi", "ecoregion")]

map(metric_names,plot_metric_density)


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
    theme(text = element_text(size = 15))

  ggsave(here::here("figures", paste0(metric_col, "density.png")))
})

### Plot Map
ecoregion_vect <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326")

US_boundary <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming", "Texas", "Colorado", "New Mexico")) %>%
  crop(ext(c(-130, -103.5, 18, 50)))

inc_ecoregion <- metric_stack_df %>%
  filter(cbi %in% c(1, 2)) %>%
  pull(ecoregion) %>% unique()

ecoregion_map <- ggplot() +
  geom_spatvector(data = US_boundary, color = "black", fill = "transparent", alpha = 0.5) +
  geom_spatvector(data = ecoregion_vect %>% filter(ECO_NAME %in% inc_ecoregion),
                  aes(fill = factor(ECO_NAME)),
                  color = "black", linewidth = 0.4) +
  scale_fill_manual(values = list(`Canadian Rocky Mountains` = "#CC6677", `Middle Rockies - Blue Mountains` = "#332288",
                                  `Utah-Wyoming Rocky Mountains` = "#117733", `Southern Rocky Mountains` = "#88CCEE",
                                  `Utah High Plateaus` = "#882255", `Colorado Plateau` = "#44AA99", `Arizona-New Mexico Mountains` = "#999933",
                                  `Apache Highlands` = "#AA4499"), drop = F, name = NULL) +
  theme_void()

ggsave(here::here("figures/ecoregion_map.png"), ecoregion_map)
