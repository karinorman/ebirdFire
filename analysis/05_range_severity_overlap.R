library(dplyr)
library(terra)
library(tidyterra)
library(purrr)
library(ggplot2)
library(cowplot)

cbi <- rast(here::here("raw_data/predict.high.severity.fire.draft.tif"))

high_sev <- cbi %>% filter(predict.high.severity.fire.draft == 2) %>%
  as.polygons()

load(here::here("data/species_list.rda"))

### read in range map rasters ###

# residents
resident_path <- here::here("data/species_ranges_wus/resident")
resident_df <- data.frame(file = list.files(resident_path)) %>%
  mutate(path = paste0(resident_path, "/", file),
         species_code = stringr::str_extract(file, "[^_]+")) %>%
  filter(species_code %in% species_list)

resident_rast <- rast(resident_df$path)
resident_polys <- map(names(resident_rast), ~resident_rast[[.x]] %>%
                        subst(., 0, NA) %>%
                        as.polygons())
names(resident_polys) <- names(resident_rast)

# breeding
breeding_path <- here::here("data/species_ranges_wus/breeding")
breeding_df <- data.frame(file = list.files(breeding_path)) %>%
  mutate(path = paste0(breeding_path, "/", file),
         species_code = stringr::str_extract(file, "[^_]+")) %>%
  filter(species_code %in% species_list)

breeding_rast <- rast(breeding_df$path)
breeding_polys <- map(names(breeding_rast), ~breeding_rast[[.x]] %>%
                        subst(., 0, NA) %>%
                        as.polygons())
names(breeding_polys) <- names(breeding_rast)

# nonbreeding
nonbreeding_path <- here::here("data/species_ranges_wus/non_breeding")
nonbreeding_df <- data.frame(file = list.files(nonbreeding_path)) %>%
  mutate(path = paste0(nonbreeding_path, "/", file),
         species_code = stringr::str_extract(file, "[^_]+")) %>%
  filter(species_code %in% species_list)

nonbreeding_rast <- rast(nonbreeding_df$path)
nonbreeding_polys <- map(names(nonbreeding_rast), ~nonbreeding_rast[[.x]] %>%
                           subst(., 0, NA) %>%
                           as.polygons())
names(nonbreeding_polys) <- names(nonbreeding_rast)

# Get the intersections
intersect <- map(resident_polys, ~ terra::intersect(.x, high_sev) %>% expanse()) %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "metric", values_to = "highsev_intersect")


get_intersect <- function(species_polys, cbi_poly){
  # get dataframe of intersected area for each species
  intersect_df = map(species_polys, ~ terra::intersect(.x, cbi_poly) %>% expanse()) %>%
    compact() %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "species_code", values_to = "intersect")

  area_df <- map(species_polys, expanse) %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "species_code", values_to = "area")

  area_df %>%
    left_join(intersect_df) %>%
    mutate(intersect = replace_na(intersect,0), percent = intersect/area) %>%
    select(species_code, percent)
}

percent_highsev_df <- get_intersect(compact(resident_polys), high_sev) %>%
  rename(resident_percent = percent) %>%
  full_join(get_intersect(compact(breeding_polys), high_sev) %>% rename(breeding_percent = percent)) %>%
  full_join(get_intersect(compact(nonbreeding_polys), high_sev) %>% rename(nonbreeding_percent = percent))
