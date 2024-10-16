library(dplyr)
library(terra)
library(tidyterra)
library(purrr)
library(ggplot2)
library(cowplot)

cbi <- rast(here::here("data/cbi.tif"))

high_sev <- cbi %>% filter(predict.high.severity.fire.draft == 2) %>%
  as.polygons()

boundary <- vect(here::here("data/study_boundary.shp"))

load(here::here("data/species_list.rda"))

### read in range map rasters ###

# residents
resident_path <- here::here("data/species_ranges/resident")
resident_df <- data.frame(file = list.files(resident_path)) %>%
  mutate(path = paste0(resident_path, "/", file),
         species_code = stringr::str_extract(file, "[^_]+")) %>%
  filter(species_code %in% species_list)

resident_rast <- rast(resident_df$path)
resident_rast_crop <- crop(resident_rast, boundary, mask = TRUE)

total_area <- map_dfr(names(resident_rast), ~ tibble(species_code = .x,
                                                     total_area = freq(resident_rast[[.x]], value = 1)$count))
wus_area <- map_dfr(names(resident_rast_crop), ~ tibble(species_code = .x,
                                                        wus_area = freq(resident_rast_crop[[.x]], value = 1)$count))

resident_polys <- map(names(resident_rast_crop), ~resident_rast_crop[[.x]] %>%
                        subst(., 0, NA) %>%
                        as.polygons())
names(resident_polys) <- names(resident_rast_crop)

# breeding
breeding_path <- here::here("data/species_ranges/breeding")
breeding_df <- data.frame(file = list.files(breeding_path)) %>%
  mutate(breeding_path = paste0(breeding_path, "/", file),
         species_code = stringr::str_extract(file, "[^_]+")) %>%
  filter(species_code %in% species_list)

# nonbreeding
nonbreeding_path <- here::here("data/species_ranges/nonbreeding")
nonbreeding_df <- data.frame(file = list.files(nonbreeding_path)) %>%
  mutate(nonbreeding_path = paste0(nonbreeding_path, "/", file),
         species_code = stringr::str_extract(file, "[^_]+")) %>%
  filter(species_code %in% species_list)

nonresident_df <- breeding_df %>%
  select(-file) %>%
  full_join(nonbreeding_df %>% select(-file)) %>%
  select(species_code, breeding_path, nonbreeding_path)

nonresident_polygons <- pmap(nonresident_df, function(species_code, breeding_path, nonbreeding_path){

  col = sym(species_code)

  if(is.na(breeding_path)) {
    rast_pair <- rast(nonbreeding_path) %>%
      crop(boundary, mask = TRUE)

  } else if (is.na(nonbreeding_path)){

    rast_pair <- rast(breeding_path) %>%
      crop(boundary, mask = TRUE)

  } else {
    rast_pair <- rast(breeding_path, nonbreeding_path) %>%
      crop(boundary, mask = TRUE)
  }

  sum_rast <- sum(rast_pair)

  if(!1 %in% unique(values(sum_rast))){
    return(NULL)
  }else{

    sum_rast %>%
      subst(., 0, NA) %>%
      subst(., 2, 1) %>%
      as.polygons() %>%
      rename(!!col := sum)
  }

})

names(nonresident_polygons) <- nonresident_df$species_code


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
    select(species_code, percent, area)
}

percent_range_highsev_df <- get_intersect(compact(resident_polys), high_sev) %>%
  mutate(range_type = "resident") %>%
  full_join(get_intersect(compact(nonresident_polygons), high_sev) %>% mutate(range_type = "nonresident"))

usethis::use_data(percent_range_highsev_df)

nonresident_area_df <- pmap(nonresident_df, function(species_code, breeding_path, nonbreeding_path){

  if(is.na(breeding_path)) {
    rast_pair <- rast(nonbreeding_path)
  } else if (is.na(nonbreeding_path)){
    rast_pair <- rast(breeding_path)
  } else {
    rast_pair <- rast(breeding_path, nonbreeding_path)
  }

  total_area <- terra::freq(sum(rast_pair)) %>%
    filter(value %in% c(1,2)) %>%
    pull(count) %>%
    sum()

  wus_area <- terra::freq(sum(crop(rast_pair, boundary, mask = TRUE))) %>%
    filter(value %in% c(1,2)) %>%
    pull(count) %>%
    sum()

  tibble(species_code , total_area, wus_area)
})

range_area_df <- nonresident_area_df %>%
  bind_rows(resident_area_df) %>%
  mutate(wus_percent = wus_area/total_area)

species_range_metrics <- range_area_df %>%
  rename(total_range_area = total_area, wus_range_area = wus_area) %>%
  full_join(percent_range_highsev_df) %>%
  rename(high_sev_percent = percent) %>%
  select(-area) %>%
  filter(wus_range_area != 0)

usethis::use_data(species_range_metrics)
