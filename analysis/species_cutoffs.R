### Find list of species to exclude from individual species models ###

library(dplyr)
library(terra)
library(tidyterra)

load(here::here("data/species_list.rda"))
boundary <- vect(here::here("data/study_boundary.shp"))

## 1. Get total range area for entire species range and breeding + nonbreeding

# Get total area and Western US area for each species - residents, and nonresidents (breeding + nonbreeding range)

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

resident_area_df <- total_area %>% left_join(wus_area)

# breeding
breeding_path <- here::here("data/species_ranges/breeding")
breeding_df <- data.frame(file = list.files(breeding_path)) %>%
  mutate(breeding_path = paste0(breeding_path, "/", file),
         species_code = stringr::str_extract(file, "[^_]+")) %>%
  filter(species_code %in% species_list) %>%
  select(-file)

# nonbreeding
nonbreeding_path <- here::here("data/species_ranges/nonbreeding")
nonbreeding_df <- data.frame(file = list.files(nonbreeding_path)) %>%
  mutate(nonbreeding_path = paste0(nonbreeding_path, "/", file),
         species_code = stringr::str_extract(file, "[^_]+")) %>%
  filter(species_code %in% species_list) %>%
  select(-file)

nonresident_df <- full_join(breeding_df, nonbreeding_df) %>%
  select(species_code, breeding_path, nonbreeding_path)

# map across data frame to get paired breeding and non breeding, add

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

usethis::use_data(range_area_df)
