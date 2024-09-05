### Find list of species to exclude from individual species models ###

library(dplyr)
library(terra)
library(tidyterra)

load(here::here("data/species_list.rda"))
boundary <- vect(here::here("data/study_boundary.shp"))

## 1. Get total range area for entire species range and breeding + nonbreeding

### read in range map rasters ###

# residents
resident_path <- here::here("data/species_ranges/resident")
resident_df <- data.frame(file = list.files(resident_path)) %>%
  mutate(path = paste0(resident_path, "/", file),
         species_code = stringr::str_extract(file, "[^_]+")) %>%
  filter(species_code %in% species_list)

resident_rast <- rast(resident_df$path)
resident_rast_crop <- crop(resident_rast, boundary, mask = TRUE)

total_area <- map_dfr(names(resident_rast)[1:3], ~ tibble(species_code = .x,
                                          total_area = freq(resident_rast[[.x]], value = 1)$count))
wus_area <- map_dfr(names(resident_rast_crop)[1:3], ~ tibble(species_code = .x,
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

nonresident_df <- full_join(breeding_df, nonbreeding_df)

# map across data frame to get paired breeding and non breeding, add


## 2. Get total area within Western US boundary

## 3.
