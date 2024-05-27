library(dplyr)
library(ebirdst)
library(rnaturalearth)
library(sf)
library(terra)
library(tidyterra)

load(here::here("data/US_ebird.rda"))

# Create a breeding bird map
# steps are: project -> mosaic -> bianarize

# resident species
resident_species <- US_ebird %>%
  filter(resident == 1) %>%
  select(species_code, path)

resident_rasts <- purrr::pmap(resident_species, function(species_code, path) {
  species_code <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    rename(!!species_code := resident) %>%
    project("epsg:4326")

} ) %>% terra::sds()

# breeding species
breeding_species <- US_ebird %>%
  filter(breeding == 1) %>%
  select(species_code, path)

breeding_rasts <- purrr::pmap(breeding_species, function(species_code, path) {
  species_code <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    .[["breeding"]] %>%
    rename(!!species_code := breeding) %>%
    project("epsg:4326")

} )  %>% terra::sds()

# nonbreeding_species
nonbreeding_species <- US_ebird %>%
  filter(nonbreeding == 1) %>%
  select(species_code, path)

nonbreeding_rasts <- purrr::pmap(nonbreeding_species, function(species_code, path) {
  species_code <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    .[["nonbreeding"]] %>%
    rename(!!species_code := nonbreeding) %>%
    project("epsg:4326")

} )  %>% terra::sds()


# bird_rast <- terra::rast("/Users/karinorman/Documents/Projects/ebirdFire/raw_data/2022/zothaw/seasonal/zothaw_abundance_seasonal_mean_3km_2022.tif")
# abetow <- classify(bird_rast, matrix(c(0,1,1), nrow = 1,byrow = TRUE), include.lowest = FALSE)
