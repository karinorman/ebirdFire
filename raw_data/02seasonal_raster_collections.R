library(dplyr)
library(ebirdst)
library(rnaturalearth)
library(sf)
library(terra)
library(tidyterra)
library(furrr)

load(here::here("data/US_ebird.rda"))

# Create a breeding bird map
# steps are: project -> mosaic -> bianarize
# here we extract seasonal layers, project to espg:4326, and save out

#### resident species #####
dir.create(here::here("data/species_rasters"))
dir.create((here::here("data/species_rasters/resident")))

resident_species <- US_ebird %>%
  filter(resident == 1) %>%
  select(species_code, path)

resident_rasts <- purrr::pmap(resident_species, function(species_code, path) {

  species_code_enquo <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    rename(!!species_code_enquo := resident) %>%
    project("epsg:4326")

  writeRaster(bird_rast, here::here("data/species_rasters/resident", paste0(species_code, "_resident.tif")))

  return(bird_rast)
} ) %>% terra::sds()

##### breeding species ####
dir.create((here::here("data/species_rasters/breeding")))

breeding_species <- US_ebird %>%
  filter(breeding == 1) %>%
  select(species_code, path)

plan(multisession, workers = 10)

breeding_rasts <- furrr::future_pmap(breeding_species[5:476,], function(species_code, path) {
  species_code_enquo <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    .[["breeding"]] %>%
    tidyterra::rename(!!species_code_enquo := breeding) %>%
    project("epsg:4326")

  writeRaster(bird_rast, here::here("data/species_rasters/breeding", paste0(species_code, "_breeding.tif")))

  return(bird_rast)

} )

#### nonbreeding_species ###

dir.create((here::here("data/species_rasters/nonbreeding")))

nonbreeding_species <- US_ebird %>%
  filter(nonbreeding == 1) %>%
  select(species_code, path)

plan(multisession, workers = 10)

breeding_rasts <- furrr::future_pmap(nonbreeding_species[7:435,], function(species_code, path) {
  species_code_enquo <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    .[["nonbreeding"]] %>%
    tidyterra::rename(!!species_code_enquo := nonbreeding) %>%
    project("epsg:4326")

  writeRaster(bird_rast, here::here("data/species_rasters/nonbreeding", paste0(species_code, "_nonbreeding.tif")))

  return(bird_rast)
} )

# bird_rast <- terra::rast("/Users/karinorman/Documents/Projects/ebirdFire/raw_data/2022/zothaw/seasonal/zothaw_abundance_seasonal_mean_3km_2022.tif")
# abetow <- classify(bird_rast, matrix(c(0,1,1), nrow = 1,byrow = TRUE), include.lowest = FALSE)
