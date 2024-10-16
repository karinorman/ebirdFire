library(dplyr)
library(ebirdst)
library(rnaturalearth)
library(sf)
library(terra)
library(tidyterra)
library(purrr)
library(furrr)

load(here::here("data/species_layers.rda"))

# get Western US and ecoregion boundaryies to crop
US_boundary <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming")) %>%
  crop(ext(c(-130, -105.5, 18, 50))) %>%
  aggregate()

# cbi is clipped to the easter border of Montana and Wyoming, so need to remove that section from the ecoregions
cbi_boundary <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming", "Texas", "Colorado", "New Mexico",
                     "Kansas", "Oklahoma"))

# ecoregions for wester boundary
ecoregions <- terra::vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326") %>%
  filter(ECO_NAME %in% c("Northern Great Plains Steppe", "Southern Rocky Mountains",
                         "Arizona-New Mexico Mountains", "Apache Highlands")) %>%
  aggregate() %>%
  crop(cbi_boundary)

boundary <- rbind(US_boundary, ecoregions) %>%
  aggregate() %>%
  fillHoles()

writeVector(boundary, here::here("data/study_boundary.shp"))

###########################
####### this makes ########
#### binary range maps ####
#### complete and crop ####
###########################

#### resident species #####
dir.create(here::here("data/species_ranges"))
dir.create((here::here("data/species_ranges/resident")))

resident_species <- species_layers %>%
  filter(resident == 1) %>%
  select(species_code, path)

plan(multisession, workers = 10)

resident_rasts <- furrr::future_pmap(resident_species, function(species_code, path) {

  species_code_enquo <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    tidyterra::rename(!!species_code_enquo := resident) %>%
    project("epsg:4326") %>%
    terra::classify(matrix(c(0, minmax(.)[2,1], 1), nrow = 1,byrow = TRUE), include.lowest = FALSE)

  writeRaster(bird_rast, here::here("data/species_ranges/resident", paste0(species_code, "_resident_range.tif")))

  return(bird_rast)
} )

dir.create(here::here("data/species_ranges_wus"))
dir.create((here::here("data/species_ranges_wus/resident")))

resident_path <- here::here("data/species_ranges/resident")
resident_files <- tibble(file = list.files(resident_path),
                         species_code = stringr::str_extract(file, "[^_]+"),
                         path =  paste0(resident_path, "/", file)

) %>%
  select(-file)

pmap(resident_files, function(species_code, path) {

  species_code_enquo <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    terra::crop(boundary, mask = T)

  writeRaster(bird_rast, here::here("data/species_ranges_wus/resident", paste0(species_code, "_resident_range_wus.tif")))

} )

##### breeding species ####
dir.create((here::here("data/species_ranges/breeding")))

breeding_species <- species_layers %>%
  filter(breeding == 1) %>%
  select(species_code, path)

plan(multisession, workers = 10)

breeding_rasts <- furrr::future_pmap(breeding_species, function(species_code, path) {
  species_code_enquo <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    .[["breeding"]] %>%
    tidyterra::rename(!!species_code_enquo := breeding) %>%
    project("epsg:4326") %>%
    terra::classify(matrix(c(0, minmax(.)[2,1], 1), nrow = 1,byrow = TRUE), include.lowest = FALSE)

  writeRaster(bird_rast, here::here("data/species_ranges/breeding", paste0(species_code, "_breeding_range.tif")))

  return(bird_rast)

} )

dir.create((here::here("data/species_ranges_wus/breeding")))

breeding_path <- here::here("data/species_ranges/breeding")
breeding_files <- tibble(file = list.files(breeding_path),
                         species_code = stringr::str_extract(file, "[^_]+"),
                         path =  paste0(breeding_path, "/", file)

) %>%
  select(-file)

pmap(breeding_files, function(species_code, path) {

  species_code_enquo <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    terra::crop(boundary, mask = T)

  writeRaster(bird_rast, here::here("data/species_ranges_wus/breeding", paste0(species_code, "_breeding_range_wus.tif")))

} )

#### nonbreeding_species ###

dir.create((here::here("data/species_ranges/nonbreeding")))

nonbreeding_species <- species_layers %>%
  filter(nonbreeding == 1) %>%
  select(species_code, path)

plan(multisession, workers = 10)

nonbreeding_rasts <- furrr::future_pmap(nonbreeding_species, function(species_code, path) {
  species_code_enquo <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    .[["nonbreeding"]] %>%
    tidyterra::rename(!!species_code_enquo := nonbreeding) %>%
    project("epsg:4326") %>%
    terra::classify(matrix(c(0, minmax(.)[2,1], 1), nrow = 1,byrow = TRUE), include.lowest = FALSE)

  writeRaster(bird_rast, here::here("data/species_ranges/nonbreeding", paste0(species_code, "_nonbreeding_range.tif")))

  return(bird_rast)
} )

dir.create((here::here("data/species_ranges_wus/non_breeding")))

nonbreeding_path <- here::here("data/species_ranges/nonbreeding")
nonbreeding_files <- tibble(file = list.files(nonbreeding_path),
                         species_code = stringr::str_extract(file, "[^_]+"),
                         path =  paste0(nonbreeding_path, "/", file)

) %>%
  select(-file)

pmap(nonbreeding_files, function(species_code, path) {

  species_code_enquo <- enquo(species_code)

  bird_rast <- terra::rast(path) %>%
    terra::crop(boundary, mask = T)

  writeRaster(bird_rast, here::here("data/species_ranges_wus/non_breeding", paste0(species_code, "_non_breeding_range_wus.tif")))

} )

