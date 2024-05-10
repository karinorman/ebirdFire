## Mission: get the data

library(dplyr)
library(ebirdst)
library(rnaturalearth)
library(sf)
library(terra)

# where the data lives
data_path <- here::here("raw_data")

# dataframe with runs for different years
runs2022 <- read.csv(here::here("raw_data/ebirdst_runs_2022.csv"))
runs2021 <- read.csv(here::here("raw_data/ebirdst_runs_2021.csv"))

# Let's start with the 2022 data, since we can get that fairly easily through the R package
purrr::map(runs2022$species_code, ~ebirdst_download_status(., path = data_path, pattern = "abundance_seasonal_mean"))

# this will make R cranky if you source it before doing the above
source(here::here("R/ebirdst_2021.R"))
source(here::here("R/access-key.R"))

missing_codes <- setdiff(runs2021 %>% filter(scientific_name != "Accipiter gentilis") %>% pull(species_code), runs2022$species_code)
purrr::map(missing_codes, ~ebirdst_download_status_2021(., path = data_path, pattern = "abundance_seasonal_mean"))

# let's get a list of birds for each year that occurred in the US
# start with one we know is here

# Let's play with some data
data_path <- here::here("raw_data")
US_boundary <- ne_states(iso_a2 = "US")

gos <- load_raster("norgos", product = "abundance", period = "seasonal", metric = "mean", resolution = "27km", path = data_path)
US_boundary_proj <- st_transform(US_boundary, st_crs(gos))

gos_mask <- crop(gos, US_boundary_proj) %>%
  mask(US_boundary_proj) %>%
  values(.) %>%
  as_tibble()

colSums(!is.na(gos_mask))

# dataframe of species codes and years
species_by_year <- bind_rows(tibble(species_code = runs2022$species_code, year = 2022),
                             tibble(species_code = missing_codes, year = 2021)) %>%
  mutate(path = data_path)

check_in_US <- function(species_code, year, path){
  if (year == 2022){
    bird_rast <- load_raster_flexyear(species_code, product = "abundance", period = "seasonal",
                                      metric = "mean", resolution = "27km", path = path)
  } else{
    file <- paste0(path, "/2021/", species_code, "/seasonal/", species_code, "_abundance_seasonal_mean_mr_2021.tif")

    bird_rast <- terra::rast(file)
  }

  US_boundary_proj <- st_transform(US_boundary, st_crs(bird_rast))

  bird_mask <- crop(bird_rast, US_boundary_proj) %>%
    mask(US_boundary_proj) %>%
    values(.) %>%
    as_tibble() %>%
    na.omit()

  stack(colSums(bird_mask)) %>%
    tidyr::pivot_wider(names_from = "ind", values_from = "values") %>%
    mutate(across(everything(), .fns = ~ifelse(.x == 0, 0, 1))) %>%
    mutate(species_code = species_code, year = year) %>%
    select(species_code, year, everything())
}

layer_types <- c("resident", "breeding", "nonbreeding", "prebreeding_migration", "postbreeding_migration")

species_layers_df <- purrr::pmap(species_by_year, check_in_US) %>%
  bind_rows() %>%
  mutate(US = if_any(layer_types, ~. != 0),
         US = tidyr::replace_na(US, replace = FALSE))

usethis::use_data(species_layers_df)

# Gotta deal with taxonomy - let's see what the deal is
taxonomy <- read.csv(here::here("raw_data/eBird-Clements-v2023-integrated-checklist-October-2023.csv"))

us_species <- species_layers_df %>%
  filter(US == TRUE) %>%
  pull(species_code)

# these are the US species that had a taxonomic change between the 2021 and 2022 releases
taxonomy %>%
  filter(species_code %in% us_species, Clements.v2023.change != "") %>%
  View()

US_ebird <- species_layers_df %>%
  filter(US == TRUE) %>%
  filter(!species_code %in%
           #remove two flycatchers that were lumped together for 2022
           c("pasfly", "corfly")) %>%
  select(-US)
