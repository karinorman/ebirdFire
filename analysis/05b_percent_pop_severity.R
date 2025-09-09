library(dplyr)
library(rnaturalearth)
library(sf)
library(terra)
library(tidyterra)
library(purrr)
library(furrr)

load(here::here("data/species_layers.rda"))
load(here::here("data/species_list.rda"))
boundary <- vect(here::here("data/study_boundary.shp"))
cbi <- rast(here::here("data/cbi.tif"))

species_layers <- species_layers %>%
  filter(species_code %in% species_list)

# Get resident species percentages
get_relabun_mets <- function(species_code, path, range_type){

  # get total global relative abundance
  species_rast <- rast(path) %>% select(!!sym(range_type)) %>% project("epsg:4326")
  total_pop <- sum(values(species_rast), na.rm = TRUE)

  # crop to western US
  crop_rast <- crop(species_rast, boundary, mask = TRUE)
  wus_pop <- sum(values(crop_rast), na.rm = TRUE)

  # get only forest areas
  cbi_resamp <- resample(cbi, crop_rast, method = "near")

  species_forest <- crop_rast %>%
    crop(cbi_resamp, mask = TRUE) %>%
    c(., cbi_resamp) %>%
    filter(predict.high.severity.fire.final %in% c(1,2))

  forest_pop <- sum(values(species_forest %>%
                             select(!!sym(range_type))), na.rm = TRUE)

  # get high severity
  species_cbi <- species_forest %>%
    filter(predict.high.severity.fire.final == 2)

  sev_pop <- sum(values(species_cbi %>% select(!!sym(range_type))), na.rm = TRUE)

  data.frame(species_code = species_code, total_pop = total_pop, wus_pop = wus_pop,
             forest_pop = forest_pop, sev_pop = sev_pop)
}

resident_pop_mets <- purrr::pmap(species_layers %>%
                                   filter(resident == 1) %>%
                                   select(species_code, path), get_relabun_mets,
                                 range_type = "resident") %>%
  bind_rows()

breeding_pop_mets <- purrr::pmap(species_layers %>%
                                   filter(breeding == 1) %>%
                                   select(species_code, path), get_relabun_mets,
                                 range_type = "breeding") %>%
  bind_rows()

nonbreeding_pop_mets <- purrr::pmap(species_layers %>%
                                   filter(nonbreeding == 1) %>%
                                   select(species_code, path), get_relabun_mets,
                                   range_type = "nonbreeding") %>%
  bind_rows()

usethis::use_data(resident_pop_mets)
usethis::use_data(breeding_pop_mets)
usethis::use_data(nonbreeding_pop_mets)

pop_mets <- breeding_pop_mets %>%
  mutate(type = "breeding") %>%
  bind_rows(nonbreeding_pop_mets %>%
              mutate(type = "nonbreeding")) %>%
  # group_by(species_code) %>%
  # summarize(across(everything(), sum)) %>%
  # ungroup() %>%
  # mutate(type = "nonresident") %>%
  bind_rows(resident_pop_mets %>%
              mutate(type = "resident")) %>%
  mutate(total_percent = sev_pop/total_pop,
         forest_percent = forest_pop/total_pop,
        sev_forest_percent = sev_pop/forest_pop)

usethis::use_data(pop_mets)

#############################################
### Let's find just SW species real quick ###
#############################################

sw_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Nevada", "Arizona", "Utah", "Texas", "New Mexico")) %>%
  aggregate()

species_df <- bind_rows(data.frame(path = list.files(here::here("data/species_ranges_wus/breeding/"), full.names = TRUE),
                                   file = list.files(here::here("data/species_ranges_wus/breeding/")),
           range_type = "breeding"),
           data.frame(path = list.files(here::here("data/species_ranges_wus/non_breeding/"), full.names = TRUE),
                      file = list.files(here::here("data/species_ranges_wus/non_breeding/")),
          range_type = "nonbreeding"),
          data.frame(path = list.files(here::here("data/species_ranges_wus/resident/"), full.names = TRUE),
                     file = list.files(here::here("data/species_ranges_wus/resident/")),
                     range_type = "resident")) %>%
  tidyr::separate(file, c("species_code"))


sw_occ_df <- pmap_dfr(species_df, function(path, species_code, range_type){

 # browser()
  rast_crop <- rast(path) %>%
    crop(sw_states) %>%
    as.data.frame()

  rast_perc <- sum(rast_crop) / dim(rast_crop)[1]

  return(data.frame(species_code = species_code, sw_perc = rast_perc,
                    range_type = range_type))
})

sw_species <- sw_occ_df %>%
  filter(sw_perc > 0.15)

usethis::use_data(sw_species)
