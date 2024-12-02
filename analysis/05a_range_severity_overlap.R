library(dplyr)
library(terra)
library(tidyterra)
library(purrr)
library(ggplot2)
library(cowplot)

cbi <- rast(here::here("data/cbi.tif"))

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

nonresident_rast <- pmap(nonresident_df, function(species_code, breeding_path, nonbreeding_path){

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

    return(sum_rast %>% rename(!!col := sum))
    # sum_rast %>%
    #   subst(., 0, NA) %>%
    #   subst(., 2, 1) %>%
    #   as.polygons() %>%
    #   rename(!!col := sum)
  }

})

#names(nonresident_polygons) <- nonresident_df$species_code
species_rasts_crop  <- c(resident_rast_crop, rast(unlist(nonresident_rast)))

cbi_resamp <- resample(cbi, species_rasts_crop, method = "near")

high_sev_int_rasts <- species_rasts_crop %>%
  crop(cbi_resamp, mask = TRUE) %>%
  c(., cbi_resamp) %>%
  filter(predict.high.severity.fire.final == 2) %>%
  select(-predict.high.severity.fire.final)

sev_area <- map_dfr(names(high_sev_int_rasts), ~ tibble(species_code = .x,
                                                        sev_area = freq(high_sev_int_rasts[[.x]], value = 1)$count))

wus_area <- map_dfr(names(species_rasts_crop), ~ tibble(species_code = .x,
                                                        wus_area = freq(species_rasts_crop[[.x]], value = 1)$count))

# Get total range areas for resident and nonresident species
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
#
#   wus_area <- terra::freq(sum(crop(rast_pair, boundary, mask = TRUE))) %>%
#     filter(value %in% c(1,2)) %>%
#     pull(count) %>%
#     sum()

  tibble(species_code, total_area)
}) %>% bind_rows()

total_area <- map_dfr(names(resident_rast), ~ tibble(species_code = .x,
                                                     total_area = freq(resident_rast[[.x]], value = 1)$count)) %>%
  bind_rows(nonresident_area_df)


species_range_metrics <- total_area %>%
  left_join(wus_area) %>%
  left_join(sev_area) %>%
  mutate(wus_percent = wus_area/total_area,
         high_sev_percent = sev_area/wus_area) %>%
  filter(wus_area != 0)

usethis::use_data(species_range_metrics)
