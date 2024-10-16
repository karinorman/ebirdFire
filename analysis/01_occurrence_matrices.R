library(dplyr)
library(terra)
library(purrr)
library(tidyr)
library(tidyterra)

#### get list of files for each species type ####

# residents
resident_path <- here::here("data/species_ranges_wus/resident")
resident_files <- paste0(resident_path, "/", list.files(resident_path))

# breeding
breeding_path <- here::here("data/species_ranges_wus/breeding")
breeding_files <- paste0(breeding_path, "/", list.files(breeding_path))

# nonbreeding
nonbreeding_path <- here::here("data/species_ranges_wus/non_breeding")
nonbreeding_files <- paste0(nonbreeding_path, "/", list.files(nonbreeding_path))


# function to read in a raster and convert to a df
raster_to_df <- function(raster_path){

  raster <- rast(raster_path)
  species <- names(raster)

  df <- terra::as.data.frame(raster, xy = TRUE, cells = TRUE) %>%
    filter(!! rlang::sym(species) == 1)

  if (dim(df)[1] == 0){
    return(NULL)
  } else{
    return(df)
  }
}

# Get list of species rasters as df, remove empty df (species don't occur in the extent)

breeding_df_list <- c(resident_files, breeding_files) %>%
  map(., raster_to_df) %>% plyr::compact()

nonbreeding_df_list <- c(resident_files, nonbreeding_files) %>%
  map(., raster_to_df) %>% compact()

# join by the cell
breeding_df <- plyr::join_all(breeding_df_list, type = "full")
nonbreeding_df <- plyr::join_all(nonbreeding_df_list, type = "full")

# save
usethis::use_data(breeding_df)
usethis::use_data(nonbreeding_df)

# get ecoregion labels for each cell
# get template raster to compare ecoregions to
ecoregions <- vect(here::here("raw_data/ecoregions/ecoregions_edc.shp")) %>%
  project("epsg:4326")

cell_coords <- breeding_df %>% select(cell, x, y)
cell_coords_vec <- vect(cell_coords, geom = c("x", "y"), crs = "epsg:4326")

ecoregion_int <- intersect(ecoregions, cell_coords_vec)
ecoregion_int_df <- terra::as.data.frame(ecoregion_int %>% select(ECO_NAME, cell))

## Create occurrence matrices, include cell id for joining later
breeding_mat <- breeding_df %>%
  mutate(across(everything(), ~tidyr::replace_na(., 0))) %>%
  select(-c(x, y)) %>%
  left_join(ecoregion_int_df)

nonbreeding_mat <- nonbreeding_df %>%
  mutate(across(everything(), ~tidyr::replace_na(., 0))) %>%
  select(-c(x, y)) %>%
  left_join(ecoregion_int_df)

# write out matrices
write.table(breeding_mat, here::here("data/breeding_occ_mat.csv"), row.names = FALSE, sep=",")
write.table(nonbreeding_mat, here::here("data/nonbreeding_occ_mat.csv"), row.names = FALSE, sep=",")

# save out list of species
species_list <- unique(c(colnames(breeding_mat %>% select(-cell )),
                         colnames(nonbreeding_mat %>% select(-cell))))
usethis::use_data(species_list)

# write out raster coords for each cell ID
write.csv(cell_coords, here::here("data/raster_coords.csv"), row.names = FALSE)
