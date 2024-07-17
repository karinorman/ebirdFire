library(dplyr)
library(terra)
library(purrr)
library(tidyr)
library(tidyterra)

boundary_shp <- vect(here::here("data/study_boundary.shp"))

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

## Create occurrence matrices, include cell id for joining later
breeding_mat <- breeding_df %>%
  select(-c(x, y)) %>%
  mutate(across(everything(), ~tidyr::replace_na(., 0)))

nonbreeding_mat <- nonbreeding_df %>%
  select(-c(x, y)) %>%
  mutate(across(everything(), ~tidyr::replace_na(., 0)))

# write out matrices
write.table(breeding_mat, here::here("data/breeding_occ_mat.csv"), row.names = FALSE, sep=",")
write.table(nonbreeding_mat, here::here("data/nonbreeding_occ_mat.csv"), row.names = FALSE, sep=",")

# write out raster coords for each cell ID
write.csv(breeding_df %>% select(cell, x, y), here::here("data/raster_coords.csv"), row.names = FALSE)
