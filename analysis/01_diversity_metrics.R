library(dplyr)
library(terra)
library(purrr)
library(tidyr)

# residents
resident_path <- here::here("data/species_ranges_wus/resident")
resident_files <- paste0(resident_path, "/", list.files(resident_path))

# breeding
breeding_path <- here::here("data/species_ranges_wus/breeding")
breeding_files <- paste0(breeding_path, "/", list.files(breeding_path))

# nonbreeding
nonbreeding_path <- here::here("data/species_ranges_wus/nonbreeding")
nonbreeding_files <- paste0(nonbreeding_path, "/", list.files(nonbreeding_path))

raster_to_df <- function(raster_path){

  browser()

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
  map(., raster_to_df) %>% compact()

nonbreeding_df_list <- c(resident_files, nonbreeding_files) %>%
  map(., raster_to_df) %>% compact()

# join by the cell
breeding_df <- plyr::join_all(breeding_df_list)
nonbreeding_df <- plyr::join_all(nonbreeding_df_list)

usethis::use_data(breeding_df)
usethis::use_data(nonbreeding_df)


comm_mat <- breeding_df %>%
  select(-c(cell, x, y)) %>%
  mutate(across(everything(), ~tidyr::replace_na(., 0)))

breeding_beta <- adespatial::beta.div(comm_mat, samp = FALSE, nperm = 10)
usethis::use_data(breeding_beta)

# get LCBD as raster
# read in raster to use as template
temp_rast <- c(rast(breeding_files[1]), rast(breeding_files[2]),
               rast(breeding_files[3]))
temp_rast_df <- terra::as.data.frame(temp_rast, cells = TRUE, xy = TRUE)

breeding_lcbd <- tibble(LCBD = breeding_beta$LCBD,
                        p.LCBD = breeding_beta$p.LCBD,
                        p.adj = breeding_beta$p.adj) %>%
  bind_cols(breeding_df %>% select(cell, x, y), .) %>%
  left_join(temp_rast_df %>% select(cell, x, y), .)

lcbd_rast <- setValues(temp_rast, breeding_lcbd)



#
# ##### Get single raster for breeding + residents #####
#
# dir.create(here::here("data/wus_breeding_csv/"))
#
# c(resident_files, breeding_files) %>%
#   map_dfr(., function(raster_path){
#     raster <- rast(raster_path)
#     print(names(raster))
#     terra::as.data.frame(raster, xy = TRUE) %>%
#       mutate(species_code = names(raster), ncells = dim(.)[1]) %>%
#       write.csv(paste0(here::here("data/wus_breeding_csv/"),
#                        names(raster), "_breeding.csv"))
#   })
#
# csv_files <- paste0(here::here("data/wus_breeding_csv/"), list.files(here::here("data/wus_breeding_csv/")))
#
# breeding_df <- map_dfr(csv_files, read.csv)
#
# #%>%
# #crop(ext(c(-130, -50, 18, 50)))
# #
# # ### Analysis for only Western US ###
# # wus_breeding <- breeding_rast %>%
# #   crop(US_boundary) %>%
# #   crop(ext(c(-125, -95, 24, 50))) %>%
# #   mask(US_boundary)
# #
# # writeRaster(wus_breeding, here::here("data/wus_breeding.tif"))
#
# #create tiles
# dir.create(here::here("data/wus_breeding_richness/"))
# tile_paths <- makeTiles(breeding_rast, 30, filename = here::here("data/wus_breeding_richness/breeding_tile_.tif"))
#
# breeding_df <- purrr::map_dfr(tile_paths, function(path){
#   tile <- rast(path)
#   as.data.frame(tile, xy = TRUE)
# })
#
# #### let's try the BAT package
# library(BAT)
# breeding_mets <-  breeding_rast %>%
#   crop(ext(c(-130, -50, 18, 50))) %>%
#   raster.beta()
#
# # Get community composition matrix from raster
# library(letsR)
#
# # get grid
# grid <- breeding_richness_rast %>%
#   crop(US_boundary) %>%
#   as.polygons(aggregate = FALSE)
#
# breeding_presab <- lets.presab.grid(breeding_rast, grid, "ID")
