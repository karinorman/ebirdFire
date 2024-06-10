library(dplyr)
library(terra)
library(purrr)

# get US boundary to crop
US_boundary <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326")

# read in range map rasters

# residents
resident_path <- here::here("data/species_ranges/resident")
resident_files <- paste0(resident_path, "/", list.files(resident_path))

# breeding
breeding_path <- here::here("data/species_ranges/breeding")
breeding_files <- paste0(breeding_path, "/", list.files(breeding_path))

# nonbreeding
nonbreeding_path <- here::here("data/species_ranges/nonbreeding")
nonbreeding_files <- paste0(nonbreeding_path, "/", list.files(nonbreeding_path))


##### Get single raster for breeding + residents #####

breeding_rast <- rast(c(resident_files, breeding_files)) #%>%
  #crop(ext(c(-130, -50, 18, 50)))

breeding_richness_rast <- app(breeding_rast, sum, na.rm = TRUE)

writeRaster(breeding_richness_rast, here::here("data/breeding_richness.tif"))
# this reaches memory limit
# breeding_df <- terra::as.data.frame(breeding_rast, xy = TRUE)


##### Get single raster for nonbreeding + residents #####
nonbreeding_rast <- rast(c(resident_files, nonbreeding_files)) #%>%
 # crop(ext(c(-130, -50, 18, 50)))

nonbreeding_richness_rast <- app(nonbreeding_rast, sum, na.rm = TRUE)
writeRaster(nonbreeding_richness_rast, here::here("data/nonbreeding_richness.tif"))

# this reaches memory limit
# nonbreeding_df <- terra::as.data.frame(nonbreeding_rast, xy = TRUE)

