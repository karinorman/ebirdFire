library(dplyr)
library(rnaturalearth)
library(sf)
library(terra)
library(ebirdst)
library(ggplot2)
library(tidyr) 


# Let's play with some data
data_path <- here::here("raw_data")
gos <- load_raster("norgos", product = "abundance", period = "seasonal", metric = "mean", resolution = "27km", path = data_path)

gos_breed <- gos[["breeding"]]

US_boundary <- ne_states(iso_a2 = "US")
US_boundary_proj <- st_transform(US_boundary, st_crs(gos_breed))

gos_breed_mask <- crop(gos_breed, US_boundary_proj) %>%
  mask(US_boundary_proj)

plot(gos_breed_mask, axes = FALSE)
