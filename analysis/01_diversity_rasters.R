library(dplyr)
library(terra)
library(tidyterra)

breeding_lcbd <- read.csv(here::here("data/breeding_lcbd.csv")) %>%
  rename(breeding_lcbd = LCBD)
nonbreeding_lcbd <- read.csv(here::here("data/nonbreeding_lcbd.csv")) %>%
  rename(nonbreeding_lcbd = LCBD)

coords <- read.csv(here::here("data/raster_coords.csv"))


lcbd <- full_join(breeding_lcbd, nonbreeding_lcbd) %>%
  full_join(coords) %>%
  select(x, y, breeding_lcbd, nonbreeding_lcbd)


lcbd_rast <- rast(lcbd, type="xyz", crs = "epsg:4326")
