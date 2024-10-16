###############################################################################
## Clip the CBI to montana and wyoming boundary per Sean Parks' instructions ##
###############################################################################

library(dplyr)
library(terra)
library(tidyterra)

# get a US boundary base map
US_boundary_states <- rnaturalearth::ne_states(iso_a2 = "US") %>%
  vect() %>%
  project("epsg:4326") %>%
  filter(name %in% c("Washington", "Oregon", "California", "Idaho", "Nevada",
                     "Montana", "Arizona", "Utah", "Wyoming", "Texas", "Colorado", "New Mexico",
                     "Kansas", "Oklahoma")) %>%
  crop(ext(c(-130, -98, 18, 50)))

US_boundary <- US_boundary_states %>% aggregate()

cbi <- rast(here::here("raw_data/predict.high.severity.fire.final.tif")) %>%
  crop(US_boundary, mask = TRUE)

writeRaster(cbi, here::here("data/cbi.tif"))
