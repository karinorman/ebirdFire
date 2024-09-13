library(dplyr)
library(terra)
library(tidyterra)
library(JuliaCall)

## Let's use JuliaCall to run the lcbd computation directly from R

####################################
### Get the data pieces we need ####
####################################

# read in the community occurrence matrices
breeding_mat <- read.csv(here::here("data/breeding_occ_mat.csv"))
nonbreeding_mat <- read.csv(here::here("data/nonbreeding_occ_mat.csv"))

# coordinates for the cell id's
coords <- read.csv(here::here("data/raster_coords.csv"))

####################################
###### set up Julia session ########
####################################

# get julia running and make sure we have the packages loaded
julia_setup(JULIA_HOME  = "/Users/karinorman/.julia/juliaup/julia-1.10.4+0.aarch64.apple.darwin14/bin")

julia_install_package_if_needed("Tables")
julia_install_package_if_needed("Statistics")
julia_install_package_if_needed("DataFrames")

julia_library("Tables")
julia_library("Statistics")
julia_library("DataFrames")

# create the two necessary functions
julia_command("function betadiv(Y::Matrix)
    # S -> squared deviations from column mean
    S = (Y .- mean(Y; dims=1)) .^ 2.0
    # SStotal -> total sum of squares
    SStotal = sum(S)
    # BDtotal -> index of beta diversity, unbiased & comparable estimator of Var(Y)
    BDtotal = SStotal / (size(Y, 1) - 1)
    # SSj -> sum of squares for species j
    SSj = sum(S; dims=1)
    # SCBDj -> species contribution to beta diversity (species j, relative)
    SCBDj = SSj ./ SStotal
    # SSi -> sum of squares for site i
    SSi = sum(S; dims=2)
    # LCBD -> local contribution to beta diversity (site i, relative)
    LCBDi = SSi ./ SStotal
    # Combine results in tuple
    res = (
        S=S, SStotal=SStotal, BDtotal=BDtotal, SSj=SSj, SCBDj=SCBDj, SSi=SSi, LCBDi=LCBDi
    )
    return res
end")

julia_command("function scaledLCBD(Y::DataFrame)
    #figure out how to propogate and identifying label as a column
    commMat = Y[:, Not([:cell])] |> Tables.matrix

    betaTuple = betadiv(commMat)
    lcbdVals = betaTuple.LCBDi
    lcbdVals = lcbdVals ./ maximum(lcbdVals)

    lcbdDF = DataFrame(cell = Y[:, :cell], LCBD = vec(lcbdVals))

    return lcbdDF
end")

# function to take df from R, pass it to julia, and do the computation
julia_lcbd <- function(df){
  ## df MUST be a community occ matrix (sites on rows species on columns) with a cell identifier column
  julia_assign("comm_df", df)

  lcbd_df = julia_eval("df = scaledLCBD(comm_df)")
}

####################################
####### Get the df's of lcbd #######
####################################

# get lcbd relative to the entire landscape
breeding_lcbd <- julia_lcbd(breeding_mat %>% select(-ECO_NAME)) %>%
  rename(breeding_lcbd = LCBD)
nonbreeding_lcbd <- julia_lcbd(nonbreeding_mat %>% select(-ECO_NAME)) %>%
  rename(nonbreeding_lcbd = LCBD)

# get lcbd within each ecoregion
ecoregion_breeding_lcbd <- breeding_mat %>%
  filter(!is.na(ECO_NAME)) %>%
  group_by(ECO_NAME) %>%
  group_map(~ julia_lcbd(.x %>% select(-ECO_NAME)) %>%
              mutate(ECO_NAME = unique(.x$ECO_NAME)), .keep = TRUE) %>%
  bind_rows() %>%
  rename(ecoregion_breeding_lcbd = LCBD) %>%
  # let's do quantile calculation for hotspots now while we've got the ecoregion labels
  group_by(ECO_NAME) %>%
  mutate(lcbd_breeding_quantile_cutoff = quantile(ecoregion_breeding_lcbd, probs = 0.95, na.rm = TRUE))

ecoregion_nonbreeding_lcbd <- nonbreeding_mat %>%
  filter(!is.na(ECO_NAME)) %>%
  group_by(ECO_NAME) %>%
  group_map(~ julia_lcbd(.x %>% select(-ECO_NAME)) %>%
              mutate(ECO_NAME = unique(.x$ECO_NAME)), .keep = TRUE) %>%
  bind_rows() %>%
  rename(ecoregion_nonbreeding_lcbd = LCBD) %>%
  # let's do quantile calculation for hotspots now while we've got the ecoregion labels
  group_by(ECO_NAME) %>%
  mutate(lcbd_nonbreeding_quantile_cutoff = quantile(ecoregion_nonbreeding_lcbd, probs = 0.95, na.rm = TRUE))

####################################
####### Get species richness #######
####################################
breeding_richness <- breeding_mat %>%
  select(cell, ECO_NAME, everything()) %>%
  #filter(!is.na(ECO_NAME)) %>%
  mutate(breeding_richness = rowSums(across(3:last_col()))) %>%
  select(cell, ECO_NAME, breeding_richness) %>%
  group_by(ECO_NAME) %>%
  mutate(breeding_quantile_cutoff = quantile(breeding_richness, probs = 0.95, na.rm = TRUE)) %>%
  ungroup()

nonbreeding_richness <- nonbreeding_mat %>%
  select(cell, ECO_NAME, everything()) %>%
  #filter(!is.na(ECO_NAME)) %>%
  mutate(nonbreeding_richness = rowSums(across(3:last_col()))) %>%
  select(cell, ECO_NAME, nonbreeding_richness) %>%
  group_by(ECO_NAME) %>%
  mutate(nonbreeding_quantile_cutoff = quantile(nonbreeding_richness, probs = 0.95, na.rm = TRUE)) %>%
  ungroup()

# join, column for each metric
metrics <- full_join(breeding_lcbd, nonbreeding_lcbd) %>%
  full_join(ecoregion_breeding_lcbd) %>%
  full_join(ecoregion_nonbreeding_lcbd) %>%
  full_join(breeding_richness) %>%
  full_join(nonbreeding_richness) %>%
  full_join(coords) %>%
  select(x, y, everything()) %>%
  select(-cell, -ECO_NAME)

# dataframe to raster
metric_rast <- rast(metrics, type="xyz", crs = "epsg:4326")
#save out
writeRaster(metric_rast, here::here("data/metric_rast.tiff"))

########################################
####### Get functional diversity #######
########################################

library(fundiversity)

load(here::here("data/avonet_ebird_matched.rda"))

breeding_cells <- select(breeding_mat, cell, ECO_NAME)
nonbreeding_cells <- select(nonbreeding_mat, cell, ECO_NAME)

# order species just to be sure
avonet_ebird_matched <- avonet_ebird_matched %>%
  arrange(species_code)

# create trait matrix using all available avonet traits
trait_mat <- avonet_ebird_matched %>%
  select(-c(species_code, taxon_concept_id, avonet_sciname, avonet_taxon_concept_id, avonet_family, avonet_order)) %>%
  mutate(across(where(is.numeric), ~scale(.x, center = TRUE, scale = TRUE))) %>%
  data.matrix()

# get pairwise distance matrix
dist_mat <- StatMatch::gower.dist(trait_mat) %>% dist()

# PCoA for dimension reduction
pcoa <- ade4::dudi.pco(dist_mat, scannf = FALSE, full = TRUE)
trait_axes <- pcoa$li[,1:4]
rownames(trait_axes) <- avonet_ebird_matched$species_code

# site by species matrix
breeding_mat <- breeding_mat %>%
  select(-cell) %>%
  select(order(colnames(.)))

nonbreeding_mat <- nonbreeding_mat %>%
  select(-cell) %>%
  select(order(colnames(.)))

# compute FRic
future::plan(future::multisession, workers = 4)
options(future.globals.maxSize = 400e7)

fric_df_breeding <- fd_fric(trait_axes, breeding_mat)
feve_df_breeding <- fd_feve(trait_axes, breeding_mat)
fdiv_df_breeding <- fd_fdiv(trait_axes, breeding_mat)
#fdis_df_breeding <- fd_fdis(trait_axes, breeding_mat)

fric_df_nonbreeding <- fd_fric(trait_axes, nonbreeding_mat)
feve_df_nonbreeding <- fd_feve(trait_axes, nonbreeding_mat)
fdiv_df_nonbreeding <- fd_fdiv(trait_axes, nonbreeding_mat)
#fdis_df_breeding <- fd_fdis(trait_axes, breeding_mat)

fd_metric_breeding <- bind_cols(breeding_cells, fric_df_breeding) %>%
  left_join(feve_df_breeding) %>%
  left_join(fdiv_df_breeding) %>%
  #left_join(fdis_df) %>%
  select(-site) %>%
  full_join(coords) %>%
  group_by(ECO_NAME) %>%
  mutate(across(c(FRic, FEve, FDiv), ~ quantile(.x, probs = 0.95, na.rm = TRUE), .names = "{col}_quantile")) %>%
  ungroup() %>%
  select(-cell, -ECO_NAME) %>%
  select(x, y, everything()) %>%
  rename_with(.cols = -c(x,y), ~paste0(.x, "_breeding"))

fd_metric_nonbreeding <- bind_cols(nonbreeding_cells, fric_df_nonbreeding) %>%
  left_join(feve_df_nonbreeding) %>%
  left_join(fdiv_df_nonbreeding) %>%
  #left_join(fdis_df) %>%
  select(-site) %>%
  full_join(coords) %>%
  group_by(ECO_NAME) %>%
  mutate(across(c(FRic, FEve, FDiv), ~ quantile(.x, probs = 0.95, na.rm = TRUE), .names = "{col}_quantile")) %>%
  ungroup() %>%
  select(-cell, -ECO_NAME) %>%
  select(x, y, everything()) %>%
  rename_with(.cols = -c(x,y), ~paste0(.x, "_nonbreeding"))

usethis::use_data(fd_metric_breeding)
usethis::use_data(fd_metric_nonbreeding)

# dataframe to raster
fd_rast <- c(rast(fd_metric_breeding, type="xyz", crs = "epsg:4326"),
             rast(fd_metric_nonbreeding, type="xyz", crs = "epsg:4326"))
#save out
writeRaster(fd_rast, here::here("data/fd_rast.tiff"))
