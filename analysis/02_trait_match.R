# Script to pair AVONET traits with ebird species

# get species list
load(here::here("data/species_list.rda"))

# and taxonomy table
tax2022 <- read.csv(here::here("raw_data/eBird_Clements_v2022.csv"))
tax2023 <- read.csv(here::here("raw_data/eBird-Clements-v2023-integrated-checklist-October-2023.csv")) #%>%
  #filter(species_code %in% species_list)

species_list <- data.frame(species_code = species_list) %>% left_join(tax2023) %>%
  filter(species_code != "ECO_NAME")

# need to map 2023 and 2022 ebird taxonomy to 2021
# get trait data
avonet <- read.csv(here::here("raw_data/AVONET2_eBird.csv"))



avonet_ebird <- species_list %>%
  select(scientific.name, taxon_concept_id) %>%
  mutate(taxon_concept_id = toupper(taxon_concept_id)) %>%
  left_join(avonet, by = c("taxon_concept_id" = "Avibase.ID2")) %>%
  mutate(Avibase.ID2 = ifelse(!is.na(Species2), taxon_concept_id, NA))

missing_sci_name <- avonet_ebird %>% filter(is.na(Species2))

match_sci_name <- missing_sci_name %>%
  select(scientific.name, taxon_concept_id) %>%
  left_join(avonet, by = c("scientific.name" = "Species2")) %>%
  mutate(Species2 = ifelse(!is.na(Avibase.ID2), scientific.name, NA))

missing_sci_name_match <- match_sci_name %>% filter(is.na(Avibase.ID2))

# two more missing that we'll match by hand based on the ebird documentation of
# differences between 2021 and 2023 taxonomy
# https://www.birds.cornell.edu/clementschecklist/introduction/updateindex/october-2022/updates-corrections-october-2022/
missing_sci_name_match %>%
  select(scientific.name, taxon_concept_id) %>%
  mutate(Species2 = case_when(
    scientific.name == "Accipiter atricapillus" ~ "Accipiter gentilis",
    scientific.name == "Sturnella lilianae" ~ "Sturnella magna",
    .default = NA
  )) %>%
  left_join(avonet, by = "Species2")
