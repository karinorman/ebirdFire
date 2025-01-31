library(dplyr)

# list of T&E species from fish and wildlife
# https://ecos.fws.gov/ecp/report/species-listings-by-tax-group?statusCategory=Listed&groupName=All%20Animals

te_list <- read.csv(here::here("raw_data/species-listings-by-tax-group-report-1-25.csv"),
                    col.names = c("sci_name", "common_name", "where_listed", "region", "status", "group")) %>%
  filter(group == "Birds",
         region %in% c(1,8,6,2))

# get ebird taxonomy
tax2023 <- read.csv(here::here("raw_data/eBird-Clements-v2023-integrated-checklist-October-2023.csv"))
load(here::here("data/species_list.rda"))


ebird_join <- te_list %>%
  select(sci_name) %>%
  distinct() %>%
  left_join(tax2023, by = c("sci_name" = "scientific.name")) %>%
  mutate(species_code = case_when(
    sci_name == "Charadrius nivosus nivosus" ~ "snoplo3",
    sci_name == "Melozone crissalis eremophilus" ~ "caltow5",
    sci_name == "Gallicolumba stairi" ~ "frgdov1",
    sci_name == "Acrocephalus luscinia" ~ "nigrew1",
    sci_name == "Branta (=Nesochen) sandvicensis" ~ "hawgoo",
    .default = species_code
  )) %>% select(sci_name, species_code)

runs2022 <- read.csv(here::here("raw_data/ebirdst_runs_2022.csv")) %>%
  # remove example raster
  filter(species_code != "yebsap-example")

runs2021 <- read.csv(here::here("raw_data/ebirdst_runs_2021.csv"))


runs2021 %>% filter(species_code %in% ebird_join$species_code)
runs2022 %>% filter(species_code %in% ebird_join$species_code)



#
# ebird_join <- te_list %>%
#   select(sci_name) %>%
#   distinct() %>%
#   left_join(avonet_ebird_matched %>%
#               select(species_code, scientific.name), by = c("sci_name" = "scientific.name")) %>%
#   select(sci_name, species_code) %>%
#   # some birds clearly don't occur in continental US, mostly Hawaiian
#   filter(!sci_name %in% c("Corvus hawaiiensis", "Palmeria dolei", "Hemignathus wilsoni",
#                           "Gallirallus owstoni", "Gallicolumba stairi", "Telespiza ultima",
#                           "Corvus kubaryi", "Telespiza cantans", "Chasiempis ibidis",
#                           "Anas laysanensis", "Myadestes lanaiensis rutha", "Loxops coccineus",
#                           "Loxioides bailleui", "Acrocephalus luscinia", "Aerodramus bartschi",
#                           "Psittirostra psittacea", "Fulica alai", "Drepanis coccinea", "Himantopus mexicanus knudseni",
#                           "Paroreomyza maculata", "Gallinula galeata sandvicensis", "Loxops mana", "Gymnomyza samoensis",
#
#
#                           )) %>%
#   # have to hand match the rest
#   mutate(species_code = case_when(
#     sci_name == "Falco femoralis septentrionalis" ~ "aplfal",
#     sci_name == "Glaucidium brasilianum cactorum" ~ "fepowl",
#     sci_name == "Lanius ludovicianus mearnsi" ~ "logshr",
#     sci_name == "Setophaga chrysoparia" ~ "gchwar",
#     sci_name == "Tympanuchus cupido attwateri" ~ "grpchi",
#     sci_name == "Grus americana" ~ "whocra",
#     sci_name == "Rallus obsoletus obsoletus" ~ "ridrai1",
#     sci_name == "Polioptila californica californica" ~ "calgna",
#     # not sure if this one is specific enough
#     sci_name == "Melozone crissalis eremophilus" ~ "caltow5",
#     # also this one
#     sci_name == "Sternula antillarum browni" ~ "leater5",
#
#
#
#   ))
#
#
# avonet_ebird_matched %>% filter(scientific.name %in% c(avonet_ebird_matched$scientific.name, avonet_ebird_matched$avonet_sciname))
