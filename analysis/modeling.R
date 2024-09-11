# Models for assessing species traits relationship to high severity range ##

library(dplyr)
library(betareg)
library(broom)

load(here::here("data/avonet_ebird_matched.rda"))
load(here::here("data/species_range_metrics.rda"))

traits_df <- species_range_metrics %>%
  filter(wus_percent > 0.1) %>%
  left_join(avonet_ebird_matched) %>%
  select(-c(scientific.name, taxon_concept_id, avonet_sciname,
            avonet_taxon_concept_id, avonet_family, avonet_order,
            total_range_area, wus_range_area, wus_percent, range_type))

predictors <- colnames(traits_df)[!colnames(traits_df) %in% c("species_code", "resident_percent",
                                                              "breeding_percent", "nonbreeding_percent", "Mass")] %>%
  c("log(Mass)",.)

#beta_fit <- betareg(reformulate(termlabels = predictors[!predictors %in% c("Trophic.Niche")], response = "high_sev_percent"), data = traits_df, link = "logit")

morphology_fit <- betareg(high_sev_percent ~ Beak.Length_Culmen + Beak.Depth + Tarsus.Length+
                            Kipps.Distance + Secondary1 + Hand.Wing.Index + Tail.Length + log(Mass),
                          data = traits_df, link = "logit")

#mass_fit <- betareg(high_sev_percent ~ log(Mass), data = traits_df, link = "logit")

habitat_fit <- betareg(high_sev_percent ~ Habitat.Density + Habitat, data = traits_df, link = "logit")
habitat_contrasts <- emmeans(habitat_fit, pairwise ~ Habitat, data = traits_df)$contrasts %>% as.data.frame()


#lifestyle_fit <- betareg(high_sev_percent ~ as.factor(Migration) + Trophic.Niche + Primary.Lifestyle, data = traits_df, link = "logit")
lifestyle_fit <- betareg(high_sev_percent ~ Migration + Trophic.Niche + Primary.Lifestyle, data = traits_df, link = "logit")
trophic_niche_contrasts <- emmeans(lifestyle_fit, pairwise ~ Trophic.Niche, data = traits_df)$contrasts %>% as.data.frame()
lifestyle_contrasts <- emmeans(lifestyle_fit, pairwise ~ Primary.Lifestyle, data = traits_df)$contrasts %>% as.data.frame()
#migration_contrasts <- emmeans(lifestyle_fit, pairwise ~ as.factor(Migration), data = traits_df)$contrasts %>% as.data.frame()


