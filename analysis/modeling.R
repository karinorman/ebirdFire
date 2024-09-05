# Models for assessing species traits relationship to high severity range ##

library(dplyr)
library(betareg)
library(broom)

load(here::here("data/avonet_ebird_matched.rda"))
load(here::here("data/percent_range_highsev_df.rda"))

traits_df <- percent_range_highsev_df %>%
  left_join(avonet_ebird_matched) %>%
  select(-c(scientific.name, taxon_concept_id, avonet_sciname,
            avonet_taxon_concept_id, avonet_family, avonet_order))

predictors <- colnames(traits_df)[!colnames(traits_df) %in% c("species_code", "resident_percent",
                                                              "breeding_percent", "nonbreeding_percent")]

fit <- glm(formula = reformulate(termlabels = predictors, response = "resident_percent"), data = traits_df)
summary(fit)

# fit <- glm(formula = "breeding_percent ~ Mass", data = traits_df)
# summary(fit)

beta_fit <- betareg(reformulate(termlabels = predictors, response = "resident_percent"), data = traits_df, link = "logit")
