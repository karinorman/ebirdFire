# Models for assessing species traits relationship to high severity range ##

library(dplyr)
library(betareg)
library(broom)
library(emmeans)

load(here::here("data/avonet_ebird_matched.rda"))
load(here::here("data/species_range_metrics.rda"))
load(here::here("data/pop_mets.rda"))

traits_df <- species_range_metrics %>%
  filter(wus_percent > 0.1) %>%
  left_join(avonet_ebird_matched) %>%
  select(-c(scientific.name, taxon_concept_id, avonet_sciname,
            avonet_taxon_concept_id, avonet_family, avonet_order,
            total_area, wus_area, wus_percent, sev_area)) %>%
  left_join(pop_mets %>% select(species_code, pop_percent = percent)) %>%
  rename(range = high_sev_percent, population = pop_percent) %>%
  pivot_longer(cols = c(range, population), names_to = "percent_type", values_to = "percent")


fit_beta <- function(percent_type, model_id, formula){

  fit <- betareg(formula = formula,
                 data = traits_df %>% filter(percent_type == !!percent_type), link = "logit")

  df <- broom::tidy(fit) %>%
    mutate(percent_type = percent_type,
           model_id = model_id)

  return(list(df, fit))
}


model_df <- data.frame(percent_type = rep(c("range", "population"), 3),
                       model_id = rep(c("morphology", "habitat", "llifestyle"), each = 2),
                       formula = rep(c("percent ~ Beak.Length_Culmen + Beak.Depth + Tarsus.Length + Kipps.Distance + Secondary1 + Hand.Wing.Index + Tail.Length + log(Mass)",
                                       "percent ~ Habitat.Density + Habitat",
                                       "percent ~ Migration + Trophic.Niche + Primary.Lifestyle"), each = 2))

model_output <- purrr::pmap(model_df, fit_beta)

model_output_df <- map(model_output, 1) %>% bind_rows()
model_fits <- map(model_output_df, 2)

model_output_df %>% filter(p.value < 0.05, term != "(Intercept)", component == "mu") %>% View()



habitat_contrasts <- emmeans(habitat_fit, pairwise ~ Habitat, data = traits_df)$contrasts %>% as.data.frame()


trophic_niche_contrasts <- emmeans(lifestyle_fit, pairwise ~ Trophic.Niche, data = traits_df)$contrasts %>% as.data.frame()
lifestyle_contrasts <- emmeans(lifestyle_fit, pairwise ~ Primary.Lifestyle, data = traits_df)$contrasts %>% as.data.frame()


