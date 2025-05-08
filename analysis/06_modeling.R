# Models for assessing species traits relationship to high severity range ##

library(dplyr)
library(betareg)
library(broom)
library(emmeans)
library(tidyr)
library(purrr)

load(here::here("data/avonet_ebird_matched.rda"))
load(here::here("data/species_range_metrics.rda"))
load(here::here("data/pop_mets.rda"))

traits_df <- species_range_metrics %>%
  #filter(wus_percent > 0.1) %>%
  left_join(avonet_ebird_matched) %>%
  select(-c(scientific.name, taxon_concept_id, avonet_sciname,
            avonet_taxon_concept_id, avonet_family, avonet_order,
            total_area, wus_area, wus_percent, sev_area)) %>%
  left_join(pop_mets %>% select(species_code, ends_with ("_percent"), type))

model_data <- traits_df %>%
  filter(forest_percent > 0.1) %>%
  rename(range = high_sev_percent, population = sev_forest_percent) %>%
  pivot_longer(cols = c(range, population), names_to = "percent_type", values_to = "percent")


fit_beta <- function(percent_type, model_id, formula, type){

  if (type == "breeding"){
    df <- model_data %>% filter(type %in% c("breeding", "resident"))
  } else {
    df <- model_data %>% filter(type %in% c("nonbreeding", "resident"))
  }

  fit <- betareg(formula = formula,
                 data = df %>% filter(percent_type == !!percent_type), link = "logit")

  output_df <- broom::tidy(fit) %>%
    mutate(percent_type = percent_type,
           model_id = model_id,
           type = type)

  return(list(output_df, fit))
}


model_df <- data.frame(percent_type = rep(c("range", "population"), 3),
                       model_id = rep(c("morphology", "habitat", "llifestyle"), each = 2),
                       formula = rep(c("percent ~ Beak.Length_Culmen + Beak.Depth + Tarsus.Length + Kipps.Distance + Secondary1 + Hand.Wing.Index + Tail.Length + log(Mass)",
                                       "percent ~ Habitat.Density + Habitat",
                                       "percent ~ Migration + Trophic.Niche + Primary.Lifestyle"), each = 2))

model_df <- bind_rows(model_df %>% mutate(type = "breeding"),
                      model_df %>% mutate(type = "nonbreeding"))

model_output <- purrr::pmap(model_df[1:10,], fit_beta)

model_output_df <- map(model_output, 1) %>% bind_rows()
model_fits <- map(model_output_df, 2)

# shorter beaks are more likely to be in high severity, and birds that like high density
model_output_df %>% filter(p.value < 0.05, term != "(Intercept)", component %in% c("mu", "mean"), percent_type == "population") %>% View()

model_output_supp <- model_output_df %>%
  filter(percent_type == "population") %>%
  select(-percent_type) %>%
  mutate(across(where(~is.numeric(.x)), round, 2))

readr::write_csv(model_output_supp, here::here("figures/model_output_supp.csv"))

# habitat_contrasts <- emmeans(habitat_fit, pairwise ~ Habitat, data = traits_df)$contrasts %>% as.data.frame()
#
#
# trophic_niche_contrasts <- emmeans(lifestyle_fit, pairwise ~ Trophic.Niche, data = traits_df)$contrasts %>% as.data.frame()
# lifestyle_contrasts <- emmeans(lifestyle_fit, pairwise ~ Primary.Lifestyle, data = traits_df)$contrasts %>% as.data.frame()


#######################
###### Histogram ######
#######################
library(ggplot2)

percent_hist <- traits_df %>%
  mutate(percent = total_percent * 100) %>%
  ggplot(aes(percent)) +
  geom_histogram(color = "#000000", fill = "#82A6B1") +
  theme_classic() +
  xlab("Percent of Global Population") +
  ylab("Species Count") +
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = NA, color = NA),
        legend.box.background = element_rect(fill = NA, color = NA),
        legend.key = element_rect(fill = "transparent"),
        legend.box = element_blank()
  )

ggsave(here::here("figures/percent_pop_hist.png"), percent_hist, width = 17, height = 10, bg = "transparent")
