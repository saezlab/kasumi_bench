source("utils.R")

# MIBI  ----
ct <- read_rds("rocs/dcis.ct.rds")
expr <- read_rds("rocs/dcis.expr.rds")
ws <- read_rds("rocs/dcis.ws.rds")

ggroc(ct, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("MIBI DCIS") + theme_classic()

ggsave("roc.mibi.pdf")

ggroc(list(sm.ct = ct$sm, ws.ct = ws$ct, sm.expr = expr, ws.expr = ws$expr), legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("MIBI DCIS") + theme_classic()

ggsave("roc.misty.mibi.pdf")

misty.results <- read_rds("DCISct.rds")
sm.repr <- sm_labels(misty.results, 0.3, 0.6)

repr.ids <- sm.repr %>% map_chr(~ .x$id[1])

resp <- read_csv("data/DCISMIBI/Tissue Feature Data/Table_S1_Patient_Feature_Table.csv")

freq.sm <- sm.repr %>%
  map_dfr(~ .x %>%
            select(-c(id, x, y)) %>%
            freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.))))) %>%
  add_column(id = repr.ids) %>%
  left_join(resp %>% select(PointNumber, Status) %>% mutate(PointNumber = as.character(PointNumber)), by = c("id" = "PointNumber")) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  select(-id)

model_reliance(freq.sm)
  

misty.results <- read_rds("DCISexpr.rds")

sm.repr <- sm_labels(misty.results, 0.3, 0.8)

repr.ids <- sm.repr %>% map_chr(~ .x$id[1])

freq.sm <- sm.repr %>%
  map_dfr(~ .x %>%
            select(-c(id, x, y)) %>%
            freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.))))) %>%
  add_column(id = repr.ids) %>%
  left_join(resp %>% select(PointNumber, Status) %>% mutate(PointNumber = as.character(PointNumber)), by = c("id" = "PointNumber")) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  select(-id)

repr_model_stats(freq.sm)


# CODEX ----

ct <- read_rds("rocs/ctcl.ct.rds")
expr <- read_rds("rocs/ctcl.expr.rds")
ws <- read_rds("rocs/ctcl.ws.rds")

lymph <- read_csv("data/LymphomaCODEX/single_cells.csv")

outcome <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  group_by(Spots, Patients) %>%
  summarize(Groups = Groups[1], .groups = "drop") %>%
  mutate(Groups = ifelse(Groups == 1, "responder", "nonresponder"))

ggroc(ct, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("CODEX CTCL") + theme_classic()

ggsave("roc.codex.pdf")

ggroc(list(sm.ct = ct$sm, ws.ct = ws$ct, sm.expr = expr, ws.expr = ws$expr), legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("CODEX CTCL") + theme_classic()

ggsave("roc.misty.codex.pdf")


misty.results  <- read_rds("CTCLct.rds")
sm.repr <- sm_labels(misty.results, 0.6, 0.6)

repr.ids <- sm.repr %>% map_chr(~ .x$id[1])


freq.sm <- sm.repr %>%
  map_dfr(~ .x %>%
            select(-c(id, x, y)) %>%
            freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.))))) %>%
  add_column(id = repr.ids) %>%
  left_join(outcome %>% mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(target)) %>%
  select(-id, -Patients)

model_reliance(freq.sm)

# IMC ----

ct <- read_rds("rocs/bc.ct.rds")
expr <- read_rds("rocs/bc.expr.rds")
ws <- read_rds("rocs/bc.ws.rds")

bmeta <- read_csv("data/BCIMC/Basel_PatientMetadata.csv")
with_seed(
  1,
  cores <- bmeta %>%
    filter(
      diseasestatus == "tumor",
      response %in% c("Sensitive", "Resistant"),
      clinical_type == "HR+HER2-", Subtype == "PR+ER+"
    ) %>%
    group_by(response) %>%
    slice_sample(n = 15) %>%
    pull(core)
)
resp <- bmeta %>% filter(core %in% cores) %>% select(core, response)

ggroc(ct, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("IMC BC") + theme_classic()

ggsave("roc.imc.pdf")

ggroc(list(sm.ct = ct$sm, ws.ct = ws$ct, sm.expr = expr, ws.expr = ws$expr), legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("IMC BC") + theme_classic()

ggsave("roc.misty.imc.pdf")

misty.results <- read_rds("BCct.rds")

sm.repr <- sm_labels(misty.results, 0.7, 0.5)

repr.ids <- sm.repr %>% map_chr(~ .x$id[1])

freq.sm <- sm.repr %>%
  map_dfr(~ .x %>%
            select(-c(id, x, y)) %>%
            freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.))))) %>%
  add_column(id = repr.ids) %>%
  left_join(resp, by = c("id" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id)

model_reliance(freq.sm)

# Sensitivity ----

sens <- read_csv("sensitivity.csv") 
  
sens.long <- sens %>% pivot_longer(-c(cut,res), names_to = "Data", values_to = "AUC") %>% 
  mutate(across(c(cut,res), as.factor)) %>% rename(Threshold = cut, Resolution = res)

ggplot(sens.long, aes(y = Resolution, x = Threshold, fill = AUC)) + 
  geom_tile() + 
  scale_fill_gradient(low = "gray50", high = "tomato3", limits = c(0.5,1)) +
  facet_wrap("Data", ncol = 2) +
  coord_equal() +
  theme_classic()

ggsave("sensitivity.pdf")


# number of clusters in the sm representation

# number of singnificant clsters after regression