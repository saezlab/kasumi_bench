source("utils.R")

# MIBI  ----

## DCIS progression ----

resp <- read_csv("data/DCISMIBI/Tissue Feature Data/Table_S1_Patient_Feature_Table.csv")

misty.results.ct <- read_rds("DCISct_ws.rds")
misty.results.expr <- read_rds("DCISexpr_ws.rds")

ws.repr.ct <- misty_labels(misty.results.ct) %>%
  left_join(
    resp %>%
      select(PointNumber, Status) %>%
      mutate(PointNumber = as.character(PointNumber)),
    by = c("sample" = "PointNumber")
  ) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  select(-sample)

ws.repr.expr <- misty_labels(misty.results.expr) %>%
  left_join(
    resp %>%
      select(PointNumber, Status) %>%
      mutate(PointNumber = as.character(PointNumber)),
    by = c("sample" = "PointNumber")
  ) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  select(-sample)



roc.ws.dcis.ct <- classify(ws.repr.ct)
roc.ws.dcis.expr <- classify(ws.repr.expr)

write_rds(list(ct = roc.ws.dcis.ct, expr = roc.ws.dcis.expr), "rocs/dcis.ws.rds")


# CODEX ----

## CTCL responders  ----

lymph <- read_csv("data/LymphomaCODEX/single_cells.csv")

outcome <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  group_by(Spots, Patients) %>%
  summarize(Groups = Groups[1], .groups = "drop")


misty.results.ct <- read_rds("CTCLct_ws.rds")
misty.results.expr <- read_rds("CTCLexpr_ws.rds")

ws.repr.ct <- misty_labels(misty.results.ct) %>% 
  left_join(outcome %>% mutate(Spots = as.character(Spots)), by = c("sample" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-sample, -Patients)

ws.repr.expr <- misty_labels(misty.results.expr) %>% 
  left_join(outcome %>% mutate(Spots = as.character(Spots)), by = c("sample" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-sample, -Patients)

roc.ws.ctcl.ct <- classify(ws.repr.ct)
roc.ws.ctcl.expr <- classify(ws.repr.expr)

write_rds(list(ct = roc.ws.ctcl.ct, expr = roc.ws.ctcl.expr), "rocs/ctcl.ws.rds")


# IMC ----

## BC responders ----

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
    arrange(core) %>%
    pull(core)
)

resp <- bmeta %>%
  filter(core %in% cores) %>%
  select(core, response)

misty.results.ct <- read_rds("BCct_ws.rds")
misty.results.expr <- read_rds("BCexpr_ws.rds")

ws.repr.ct <- misty_labels(misty.results.ct) %>% 
  left_join(resp, by = c("sample" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-sample)

ws.repr.expr <- misty_labels(misty.results.expr) %>% 
  left_join(resp, by = c("sample" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-sample)

roc.ws.bc.ct <- classify(ws.repr.ct)
roc.ws.bc.expr <- classify(ws.repr.expr)

write_rds(list(ct = roc.ws.bc.ct, expr = roc.ws.bc.expr), "rocs/bc.ws.rds")

