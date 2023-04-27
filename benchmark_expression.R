source("utils.R")

options(future.globals.maxSize = 2024^3)

# MIBI  ----

## DCIS progression ----

dcis <- read_csv("data/DCISMIBI/Single Cell Data/Single_Cell_Data.csv")

resp <- read_csv("data/DCISMIBI/Tissue Feature Data/Table_S1_Patient_Feature_Table.csv")
panel <- read_csv("data/DCISMIBI/panel.csv") %>% colnames()

points <- dcis %>%
  filter(Tissue_Type == "DCIS") %>%
  pull(Point_Num) %>%
  unique() %>%
  intersect(resp %>% filter(Status %in% c("progressor", "nonprogressor")) %>% pull(PointNumber))

all.cells.dcis <- points %>% map(\(id){
  dcis %>%
    filter(Point_Num == id) %>%
    select(all_of(panel))
})

names(all.cells.dcis) <- points

all.positions.dcis <- points %>% map(\(id){
  dcis %>%
    filter(Point_Num == id) %>%
    select(label) %>%
    left_join(read_csv(paste0("data/DCISMIBI/Image Data/Segmetation_Outlines_and_Labels_Mendeley/", id, ".csv"), show_col_types = FALSE), by = "label") %>%
    select("centroid-0", "centroid-1") %>%
    `colnames<-`(c("x", "y"))
})

## SM DCIS ----

plan(multisession, workers = 9)
misty.results <- sm_train(all.cells.dcis, all.positions.dcis, 100, 200, 20, 2, "DCISexpr", "gaussian")

plan(multisession, workers = 9)
param.opt <- optimal_smclust(misty.results, resp %>% select(PointNumber, Status) %>%
  rename(id = PointNumber, target = Status) %>%
  filter(target %in% c("progressor", "nonprogressor")) %>%
  mutate(id = as.character(id), target = as.factor(target)))

sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])

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


roc.sm.dcis <- classify(freq.sm)

write_rds(roc.sm.dcis, "rocs/dcis.expr.rds")

## WS DCIS ----

misty.results <- misty_train(all.cells.dcis, all.positions.dcis, 100, "DCISexpr", "gaussian")

# CODEX ----

## CTCL responders  ----

lymph <- read_csv("data/LymphomaCODEX/single_cells.csv")
spots <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  pull(Spots) %>%
  unique()
outcome <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  group_by(Spots, Patients) %>%
  summarize(Groups = Groups[1], .groups = "drop")

panel <- read_csv("data/LymphomaCODEX/panel.csv") %>% colnames()

all.cells.lymph <- spots %>% map(\(id){
  lymph %>%
    filter(Spots == id) %>%
    select(all_of(panel)) %>%
    rename_with(make.names)
})

names(all.cells.lymph) <- spots

all.positions.lymph <- spots %>% map(\(id){
  lymph %>%
    filter(Spots == id) %>%
    select(X, Y) %>%
    `colnames<-`(c("x", "y"))
})

## SM CTCL ----

# 10 neighbors as in publication, 200px = 75um window
plan(multisession, workers = 9)
misty.results <- sm_train(all.cells.lymph, all.positions.lymph, 100, 200, 20, 2, "CTCLexpr", "gaussian")

plan(multisession, workers = 9)
param.opt <- optimal_smclust(misty.results, outcome %>% select(-Patients) %>%
  rename(id = Spots, target = Groups) %>%
  mutate(id = as.character(id), target = as.factor(make.names(target))))

sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])

repr.ids <- sm.repr %>% map_chr(~ .x$id[1])

freq.sm <- sm.repr %>%
  map_dfr(~ .x %>%
    select(-c(id, x, y)) %>%
    freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.))))) %>%
  add_column(id = repr.ids) %>%
  left_join(outcome %>% mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)


roc.sm.ctcl <- classify(freq.sm)

write_rds(roc.sm.ctcl, "rocs/ctcl.expr.rds")

## WS CTCL ----

misty.results <- misty_train(all.cells.lymph, all.positions.lymph, 100, "CTCLexpr", "gaussian")

# IMC ----

## BC responders ----
bmeta <- read_csv("data/BCIMC/Basel_PatientMetadata.csv")

panel <- read_csv("data/BCIMC/Basel_Zuri_StainingPanel.csv") %>%
  select(Target, FullStack) %>%
  filter(complete.cases(.)) %>%
  distinct()

wi <- read_csv("data/BCIMC/Basel_Zuri_WholeImage.csv")

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


obj.ids <- bmeta %>%
  filter(core %in% cores) %>%
  select(core, FileName_FullStack) %>%
  mutate(fp = paste0(
    str_extract(FileName_FullStack, "(.+?_){7}"), "[0-9]+_",
    str_extract(core, "X.+")
  )) %>%
  arrange(core) %>%
  pull(fp) %>%
  map_int(~ (wi %>% pluck("ImageNumber", str_which(wi %>% pull(FileName_FullStack), .x))))


all.objs <- obj.ids %>% map(\(id){
  path <- paste0("data/BCIMC/single_images/ ", id, ".csv")
  data <- read_csv(path, show_col_types = FALSE)

  intensities <- data %>% select(contains("Intensity_MeanIntensity_FullStack"))

  markers <- tibble(channel = colnames(intensities) %>% str_extract("[0-9]+") %>% as.numeric()) %>%
    left_join(panel, by = c("channel" = "FullStack"))

  to.remove <- which(markers$channel < 9 | is.na(markers$Target) |
    markers$channel > 47 | markers$channel %in% c(26, 32, 36, 42))

  pos <- data %>% select(Location_Center_X, Location_Center_Y)

  pos.complete <- which(complete.cases(pos))

  expr <- intensities %>%
    select(-all_of(to.remove)) %>%
    slice(pos.complete)
  colnames(expr) <- markers %>%
    slice(-to.remove) %>%
    pull(Target) %>%
    make.names()
  pos <- pos %>% slice(pos.complete)

  list(expr = expr, pos = pos)
})

all.cells.bc <- all.objs %>% map(~ .x[["expr"]])
names(all.cells.bc) <- cores
all.positions.bc <- all.objs %>% map(~ .x[["pos"]])

## SM BC ----

misty.results <- sm_train(all.cells.bc, all.positions.bc, 50, 100, 20, 2, "BCexpr", "gaussian")

resp <- bmeta %>%
  filter(core %in% cores) %>%
  select(core, response)

param.opt <- optimal_smclust(misty.results, resp %>%
  rename(id = core, target = response) %>%
  mutate(id = as.character(id), target = as.factor(make.names(target))))

sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])

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

roc.sm.bc <- classify(freq.sm)

write_rds(roc.sm.bc, "rocs/bc.expr.rds")

## WS BC ----
misty.results <- misty_train(all.cells.bc, all.positions.bc, 50, "BCexpr", "gaussian")
