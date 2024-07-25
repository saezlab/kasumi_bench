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
    left_join(read_csv(paste0("data/DCISMIBI/Image Data/Segmetation_Outlines_and_Labels_Mendeley/", id, ".csv"),
      show_col_types = FALSE
    ), by = "label") %>%
    select("centroid-0", "centroid-1") %>%
    `colnames<-`(c("x", "y"))
})

## SM DCIS ----

repr.ids <- roc.sm.dcis <- NULL
max.auc <- 0

rocs <- c(100, 200, 300, 400, 500) %>% map_dbl(\(ws){
  plan(multisession, workers = 9)
  misty.results <- sm_train(all.cells.dcis, all.positions.dcis, 100, ws, 20, paste0("DCISexpr", ws, ".sqm"), "gaussian")

  plan(multisession, workers = 9)
  param.opt <- optimal_smclust(misty.results, resp %>% select(PointNumber, Status) %>%
    rename(id = PointNumber, target = Status) %>%
    filter(target %in% c("progressor", "nonprogressor")) %>%
    mutate(id = as.character(id), target = as.factor(target)))

  sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])


  freq.sm <- sm.repr %>%
    left_join(resp %>% select(PointNumber, Status) %>%
      mutate(PointNumber = as.character(PointNumber)), by = c("id" = "PointNumber")) %>%
    rename(target = Status) %>%
    mutate(target = as.factor(target)) %>%
    select(-id)


  roc.sm.dcis <- classify(freq.sm)

  if (roc.sm.dcis$auc > max.auc) {
    roc.sm.dcis <<- roc.sm.dcis
    repr.ids <<- repr.ids
    max.auc <<- roc.sm.dcis$auc
  }

  roc.sm.dcis$auc
})

write_rds(roc.sm.dcis, "rocs/dcis.expr.rds")

## WS DCIS ----

misty.results <- misty_train(all.cells.dcis, all.positions.dcis, 100, "DCISexpr_ws.sqm", "gaussian")

## Baseline ----

freq.dcis <- all.cells.dcis %>%
  map(freq_repr) %>%
  list_transpose() %>%
  as_tibble() %>%
  mutate(PointNumber = points) %>%
  left_join(resp %>% select(PointNumber, Status), by = "PointNumber") %>%
  select(-PointNumber) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target))

roc.fr <- classify(freq.dcis)

write_rds(roc.fr, "rocs/dcis.expr.base.rds")

## Alternatives ----

freq.banksy.ct <- banksy_labels(all.cells.dcis, all.positions.dcis, 10, 0.2)


rocs.banksy.ct <- freq.banksy.ct %>% map(\(fbc){
  freq.banksy <- fbc %>%
    left_join(
      resp %>%
        select(PointNumber, Status) %>%
        mutate(PointNumber = as.character(PointNumber)),
      by = join_by(id == PointNumber)
    ) %>%
    select(-id) %>%
    rename(target = Status) %>%
    mutate(target = as.factor(target))

  classify(freq.banksy)
})
roc.banksy.ct <- rocs.banksy.ct[[which.max(rocs.banksy.ct %>% map_dbl(~ .x$auc))]]

freq.banksy.niche <- banksy_labels(all.cells.dcis, all.positions.dcis, 10, 0.8)

rocs.banksy.niche <- freq.banksy.niche %>% map(\(fbn){
  freq.banksy <- fbn %>%
    left_join(
      resp %>%
        select(PointNumber, Status) %>%
        mutate(PointNumber = as.character(PointNumber)),
      by = join_by(id == PointNumber)
    ) %>%
    select(-id) %>%
    rename(target = Status) %>%
    mutate(target = as.factor(target))

  classify(freq.banksy)
})
roc.banksy.niche <- rocs.banksy.niche[[which.max(rocs.banksy.niche %>% map_dbl(~ .x$auc))]]

export_anndata(all.cells.dcis, all.posistions.dcis, "dcis.r5ad")

freq.cc <- read_csv("dcis_l1_k3.csv")

roc.cc.l1 <- freq.cc %>%
  left_join(
    resp %>%
      select(PointNumber, Status),
    by = join_by(sample == PointNumber)
  ) %>%
  select(-sample) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  classify()

freq.cc <- read_csv("dcis_l3_k3.csv")

roc.cc.l3 <- freq.cc %>%
  left_join(
    resp %>%
      select(PointNumber, Status),
    by = join_by(sample == PointNumber)
  ) %>%
  select(-sample) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  classify()

write_rds(
  list(
    banksy.ct = test$banksy.ct, banksy.niche = test$banksy.niche,
    cc.l1 = roc.cc.l1, cc.l3 = roc.cc.l3
  ),
  "rocs/dcis.expr.alts.rds"
)


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

# 10 neighbors as in publication, 400px = 150um window
plan(multisession, workers = 9)
misty.results <- sm_train(all.cells.lymph, all.positions.lymph, 100, 400, 20, "CTCLexpr400.sqm", "gaussian")

plan(multisession, workers = 9)
param.opt <- optimal_smclust(misty.results, outcome %>% select(-Patients) %>%
  rename(id = Spots, target = Groups) %>%
  mutate(id = as.character(id), target = as.factor(make.names(target))))

sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])

freq.sm <- sm.repr %>%
  left_join(outcome %>% mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)


roc.sm.ctcl <- classify(freq.sm)

write_rds(roc.sm.ctcl, "rocs/ctcl.expr.rds")

## WS CTCL ----

misty.results <- misty_train(all.cells.lymph, all.positions.lymph, 100, "CTCLexpr_ws.sqm", "gaussian")

## Baseline ----

freq.lymph <- all.cells.lymph %>%
  map(freq_repr) %>%
  list_transpose() %>%
  as_tibble() %>%
  mutate(Spots = spots) %>%
  left_join(outcome, by = "Spots") %>%
  select(-Spots) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target)))


roc.fr <- classify(freq.lymph)

write_rds(roc.fr, "rocs/ctcl.expr.base.rds")

## Alternatives ----

freq.banksy.ct <- banksy_labels(all.cells.lymph, all.positions.lymph, 10, 0.2)

rocs.banksy.ct <- freq.banksy.ct %>% map(\(fbc){
  freq.banksy <- fbc %>%
    left_join(
      outcome %>% mutate(Spots = as.character(Spots)),
      by = join_by(id == Spots)
    ) %>%
    select(-Patients) %>%
    rename(target = Groups) %>%
    mutate(target = as.factor(make.names(target)))

  classify(freq.banksy)
})

roc.banksy.ct <- rocs.banksy.ct[[which.max(rocs.banksy.ct %>% map_dbl(~ .x$auc))]]

freq.banksy.niche <- banksy_labels(all.cells.lymph, all.positions.lymph, 10, 0.8)

rocs.banksy.niche <- freq.banksy.niche %>% map(\(fbn){
  freq.banksy <- fbn %>%
    left_join(
      outcome %>% mutate(Spots = as.character(Spots)),
      by = join_by(id == Spots)
    ) %>%
    select(-Patients) %>%
    rename(target = Groups) %>%
    mutate(target = as.factor(make.names(target)))

  classify(freq.banksy)
})

roc.banksy.niche <- rocs.banksy.niche[[which.max(rocs.banksy.niche %>% map_dbl(~ .x$auc))]]

export_anndata(all.cells.lymph, all.posistions.lymph, "lymph.r5ad")

freq.cc <- read_csv("lymph_l1_k5.csv")

roc.cc.l1 <- freq.cc %>%
  left_join(
    outcome,
    by = join_by(sample == Spots)
  ) %>%
  select(-Patients) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  classify()

freq.cc <- read_csv("lymph_l3_k4.csv")

roc.cc.l3 <- freq.cc %>%
  left_join(
    outcome,
    by = join_by(sample == Spots)
  ) %>%
  select(-Patients) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  classify()

write_rds(
  list(
    banksy.ct = test$banksy.ct, banksy.niche = test$banksy.niche,
    cc.l1 = roc.cc.l1, cc.l3 = roc.cc.l3
  ),
  "rocs/ctcl.expr.alts.rds"
)


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
plan(multisession, workers = 9)
misty.results <- sm_train(all.cells.bc, all.positions.bc, 50, 200, 20, "BCexpr200.sqm", "gaussian")

resp <- bmeta %>%
  filter(core %in% cores) %>%
  select(core, response)
plan(multisession, workers = 9)
param.opt <- optimal_smclust(misty.results, resp %>%
  rename(id = core, target = response) %>%
  mutate(id = as.character(id), target = as.factor(make.names(target))))

sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])

freq.sm <- sm.repr %>%
  left_join(resp, by = c("id" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id)

roc.sm.bc <- classify(freq.sm)

write_rds(roc.sm.bc, "rocs/bc.expr.rds")

## WS BC ----
misty.results <- misty_train(all.cells.bc, all.positions.bc, 50, "BCexpr_ws.sqm", "gaussian")

## Baseline ----

freq.bc <- all.cells.bc %>%
  map(freq_repr) %>%
  list_transpose() %>%
  as_tibble() %>%
  mutate(core = cores) %>%
  left_join(resp, by = "core") %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-core)

roc.fr <- classify(freq.bc)
write_rds(roc.fr, "rocs/bc.expr.base.rds")

## Alternatives ----

freq.banksy.ct <- banksy_labels(all.cells.bc, all.positions.bc, 10, 0.2)

rocs.banksy.ct <- freq.banksy.ct %>% map(\(fbc){
  freq.banksy <- fbc %>%
    left_join(
      resp,
      by = join_by(id == core)
    ) %>%
    rename(target = response) %>%
    mutate(target = as.factor(make.names(target)))

  classify(freq.banksy)
})

roc.banksy.ct <- rocs.banksy.ct[[which.max(rocs.banksy.ct %>% map_dbl(~ .x$auc))]]

freq.banksy.niche <- banksy_labels(all.cells.lymph, all.positions.lymph, 10, 0.8)

rocs.banksy.niche <- freq.banksy.niche %>% map(\(fbn){
  freq.banksy <- fbn %>%
    left_join(
      resp,
      by = join_by(id == core)
    ) %>%
    rename(target = response) %>%
    mutate(target = as.factor(make.names(target)))

  classify(freq.banksy)
})

roc.banksy.niche <- rocs.banksy.niche[[which.max(rocs.banksy.niche %>% map_dbl(~ .x$auc))]]

export_anndata(all.cells.bc, all.posistions.bc, "bc.r5ad")

freq.cc <- read_csv("bc_l1_k4.csv")

roc.cc.l1 <- freq.cc %>%
  left_join(
    resp,
    by = join_by(sample == core)
  ) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  classify()

freq.cc <- read_csv("bc_l3_k2.csv")

roc.cc.l3 <- freq.cc %>%
  left_join(
    resp,
    by = join_by(sample == core)
  ) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  classify()

write_rds(
  list(
    banksy.ct = test$banksy.ct, banksy.niche = test$banksy.niche,
    cc.l1 = roc.cc.l1, cc.l3 = roc.cc.l3
  ),
  "rocs/bc.expr.alts.rds"
)
