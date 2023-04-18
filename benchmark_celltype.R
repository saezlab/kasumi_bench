source("utils.R")


plan(multisession, workers = 5)


# MIBI  ----

## DCIS progression ----

dcis <- read_csv("data/DCISMIBI/Single Cell Data/Single_Cell_Data.csv")

# can also be dcis$manual_gating_subineage %>% unique()
cts <- dcis$sublineage %>% unique()
resp <- read_csv("data/DCISMIBI/Tissue Feature Data/Table_S1_Patient_Feature_Table.csv")

points <- dcis %>%
  filter(Tissue_Type == "DCIS") %>%
  pull(Point_Num) %>%
  unique() %>%
  intersect(resp %>% filter(Status %in% c("progressor", "nonprogressor")) %>% pull(PointNumber))

all.cells.dcis <- points %>% map(\(id){
  dcis %>%
    filter(Point_Num == id) %>%
    pull(sublineage) %>%
    map(~ .x == cts) %>%
    rlist::list.rbind() %>%
    `colnames<-`(make.names(cts)) %>%
    as_tibble()
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

misty.results <- sm_train(all.cells.dcis, all.positions.dcis, 10, 150, 20, 2, "DCISct")

param.opt <- optimal_smclust(misty.results, resp %>% select(PointNumber, Status) %>%
  rename(id = PointNumber, target = Status) %>%
  filter(target %in% c("progressor", "nonprogressor")) %>%
  mutate(id = as.character(id), target = as.factor(target)))

sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])
# sm.repr <- sm_labels(misty.results, cuts = 0.5, res = 0.6)

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

# 0.779
roc.sm <- classify(freq.sm)

## WS DCIS ----

misty.results <- misty_train(all.cells.dcis, all.positions.dcis, 10, "DCISct")

## Alternatives DCIS ----

cn.results <- cn_train(all.cells.dcis, all.positions.dcis)
cn.repr <- cn_labels(cn.results, k = 17)

repr.ids.cn <- cn.repr %>% map_chr(~ .x$id[1])

freq.cn <- cn.repr %>%
  map_dfr(~ .x %>%
    select(-c(id, x, y)) %>%
    freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3))) %>%
  add_column(id = repr.ids.cn) %>%
  left_join(resp %>% select(PointNumber, Status) %>% mutate(PointNumber = as.character(PointNumber)), by = c("id" = "PointNumber")) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  select(-id)

# 0.502
roc.cn <- classify(freq.cn)

freq.dcis <- all.cells.dcis %>%
  map(freq_repr) %>%
  list_transpose() %>%
  as_tibble() %>%
  mutate(PointNumber = points) %>%
  left_join(resp %>% select(PointNumber, Status), by = "PointNumber") %>%
  filter(PointNumber %in% repr.ids) %>%
  select(-PointNumber) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target))

# 0.448
roc.fr <- classify(freq.dcis)

markov.dcis <- all.cells.dcis %>%
  map2(all.positions.dcis, markov_repr) %>%
  list_transpose() %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(PointNumber = points) %>%
  left_join(resp %>% select(PointNumber, Status), by = "PointNumber") %>%
  filter(PointNumber %in% repr.ids) %>%
  select(-PointNumber) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target))

# 0.516
roc.pa <- classify(markov.dcis)

write_rds(list(sm = roc.sm, cn = roc.cn, fr = roc.fr, pa = roc.pa), "rocs/dcis.ct.rds")

ggroc(list(sm = roc.sm, cn = roc.cn, fr = roc.fr, pa = roc.pa), legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("MIBI DCIS") + theme_classic()

ggsave("roc.mibi.pdf")


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

cts <- lymph %>%
  pull(ClusterName) %>%
  unique()

all.cells.lymph <- spots %>% map(\(id){
  lymph %>%
    filter(Spots == id) %>%
    pull(ClusterName) %>%
    map(~ .x == cts) %>%
    rlist::list.rbind() %>%
    `colnames<-`(make.names(cts)) %>%
    as_tibble()
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
misty.results <- sm_train(all.cells.lymph, all.positions.lymph, 10, 200, 20, 2, "CTCLct")

param.opt <- optimal_smclust(misty.results, outcome %>% select(-Patients) %>%
  rename(id = Spots, target = Groups) %>%
  mutate(id = as.character(id), target = as.factor(make.names(target))), minsamp = 0.2)


sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])
# sm.repr <- sm_labels(misty.results, cuts = 0.6, res = 0.9)

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

# 0.932
roc.sm <- classify(freq.sm)

## WS CTCL ----

misty.results <- misty_train(all.cells.lymph, all.positions.lymph, 10, "CTCLct")

## Alternatives CTCL ----

cn.results <- cn_train(all.cells.lymph, all.positions.lymph)
cn.repr <- cn_labels(cn.results, k = 10)

repr.ids.cn <- cn.repr %>% map_chr(~ .x$id[1])

freq.cn <- cn.repr %>%
  map_dfr(~ .x %>%
    select(-c(id, x, y)) %>%
    freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3))) %>%
  add_column(id = repr.ids.cn) %>%
  filter(id %in% repr.ids) %>%
  left_join(outcome %>% mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)

# 0.625
roc.cn <- classify(freq.cn)


freq.lymph <- all.cells.lymph %>%
  map(freq_repr) %>%
  list_transpose() %>%
  as_tibble() %>%
  mutate(Spots = spots) %>%
  left_join(outcome, by = "Spots") %>%
  filter(Spots %in% repr.ids) %>%
  select(-Spots) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target)))


# 0.4875
roc.fr <- classify(freq.lymph)

markov.lymph <- all.cells.lymph %>%
  map2(all.positions.lymph, markov_repr) %>%
  list_transpose() %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(Spots = spots) %>%
  left_join(outcome, by = "Spots") %>%
  filter(Spots %in% repr.ids) %>%
  select(-Spots) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target)))

# 0.658
roc.pa <- classify(markov.lymph)

write_rds(list(sm = roc.sm, cn = roc.cn, fr = roc.fr, pa = roc.pa), "rocs/ctcl.ct.rds")

ggroc(list(sm = roc.sm, cn = roc.cn, fr = roc.fr, pa = roc.pa), legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("CODEX CTCL") + theme_classic()

ggsave("roc.codex.pdf")

# IMC ----

## BC responders ----
bmeta <- read_csv("data/BCIMC/Basel_PatientMetadata.csv")
cores <- bmeta %>%
  filter(
    diseasestatus == "tumor",
    response %in% c("Sensitive", "Resistant"),
    clinical_type == "HR+HER2-", Subtype == "PR+ER+"
  ) %>%
  pull(core)
bclusters <- read_csv("data/BCIMC/Cluster_labels/Basel_metaclusters.csv") %>%
  mutate(core = str_remove(id, "_\\d+$")) %>%
  filter(core %in% cores)
scloc <- read_csv("data/BCIMC/singlecell_locations/Basel_SC_locations.csv") %>%
  filter(core %in% cores)

cts <- bclusters %>%
  pull(cluster) %>%
  unique()

all.cells.bc <- cores %>% map(\(clus){
  bclusters %>%
    filter(core == clus) %>%
    pull(cluster) %>%
    map(~ .x == cts) %>%
    rlist::list.rbind() %>%
    `colnames<-`(make.names(cts)) %>%
    as_tibble()
})

names(all.cells.bc) <- cores

all.positions.bc <- cores %>% map(\(clus){
  scloc %>%
    filter(core == clus) %>%
    select(Location_Center_X, Location_Center_Y) %>%
    `colnames<-`(c("x", "y"))
})


## SM BC ----

misty.results <- sm_train(all.cells.bc, all.positions.bc, 10, 100, 20, 2, "BCct")

resp <- bmeta %>%
  filter(
    diseasestatus == "tumor",
    response %in% c("Sensitive", "Resistant"),
    clinical_type == "HR+HER2-", Subtype == "PR+ER+"
  ) %>%
  select(core, response)

param.opt <- optimal_smclust(misty.results, resp %>%
  rename(id = core, target = response) %>%
  mutate(id = as.character(id), target = as.factor(make.names(target))))

sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])
# sm.repr <- sm_labels(misty.results, cuts = 0.5, res = 0.8)

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

roc.sm <- classify(freq.sm)

## WS BC ----

misty.results <- misty_train(all.cells.bc, all.positions.bc, 10, "BCct")


## Alternatives BC ----

cn.results <- cn_train(all.cells.bc, all.positions.bc)
cn.repr <- cn_labels(cn.results, k = 6)

repr.ids.cn <- cn.repr %>% map_chr(~ .x$id[1])

freq.cn <- cn.repr %>%
  map_dfr(~ .x %>%
    select(-c(id, x, y)) %>%
    freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3))) %>%
  add_column(id = repr.ids.cn) %>%
  filter(id %in% repr.ids) %>%
  left_join(resp, by = c("id" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id)

# 0.4475
roc.cn <- classify(freq.cn)

freq.bc <- all.cells.bc %>%
  map(freq_repr) %>%
  list_transpose() %>%
  as_tibble() %>%
  mutate(core = cores) %>%
  filter(core %in% repr.ids) %>%
  left_join(resp, by = "core") %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-core)

# 0.509
roc.fr <- classify(freq.bc)

markov.bc <- all.cells.bc %>%
  map2(all.positions.bc, markov_repr) %>%
  list_transpose() %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(core = cores) %>%
  filter(core %in% repr.ids) %>%
  left_join(resp, by = "core") %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-core)

# 0.658
roc.pa <- classify(markov.bc)

write_rds(list(sm = roc.sm, cn = roc.cn, fr = roc.fr, pa = roc.pa), "rocs/bc.ct.rds")

ggroc(list(sm = roc.sm, cn = roc.cn, fr = roc.fr, pa = roc.pa), legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("IMC BC") + theme_classic()

ggsave("roc.imc.pdf")
