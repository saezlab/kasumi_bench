source("utils.R")

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

repr.ids <- roc.sm <- NULL
max.auc <- 0

rocs <- c(100, 200, 300, 400, 500) %>% map_dbl(\(ws){
  plan(multisession, workers = 8)

  misty.results <- sm_train(all.cells.dcis, all.positions.dcis, 10, ws, 20, paste0("DCISct", ws, ".sqm"))

  plan(multisession, workers = 8)

  param.opt <- optimal_smclust(misty.results, resp %>% select(PointNumber, Status) %>%
    rename(id = PointNumber, target = Status) %>%
    filter(target %in% c("progressor", "nonprogressor")) %>%
    mutate(id = as.character(id), target = as.factor(target)))

  sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])

  repr.ids <- sm.repr %>% pull(id)

  freq.sm <- sm.repr %>%
    left_join(resp %>% select(PointNumber, Status) %>%
      mutate(PointNumber = as.character(PointNumber)), by = c("id" = "PointNumber")) %>%
    rename(target = Status) %>%
    mutate(target = as.factor(target)) %>%
    select(-id)

  roc.sm <- classify(freq.sm)

  if (roc.sm$auc > max.auc) {
    roc.sm <<- roc.sm
    repr.ids <<- repr.ids
    max.auc <<- roc.sm$auc
  }

  roc.sm$auc
})

write_rds(rocs, "rocs/wrocs.dcis.rds")

## WS DCIS ----

misty.results <- misty_train(all.cells.dcis, all.positions.dcis, 10, "DCISct_ws.sqm")

## Alternatives DCIS ----

cn.results <- cn_train(all.cells.dcis, all.positions.dcis, 10)
cn.repr <- cn_labels(cn.results, k = 17)

freq.cn <- cn.repr %>%
  filter(id %in% repr.ids) %>%
  left_join(resp %>% select(PointNumber, Status) %>%
    mutate(PointNumber = as.character(PointNumber)), by = c("id" = "PointNumber")) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  select(-id)

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

roc.fr <- classify(freq.dcis)

pa.dcis <- all.cells.dcis %>%
  map2(all.positions.dcis, pa_repr) %>%
  list_transpose() %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(PointNumber = points) %>%
  left_join(resp %>% select(PointNumber, Status), by = "PointNumber") %>%
  filter(PointNumber %in% repr.ids) %>%
  select(-PointNumber) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target))

roc.pa <- classify(pa.dcis)

csea.dcis <- all.cells.dcis %>%
  map2(all.positions.dcis, csea_repr) %>%
  list_transpose() %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(PointNumber = points) %>%
  left_join(resp %>% select(PointNumber, Status), by = "PointNumber") %>%
  filter(PointNumber %in% repr.ids) %>%
  select(-PointNumber) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target))

roc.csea <- classify(csea.dcis)

wc.repr <- bin_count_cluster(all.cells.dcis, all.positions.dcis, 200, 50, 33)
freq.wc <- wc.repr %>%
  left_join(resp %>% select(PointNumber, Status) %>%
    mutate(PointNumber = as.character(PointNumber)), by = c("id" = "PointNumber")) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  select(-id)

roc.wc <- classify(freq.wc)

write_rds(list(sm = roc.sm, cn = roc.cn, fr = roc.fr, pa = roc.pa, csea = roc.csea, wc = roc.wc), "rocs/dcis.ct.rds")

## ARI ----

misty.results <- read_rds("DCISct200.rds")
sm.repr <- sm_repr(misty.results, cuts = 0.9, res = 0.6)
repr.ids <- sm.repr %>% map_chr(~ .x$id[1])

sm.freq <- sm.repr %>%
  map_dfr(~ .x %>%
            select(-c(id, x, y)) %>%
            freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.)))))

persistent <- sm.freq %>% 
  colnames() %>% str_remove_all("\\.")


all.sm <- sm.repr %>% reduce(~rbind(.x, .y), .init = tibble())
clusters <- all.sm %>% select(-c(id,x,y)) %>% apply(1,which)

all.sm.clusters <- all.sm %>% select(c(id,x,y)) %>% cbind(cluster = clusters)

wc.repr <- bin_count_cluster(all.cells.dcis, all.positions.dcis, 200, 50, 33, freq=FALSE)

sm.wc <- all.sm.clusters %>% left_join(wc.repr, by = join_by(id,x,y)) %>% mutate(cluster.x = ifelse(cluster.x %in% persistent, cluster.x, 0))
sm.wc.p <- all.sm.clusters %>% left_join(wc.repr, by = join_by(id,x,y)) %>% filter(cluster.x %in% persistent)


aris <- sm.wc %>% select(-c(x,y)) %>% group_by(id) %>% 
  summarise(ARI = genieclust::adjusted_rand_score(cluster.x, cluster.y, clipped = TRUE))
aris.p <- sm.wc.p %>% select(-c(x,y)) %>% group_by(id) %>% 
  summarise(ARI.p = genieclust::adjusted_rand_score(cluster.x, cluster.y, clipped = TRUE))

write_rds(left_join(aris,aris.p, by = "id"), "ari/dcis.ct.rds")


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
    as_tibble(.name_repair = "unique")
})

names(all.cells.lymph) <- spots

all.positions.lymph <- spots %>% map(\(id){
  lymph %>%
    filter(Spots == id) %>%
    select(X, Y) %>%
    `colnames<-`(c("x", "y"))
})


## SM CTCL ----

repr.ids <- roc.sm <- NULL
max.auc <- 0

rocs <- c(100, 200, 300, 400, 500) %>% map_dbl(\(ws){
  # 10 neighbors as in publication, 200px = 75um window
  plan(multisession, workers = 8)
  misty.results <- sm_train(all.cells.lymph, all.positions.lymph, 10, ws, 20, paste0("CTCLct", ws, ".sqm"))

  plan(multisession, workers = 8)
  param.opt <- optimal_smclust(misty.results, outcome %>% select(-Patients) %>%
    rename(id = Spots, target = Groups) %>%
    mutate(id = as.character(id), target = as.factor(make.names(target))))

  sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])

  repr.ids <- sm.repr %>% pull(id)

  freq.sm <- sm.repr %>%
    left_join(outcome %>%
      mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
    rename(target = Groups) %>%
    mutate(target = as.factor(make.names(target))) %>%
    select(-id, -Patients)

  roc.sm <- classify(freq.sm)

  if (roc.sm$auc > max.auc) {
    roc.sm <<- roc.sm
    repr.ids <<- repr.ids
    max.auc <<- roc.sm$auc
  }

  roc.sm$auc
})

write_rds(rocs, "rocs/wrocs.ctcl.rds")


## WS CTCL ----

misty.results <- misty_train(all.cells.lymph, all.positions.lymph, 10, "CTCLct_ws.sqm")

## Alternatives CTCL ----

cn.results <- cn_train(all.cells.lymph, all.positions.lymph, 10)
cn.repr <- cn_labels(cn.results, k = 10)

freq.cn <- cn.repr %>%
  filter(id %in% repr.ids) %>%
  left_join(outcome %>%
    mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)

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


roc.fr <- classify(freq.lymph)

pa.lymph <- all.cells.lymph %>%
  map2(all.positions.lymph, pa_repr) %>%
  list_transpose() %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(Spots = spots) %>%
  left_join(outcome, by = "Spots") %>%
  filter(Spots %in% repr.ids) %>%
  select(-Spots) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target)))

roc.pa <- classify(pa.lymph)

csea.lymph <- all.cells.lymph %>%
  map2(all.positions.lymph, csea_repr) %>%
  list_transpose() %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(Spots = spots) %>%
  left_join(outcome, by = "Spots") %>%
  filter(Spots %in% repr.ids) %>%
  select(-Spots) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target)))

roc.csea <- classify(csea.lymph)


wc.repr <- bin_count_cluster(all.cells.lymph, all.positions.lymph, 400, 50, ncol(sm.repr) - 1)
freq.wc <- wc.repr %>%
  left_join(outcome %>%
    mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  filter(id %in% repr.ids) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)

roc.wc <- classify(freq.wc)

write_rds(list(sm = roc.sm, cn = roc.cn, fr = roc.fr, pa = roc.pa, csea = roc.csea, wc = roc.wc), "rocs/ctcl.ct.rds")

## ARI ----

misty.results <- read_rds("CTCLct400.rds")
sm.repr <- sm_repr(misty.results, cuts = 0.4, res = 0.9)
repr.ids <- sm.repr %>% map_chr(~ .x$id[1])

sm.freq <- sm.repr %>%
  map_dfr(~ .x %>%
            select(-c(id, x, y)) %>%
            freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.)))))

persistent <- sm.freq %>% 
  colnames() %>% str_remove_all("\\.")


all.sm <- sm.repr %>% reduce(~rbind(.x, .y), .init = tibble())
clusters <- all.sm %>% select(-c(id,x,y)) %>% apply(1,which)

all.sm.clusters <- all.sm %>% select(c(id,x,y)) %>% cbind(cluster = clusters)

wc.repr <- bin_count_cluster(all.cells.lymph, all.positions.lymph, 400, 50, 16, freq=FALSE)

sm.wc <- all.sm.clusters %>% left_join(wc.repr, by = join_by(id,x,y)) %>% mutate(cluster.x = ifelse(cluster.x %in% persistent, cluster.x, 0))
sm.wc.p <- all.sm.clusters %>% left_join(wc.repr, by = join_by(id,x,y)) %>% filter(cluster.x %in% persistent)


aris <- sm.wc %>% select(-c(x,y)) %>% group_by(id) %>% 
  summarise(ARI = genieclust::adjusted_rand_score(cluster.x, cluster.y, clipped = TRUE))
aris.p <- sm.wc.p %>% select(-c(x,y)) %>% group_by(id) %>% 
  summarise(ARI.p = genieclust::adjusted_rand_score(cluster.x, cluster.y, clipped = TRUE))

write_rds(left_join(aris,aris.p, by = "id"), "ari/ctcl.ct.rds")

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
    pull(core)
)
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

resp <- bmeta %>%
  filter(core %in% cores) %>%
  select(core, response)

repr.ids <- roc.sm <- NULL
max.auc <- 0

rocs <- c(100, 200, 300, 400, 500) %>% map_dbl(\(ws){
  plan(multisession, workers = 8)
  misty.results <- sm_train(all.cells.bc, all.positions.bc, 10, ws, 20, paste0("BCct", ws, ".sqm"))


  plan(multisession, workers = 8)
  param.opt <- optimal_smclust(misty.results, resp %>%
    rename(id = core, target = response) %>%
    mutate(id = as.character(id), target = as.factor(make.names(target))))

  sm.repr <- sm_labels(misty.results, cuts = param.opt["cut"], res = param.opt["res"])

  repr.ids <- sm.repr %>% pull(id)

  freq.sm <- sm.repr %>%
    left_join(resp, by = c("id" = "core")) %>%
    rename(target = response) %>%
    mutate(target = as.factor(make.names(target))) %>%
    select(-id)

  roc.sm <- classify(freq.sm)

  if (roc.sm$auc > max.auc) {
    roc.sm <<- roc.sm
    repr.ids <<- repr.ids
    max.auc <<- roc.sm$auc
  }

  roc.sm$auc
})

write_rds(rocs, "rocs/wrocs.bc.rds")

## WS BC ----

misty.results <- misty_train(all.cells.bc, all.positions.bc, 10, "BCct_ws.sqm")


## Alternatives BC ----

cn.results <- cn_train(all.cells.bc, all.positions.bc, 10)
cn.repr <- cn_labels(cn.results, k = 6)

freq.cn <- cn.repr %>%
  filter(id %in% repr.ids) %>%
  left_join(resp, by = c("id" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id)

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

roc.fr <- classify(freq.bc)

pa.bc <- all.cells.bc %>%
  map2(all.positions.bc, pa_repr) %>%
  list_transpose() %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(core = cores) %>%
  filter(core %in% repr.ids) %>%
  left_join(resp, by = "core") %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-core)

roc.pa <- classify(pa.bc)

csea.bc <- all.cells.bc %>%
  map2(all.positions.bc, csea_repr) %>%
  list_transpose() %>%
  as_tibble(.name_repair = "unique") %>%
  mutate(core = cores) %>%
  filter(core %in% repr.ids) %>%
  left_join(resp, by = "core") %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-core)

roc.csea <- classify(csea.bc)

wc.repr <- bin_count_cluster(all.cells.bc, all.positions.bc, 200, 50, 9)

freq.wc <- wc.repr %>%
  filter(id %in% repr.ids) %>%
  left_join(resp, by = c("id" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id)

roc.wc <- classify(freq.wc)

write_rds(list(sm = roc.sm, cn = roc.cn, fr = roc.fr, pa = roc.pa, csea = roc.csea, wc = roc.wc), "rocs/bc.ct.rds")

## ARI ----

misty.results <- read_rds("BCct200.rds")
sm.repr <- sm_repr(misty.results, cuts = 0.3, res = 0.8)
repr.ids <- sm.repr %>% map_chr(~ .x$id[1])

sm.freq <- sm.repr %>%
  map_dfr(~ .x %>%
            select(-c(id, x, y)) %>%
            freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.)))))

persistent <- sm.freq %>% 
  colnames() %>% str_remove_all("\\.")


all.sm <- sm.repr %>% reduce(~rbind(.x, .y), .init = tibble())
clusters <- all.sm %>% select(-c(id,x,y)) %>% apply(1,which)

all.sm.clusters <- all.sm %>% select(c(id,x,y)) %>% cbind(cluster = clusters)

wc.repr <- bin_count_cluster(all.cells.bc, all.positions.bc, 200, 50, 9, freq=FALSE)

sm.wc <- all.sm.clusters %>% left_join(wc.repr, by = join_by(id,x,y)) %>% mutate(cluster.x = ifelse(cluster.x %in% persistent, cluster.x, 0))
sm.wc.p <- all.sm.clusters %>% left_join(wc.repr, by = join_by(id,x,y)) %>% filter(cluster.x %in% persistent)


aris <- sm.wc %>% select(-c(x,y)) %>% group_by(id) %>% 
  summarise(ARI = genieclust::adjusted_rand_score(cluster.x, cluster.y, clipped = TRUE))
aris.p <- sm.wc.p %>% select(-c(x,y)) %>% group_by(id) %>% 
  summarise(ARI.p = genieclust::adjusted_rand_score(cluster.x, cluster.y, clipped = TRUE))

write_rds(left_join(aris,aris.p, by = "id"), "ari/bc.ct.rds")
