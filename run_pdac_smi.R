library(anndata)
library(uwot)
source("utils.R")

plan(multisession, workers = 8)
options(future.globals.maxSize = 2024^3)

# Read data ----
all <- read_h5ad("data/PDACSMI/raw_meta_data_final.h5ad")

all.feat <- all$obs %>%
  filter(treatment_status %in% c("CRT", "Untreated")) %>%
  dplyr::select(
    sample_id, fov, donor, treatment_status,
    CenterX_local_px, CenterY_local_px, annotation_majortypes, annotation_subtypes
  )


# Collect condition information by sample ----

ids <- all.feat %>%
  group_by(sample_id, fov) %>%
  summarise(id = paste0(sample_id[1], "_", fov[1]), .groups = "drop") %>%
  pull(id)

conditions <- all.feat %>%
  group_by(sample_id, fov) %>%
  summarise(
    id = paste0(sample_id[1], "_", fov[1]),
    treatment_status = as.character(treatment_status[1]), .groups = "drop"
  )


# Collect position and cell type information ----

# fov size 0.9 x 0.7 mm^2
# SMI CosMx
# 180nm per px
all.positions <- all.feat %>%
  group_by(sample_id, fov) %>%
  dplyr::select(CenterX_local_px, CenterY_local_px) %>%
  group_split(.keep = FALSE)

names(all.positions) <- ids

# Cell type information
cts <- levels(all.feat$annotation_subtypes)

all.cells <- all.feat %>%
  group_by(sample_id, fov) %>%
  dplyr::select(annotation_subtypes) %>%
  group_split(.keep = FALSE) %>%
  map(\(samp) samp %>%
    pull(annotation_subtypes) %>%
    map(~ .x == cts) %>%
    rlist::list.rbind() %>%
    `colnames<-`(make.names(cts)) %>%
    as_tibble())

names(all.cells) <- ids

# Run Kasumi from benchmark implementation utils ----

paraview.neighbors <- 10
window.size <- 1200
minimum.cells <- 30

kasumi.results <- sm_train(
  all.cells, all.positions, paraview.neighbors,
  window.size, minimum.cells, "PDACct1200.sqm"
)


# Find an optimal parameter combination given a downstream task
# Takes a long time to compute and there is no need since most of the parameter
# combinations will lead to good separation of samples. Skip to precalculated values.

# plan(multisession, workers = 8)
# param.opt <- optimal_smclust(kasumi.results, conditions %>%
#   select(id, treatment_status) %>%
#   mutate(
#     target = treatment_status
#   ) %>%
#   mutate(
#     target = as.factor(make.names(target))
#   ))


# Generate representation ----
sm.repr <- sm_labels(kasumi.results, cuts = 0.6, res = 0.6)


# Classifiction task ----
freq.sm <- sm.repr %>%
  left_join(conditions %>% select(id, treatment_status), by = "id") %>%
  rename(target = treatment_status) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id)

roc.sm <- classify(freq.sm)

roc.sm$auc

# Model reliance ----

mr <- model_reliance(freq.sm)
ggsave("PDACreliance.pdf", height = 6, width = 5)

# Kasumi representation UMAP ----

kasumi.umap <- freq.sm %>%
  select(-target) %>%
  umap(n_neighbors = 10, seed = 1)

data.frame(kasumi.umap[, 1:2]) %>%
  mutate(Treatment = freq.sm$target) %>%
  rename(UMAP1 = X1, UMAP2 = X2) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Treatment)) +
  geom_point() +
  theme_classic()

ggsave("PDACumap.pdf", height = 5, width = 6)

# Heatmaps ----

sm.repr.ext <- sm_labels(kasumi.results, 0.6, 0.6, freq = FALSE)

kasumi.cluster.13 <- describe_cluster(sm.repr.ext, 13, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.13, trim = 1)
ggsave("clusterfigs/pdac13.pdf", width = 4, height = 4)
plot_interaction_heatmap(kasumi.cluster.13, "paraview.10", trim = 1, clean = TRUE)
ggsave("clusterfigs/pdac13h.pdf", width = 5, height = 4)
plot_interaction_heatmap(kasumi.cluster.13, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/pdac13c.pdf", width = 5, height = 4)

kasumi.cluster.1 <- describe_cluster(sm.repr.ext, 1, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.1, trim = 1)
ggsave("clusterfigs/pdac1.pdf", width = 4, height = 4)
plot_interaction_heatmap(kasumi.cluster.1, "paraview.10", trim = 1, clean = TRUE)
ggsave("clusterfigs/pdac1h.pdf", width = 5, height = 3)
plot_interaction_heatmap(kasumi.cluster.1, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/pdac1c.pdf", width = 5, height = 3)

kasumi.cluster.3 <- describe_cluster(sm.repr.ext, 3, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.3, trim = 1)
ggsave("clusterfigs/pdac3.pdf", width = 4, height = 4)
plot_interaction_heatmap(kasumi.cluster.3, "paraview.10", trim = 1, clean = TRUE)
ggsave("clusterfigs/pdac3h.pdf", width = 5, height = 3)
plot_interaction_heatmap(kasumi.cluster.3, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/pdac3c.pdf", width = 5, height = 3)



kasumi.cluster.30 <- describe_cluster(sm.repr.ext, 30, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.30, trim = 1)
ggsave("clusterfigs/pdac30.pdf", width = 4, height = 4)
plot_interaction_heatmap(kasumi.cluster.30, "paraview.10", trim = 1, clean = TRUE)
ggsave("clusterfigs/pdac30h.pdf", width = 5, height = 3)
plot_interaction_heatmap(kasumi.cluster.30, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/pdac30c.pdf", width = 5, height = 3)

kasumi.cluster.116 <- describe_cluster(sm.repr.ext, 116, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.116, trim = 1)
ggsave("clusterfigs/pdac116.pdf", width = 4, height = 4)
plot_interaction_heatmap(kasumi.cluster.116, "paraview.10", trim = 1, clean = TRUE)
ggsave("clusterfigs/pdac116h.pdf", width = 5, height = 4)
plot_interaction_heatmap(kasumi.cluster.116, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/pdac116c.pdf", width = 5, height = 4)

kasumi.cluster.18 <- describe_cluster(sm.repr.ext, 18, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.18, trim = 1)
ggsave("clusterfigs/pdac18.pdf", width = 4, height = 4)
plot_interaction_heatmap(kasumi.cluster.18, "paraview.10", trim = 1, clean = TRUE)
ggsave("clusterfigs/pdac18h.pdf", width = 5, height = 4)
plot_interaction_heatmap(kasumi.cluster.18, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/pdac18c.pdf", width = 5, height = 4)



# Baseline

baseline <- all.cells[sm.repr$id] %>%
  map_dfr(~ colSums(.x) / sum(.x)) %>%
  mutate(id = sm.repr$id, .before = 1)

baseline.umap <- umap(baseline %>% select(-id),n_neighbors = 10, seed = 1)
toplot <- cbind(data.frame(baseline.umap), id = baseline$id) %>%
  mutate(id = str_extract(id, "(T|U)")) %>%
  rename(UMAP1 = X1, UMAP2 = X2, Treatment = id) %>%
  mutate(Treatment = ifelse(Treatment == "T", "CRT", "Untreated"))

ggplot(toplot, aes(x = UMAP1, y = UMAP2, color = Treatment)) +
  geom_point() +
  theme_classic()

ggsave("PDACbaseumap.pdf", height = 5, width = 6)

# 0.9637
baseline.roc <- classify(baseline %>% mutate(target = str_extract(baseline$id, "(T|U)")) %>% select(-id))
baseline.roc$auc

baseline.mr <- model_reliance(baseline %>% mutate(target = str_extract(baseline$id, "(T|U)")) %>%
  mutate(target = as.factor(ifelse(target == "T", "CRT", "Untreated"))) %>%
  select(-id))

ggsave("PDACbasereliance.pdf", height = 6, width = 5)
 