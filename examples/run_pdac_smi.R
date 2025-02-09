library(anndata)
library(uwot)
source("utils.R")

plan(multisession, workers = 8)
options(future.globals.maxSize = 2024^3)

# Download data ----
if (!file.exists("raw_meta_data_final.h5ad")) {
  download.file("https://data.mendeley.com/public-files/datasets/kx6b69n3cb/files/26999a73-4b05-4f47-9cab-023cd0de80a9/file_downloaded", "raw_meta_data_final.h5ad", method = "libcurl")
}
all <- read_h5ad("raw_meta_data_final.h5ad")

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

plan(multisession, workers = 8)
param.opt <- optimal_smclust(kasumi.results, conditions %>%
   select(id, treatment_status) %>%
   mutate(
     target = treatment_status) %>%
    mutate(
     target = as.factor(make.names(target))
   ))


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


# Kasumi representation UMAP ----

kasumi.umap <- freq.sm %>%
  select(-target) %>%
  umap(seed = 1)

data.frame(kasumi.umap[, 1:2]) %>%
  mutate(Treatment = freq.sm$target) %>%
  ggplot(aes(x = X1, y = X2, color = Treatment)) +
  geom_point() +
  theme_classic()


# Plots ----

sm.repr.ext <- sm_labels(kasumi.results, 0.6, 0.6, freq = FALSE)

kasumi.cluster.13 <- describe_cluster(sm.repr.ext, 13, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.13, trim = 1)
plot_interaction_heatmap(kasumi.cluster.13, "paraview.10", trim = 1, clean = TRUE)
plot_interaction_heatmap(kasumi.cluster.13, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)

kasumi.cluster.1 <- describe_cluster(sm.repr.ext, 1, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.1, trim = 1)
plot_interaction_heatmap(kasumi.cluster.1, "paraview.10", trim = 1, clean = TRUE)
plot_interaction_heatmap(kasumi.cluster.1, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)

kasumi.cluster.30 <- describe_cluster(sm.repr.ext, 30, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.30, trim = 1)
plot_interaction_heatmap(kasumi.cluster.30, "paraview.10", trim = 1, clean = TRUE)
plot_interaction_heatmap(kasumi.cluster.30, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)

kasumi.cluster.116 <- describe_cluster(sm.repr.ext, 116, "PDACct1200.sqm")

plot_improvement_stats(kasumi.cluster.116, trim = 1)
plot_interaction_heatmap(kasumi.cluster.116, "paraview.10", trim = 1, clean = TRUE)
plot_interaction_heatmap(kasumi.cluster.116, "paraview.10", trim = 1, clean = TRUE, correlation = TRUE)


# Baseline

baseline <- all.cells[sm.repr$id] %>%
  map_dfr(~ colSums(.x)/sum(.x)) %>%
  mutate(id = sm.repr$id, .before = 1)

baseline.umap <- umap(baseline %>% select(-id), seed = 1)
toplot <- cbind(data.frame(baseline.umap), id = baseline$id) %>% mutate(id = str_extract(id, "(T|U)"))
ggplot(toplot, aes(x = X1, y = X2, color = id)) +
  geom_point() +
  theme_classic()

baseline.roc <- classify(baseline %>% mutate(target = str_extract(baseline$id, "(T|U)")) %>% select(-id))
baseline.roc$auc

baseline.mr <- model_reliance(baseline %>% mutate(target = as.factor(str_extract(baseline$id, "(T|U)"))) %>% select(-id))
