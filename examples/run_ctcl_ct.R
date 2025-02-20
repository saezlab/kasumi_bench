source("../utils.R")
library(readxl)

# Download data ----
if(!file.exists("CTCL.xlsx"))
  download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-26974-6/MediaObjects/41467_2021_26974_MOESM3_ESM.xlsx", "CTCL.xlsx")
lymph <- read_xlsx("CTCL.xlsx", skip=2)


# Collect condition information by sample ----
spots <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  pull(Spots) %>%
  unique()
outcome <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  group_by(Spots, Patients) %>%
  summarize(Groups = Groups[1], .groups = "drop")


# Collect position and cell type information ----

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


# Run Kasumi from benchmark implementation utils ----

plan(multisession, workers = 8)
kasumi.results <- sm_train(all.cells.lymph, all.positions.lymph, 10, 400, 20, "CTCLct400.sqm")


# Find an optimal parameter combination given a downstream task ----

# You can skip this step and use precalculated values for the next step

# sm.repr <- sm_labels(kasumi.results, cuts = 0.4, res = 0.9)

plan(multisession, workers = 8)
param.opt <- optimal_smclust(kasumi.results, outcome %>% select(-Patients) %>%
                                 rename(id = Spots, target = Groups) %>%
                                 mutate(id = as.character(id), target = as.factor(make.names(target))))


# Generate representation ----
sm.repr <- sm_labels(kasumi.results, cuts = param.opt["cut"], res = param.opt["res"])


# Classifiction task ----

freq.sm <- sm.repr %>%
  left_join(outcome %>%
              mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)

roc.sm <- classify(freq.sm)

roc.sm$auc


# Model reliance ----

mr <- model_reliance(freq.sm)


# Plotting ----
sm.repr.ext <- sm_labels(kasumi.results, 0.4, 0.9, freq = FALSE)

kasumi.cluster.23 <- describe_cluster(sm.repr.ext, 23, "CTCLct400.sqm")

plot_improvement_stats(kasumi.cluster.23, trim = 1)
plot_interaction_heatmap(kasumi.cluster.23, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE)
plot_interaction_heatmap(kasumi.cluster.23, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)


misty.cluster.16 <- describe_cluster(sm.repr.ext, 16, "CTCLct400.sqm")

plot_improvement_stats(kasumi.cluster.16, trim = 1)
plot_interaction_heatmap(kasumi.cluster.16, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE)
plot_interaction_heatmap(kasumi.cluster.16, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)


misty.cluster.4 <- describe_cluster(sm.repr.ext, 4, "CTCLct400.sqm")

plot_improvement_stats(kasumi.cluster.4, trim = 1)
plot_interaction_heatmap(kasumi.cluster.4, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE)
plot_interaction_heatmap(kasumi.cluster.4, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)


kasumi.cluster.13 <- describe_cluster(sm.repr.ext, 13, "CTCLct400.sqm")

plot_improvement_stats(kasumi.cluster.13, trim = 1)
plot_interaction_heatmap(kasumi.cluster.13, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE)
plot_interaction_heatmap(kasumi.cluster.13, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)

