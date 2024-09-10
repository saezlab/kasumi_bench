source("utils.R")

# MIBI  ----
ct <- read_rds("rocs/dcis.ct.rds")
expr <- read_rds("rocs/dcis.expr.rds")
ws <- read_rds("rocs/dcis.ws.rds")
expr.base <- read_rds("rocs/dcis.expr.base.rds")
expr.alts <- read_rds("rocs/dcis.expr.alts.rds")

ggroc(list.append(list.remove(ct, "pa"), ws.ct = ws$ct), legacy.axes = TRUE, linewidth = pi / 10) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2 / 3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("MIBI DCIS") + theme_light()

ggsave("roc.mibi.pdf", width = 4, height = 3)

ggroc(
  list(
    sm.ct = ct$sm, banksy = expr.alts$banksy.ct, base = expr.base, cc = expr.alts$cc.l1, wc = expr.alts$wc, ws.expr = ws$expr, sm.expr = expr, utag = expr.alts$utag
  ),
  legacy.axes = TRUE, linewidth = pi / 10
) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2 / 3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("MIBI DCIS") + theme_light()

ggsave("roc.misty.mibi.pdf", width = 4, height = 3)


# CODEX ----

ct <- read_rds("rocs/ctcl.ct.rds")
expr <- read_rds("rocs/ctcl.expr.rds")
ws <- read_rds("rocs/ctcl.ws.rds")
expr.base <- read_rds("rocs/ctcl.expr.base.rds")
expr.alts <- read_rds("rocs/ctcl.expr.alts.rds")

lymph <- read_csv("data/LymphomaCODEX/single_cells.csv")

spots <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  pull(Spots) %>%
  unique()

outcome <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  group_by(Spots, Patients) %>%
  summarize(Groups = Groups[1], .groups = "drop")

ggroc(list.append(list.remove(ct, "pa"), ws.ct = ws$ct), legacy.axes = TRUE, linewidth = pi / 10) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2 / 3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("CODEX CTCL") + theme_light()

ggsave("roc.codex.pdf", width = 4, height = 3)

ggroc(
  list(
    sm.ct = ct$sm, banksy = expr.alts$banksy.ct, base = expr.base, cc = expr.alts$cc.l1, wc = expr.alts$wc, ws.expr = ws$expr, sm.expr = expr, utag = expr.alts$utag
  ),
  legacy.axes = TRUE, linewidth = pi / 10
) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2 / 3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("CODEX CTCL") + theme_light()

ggsave("roc.misty.codex.pdf", width = 4, height = 3)


misty.results <- read_rds("CTCLct400.rds")
sm.repr <- sm_labels(misty.results, 0.4, 0.9)

# X1 responders, X2 non-responders
freq.sm <- sm.repr %>%
  left_join(outcome %>%
    mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)

mr <- model_reliance(freq.sm)

ggsave("CTCLreliance.pdf", height = 6, width = 4)

sm.repr.ext <- sm_labels(misty.results, 0.4, 0.9, freq = FALSE)

misty.cluster.23 <- describe_cluster(sm.repr.ext, 23, "CTCLct400.sqm")

plot_improvement_stats(misty.cluster.23, trim = 1)
ggsave("clusterfigs/ctcl23.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.23, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE)
ggsave("clusterfigs/ctcl23h.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.23, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/ctcl23c.pdf", width = 4, height = 4)

misty.cluster.23$importances.aggregated <- misty.cluster.23$importances.aggregated %>%
  mutate(Predictor = str_remove(Predictor, "^p\\."))

plot_interaction_communities(misty.cluster.23,
  view = "paraview.10",
  trim = 1, cutoff = 0.5, path = "clusterfigs/ctcl23.graphml"
)

misty.cluster.16 <- describe_cluster(sm.repr.ext, 16, "CTCLct400.sqm")

plot_improvement_stats(misty.cluster.16, trim = 1)
ggsave("clusterfigs/ctcl16.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.16, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE)
ggsave("clusterfigs/ctcl16h.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.16, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/ctcl16c.pdf", width = 4, height = 4)

misty.cluster.16$importances.aggregated <- misty.cluster.16$importances.aggregated %>%
  mutate(Predictor = str_remove(Predictor, "^p\\."))

plot_interaction_communities(misty.cluster.16,
  view = "paraview.10",
  trim = 1, cutoff = 0.5, path = "clusterfigs/ctcl16.graphml"
)


misty.cluster.4 <- describe_cluster(sm.repr.ext, 4, "CTCLct400.sqm")

plot_improvement_stats(misty.cluster.4, trim = 1)
ggsave("clusterfigs/ctcl4.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.4, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE)
ggsave("clusterfigs/ctcl4h.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.4, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/ctcl4c.pdf", width = 4, height = 4)

misty.cluster.4$importances.aggregated <- misty.cluster.4$importances.aggregated %>%
  mutate(Predictor = str_remove(Predictor, "^p\\."))

plot_interaction_communities(misty.cluster.4,
  view = "paraview.10",
  trim = 1, cutoff = 0.5, path = "clusterfigs/ctcl4.graphml"
)


misty.cluster.13 <- describe_cluster(sm.repr.ext, 13, "CTCLct400.sqm")

plot_improvement_stats(misty.cluster.13, trim = 1)
ggsave("clusterfigs/ctcl13.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.13, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE)
ggsave("clusterfigs/ctcl13h.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.13, "paraview.10", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/ctcl13c.pdf", width = 4, height = 4)

misty.cluster.13$importances.aggregated <- misty.cluster.13$importances.aggregated %>%
  mutate(Predictor = str_remove(Predictor, "^p\\."))

plot_interaction_communities(misty.cluster.13,
  view = "paraview.10",
  trim = 1, cutoff = 0.5, path = "clusterfigs/ctcl13.graphml"
)


# IMC ----

ct <- read_rds("rocs/bc.ct.rds")
expr <- read_rds("rocs/bc.expr.rds")
ws <- read_rds("rocs/bc.ws.rds")
expr.base <- read_rds("rocs/bc.expr.base.rds")
expr.alts <- read_rds("rocs/bc.expr.alts.rds")

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
resp <- bmeta %>%
  filter(core %in% cores) %>%
  select(core, response)

ggroc(list.append(list.remove(ct, "pa"), ws.ct = ws$ct), legacy.axes = TRUE, linewidth = pi / 10) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2 / 3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("IMC BC") + theme_light()

ggsave("roc.imc.pdf", width = 4, height = 3)

ggroc(
  list(
    sm.ct = ct$sm, banksy = expr.alts$banksy.ct, base = expr.base, cc = expr.alts$cc.l1, wc = expr.alts$wc, ws.expr = ws$expr, sm.expr = expr, utag = expr.alts$utag
  ),
  legacy.axes = TRUE, linewidth = pi / 10
) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2 / 3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("IMC BC") + theme_light()

ggsave("roc.misty.imc.pdf", width = 4, height = 3)


misty.results <- read_rds("BCexpr200.rds")

sm.repr <- sm_labels(misty.results, 0.3, 0.8)

freq.sm <- sm.repr %>%
  left_join(resp, by = c("id" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id)

mr <- model_reliance(freq.sm)
ggsave("BCreliance_expr.pdf", height = 5, width = 4)

sm.repr.ext <- sm_labels(misty.results, 0.3, 0.8, freq = FALSE)

misty.cluster.37 <- describe_cluster(sm.repr.ext, 37, "BCexpr200.sqm")

plot_improvement_stats(misty.cluster.37, trim = 1)
ggsave("clusterfigs/bc37.pdf", width = 5, height = 4)
plot_interaction_heatmap(misty.cluster.37, "intraview", trim = 1, cutoff = 1, clean = TRUE)
ggsave("clusterfigs/bc37h_intra.pdf", width = 5, height = 5)
plot_interaction_heatmap(misty.cluster.37, "intraview", trim = 1, cutoff = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/bc37c_intra.pdf", width = 5, height = 5)
plot_interaction_heatmap(misty.cluster.37, "paraview.50", trim = 1, cutoff = 0.5, clean = TRUE)
ggsave("clusterfigs/bc37h.pdf", width = 5, height = 4)
plot_interaction_heatmap(misty.cluster.37, "paraview.50", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/bc37c.pdf", width = 5, height = 4)


misty.cluster.5 <- describe_cluster(sm.repr.ext, 5, "BCexpr200.sqm")

plot_improvement_stats(misty.cluster.5, trim = 1)
ggsave("clusterfigs/bc5.pdf", width = 5, height = 4)
plot_interaction_heatmap(misty.cluster.5, "intraview", trim = 1, cutoff = 1, clean = TRUE)
ggsave("clusterfigs/bc5h_intra.pdf", width = 5, height = 5)
plot_interaction_heatmap(misty.cluster.5, "intraview", trim = 1, cutoff = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/bc5c_intra.pdf", width = 5, height = 5)
plot_interaction_heatmap(misty.cluster.5, "paraview.50", trim = 1, cutoff = 0.5, clean = TRUE)
ggsave("clusterfigs/bc5h.pdf", width = 5, height = 4)
plot_interaction_heatmap(misty.cluster.5, "paraview.50", trim = 1, cutoff = 0.5, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/bc5c.pdf", width = 5, height = 4)



# Sensitivity ----

ws.dcis <- read_rds("rocs/wrocs.dcis.rds")
ws.ctcl <- read_rds("rocs/wrocs.ctcl.rds")
ws.bc <- read_rds("rocs/wrocs.bc.rds")

ws <- tibble(Window = seq(100, 500, by = 100), DCIS = ws.dcis, CTCL = ws.ctcl, BC = ws.bc) %>%
  pivot_longer(-Window, names_to = "Data", values_to = "AUROC")

ggplot(ws, aes(x = Window, y = AUROC, color = Data)) +
  geom_point() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  ylim(c(0.5, 1)) +
  theme_light()

ggsave("ws.pdf", width = 4, height = 3)


sens <- read_csv("sensitivity.csv")

sens.long <- sens %>%
  pivot_longer(-c(cut, res), names_to = "Data", values_to = "AUC") %>%
  mutate(across(c(cut, res), as.factor)) %>%
  rename(Threshold = cut, Resolution = res)

ggplot(sens.long, aes(y = Resolution, x = Threshold, fill = AUC)) +
  geom_tile() +
  scale_fill_gradient(low = "gray50", high = "tomato3", limits = c(0.5, 1)) +
  facet_wrap("Data", ncol = 2) +
  coord_equal() +
  theme_classic()

ggsave("sensitivity.pdf")




# Kasumi sees CTCL ct ----

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


test <- lymph %>%
  filter(Spots %in% c(2, 39, 57, 22)) %>%
  select(Spots, ClusterName, X, Y) %>%
  mutate(Spots = factor(Spots, levels = c(2, 39, 57, 22)))
ggplot(test, aes(x = X, y = Y, color = ClusterName)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(vars(Spots)) +
  scale_color_manual(values = as.vector(pals::cols25(n = 20))) +
  coord_equal() +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("ct_distro.pdf", width = 8, height = 6)


all.lymph.g1 <- lymph %>%
  filter(Spots %in% (outcome %>% filter(Groups == 1) %>% pull(Spots))) %>%
  select(Spots, ClusterName, X, Y)
ggplot(all.lymph.g1, aes(x = X, y = Y, color = ClusterName)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(vars(Spots), ncol = 6) +
  scale_color_manual(values = as.vector(pals::cols25(n = 20))) +
  coord_equal() +
  labs(color = "Cell type") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("ct_distro_g1_lymph.pdf", width = 16, height = 8)

all.lymph.g2 <- lymph %>%
  filter(Spots %in% (outcome %>% filter(Groups == 2) %>% pull(Spots))) %>%
  select(Spots, ClusterName, X, Y)
ggplot(all.lymph.g2, aes(x = X, y = Y, color = ClusterName)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(vars(Spots), ncol = 6) +
  scale_color_manual(values = append(pals::cols25(n = 20), values = "#000000", after = 6)) +
  coord_equal() +
  labs(color = "Cell type") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("ct_distro_g2_lymph.pdf", width = 16, height = 8)



misty.results <- read_rds("CTCLct400.rds")
sm.repr <- sm_repr(misty.results, cuts = 0.4, res = 0.9)
repr.ids <- sm.repr %>% map_chr(~ .x$id[1])

sm.freq <- sm.repr %>%
  map_dfr(~ .x %>%
    select(-c(id, x, y)) %>%
    freq_repr()) %>%
  select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.)))))

persistent <- sm.freq %>%
  colnames() %>%
  str_remove_all("\\.")


test <- rbind(sm.repr[[7]], sm.repr[[17]], sm.repr[[26]], sm.repr[[9]])
clusters <- test %>%
  select(-c(id, x, y)) %>%
  apply(1, which)

kasumirep <- test %>%
  select(c(id, x, y)) %>%
  cbind(cluster = clusters) %>%
  filter(cluster %in% persistent) %>%
  mutate(cluster = as.factor(cluster), id = factor(id, levels = c(2, 39, 57, 22)))

ggplot(kasumirep, aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 200, width = 200, alpha = 0.7) +
  facet_wrap(id ~ .) +
  scale_fill_manual(values = as.vector(pals::cols25(n = 15))) +
  coord_equal() +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("kasumi_sees_tiles.pdf", width = 8, height = 6)


ids <- outcome %>%
  filter(Groups == 1) %>%
  pull(Spots)
repr.map <- as.character(ids) %>% map_int(~ which(repr.ids == .x))

# color scheme inconsistency
orig <- pals::cols25(n = 13)
plus <- pals::cols25(n = 16)
ext <- append(append(orig, plus[c(14, 15)], after = 6), plus[16], after = 12)

all.lymph.g1 <- repr.map %>% reduce(~ rbind(.x, sm.repr[[.y]]), .init = tibble())
clusters <- all.lymph.g1 %>%
  select(-c(id, x, y)) %>%
  apply(1, which)

kasumirep <- all.lymph.g1 %>%
  select(c(id, x, y)) %>%
  cbind(cluster = clusters) %>%
  filter(cluster %in% persistent) %>%
  mutate(cluster = as.factor(cluster), id = factor(id, levels = ids))

ggplot(kasumirep, aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 200, width = 200, alpha = 0.7) +
  facet_wrap(id ~ ., ncol = 6) +
  scale_fill_manual(values = ext) +
  coord_equal() +
  labs(x = "X", y = "Y", fill = "Kasumi\ncluster") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("kasumi_sees_g1_lymph.pdf", width = 16, height = 8)

wc.repr <- bin_count_cluster(all.cells.lymph, all.positions.lymph, 400, 20, 50, 16, freq = FALSE)

ggplot(wc.repr %>% filter(id %in% ids) %>% mutate(cluster = as.factor(cluster), id = as.integer(id)), aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 200, width = 200, alpha = 0.7) +
  facet_wrap(id ~ ., ncol = 6) +
  scale_fill_manual(values = ext) +
  coord_equal() +
  labs(x = "X", y = "Y", fill = "Window\ncluster") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("windows_g1_lymph.pdf", width = 16, height = 8)


ids <- outcome %>%
  filter(Groups == 2) %>%
  pull(Spots)
repr.map <- as.character(ids) %>% map_int(~ which(repr.ids == .x))

all.lymph.g2 <- repr.map %>% reduce(~ rbind(.x, sm.repr[[.y]]), .init = tibble())
clusters <- all.lymph.g2 %>%
  select(-c(id, x, y)) %>%
  apply(1, which)

kasumirep <- all.lymph.g2 %>%
  select(c(id, x, y)) %>%
  cbind(cluster = clusters) %>%
  filter(cluster %in% persistent) %>%
  mutate(cluster = as.factor(cluster), id = factor(id, levels = ids))

ggplot(kasumirep, aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 200, width = 200, alpha = 0.7) +
  facet_wrap(id ~ ., ncol = 6) +
  scale_fill_manual(values = ext) +
  coord_equal() +
  labs(x = "X", y = "Y", fill = "Kasumi\ncluster") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("kasumi_sees_g2_lymph.pdf", width = 16, height = 8)


ggplot(wc.repr %>% filter(id %in% ids) %>% mutate(cluster = as.factor(cluster), id = as.integer(id)), aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 200, width = 200, alpha = 0.7) +
  facet_wrap(id ~ ., ncol = 6) +
  scale_fill_manual(values = ext) +
  coord_equal() +
  labs(x = "X", y = "Y", fill = "Window\ncluster") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("windows_g2_lymph.pdf", width = 16, height = 8)


# Kasumi sees BC expr ----

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

resp <- bmeta %>%
  filter(core %in% cores) %>%
  select(core, response)


misty.results <- read_rds("BCexpr200.rds")
sm.repr <- sm_repr(misty.results, 0.3, 0.8)
repr.ids <- sm.repr %>% map_chr(~ .x$id[1])
sm.freq <- sm_labels(misty.results, 0.3, 0.8)

persistent <- sm.freq %>%
  select(-id) %>%
  colnames() %>%
  str_remove_all("\\.")

test <- rbind(sm.repr[[8]], sm.repr[[29]], sm.repr[[7]], sm.repr[[19]])
clusters <- test %>%
  select(-c(id, x, y)) %>%
  apply(1, which)

kasumirep <- test %>%
  select(c(id, x, y)) %>%
  cbind(cluster = clusters) %>%
  filter(cluster %in% persistent) %>%
  mutate(cluster = as.factor(cluster), id = factor(id, levels = c(repr.ids[8], repr.ids[29], repr.ids[7], repr.ids[19])))

ggplot(kasumirep, aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 100, width = 100, alpha = 0.7) +
  facet_wrap(id ~ .) +
  scale_fill_manual(values = as.vector(pals::cols25(n = 15))) +
  coord_equal() +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("kasumi_sees_tiles_expr.pdf", width = 8, height = 6)

ids <- resp %>%
  filter(response == "Sensitive") %>%
  pull(core)
repr.map <- ids %>% map_int(~ which(repr.ids == .x))

ext <- pals::cols25(n = 15)


all.bc.sens <- repr.map %>% reduce(~ rbind(.x, sm.repr[[.y]]), .init = tibble())
clusters <- all.bc.sens %>%
  select(-c(id, x, y)) %>%
  apply(1, which)

kasumirep <- all.bc.sens %>%
  select(c(id, x, y)) %>%
  cbind(cluster = clusters) %>%
  filter(cluster %in% persistent) %>%
  mutate(cluster = as.factor(cluster), id = factor(id, levels = ids))

ggplot(kasumirep, aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 100, width = 100, alpha = 0.7) +
  facet_wrap(id ~ ., ncol = 6) +
  scale_fill_manual(values = ext) +
  coord_equal() +
  labs(x = "X", y = "Y", fill = "Kasumi\ncluster") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("kasumi_sees_sens_bc.pdf", width = 16, height = 8)

wc.repr <- bin_count_cluster(all.cells.bc, all.positions.bc, 200, 20, 50, 15, freq = FALSE)

ggplot(wc.repr %>% filter(id %in% ids) %>% mutate(cluster = as.factor(cluster), id = factor(id, levels = ids)), aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 100, width = 100, alpha = 0.7) +
  facet_wrap(id ~ ., ncol = 6) +
  scale_fill_manual(values = ext) +
  coord_equal() +
  labs(x = "X", y = "Y", fill = "Window\ncluster") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("windows_sens_bc.pdf", width = 16, height = 8)


banksy.single <- banksy_single(all.cells.bc, all.positions.bc, 10, 0.2, 0.5)

ggplot(banksy.single %>% filter(id %in% ids) %>% mutate(id = factor(id, levels = ids)), aes(x = x, y = y, color = cluster)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(vars(id), ncol = 6) +
  scale_color_manual(values = as.vector(pals::cols25(n = 10))) +
  coord_equal() +
  labs(x = "X", y = "Y", color = "Banksy\nclusters") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("banksy_sens_bc.pdf", width = 16, height = 8)

ids <- resp %>%
  filter(response == "Resistant") %>%
  pull(core)
repr.map <- ids[-9] %>% map_int(~ which(repr.ids == .x))

all.bc.res <- repr.map %>% purrr::reduce(~ rbind(.x, sm.repr[[.y]]), .init = tibble())
clusters <- all.bc.res %>%
  select(-c(id, x, y)) %>%
  apply(1, which)

kasumirep <- all.bc.res %>%
  select(c(id, x, y)) %>%
  cbind(cluster = clusters) %>%
  filter(cluster %in% persistent) %>%
  mutate(cluster = as.factor(cluster), id = factor(id, levels = ids))

ggplot(kasumirep, aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 100, width = 100, alpha = 0.7) +
  facet_wrap(id ~ ., ncol = 6) +
  scale_fill_manual(values = ext) +
  coord_equal() +
  labs(x = "X", y = "Y", fill = "Kasumi\ncluster") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("kasumi_sees_res_bc.pdf", width = 16, height = 8)

ids <- ids[-4]
ggplot(wc.repr %>% filter(id %in% ids) %>% mutate(cluster = as.factor(cluster), id = factor(id, levels = ids)), aes(x = x, y = y, fill = cluster)) +
  geom_tile(height = 100, width = 100, alpha = 0.7) +
  facet_wrap(id ~ ., ncol = 6) +
  scale_fill_manual(values = ext) +
  coord_equal() +
  labs(x = "X", y = "Y", fill = "Window\ncluster") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("windows_res_bc.pdf", width = 16, height = 8)

ids <- ids[-8]
ggplot(banksy.single %>% filter(id %in% ids) %>% mutate(id = factor(id, levels = ids)), aes(x = x, y = y, color = cluster)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(vars(id), ncol = 6) +
  scale_color_manual(values = as.vector(pals::cols25(n = 10))) +
  coord_equal() +
  labs(x = "X", y = "Y", color = "Banksy\nclusters") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("banksy_res_bc.pdf", width = 16, height = 8)

# Complexity ----

comp <- read_csv("comp.csv") %>%
  add_column(type = c(rep("Cell type", 117), rep("Marker", 117))) %>%
  mutate(Features = as.factor(Features))

ggplot(comp, aes(x = Cells, y = `t(sec)`, color = Features)) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~ type + Features) +
  theme_classic()

ggsave("complexity.pdf")

# ARI ----

dcis.ct <- read_rds("ari/dcis.ct.rds") %>% add_column(data = "DCIS")
ctcl.ct <- read_rds("ari/ctcl.ct.rds") %>% add_column(data = "CTCL")
bc.ct <- read_rds("ari/bc.ct.rds") %>% add_column(data = "BC")

ggplot(rbind(dcis.ct, ctcl.ct, bc.ct) %>% mutate(data = factor(data, levels = c("DCIS", "CTCL", "BC"))), aes(x = data, y = ARI.p, fill = data)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "ARI") +
  ggtitle("Cell type") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("ari_ct.pdf", width = 6, height = 6)



dcis.expr <- read_rds("ari/dcis.expr.rds") %>% add_column(data = "DCIS")
ctcl.expr <- read_rds("ari/ctcl.expr.rds") %>% add_column(data = "CTCL")
bc.expr <- read_rds("ari/bc.expr.rds") %>% add_column(data = "BC")

ggplot(rbind(dcis.expr, ctcl.expr, bc.expr) %>% mutate(data = factor(data, levels = c("DCIS", "CTCL", "BC"))), aes(x = data, y = ARI.p, fill = data)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "ARI") +
  ggtitle("Marker abundance") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("ari_expr.pdf", width = 6, height = 6)
