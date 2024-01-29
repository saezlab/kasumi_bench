source("utils.R")

# MIBI  ----
ct <- read_rds("rocs/dcis.ct.rds")
expr <- read_rds("rocs/dcis.expr.rds")
ws <- read_rds("rocs/dcis.ws.rds")
expr.base  <- read_rds("rocs/dcis.expr.base.rds")

resp <- read_csv("data/DCISMIBI/Tissue Feature Data/Table_S1_Patient_Feature_Table.csv")

ggroc(list.append(list.remove(ct, "pa"), ws.ct = ws$ct), legacy.axes = TRUE, linewidth = pi/10) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2/3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("MIBI DCIS") + theme_light()

ggsave("roc.mibi.pdf", width = 4, height = 3)

ggroc(list(sm.ct = ct$sm, ws.ct = ws$ct, sm.expr = expr, ws.expr = ws$expr,  base = expr.base), 
      legacy.axes = TRUE, linewidth = pi/10) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2/3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("MIBI DCIS") + theme_light()

ggsave("roc.misty.mibi.pdf", width = 4, height = 3)

misty.results <- read_rds("DCISct200.rds")
sm.repr <- sm_labels(misty.results, 0.8, 0.5)

freq.sm <- sm.repr %>%
  left_join(resp %>% select(PointNumber, Status) %>%
    mutate(PointNumber = as.character(PointNumber)), by = c("id" = "PointNumber")) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  select(-id)

mr <- model_reliance(freq.sm)

ggsave("DCISreliance.pdf", height = 5, width = 4)

sm.repr.ext <- sm_labels(misty.results, 0.8, 0.5, freq = FALSE)

mr.clusters <- mr %>% filter(abs(sMR) >= 1) %>% arrange(-sMR) %>% pull(Cluster)

cl.list <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "DCISct200.sqm")
  plot_interaction_heatmap(misty.cluster, "paraview.10", trim = 1, cutoff = 1, clean = TRUE)
  last_plot()
})

plot_grid(plotlist = cl.list, labels = mr.clusters)

ggsave("clusterfigs/DCISct200.pdf", width = 20, height = 20)


misty.results <- read_rds("DCISexpr200.rds")
 
sm.repr <- sm_labels(misty.results, 0.6, 0.8)

freq.sm <- sm.repr %>%
  left_join(resp %>% select(PointNumber, Status) %>%
    mutate(PointNumber = as.character(PointNumber)), by = c("id" = "PointNumber")) %>%
  rename(target = Status) %>%
  mutate(target = as.factor(target)) %>%
  select(-id)

mr <- model_reliance(freq.sm)
ggsave("DCISreliance_expr.pdf", height = 5, width = 4)

sm.repr.ext <- sm_labels(misty.results, 0.6, 0.8, freq = FALSE)
mr.clusters <- mr %>% filter(abs(sMR) >= 1) %>% arrange(-sMR) %>% pull(Cluster)

cl.list <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "DCISexpr200.sqm")
  plot_interaction_heatmap(misty.cluster, "paraview.100", trim = 1, cutoff = 0.5, clean = TRUE)
  last_plot()
}) 

plot_grid(plotlist = cl.list, labels = mr.clusters)

ggsave("clusterfigs/DCISexpr200.pdf", width = 10, height = 8)

cl.list <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "DCISexpr200.sqm")
  plot_interaction_heatmap(misty.cluster, "intraview", trim = 1, cutoff = 0.5, clean = TRUE)
  last_plot()
}) 

plot_grid(plotlist = cl.list, labels = mr.clusters)

ggsave("clusterfigs/DCISexpr200.intra.pdf", width = 10, height = 8)


# CODEX ----

ct <- read_rds("rocs/ctcl.ct.rds")
expr <- read_rds("rocs/ctcl.expr.rds")
ws <- read_rds("rocs/ctcl.ws.rds")
expr.base  <- read_rds("rocs/ctcl.expr.base.rds")

lymph <- read_csv("data/LymphomaCODEX/single_cells.csv")

outcome <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  group_by(Spots, Patients) %>%
  summarize(Groups = Groups[1], .groups = "drop")

ggroc(list.append(list.remove(ct, "pa"), ws.ct = ws$ct), legacy.axes = TRUE, linewidth = pi/10) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2/3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("CODEX CTCL") + theme_light()

ggsave("roc.codex.pdf", width = 4, height = 3)

ggroc(list(sm.ct = ct$sm, ws.ct = ws$ct, sm.expr = expr, ws.expr = ws$expr, base = expr.base), 
      legacy.axes = TRUE, linewidth = pi/10) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2/3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("CODEX CTCL") + theme_light()

ggsave("roc.misty.codex.pdf", width = 4, height = 3)


misty.results <- read_rds("CTCLct400.rds")
sm.repr <- sm_labels(misty.results, 0.4, 0.9)

freq.sm <- sm.repr %>%
  left_join(outcome %>%
              mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)

mr <- model_reliance(freq.sm)

ggsave("CTCLreliance.pdf", height = 6, width = 4)

sm.repr.ext <- sm_labels(misty.results, 0.4, 0.9, freq = FALSE)

mr.clusters <- mr %>% filter(abs(sMR) >= 1) %>% arrange(-sMR) %>% pull(Cluster)

cl.list <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "CTCLct400.sqm")
  plot_interaction_heatmap(misty.cluster, "paraview.10", trim = 1, cutoff = 1, clean = TRUE)
  last_plot()
})

plot_grid(plotlist = cl.list, labels = mr.clusters)

ggsave("clusterfigs/CTCLct400.pdf", width = 16, height = 12)

cl.list.cor <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "CTCLct400.sqm")
  plot_interaction_heatmap(misty.cluster, "paraview.10", trim = 1, cutoff = 1, 
                           clean = TRUE, correlation = TRUE)
  last_plot()
})

plot_grid(plotlist = cl.list.cor, labels = mr.clusters)

ggsave("clusterfigs/CTCLct400_cor.pdf", width = 16, height = 12)

misty.cluster.23 <- describe_cluster(sm.repr.ext, 23, "CTCLct400.sqm")

plot_improvement_stats(misty.cluster.23, trim = 1)
ggsave("clusterfigs/ctcl23.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.23, "paraview.10", trim = 1, cutoff = 1, clean = TRUE)
ggsave("clusterfigs/ctcl23h.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.23, "paraview.10", trim = 1, cutoff = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/ctcl23c.pdf", width = 4, height = 4)

misty.cluster.23$importances.aggregated <- misty.cluster.23$importances.aggregated %>% 
  mutate(Predictor = str_remove(Predictor, "^p\\."))

plot_interaction_communities(misty.cluster.23, view = "paraview.10", 
                             trim = 1, cutoff = 1, path = "clusterfigs/ctcl23.graphml")

misty.cluster.8 <- describe_cluster(sm.repr.ext, 8, "CTCLct400.sqm")

plot_improvement_stats(misty.cluster.8, trim = 1)
ggsave("clusterfigs/ctcl8.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.8, "paraview.10", trim = 1, cutoff = 1, clean = TRUE)
ggsave("clusterfigs/ctcl8h.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.8, "paraview.10", trim = 1, cutoff = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/ctcl8c.pdf", width = 4, height = 4)

misty.cluster.8$importances.aggregated <- misty.cluster.8$importances.aggregated %>% 
  mutate(Predictor = str_remove(Predictor, "^p\\."))

plot_interaction_communities(misty.cluster.8, view = "paraview.10", 
                             trim = 1, cutoff = 1, path = "clusterfigs/ctcl8.graphml")


misty.cluster.4 <- describe_cluster(sm.repr.ext, 4, "CTCLct400.sqm")

plot_improvement_stats(misty.cluster.4, trim =1)
ggsave("clusterfigs/ctcl4.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.4, "paraview.10", trim = 1, cutoff = 1, clean = TRUE)
ggsave("clusterfigs/ctcl4h.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.4, "paraview.10", trim = 1, cutoff = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/ctcl4c.pdf", width = 4, height = 4)

misty.cluster.4$importances.aggregated <- misty.cluster.4$importances.aggregated %>% 
  mutate(Predictor = str_remove(Predictor, "^p\\."))

plot_interaction_communities(misty.cluster.4, view = "paraview.10", 
                             trim = 1, cutoff = 1, path = "clusterfigs/ctcl4.graphml")


misty.cluster.5 <- describe_cluster(sm.repr.ext, 5, "CTCLct400.sqm")

plot_improvement_stats(misty.cluster.5, trim =1)
ggsave("clusterfigs/ctcl5.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.5, "paraview.10", trim = 1, cutoff = 1, clean = TRUE)
ggsave("clusterfigs/ctcl5h.pdf", width = 4, height = 4)
plot_interaction_heatmap(misty.cluster.5, "paraview.10", trim = 1, cutoff = 1, clean = TRUE, correlation = TRUE)
ggsave("clusterfigs/ctcl5c.pdf", width = 4, height = 4)

misty.cluster.5$importances.aggregated <- misty.cluster.5$importances.aggregated %>% 
  mutate(Predictor = str_remove(Predictor, "^p\\."))

plot_interaction_communities(misty.cluster.5, view = "paraview.10", 
                             trim = 1, cutoff = 1, path = "clusterfigs/ctcl5.graphml")


misty.results <- read_rds("CTCLexpr400.rds")
sm.repr <- sm_labels(misty.results, 0.3, 0.5)

freq.sm <- sm.repr %>%
  left_join(outcome %>%
              mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)

mr <- model_reliance(freq.sm)

ggsave("CTCLreliance_expr.pdf", height = 6, width = 4)


sm.repr.ext <- sm_labels(misty.results, 0.3, 0.5, freq = FALSE)

mr.clusters <- mr %>% filter(abs(sMR) >= 1) %>% arrange(-sMR) %>% pull(Cluster)

cl.list <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "CTCLexpr400.sqm")
  plot_interaction_heatmap(misty.cluster, "paraview.100", trim = 1, cutoff = 0.5, clean = TRUE)
  last_plot()
})

plot_grid(plotlist = cl.list, labels = mr.clusters)

ggsave("clusterfigs/CTCLexpr400.pdf", width = 15, height = 8)

cl.list <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "CTCLexpr400.sqm")
  plot_interaction_heatmap(misty.cluster, "intraview", trim = 1, cutoff = 1, clean = TRUE)
  last_plot()
})

plot_grid(plotlist = cl.list, labels = mr.clusters)

ggsave("clusterfigs/CTCLexpr400.intra.pdf", width = 15, height = 8)


# IMC ----

ct <- read_rds("rocs/bc.ct.rds")
expr <- read_rds("rocs/bc.expr.rds")
ws <- read_rds("rocs/bc.ws.rds")
expr.base  <- read_rds("rocs/bc.expr.base.rds")

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

ggroc(list.append(list.remove(ct, "pa"), ws.ct = ws$ct), legacy.axes = TRUE, linewidth = pi/10) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2/3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("IMC BC") + theme_light()

ggsave("roc.imc.pdf", width = 4, height = 3)

ggroc(list(sm.ct = ct$sm, ws.ct = ws$ct, sm.expr = expr, ws.expr = ws$expr, base = expr.base), 
      legacy.axes = TRUE, linewidth = pi/10) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", linetype = "dotted", linewidth = 2/3) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL) +
  ggtitle("IMC BC") + theme_light()

ggsave("roc.misty.imc.pdf", width = 4, height = 3)

misty.results <- read_rds("BCct200.rds")

sm.repr <- sm_labels(misty.results, 0.3, 0.8)

freq.sm <- sm.repr %>%
  left_join(resp, by = c("id" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id)

mr <- model_reliance(freq.sm)
ggsave("BCreliance.pdf", height = 5, width = 4)

sm.repr.ext <- sm_labels(misty.results, 0.3, 0.8, freq = FALSE)

mr.clusters <- mr %>% filter(abs(sMR) >= 1) %>% arrange(-sMR) %>% pull(Cluster)

cl.list <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "BCct200.sqm")
  plot_interaction_heatmap(misty.cluster, "paraview.10", trim = 1, cutoff = 1, clean = TRUE)
  last_plot()
})

plot_grid(plotlist = cl.list, labels = mr.clusters)

ggsave("clusterfigs/BCct200.pdf", width = 12, height = 12)



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

mr.clusters <- mr %>% filter(abs(sMR) >= 1) %>% arrange(-sMR) %>% pull(Cluster)

cl.list <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "BCexpr200.sqm")
  plot_interaction_heatmap(misty.cluster, "paraview.50", trim = 1, cutoff = 0.5, clean = TRUE)
  last_plot()
})

plot_grid(plotlist = cl.list, labels = mr.clusters)

ggsave("clusterfigs/BCexpr200.pdf", width = 20, height = 16)

cl.list <- mr.clusters %>% map(\(cl){
  misty.cluster <- describe_cluster(sm.repr.ext, cl, "BCexpr200.sqm")
  plot_interaction_heatmap(misty.cluster, "intraview", trim = 1, cutoff = 1, clean = TRUE)
  last_plot()
})

plot_grid(plotlist = cl.list, labels = mr.clusters)

ggsave("clusterfigs/BCexpr200.intra.pdf", width = 20, height = 16)



misty.cluster.5 <- describe_cluster(sm.repr.ext, 5, "BCexpr200.sqm")

plot_improvement_stats(misty.cluster.5, trim = 1)
ggsave("clusterfigs/bc5.pdf", width = 5, height = 4)
plot_interaction_heatmap(misty.cluster.5, "intraview", trim = 1, cutoff = 1, clean = TRUE)
ggsave("clusterfigs/bc5h_intra.pdf", width = 5, height = 5)
plot_interaction_heatmap(misty.cluster.5, "paraview.50", trim = 1, cutoff = 0.5, clean = TRUE)
ggsave("clusterfigs/bc5h.pdf", width = 5, height = 4)


misty.cluster.37 <- describe_cluster(sm.repr.ext, 37, "BCexpr200.sqm")

plot_improvement_stats(misty.cluster.37, trim = 1)
ggsave("clusterfigs/bc37.pdf", width = 5, height = 4)
plot_interaction_heatmap(misty.cluster.37, "intraview", trim = 1, cutoff = 1, clean = TRUE)
ggsave("clusterfigs/bc37h_intra.pdf", width = 5, height = 4)
plot_interaction_heatmap(misty.cluster.37, "paraview.50", trim = 1, cutoff = 0.5, clean = TRUE)
ggsave("clusterfigs/bc37h.pdf", width = 5, height = 4)

# Sensitivity ----

ws.dcis <- read_rds("rocs/wrocs.dcis.rds")
ws.ctcl <- read_rds("rocs/wrocs.ctcl.rds")
ws.bc <- read_rds("rocs/wrocs.bc.rds")

ws <- tibble(Window = seq(100,500, by = 100), DCIS = ws.dcis, CTCL = ws.ctcl, BC = ws.bc) %>%
  pivot_longer(-Window, names_to = "Data", values_to = "AUROC")

ggplot(ws, aes(x = Window, y = AUROC, color = Data)) + geom_point() + geom_line() + 
  scale_color_brewer(palette = "Set1") + ylim(c(0.5,1)) + theme_light()

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




# Kasumi sees ----
test <- lymph %>% filter(Spots %in% c(1,13)) %>% select(Spots,ClusterName,X,Y)
ggplot(test, aes(x = X, y=Y, color = ClusterName)) + geom_point(size=0.5) + facet_wrap(vars(Spots)) + 
  scale_color_manual(values=as.vector(pals::glasbey(n=19))) + coord_equal() + 
  theme_classic() + theme(legend.position = "bottom")

ggsave("ct_distro.pdf", width=8, height=6)

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


test <- rbind(sm.repr[[1]],sm.repr[[3]])
clusters <- test %>% select(-c(id,x,y)) %>% apply(1,which)

kasumirep <- test %>% select(c(id,x,y)) %>% cbind(cluster = clusters) %>% filter(cluster %in% persistent) %>% mutate(cluster = as.factor(cluster))

ggplot(kasumirep, aes(x = x, y=y, color = cluster)) + geom_point(size = 2.5) + #geom_tile(height = 200, width = 200, alpha = 0.8) + 
  scale_color_manual(values=as.vector(pals::glasbey(n=11))) + 
  facet_wrap(vars(id)) + coord_equal()  + theme_classic() + theme(legend.position = "bottom")

ggsave("kasumi_sees_tiles_dots.pdf", width = 8, height = 6)

# kasumi.freqs <- c(3,13) %>% map_dfr(~
# sm.repr[[.x]] %>% select(-c(id,x,y)) %>% apply(2,sum) %>% t() %>% as_tibble() %>% pivot_longer(everything()) %>% 
#   mutate(name = str_remove_all(name, "\\.")) %>% filter(value != 0, name %in% persistent) %>% add_column(id = sm.repr[[.x]]$id[1])
# )
# 
# ggplot(kasumi.freqs, aes(x = name, y = value, fill =id)) + 
#   geom_bar(stat="identity", position=position_dodge(1), width = 0.5) + 
#   scale_fill_brewer(palette = "Set2") +
#   theme_classic()
# 
# ct.freqs <- c(3,13) %>% map_dfr(~
#                                       all.cells.lymph[[.x]] %>% colSums() %>% t() %>% as_tibble() %>% pivot_longer(everything()) %>% 
#                                       add_column(id = as.character(.x))
# )
# 
# ggplot(ct.freqs, aes(x = name, y = value, fill = id)) + 
#   geom_bar(stat="identity", position=position_dodge(1), width = 0.5) + 
#   scale_fill_brewer(palette = "Set2") +
#   theme_classic()
# kasumi.freqs <- c(3,13) %>% map_dfr(~
# sm.repr[[.x]] %>% select(-c(id,x,y)) %>% apply(2,sum) %>% t() %>% as_tibble() %>% pivot_longer(everything()) %>% 
#   mutate(name = str_remove_all(name, "\\.")) %>% filter(value != 0, name %in% persistent) %>% add_column(id = sm.repr[[.x]]$id[1])
# )
# 
# ggplot(kasumi.freqs, aes(x = name, y = value, fill =id)) + 
#   geom_bar(stat="identity", position=position_dodge(1), width = 0.5) + 
#   scale_fill_brewer(palette = "Set2") +
#   theme_classic()
# 
# ct.freqs <- c(3,13) %>% map_dfr(~
#                                       all.cells.lymph[[.x]] %>% colSums() %>% t() %>% as_tibble() %>% pivot_longer(everything()) %>% 
#                                       add_column(id = as.character(.x))
# )
# 
# ggplot(ct.freqs, aes(x = name, y = value, fill = id)) + 
#   geom_bar(stat="identity", position=position_dodge(1), width = 0.5) + 
#   scale_fill_brewer(palette = "Set2") +
#   theme_classic()


