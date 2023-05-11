library(mistyR) # requires >= 1.99.2
library(future)
library(tidyverse)
library(furrr)
library(knn.covertree)
library(proxy)
library(igraph)
library(readxl)
library(caret)
library(pROC)
library(ClusterR)
library(withr)


# these two functions for highest level first and second order representation per slide
# should work in general for any labeling (binary) of the samples
# labels can be: cell types, cell neighborhoods, sliding misty labels
freq_repr <- function(labels) {
  colSums(labels) / sum(labels)
}

pa_repr <- function(labels, positions) {
  dist1 <- find_knn(positions, 1)$dist
  threshold <- 2 * sd(dist1) + mean(dist1)
  misty.views <- create_initial_view(labels) %>%
    add_juxtaview(positions, neighbor.thr = threshold)
  suppressMessages(
    neighb <- misty.views[[paste0("juxtaview.", threshold)]]$data %>%
      mutate(id = apply(labels, 1, which)) %>%
      add_row(id = seq_len(ncol(labels))) %>%
      replace(is.na(.), 0) %>%
      group_by(id) %>%
      group_modify(~ colSums(.x) %>% as_tibble_row()) %>%
      ungroup() %>%
      select(-id)
  )

  (neighb / colSums(labels)) %>%
    replace(is.na(.), 0) %>%
    as.matrix() %>%
    as.numeric()
}

leiden_onsim <- function(representation, minsim = 0.8, resolution = 0.8, measure = "cosine") {
  sim <- simil(representation, measure)

  sim[sim < minsim] <- 0
  with_seed(
    1,
    groups <-
      graph.adjacency(sim %>% as.matrix(), mode = "undirected", weighted = TRUE) %>%
      cluster_leiden(resolution_parameter = resolution, n_iterations = -1)
  )

  return(groups$membership)
}

kmeans_ondist <- function(representation, k = 10) {
  with_seed(
    1,
    clust <- KMeans_rcpp(representation, k)
  )
  return(clust$clusters)
}

cn_train <- function(all.cells, all.positions) {
  concat <- seq_along(all.cells) %>% map_dfr(\(i){
    misty.views <- create_initial_view(all.cells[[i]]) %>%
      add_paraview(all.positions[[i]], 10, family = "constant", cache = TRUE)
    misty.views[["paraview.10"]]$data %>%
      add_column(id = names(all.cells)[i]) %>%
      cbind(all.positions[[i]])
  })

  return(concat)
}


cn_labels <- function(neighborhoods, k) {
  # clusters <- leiden_onsim(neighborhoods %>% select(-c(id, x, y)), cut, res)
  clusters <- kmeans_ondist(neighborhoods %>% select(-c(id, x, y)), k)
  samps <- neighborhoods %>% select(id, x, y)

  suppressMessages(
    cn.repr <- map(clusters, ~ .x == seq(length(unique(clusters)))) %>% reduce(rbind) %>%
      as_tibble(.name_repair = "unique") %>% cbind(samps) %>% group_by(id) %>%
      group_split()
  )

  repr.ids.cn <- cn.repr %>% map_chr(~ .x$id[1])

  freq.cn <- cn.repr %>%
    map_dfr(~ .x %>%
      select(-c(id, x, y)) %>%
      freq_repr()) %>%
    select(where(~ (sd(.) > 1e-3))) %>%
    add_column(id = repr.ids.cn)

  return(freq.cn)
}


sm_train <- function(all.cells, all.positions, l, window, minu, top.folder,
                     family = "constant") {
  if (file.exists(paste0(top.folder, ".rds"))) {
    misty.results <- read_rds(paste0(top.folder, ".rds"))
  } else {
    outputs <- seq_along(all.cells) %>%
      map(\(i){
        if (dir.exists(paste0(top.folder, "/sample", names(all.cells)[i], "/"))) {
          folders <- list.dirs(paste0(top.folder, "/sample", names(all.cells)[i], "/"))[-1]
        } else {
          misty.views <- create_initial_view(all.cells[[i]]) %>%
            add_paraview(all.positions[[i]], l,
              family = family, cached = TRUE,
              prefix = ifelse(family == "constant", "p.", "")
            )

          folders <- run_sliding_misty(misty.views, all.positions[[i]], window,
            minu = minu,
            results.folder = paste0(top.folder, "/sample", names(all.cells)[i], "/"),
            bypass.intra = (family == "constant"),
            cv.strict = (family != "constant")
          )
        }

        ifelse(length(folders) >= 10, return(folders), return(NA))
      }) %>%
      unlist()

    misty.results <- collect_results(outputs[!is.na(outputs)])
    write_rds(misty.results, paste0(top.folder, ".rds"), "gz")
  }
  return(misty.results)
}

misty_train <- function(all.cells, all.positions, l, top.folder, family = "constant") {
  if (file.exists(paste0(top.folder, "_ws.rds"))) {
    misty.results <- read_rds(paste0(top.folder, "_ws.rds"))
  } else {
    outputs <- seq_along(all.cells) %>%
      map(\(i){
        misty.views <- create_initial_view(all.cells[[i]]) %>%
          add_paraview(all.positions[[i]], l,
            family = family, cached = TRUE,
            prefix = ifelse(family == "constant", "p.", "")
          ) %>%
          select_markers("intraview", where(~ sd(.) != 0))

        run_misty(misty.views,
          results.folder = paste0(top.folder, "_ws/sample", names(all.cells)[i], "/"),
          bypass.intra = (family == "constant")
        )
      }) %>%
      unlist()

    misty.results <- collect_results(outputs)
    write_rds(misty.results, paste0(top.folder, "_ws.rds"), "gz")
  }
  return(misty.results)
}


sm_labels <- function(misty.results, cuts, res, cutoff = 0, trim = 1) {
  # trimming matters
  sig <- extract_signature(misty.results, type = "i", intersect.targets = FALSE, trim = trim)
  sig[is.na(sig)] <- floor(min(sig %>% select(-sample), na.rm = TRUE))
  sig <- sig %>% mutate(across(!sample, ~ ifelse(.x <= cutoff, 0, .x)))

  keep <- which(sig %>% select(-sample, -contains("intra_")) %>% rowSums() != 0)

  samps <- sig %>%
    slice(keep) %>%
    select(sample) %>%
    mutate(
      id = str_extract(sample, "sample.*/") %>% str_remove("/") %>% str_remove("^sample"),
      box = str_extract(sample, "/[0-9].*$") %>% str_remove("/")
    ) %>%
    rowwise(id) %>%
    summarize(rebox = box %>%
      str_split("_", simplify = T) %>%
      as.numeric() %>% list(), .groups = "drop") %>%
    rowwise(id) %>%
    summarize(
      xcenter = (rebox[3] + rebox[1]) / 2,
      ycenter = (rebox[4] + rebox[2]) / 2, .groups = "drop"
    )

  # the filtering here also matters
  clean <- sig %>%
    select(-sample, -contains("_.novar")) %>%
    slice(keep) %>%
    select(where(~ sum(.) != 0))


  clusters <- leiden_onsim(clean, cuts, res)

  suppressMessages(
    sm.repr <- map(clusters, ~ .x == seq(length(unique(clusters)))) %>% reduce(rbind) %>%
      as_tibble(.name_repair = "unique") %>% cbind(samps) %>% group_by(id) %>%
      rename(x = xcenter, y = ycenter) %>%
      group_split()
  )

  repr.ids <- sm.repr %>% map_chr(~ .x$id[1])

  # technically we can also use pa_repr instead of freq_repr or combine both
  freq.sm <- sm.repr %>%
    map_dfr(~ .x %>%
      select(-c(id, x, y)) %>%
      freq_repr()) %>%
    select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.))))) %>%
    add_column(id = repr.ids, .before = 1)

  return(freq.sm)
}

misty_labels <- function(misty.results, cutoff = 0, trim = 1) {
  sig <- extract_signature(misty.results, type = "i", intersect.targets = FALSE, trim = trim)
  sig[is.na(sig)] <- floor(min(sig %>% select(-sample), na.rm = TRUE))
  sig <- sig %>% mutate(across(!sample, ~ ifelse(.x <= cutoff, 0, .x)))

  keep <- which(sig %>% select(-sample, -contains("intra_")) %>% rowSums() != 0)

  ws.repr <- sig %>%
    mutate(sample = str_extract(sample, "sample.*") %>%
      str_remove("^sample")) %>%
    select(-contains("_.novar")) %>%
    slice(keep) %>%
    column_to_rownames("sample") %>%
    select(where(~ sum(.) != 0)) %>%
    select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.))))) %>%
    rownames_to_column("sample")

  return(ws.repr)
}


# the column target in the representation table is the ground truth
# returns ROC based on 10-fold cv predictions
classify <- function(representation) {
  with_seed(
    1,
    suppressWarnings(
      model <- train(target ~ ., representation,
        method = "glm", metric = "ROC",
        trControl = trainControl(
          method = "cv", number = 10,
          classProbs = TRUE,
          summaryFunction = twoClassSummary,
          savePredictions = TRUE
        )
      )
    )
  )

  roc(model$pred$obs, model$pred[, 3], quiet = TRUE)
}

optimal_smclust <- function(misty.results, true.labels) {
  grid <- seq(0.1, 0.9, 0.1) %>% future_map_dfr(\(cuts){
    seq(0.5, 0.9, 0.1) %>% map_dfr(\(res){
      freq.sm <- sm_labels(misty.results, cuts, res)

      perf <- try(
        freq.sm %>%
          add_column(id = repr.ids) %>% left_join(true.labels, by = "id") %>%
          drop_na() %>% select(-id) %>% classify()
      )

      auc <- ifelse(class(perf) == "try-error", 0, perf$auc)
      print(paste(cuts, res, as.numeric(auc)))
      tibble_row(cut = cuts, res = res, perf = as.numeric(auc))
    })
  }, .options = furrr_options(seed = NULL))
  grid[which.max(grid$perf), ] %>% unlist()
}


# Fisher, Rudin, Dominici, JMLR, 2019
model_reliance <- function(freq.sm) {
  model <- glm(target ~ ., freq.sm, family = "binomial")

  eorig <- classify(freq.sm)
  cat(paste0("AUC: ", eorig$auc))

  with_seed(
    1,
    splitr <- runif(nrow(freq.sm)) %>% rank()
  )

  nas <- names(which(is.na(coef(model)[-1])))

  eswitch <- freq.sm %>%
    select(-target, -nas) %>%
    colnames() %>%
    map_dbl(\(cname){
      classify(freq.sm %>% select(-nas) %>%
        mutate(!!cname := freq.sm[splitr, cname] %>% unlist()))$auc
    })

  mr <- sign(coef(model, complete = FALSE)[-1]) * (1 - eswitch) / (1 - eorig$auc)

  ggplot(tibble(Cluster = as.factor(names(mr)), sMR = mr) %>%
    mutate(Cluster = str_remove_all(Cluster, "\\.")) %>%
    mutate(Cluster = fct_reorder(Cluster, sMR)), aes(x = Cluster, y = sMR)) +
    geom_segment(aes(x = Cluster, xend = Cluster, y = 0, yend = sMR)) +
    geom_point(aes(x = Cluster, y = sMR, color = sMR)) +
    scale_color_steps2(low = "darkgreen", mid = "white", high = "blue3") +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_hline(yintercept = 1, color = "gray70", linetype = "dashed") +
    geom_hline(yintercept = -1, color = "gray70", linetype = "dashed") +
    geom_label(label = paste("\u2190", ifelse(eorig$direction == ">", eorig$levels[2], eorig$levels[1])), x = length(mr), y = -1) +
    geom_label(label = paste(ifelse(eorig$direction == ">", eorig$levels[1], eorig$levels[2]), "\u2192"), x = 1, y = 1) +
    coord_flip() +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}


describe_cluster <- function(sm.repr, cluster, folder.prefix) {
  cname <- paste0("...", cluster)
  sm.repr.all <- reduce(sm.repr, rbind)
  left <- sm.repr.all %>%
    filter(if_any(!!cname)) %>%
    select(id, x, y)

  all.folders <- list.files(folder.prefix, paste0("sample(", paste0(unique(left$id), collapse = "|"), ")"), full.names = TRUE) %>%
    map(~ list.dirs(.)[-1]) %>%
    unlist()

  right <- tibble(sample = all.folders) %>%
    mutate(
      id = str_extract(sample, "sample.*/") %>% str_remove("/") %>% str_remove("^sample"),
      box = str_extract(sample, "/[0-9].*$") %>% str_remove("/")
    ) %>%
    rowwise(id) %>%
    summarize(sample = sample, rebox = box %>%
      str_split("_", simplify = T) %>%
      as.numeric() %>% list(), .groups = "drop") %>%
    rowwise(id) %>%
    summarize(
      sample = sample,
      xcenter = (rebox[3] + rebox[1]) / 2,
      ycenter = (rebox[4] + rebox[2]) / 2, .groups = "drop"
    )

  left %>%
    left_join(right, by = c("id", "x" = "xcenter", "y" = "ycenter")) %>%
    pull(sample) %>%
    collect_results()
}
