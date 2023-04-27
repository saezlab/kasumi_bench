library(mistyR) # requires >= 1.99.0
library(future)
library(tidyverse)
library(furrr)
library(knn.covertree)
library(proxy)
library(igraph)
library(readxl)
library(ranger)
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


markov_repr <- function(labels, positions) {
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

# representation can be cns or sliding misty signatures
leiden_onsnn <- function(representation, nn = 30, minjac = 0.15, resolution = 0.8) {
  neighbors <- find_knn(representation, nn, distance = "cosine")
  with_seed(
    1,
    snn <- simil(neighbors$index, "eJaccard")
  )
  snn[snn < minjac] <- 0

  groups <-
    graph.adjacency(snn %>% as.matrix(), mode = "undirected", weighted = TRUE) %>%
    cluster_leiden(resolution_parameter = resolution, n_iterations = -1)

  return(groups$membership)
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

# cns should be calculated on all samples
# names(all.cells) used as sample ids
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


  # return one-hot representation for all all.cells and if needed all.positions
  suppressMessages(
    cn.repr <- map(clusters, ~ .x == seq(length(unique(clusters)))) %>% reduce(rbind) %>%
      as_tibble(.name_repair = "unique") %>% cbind(samps) %>% group_by(id) %>%
      group_split()
  )

  return(cn.repr)
}


sm_train <- function(all.cells, all.positions, l, window, minu, minm, top.folder,
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
            minu = minu, minm = minm,
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

  # think about gain.R2 as node weights
  clusters <- leiden_onsim(clean, cuts, res)

  suppressMessages(
    sm.repr <- map(clusters, ~ .x == seq(length(unique(clusters)))) %>% reduce(rbind) %>%
      as_tibble(.name_repair = "unique") %>% cbind(samps) %>% group_by(id) %>% rename(x = xcenter, y = ycenter) %>%
      group_split()
  )

  return(sm.repr)
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
classify <- function(representation, plot = FALSE) {
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

# the column target in the representation table is the ground truth
# returns macro F1 per class based on 10-fold cv predictions
classify_rf <- function(representation) {
  model <- ranger(target ~ ., representation,
    seed = 1,
    classification = TRUE, probability = TRUE
  )

  roc(representation$target, model$predictions[, 1], quiet = TRUE, smooth = TRUE, smooth.n = 10)
}

optimal_smclust <- function(misty.results, true.labels, funct = classify) {
    grid <- seq(0.1, 0.9, 0.1) %>% future_map_dfr(\(cuts){
      seq(0.5, 0.9, 0.1) %>% map_dfr(\(res){
        sm.repr <- sm_labels(misty.results, cuts, res)
  
        repr.ids <- sm.repr %>% map_chr(~ .x$id[1])
  
        perf <- try(
          sm.repr %>% map_dfr(~ .x %>%
            select(-c(id, x, y)) %>%
            freq_repr()) %>%
            select(where(~ (sd(.) > 1e-3) & (sum(. > 0) >= max(5, 0.1 * length(.))))) %>%
            add_column(id = repr.ids) %>% left_join(true.labels, by = "id") %>%
            drop_na() %>% select(-id) %>% funct()
        )
  
        auc <- ifelse(class(perf) == "try-error", 0, perf$auc)
        print(paste(cuts, res, as.numeric(auc)))
        tibble_row(cut = cuts, res = res, perf = as.numeric(auc))
      })
    }, .options = furrr_options(seed = NULL))
  grid[which.max(grid$perf), ] %>% unlist()
}
