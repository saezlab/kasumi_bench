library(mistyR) # requires >= 1.99.4
library(kasumi)
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
library(DBI)
library(RSQLite)
library(rlist)
library(cowplot)
library(extraDistr)
library(ggrepel)
library(spdep)
library(tictoc)
# library(Banksy)
# library(SpatialExperiment)
# library(anndata)


# these two functions for highest level first and second order representation per slide
# should work in general for any labeling (binary) of the samples
# labels can be: cell types, cell neighborhoods, sliding misty labels
freq_repr <- function(labels) {
  colSums(labels) / sum(labels)
}

pa_repr <- function(labels, positions) {
  dist1 <- find_knn(positions, 1)$dist
  threshold <- 2 * sd(dist1) + mean(dist1)

  suppressMessages({
    misty.views <- create_initial_view(labels) %>%
      add_juxtaview(positions, neighbor.thr = threshold)

    neighb <- misty.views[[paste0("juxtaview.", threshold)]] %>%
      mutate(id = apply(labels, 1, which) %>% as.numeric()) %>%
      add_row(id = seq_len(ncol(labels))) %>%
      replace(is.na(.), 0) %>%
      group_by(id) %>%
      group_modify(~ colSums(.x) %>% as_tibble_row()) %>%
      ungroup() %>%
      select(-id)
  })

  (neighb / colSums(labels)) %>%
    replace(is.na(.) | . == Inf, 0) %>%
    as.matrix() %>%
    as.numeric()
}

csea_repr <- function(labels, positions) {
  dist1 <- find_knn(positions, 1)$dist
  threshold <- 2 * sd(dist1) + mean(dist1)

  message("\nPermuting")
  suppressMessages({
    misty.views <- create_initial_view(labels) %>%
      add_juxtaview(positions, neighbor.thr = threshold)

    neighb <- misty.views[[paste0("juxtaview.", threshold)]] %>%
      mutate(id = apply(labels, 1, which)) %>%
      add_row(id = seq_len(ncol(labels))) %>%
      replace(is.na(.), 0) %>%
      group_by(id) %>%
      group_modify(~ colSums(.x) %>% as_tibble_row()) %>%
      ungroup() %>%
      select(-id)


    neighb.perm <- seq_len(100) %>% map_dfr(\(i){
      with_seed(
        i,
        perm.pos <- positions %>% sample_frac()
      )

      misty.views <- create_initial_view(labels) %>%
        add_juxtaview(perm.pos, neighbor.thr = threshold)

      misty.views[[paste0("juxtaview.", threshold)]] %>%
        mutate(id = apply(labels, 1, which)) %>%
        add_row(id = seq_len(ncol(labels))) %>%
        replace(is.na(.), 0) %>%
        group_by(id) %>%
        group_modify(~ colSums(.x) %>% as_tibble_row()) %>%
        ungroup()
    })
  })

  means <- neighb.perm %>%
    group_by(id) %>%
    group_modify(~ colMeans(.x) %>% as_tibble_row()) %>%
    ungroup() %>%
    select(-id)

  sds <- neighb.perm %>%
    group_by(id) %>%
    group_modify(~ apply(.x, 2, sd) %>% as_tibble_row()) %>%
    ungroup() %>%
    select(-id)

  # symmetric z-scores
  z <- apply(((neighb - means) / sds), 2, replace_na, 0)
  z.serial <- ((z + t(z)) / 2)[upper.tri(z, diag = T)]
  replace(z.serial, is.infinite(z.serial) | is.nan(z.serial), 0)
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

cn_train <- function(all.cells, all.positions, k) {
  concat <- seq_along(all.cells) %>% map_dfr(\(i){
    misty.views <- create_initial_view(all.cells[[i]]) %>%
      add_paraview(all.positions[[i]], k, family = "constant", cache = TRUE)
    misty.views[[paste0("paraview.", k)]] %>%
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




banksy_labels <- function(all.cells, all.positions, k, l) {
  spe_list <- map2(
    names(all.cells), all.positions,
    \(id, pos) SpatialExperiment::SpatialExperiment(
      assay = list(counts = t(all.cells[[id]])),
      spatialCoords = as.matrix(pos),
      sample_id = paste0("banksy", id)
    )
  )

  spe_list <- map(spe_list, \(spe) Banksy::computeBanksy(spe,
    assay_name = "counts",
    compute_agf = TRUE, k_geom = k
  ))
  spe_joint <- do.call(cbind, spe_list)
  rm(spe_list)

  # both cell type and niches
  # l <- c(0.2,0.8)

  use_agf <- FALSE
  spe_joint <- Banksy::runBanksyPCA(spe_joint,
    use_agf = use_agf, lambda = l,
    group = "sample_id", seed = 1000
  )
  spe_joint <- Banksy::runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = l, seed = 1000)
  spe_joint <- Banksy::clusterBanksy(spe_joint, use_agf = use_agf, lambda = l, seed = 1000, resolution = seq(0.5, 0.9, 0.1))
  spe_list <- lapply(
    paste0("banksy", names(all.cells)),
    \(id) spe_joint[, spe_joint$sample_id == id]
  )
  rm(spe_joint)

  result.list <- seq(2, 6) %>% map(\(res){
    spe_list %>%
      map_dfr(\(spe){
        labels <- table(colData(spe)[, res])
        labels / sum(labels)
      }) %>%
      add_column(id = names(all.cells))
  })
}


sm_train <- function(all.cells, all.positions, l, window, minu, db.file,
                     family = "constant") {
  tic()
  tic.clearlog()
  write_lines(db.file, "tictoclog.txt", append=TRUE)
  if (file.exists(paste0(str_remove(db.file, ".sqm"), ".rds"))) {
    misty.results <- read_rds(paste0(str_remove(db.file, ".sqm"), ".rds"))
  } else {
    outputs <- seq_along(all.cells) %>%
      walk(\(i){
        misty.views <- create_initial_view(all.cells[[i]]) %>%
          add_paraview(all.positions[[i]], l,
            family = family, cached = TRUE,
            prefix = ifelse(family == "constant", "p.", "")
          )
        tic(quiet = TRUE)
        folders <- run_kasumi(misty.views, all.positions[[i]], window,
          sample.id = paste0("sample", names(all.cells)[i]),
          results.db = db.file,
          bypass.intra = (family == "constant"),
          cv.strict = (family != "constant"),
          minu = minu
        )
        toc(log = TRUE)
      })
    write_lines(unlist(tic.log(format = TRUE)), "tictoclog.txt", append = TRUE)
    toc()
    
    misty.results <- collect_results(db.file)
    write_rds(misty.results, paste0(str_remove(db.file, ".sqm"), ".rds"), "gz")
  }
  return(misty.results)
}

misty_train <- function(all.cells, all.positions, l, db.file, family = "constant") {
  if (file.exists(paste0(str_remove(db.file, ".sqm"), ".rds"))) {
    misty.results <- read_rds(paste0(str_remove(db.file, ".sqm"), ".rds"))
  } else {
    outputs <- seq_along(all.cells) %>%
      walk(\(i){
        misty.views <- create_initial_view(all.cells[[i]]) %>%
          add_paraview(all.positions[[i]], l,
            family = family, cached = TRUE,
            prefix = ifelse(family == "constant", "p.", "")
          ) %>%
          select_markers("intraview", where(~ sd(.) != 0))
        
        run_misty(misty.views,
          sample.id = paste0("sample", names(all.cells)[i]),
          results.db = db.file,
          bypass.intra = (family == "constant")
        )
        
      })
    
    write_lines(unlist(tic.log(format = TRUE)), "tictoclog.txt", append = TRUE)

    misty.results <- collect_results(db.file)
    write_rds(misty.results, paste0(str_remove(db.file, ".sqm"), ".rds"), "gz")
  }
  return(misty.results)
}

sm_repr <- function(misty.results, cuts, res, cutoff = 0, trim = 1) {
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


  return(sm.repr)
}


sm_labels <- function(misty.results, cuts, res, cutoff = 0, trim = 1, freq = TRUE) {
  sm.repr <- sm_repr(misty.results, cuts, res, cutoff, trim)

  if (!freq) {
    return(sm.repr)
  }

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
classify <- function(representation, dir = "auto") {
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

  roc(model$pred$obs, model$pred[, 3], direction = dir, quiet = TRUE)
}

optimal_smclust <- function(misty.results, true.labels) {
  grid <- seq(0.1, 0.9, 0.1) %>% future_map_dfr(\(cuts){
    seq(0.5, 0.9, 0.1) %>% map_dfr(\(res){
      freq.sm <- sm_labels(misty.results, cuts, res)

      perf <- try(
        freq.sm %>%
          left_join(true.labels, by = "id") %>%
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
# signed Model Reliance
model_reliance <- function(freq.sm) {
  model <- glm(target ~ ., freq.sm, family = "binomial")

  eorig <- classify(freq.sm)
  dir <- eorig$direction
  dir.sign <- ifelse(dir == ">", -1, 1)
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
        mutate(!!cname := freq.sm[splitr, cname] %>% unlist()), dir)$auc
    })


  mr <- dir.sign * sign(coef(model, complete = FALSE)[-1]) * (1 - eswitch) / (1 - eorig$auc)

  toreturn <- tibble(Cluster = as.factor(names(mr)), sMR = mr) %>%
    mutate(Cluster = str_remove_all(Cluster, "\\."))

  print(ggplot(toreturn %>% mutate(Cluster = fct_reorder(Cluster, sMR)), aes(x = Cluster, y = sMR)) +
    geom_segment(aes(x = Cluster, xend = Cluster, y = 0, yend = sMR)) +
    geom_point(aes(x = Cluster, y = sMR, color = sMR)) +
    scale_color_steps2(low = "darkgreen", mid = "white", high = "blue3") +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_hline(yintercept = 1, color = "gray70", linetype = "dashed") +
    geom_hline(yintercept = -1, color = "gray70", linetype = "dashed") +
    geom_label(label = paste("\u2190", eorig$levels[2]), x = length(mr), y = -1) +
    geom_label(label = paste(eorig$levels[1], "\u2192"), x = 1, y = 1) +
    coord_flip() +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    ))

  return(toreturn)
}


describe_cluster <- function(sm.repr, cluster, db.file) {
  cname <- paste0("...", cluster)
  sm.repr.all <- reduce(sm.repr, rbind)
  left <- sm.repr.all %>%
    filter(if_any(!!cname)) %>%
    select(id, x, y)

  dbcon <- dbConnect(RSQLite::SQLite(), db.file)
  samples <- dbGetQuery(dbcon, "SELECT DISTINCT sample FROM contributions") %>% unlist()
  matching <- grep(paste0("sample(", paste0(unique(left$id), collapse = "|"), ")"), samples, value = TRUE)
  dbDisconnect(dbcon)

  right <- tibble(sample = matching) %>%
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

  pattern <- paste0("(", paste0(left %>%
    left_join(right, by = c("id", "x" = "xcenter", "y" = "ycenter")) %>%
    pull(sample), collapse = "|"), ")")

  collect_results(db.file, pattern)
}



# Shared functions for interpretation

moran_score <- function(repr, clust, stride = 100, pval = FALSE) {
  sapply(repr, function(X) moran_score_sample(X, clust, stride, pval))
}

moran_score_sample <- function(sample_rep, clust, stride, pval = FALSE) {
  clust_bool <- sample_rep[[clust]]
  if (sum(clust_bool) <= 1) {
    return(NA)
  }

  x <- sample_rep$x
  y <- sample_rep$y
  distmat <- dist(cbind(x, y))
  neighbors <- as.matrix(distmat) <= stride
  diag(neighbors) <- FALSE # Self-edges
  lw <- mat2listw(neighbors, style = "W")

  if (pval) {
    I <- moran.test(as.numeric(clust_bool), lw)
    return(c(I$estimate[1], I$estimate[2], I$p.value))
  } else {
    I <- moran(clust_bool, lw, length(lw$n), Szero(lw))[1]
    return(unlist(I))
  }
}


export_anndata <- function(all.cells, all.positions, filename = "export.h5ad") {
  concat <- reduce2(names(all.cells), all.positions, \(acc, id, pos)
  list(
    X = rbind(acc[["X"]], all.cells[[id]]),
    spatial = rbind(acc[["spatial"]], pos),
    sample = c(acc[["sample"]], rep(id, nrow(pos)))
  ),
  .init = list(
    X = data.frame(),
    spatial = data.frame(),
    sample = c()
  )
  )
  to.export <- anndata::AnnData(
    X = as.matrix(concat[["X"]]),
    obs = data.frame(sample = concat[["sample"]]),
    obsm = list(spatial = as.matrix(concat[["spatial"]]))
  )
  anndata::write_h5ad(to.export, filename)
}


bin_count_cluster <- function(all.expr, all.positions, window, overlap, k, freq = TRUE){
  all.windows <- seq_along(all.expr) %>% map_dfr(\(i){
    expr <- all.expr[[i]]
    positions <- all.positions[[i]]
    
    x <- tibble::tibble(
      xl = seq(
        min(positions[, 1]),
        max(positions[, 1]),
        window - window * overlap / 100
      ),
      xu = xl + window
    ) %>%
      dplyr::filter(xl < max(positions[, 1])) %>%
      dplyr::mutate(xu = pmin(xu, max(positions[, 1]))) %>%
      round(2)
    
    y <- tibble::tibble(
      yl = seq(
        min(positions[, 2]),
        max(positions[, 2]),
        window - window * overlap / 100
      ),
      yu = yl + window
    ) %>%
      dplyr::filter(yl < max(positions[, 2])) %>%
      dplyr::mutate(yu = pmin(yu, max(positions[, 2]))) %>%
      round(2)
    
    tiles <- tidyr::expand_grid(x, y)
    
    tiles %>%
      pmap_dfr(\(xl, xu, yl, yu){
        selected.rows <- which(
          positions[, 1] >= xl & positions[, 1] <= xu &
            positions[, 2] >= yl & positions[, 2] <= yu
        )
        
        c(expr %>%
          slice(selected.rows) %>%
          colSums(), x = (xu+xl)/2, y = (yu+yl)/2)
      }) %>%
      add_column(id = names(all.expr)[i])
  })
  
  clusters <- kmeans_ondist(all.windows %>% select(-id,x,y), k)
  
  wc.repr <- cbind(cluster = clusters, all.windows %>% select(id,x,y)) %>%
    as_tibble()
  
  if(!freq) return(wc.repr)
  
  totals <- wc.repr %>% select(-c(x,y)) %>%
    group_by(id, cluster) %>%
    tally() %>%
    ungroup() %>%
    pivot_wider(names_from = "cluster", values_from = "n")  %>%
    replace(is.na(.), 0)
  
  t(apply(totals %>% select(-id), 1, \(x) x/sum(x))) %>% as_tibble() %>% add_column(id = totals$id)
  
}


