source("utils.R")
library(Matrix)
library(heatmaply)
library(umap)
library(spatstat.geom)

top.folder <- "DCISct"
sm.results <- read_rds(paste0(top.folder, ".rds"), "gz")

#' Calculate Interaction and Localization Scores
#'
#' This function calculates the proportion of windows in which an interaction is
#' found and how clustered or spread these windows are, based on the average 
#' nearest-neighbor distance between such windows.
#'
#' @param clean A data frame containing the clean interaction data, with one
#' interaction per column.
#' @param inter The name of the column representing the interaction values.
#' @param windows.coord The coordinates and origin of each window.
#' 
#' @return A numeric vector containing the interaction and localization scores.
#'   - The interaction score represents the fraction of windows having 
#'   the interaction.
#'   - The localization score represents the average nearest neighbor distance 
#'   (nn-distance) for windows with interaction.
#'
#' @examples
#' # Example usage:
#' data <- read.csv("data.csv")
#' scores <- get_nnd_score(data, "interactions")
#' print(scores)
#'
get_nnd_score <- function(clean, inter, windows.coord) {
  inter.distances <- clean %>% 
    select(inter) %>% 
    cbind(windows.coord) %>% 
    filter(!!as.symbol(inter) > 0) %>% 
    group_by(sampleID) %>%
    mutate(nnd = nndist(x1, y1))
  
  # Skip infinite values (single window in FOV)
  inter.distances <- inter.distances %>%
    mutate(nnd = if_else(is.infinite(nnd), NA, nnd))
  
  # Return fraction of windows having the interaction
  # and localization score as average nn-distance
  res <- c(
    nrow(inter.distances)/nrow(clean), 
    mean(inter.distances$nnd, na.rm = TRUE)
          )
  return(res)
}

#' Perform Leiden Clustering on UMAP Representation
#'
#' This function performs Leiden clustering on a UMAP representation of data.
#'
#' @param representation.umap A list containing the UMAP representation of data.
#'   The list should contain the following components:
#'   - `indexes`: A data frame with row indexes representing data points and column indexes representing the nearest neighbors.
#'   - `distances`: A matrix of distances between data points and their nearest neighbors.
#' @param resolution The resolution parameter for Leiden clustering. Higher values lead to more clusters.
#'   Default is 0.5.
#' 
#' @return A numeric vector containing the cluster membership for each data point.
#'
#' @examples
#' # Example usage:
#' library(umap)
#' umap_data <- umap(raw_data)
#' cluster_memberships <- leiden_onumap(umap_data$knn, resolution = 0.6)
#' print(cluster_memberships)
#'
leiden_onumap <- function(representation.umap, resolution = 0.5) {
  sim <- Matrix(nrow = nrow(representation.umap$indexes), 
                ncol = nrow(representation.umap$indexes), 
                data = 0, sparse = TRUE)
  
  # We need to convert distances to similarity
  maxval = max(representation.umap$distances)
  for (r in row.names(representation.umap$indexes)){
    knn = representation.umap$indexes[r,]
    # We make use of the fact that each point is most similar to itself
    sim[knn[1], knn] <- maxval - representation.umap$distances[r,]
    sim[knn, knn[1]] <- maxval - representation.umap$distances[r,]
  }
  # Remove self-edges
  diag(sim) = 0
  
  # Range matters for Leiden's resolution
  sim <- sim / max(sim)
  
  with_seed(
    1,
    groups <-
      graph.adjacency(sim, mode = "undirected", weighted = TRUE) %>%
      cluster_leiden(resolution_parameter = resolution, n_iterations = -1)
  )
  
  return(groups$membership)
}

sm_interactions <- function(misty.results, resolution, 
                            cutoff = 0, trim = 1,
                            save_heatmap = FALSE) {
  # trimming matters
  sig <- extract_signature(misty.results, type = "i", 
                           intersect.targets = FALSE, trim = trim)
  sig[is.na(sig)] <- floor(min(sig %>% select(-sample), na.rm = TRUE))
  sig <- sig %>% mutate(across(!sample, ~ ifelse(.x <= cutoff, 0, .x)))
  
  keep <- which(sig %>% select(-sample, -contains("intra_")) %>% rowSums() != 0)
  
  # the filtering here also matters
  clean <- sig %>%
    select(-sample, -contains("_.novar")) %>%
    slice(keep) %>%
    select(where(~ sum(.) != 0))
    
  windows.coord <- sig %>%
    select(sample) %>%
    slice(keep) %>%
    mutate(sampleID = sapply(sample, 
                             function(x) nth(str_split(x, "/")[[1]], -2)),
           coords = str_extract(sample, 
                                "[\\d|\\.]*_[\\d|\\.]*_[\\d|\\.]*_[\\d|\\.]*$")
           ) %>%
    separate(coords, c("x1", "y1", "x2", "y2"), sep = "_")   
  
  
  interaction.scores <- t(sapply(colnames(clean), 
                                 function(inter) get_nnd_score(clean, 
                                                               inter, 
                                                               windows.coord)))
  colnames(interaction.scores) <- c("Frequency", "NND")
  # Frequency could also be computed here as colSums(clean > 0) / nrow(clean)
  interactions.stats <- data.frame(Interaction = names(clean))
  
  # Defining clusters in hard because the data is sparse
  fit <- clean %>% 
    t %>%
    scale %>%
    umap(n_components = 10, random_state = 3895)

  # For visualization only: do not cluster on UMAP embedding
  if (is_character(save_heatmap)) { 
    heatmaply(fit$layout, file = save_heatmap)
    dev.off()
  }
  
  interactions.cluster <- leiden_onumap(fit$knn, resolution)

  interactions.stats <- cbind(interactions.stats, 
                              interaction.scores, 
                              interactions.cluster)
  colnames(interactions.stats)[ncol(interactions.stats)] <- "Cluster"
  interactions.stats$Cluster <- as.factor(interactions.stats$Cluster)
  
  return(interactions.stats)
}

interactions.stats <- sm_interactions(sm.results, 0.05, 
                            save_heatmap = "DCISct_interaction_embedding.html")
gp <- ggplot(interactions.stats[interactions.stats$Frequency > 0.1, ],
             aes(x = Frequency, y = NND, color = Cluster)) +
  geom_point()
gp
