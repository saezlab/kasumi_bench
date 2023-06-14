source("utils.R")
library(Matrix)
library(heatmaply)
library(umap)
library(spatstat.geom)

top.folder <- "DCISct"
sm.results <- read_rds(paste0(top.folder, ".rds"), "gz")

#' Calculates the normalized fraction of a network formed by the largest 
#' connected component (LCC) from a table of coordinates.
#'
#' @param coordtable A data frame or tibble containing the coordinates of 
#'                   points. The table must have columns named 'x1' and 'y1'.
#' @param stride The maximum distance threshold for considering two points
#'               as neighbors. Typically the stride used between windows.
#'
#' @return A numeric value representing the normalized fraction of the largest 
#'         connected component (LCC) in the network.
#'
#' @description
#' This function calculates the normalized fraction of a network formed by the 
#' largest connected component (LCC) from the given coordinate table. 
#' It constructs a network of neighboring points based on a distance threshold 
#' ('stride'). The largest connected component is then identified in the 
#' network, and the size of the LCC is determined. The normalized fraction is 
#' computed as the ratio of the number of points in the LCC to the total number
#' of points in the coordinate table, resulting in a value between 0 and 1.
#'
#' @examples
#' # Create a coordinate table
#' coords <- data.frame(x1 = c(1, 2, 3), y1 = c(1, 2, 3))
#'
#' # Calculate the normalized fraction of the largest connected component
#' fraction <- lcc_density_from_coordinates(coords$x1, coords$x2 stride = 1.5)
#'
#' @export
lcc_density_from_coordinates <- function(x, y, stride){
  distmat <- dist(cbind(x,y))
  
  neighbors <- as.matrix(distmat) <= stride
  diag(neighbors) <- FALSE # Self-edges
  
  # Note: we could avoid computing the full distance matrix and use a k-d tree
  # to find neighbors but existing implementations have a lot of dependencies
  
  # Make network of neighboring tiles
  ngraph <- graph.adjacency(neighbors, mode = "undirected")
  lccsize <- max(clusters(ngraph)$csize)
  
  # We want a null score if no 2 tiles are connected
  # and a perfect score if they're all connected
  return((lccsize - 1)/ (length(ngraph) - 1))
}

#' Calculates the Largest Connected Component (LCC) score based for a given
#' interaction, and given coordinates where it occurs.
#'
#' @param clean A data frame containing cleaned data.
#' @param inter The column name for the interaction variable in the 'clean'
#'              data frame.
#' @param windows.coord A data frame containing coordinates for windows.
#' @param stride The maximum distance threshold for considering two windows as 
#'               neighbors (default = 100).
#'
#' @return A numeric vector with two elements:
#'         - Fraction of windows with the interaction.
#'         - Average Local Clustering Coefficient (LCC) density of windows
#'           with the interaction.
#'
#' @description
#' This function calculates the LCC score by combining the interaction data 
#' and window coordinates. It first selects the relevant rows from 'clean' 
#' based on the interaction variable. Then, it combines the selected rows with 
#' the window coordinates using a left-join operation. Next, it calculates the
#' LCC density for each group of windows based on coordinates, using the 
#' specified 'stride' parameter. Finally, it returns a numeric vector with the 
#' fraction of windows having the interaction and the average LCC density of
#' those windows.
#'
#' @examples
#' # Create a sample 'clean' data frame
#' clean <- data.frame(
#'   sampleID = c(1, 2, 2, 4, 5),
#'   interaction = c(0, 1, 1, 0, 1)
#' )
#'
#' # Create a sample 'windows.coord' data frame
#' windows.coord <- data.frame(
#'   sampleID = c(1, 2, 2, 4, 5),
#'   x1 = c(10, 20, 30, 40, 50),
#'   y1 = c(100, 200, 300, 400, 500)
#' )
#'
#' # Calculate the LCC score
#' score <- get_lcc_score(clean, "interaction", windows.coord, stride = 100)
#'
#' @export
get_lcc_score <- function(clean, inter, windows.coord, stride = 100) {
  inter.windows <- clean %>% 
    select(inter) %>% 
    cbind(windows.coord) %>% 
    filter(!!as.symbol(inter) > 0) 
  
  inter.clustering <- inter.windows %>% 
    mutate(x1 = as.integer(x1), y1 = as.integer(y1)) %>%
    group_by(sampleID) %>%
    summarize(lcc = lcc_density_from_coordinates(x1, y1, stride))
  
  # Return fraction of windows having the interaction
  # and localization score as average LCC density
  res <- c(
    nrow(inter.windows)/nrow(clean), 
    mean(inter.clustering$lcc, na.rm = TRUE)
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
  diag(sim) <- 0
  
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
                                 function(inter) get_lcc_score(clean, 
                                                               inter, 
                                                               windows.coord)))
  colnames(interaction.scores) <- c("Frequency", "LCC")
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
             aes(x = Frequency, y = LCC, color = Cluster)) +
  geom_point()
gp
