top.folder = "DCISct"
sm.results = read_rds(paste0(top.folder, ".rds"), "gz")

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
#'   - The interaction score represents the fraction of windows having the interaction.
#'   - The localization score represents the average nearest neighbor distance (nn-distance) for windows with interaction.
#'
#' @examples
#' # Example usage:
#' data <- read.csv("data.csv")
#' scores <- get_nnd_score(data, "interactions")
#' print(scores)
#'
get_nnd_score <- function(clean, inter, windows.coord){
  inter.distances <- clean %>% select(inter) %>% 
    cbind(windows.coord) %>% 
    filter(!!as.symbol(inter) > 0) %>% 
    group_by(sampleID) %>%
    mutate(nnd = nndist(x1, y1))
  
  # Skip infinite values (single window in FOV)
  inter.distances <- inter.distances %>%
    mutate(nnd = if_else(is.infinite(nnd), NA, nnd))
  
  # Return fraction of windows having the interaction
  # and localization score as average nn-distance
  res = c(
    nrow(inter.distances)/nrow(clean), 
    mean(inter.distances$nnd, na.rm = TRUE)
          )
  return(res)
}

sm_interactions <- function(misty.results, cuts, res, cutoff = 0, trim = 1) {
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
    
  windows.coord <- sig %>%
    select(sample) %>%
    slice(keep) %>%
    mutate(sampleID = sapply(sample, function(x) nth(str_split(x, "/")[[1]], -2)),
           coords = str_extract(sample, "[\\d|\\.]*_[\\d|\\.]*_[\\d|\\.]*_[\\d|\\.]*$")) %>%
    separate(coords, c("x1", "y1", "x2", "y2"), sep = "_")   
  
  
  interaction.scores <- t(sapply(colnames(clean), 
                                 function(inter) get_nnd_score(clean, inter, windows.coord)))
  colnames(interaction.scores) <- c("Frequency", "NND")
  # Frequency could also be computed here as colSums(clean > 0) / nrow(clean)
  interactions.stats <- data.frame(Interaction = names(clean))
  
  # Defining clusters in hard because the data is sparse
  
  d <- dist(as_tibble(t(clean))) # euclidean distances between the rows
  fit <- cmdscale(d,eig=TRUE, k=10) # k is the number of dim
  heatmaply(fit$points)
  
  fit = umap(t(clean), n_components = 10)
  heatmaply(fit$layout)
  interactions.cluster = leiden_onsim(as_tibble(fit$layout), 0.1, 0.3)
  # Not deterministic without fixing a seed for umap
  
  interactions.stats <- cbind(interactions.stats, interaction.scores, interactions.cluster)
  colnames(interactions.stats)[ncol(interactions.stats)] <- "Cluster"
  interactions.stats$Cluster = as.factor(interactions.stats$Cluster)
  
  gp <- ggplot(interactions.stats[interactions.stats$Frequency > 0.1, ],
               aes(x = Frequency, y = NND, color = Cluster)) +
    geom_point()
  
  print(gp)
}

sm_interactions(sm.results, 0.1, 0.5)