source("utils.R")
library(extraDistr)
library(ggrepel)
library(spdep)

# Shared functions

moran_score <- function(repr, clust, stride = 100, pval = FALSE){
    sapply(repr, function(X) moran_score_sample(X, clust, stride, pval))
    }

moran_score_sample <- function(sample_rep, clust, stride, pval = FALSE) {
    clust_bool <- sample_rep[[clust]]
    if (sum(clust_bool) <= 1){
        return(NA)
    }

    x <- sample_rep$x
    y <- sample_rep$y
    distmat <- dist(cbind(x,y))
    neighbors <- as.matrix(distmat) <= stride
    diag(neighbors) <- FALSE # Self-edges
    lw <- mat2listw(neighbors, style="W")
    
    if (pval) {
        I <- moran.test(as.numeric(clust_bool), lw)
        return(c(I$estimate[1], I$estimate[2], I$p.value))
    } else {
        I <- moran(clust_bool, lw, length(lw$n), Szero(lw))[1]
        return(unlist(I))        
    }
}

# CTCL

# Load and format Kasumi results

misty.results <- read_rds("CTCLct400.rds")
sm.repr <- sm_labels(misty.results, 0.4, 0.9)
sm.repr.ext <- sm_labels(misty.results, 0.4, 0.9, freq = FALSE)

freq_all_clusters <- sapply(sm.repr.ext, 
      function(x) x %>%
           select(-id, -x, -y) %>%
           colSums)

# Compile cluster stats

cluster_stats <- data.frame(samples = rowSums(freq_all_clusters > 0))

cluster_stats$names <- sm.repr.ext[[1]] %>%
           select(-id, -x, -y) %>%
           names

cluster_stats$windows <- rowSums(freq_all_clusters)
cluster_stats$selected <- cluster_stats$names %in% names(sm.repr)

# Number of samples, pop size per class, number of draws
expected_samples_for_x_windows <- function(x, repet = 50000){
    X = rmvhyper(repet, sapply(sm.repr.ext, nrow), x)
    X = rowSums(X > 0)
    return(quantile(x = X, probs = c(0.1,0.5,0.9)))
}

with_seed(
  42,
    expected_samples_per_windows <- sapply(min(cluster_stats$windows):max(cluster_stats$windows), 
            expected_samples_for_x_windows)
)

expected_samples_per_windows = data.frame(t(expected_samples_per_windows))
expected_samples_per_windows$windows = min(cluster_stats$windows):max(cluster_stats$windows)
expected_samples_per_windows$selected = "Expectation"

lymph <- read_csv("data/LymphomaCODEX/single_cells.csv")

spots <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  pull(Spots) %>%
  unique()

outcome <- lymph %>%
  filter(Groups %in% c(1, 2)) %>%
  group_by(Spots, Patients) %>%
  summarize(Groups = Groups[1], .groups = "drop")

freq.expr <- sm.repr %>%
  left_join(outcome %>%
              mutate(Spots = as.character(Spots)), by = c("id" = "Spots")) %>%
  rename(target = Groups) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id, -Patients)

reliance <- model_reliance(freq.expr)
ggsave("mr_docker.pdf", last_plot())

cluster_stats$names <- str_remove(string = cluster_stats$names, pattern = "...")
cluster_stats <- cluster_stats %>%
    left_join(reliance, by = join_by(names == Cluster))

# Label clusters only if they have a significant sMR
cluster_stats$label <- ifelse(abs(cluster_stats$sMR) > 1, cluster_stats$names, NA)
# and if they have an unusual number of samples compared to their number of windows
lowsamples <- (cluster_stats %>% 
    left_join(expected_samples_per_windows, by = "windows") %>%
    mutate(lowsamples = X10. > samples))$lowsamples
cluster_stats$label <- ifelse(lowsamples, cluster_stats$label, NA)

# Plot and save
gp <- ggplot(cluster_stats, aes(y = samples, x = windows, color = sMR)) +
    geom_point() +
    geom_hline(yintercept = 5, lty = 2, color = "gray50") +
    theme_classic() + 
    geom_line(data = expected_samples_per_windows,
               aes(y = X50., x = windows), color = "gray50") + 
     geom_line(data = expected_samples_per_windows,
               aes(y = X90., x = windows), color = "gray75", lty = 2) + 
     geom_line(data = expected_samples_per_windows,
               aes(y = X10., x = windows), color = "gray75", lty = 2) +
     geom_text_repel(aes(label = label)) +
     ggtitle("Frequency of window clusters") +
     xlab("Number of windows") + 
     ylab("Number of samples") + 
    scale_color_steps2(low = "#008837", mid = "white", high = "#7b3294", na.value="gray75")
ggsave("CTCLct400_cluster_frequency.pdf", gp)
ggsave("CTCLct400_cluster_frequency.png", gp)

# For all the selected clusters
# (Non-randomly spread across samples and useful for classification)
contig_scores <- cluster_stats$label %>%
    na.omit %>%
    sapply(function(x) na.omit(moran_score(sm.repr.ext, paste0("...", x), stride = 200, pval = TRUE)), simplify = FALSE) 

contig_clusters <- lapply(names(contig_scores), function(name) {
  rep(name, length(na.omit(unlist(contig_scores[[name]])))/3)
})

contiguity_df <- data.frame(clusters = rep(unlist(contig_clusters), 2))
unlist_scores <- na.omit(unlist(contig_scores))
contiguity_df$scores <- c(unlist_scores[3*1:(length(unlist_scores)/3) - 1], unlist_scores[3*1:(length(unlist_scores)/3) - 2])
contiguity_df$type <- rep(c("Expected", "Observed"), each = length(unlist(contig_clusters)))
contiguity_pval <- unlist_scores[3*1:(length(unlist_scores)/3)]

pval_df <- data.frame(pval = contiguity_pval, clusters = unlist(contig_clusters))
pval_df <- pval_df %>% group_by(clusters) %>%
    mutate(cor_p = p.adjust(pval, "bonferroni")) %>%
    summarize(mincorp = min(cor_p)) %>%
    mutate(mincorp = p.adjust(mincorp, "BH"))
pval_df$text <- paste("p =", signif(pval_df$mincorp, 3))

# Define the color mapping and the plot scales
color_mapping <- c("Expected" = "#CCCCCC", "Observed" = "#8FA1CC")
xmin <- min(contiguity_df$scores)
xmax <- max(contiguity_df$scores) + 0.2

gp <- ggplot(contiguity_df, aes(x = clusters, y = scores)) +
    geom_violin(aes(fill = type), draw_quantiles = c(0.25, 0.5, 0.75), scale = "width") +
    scale_fill_manual(values = color_mapping) +
    theme_classic() +
    xlab("Clusters") + ylab("Moran's I") +
    coord_flip() +
    ylim(xmin,xmax) +
    geom_text(aes(label = text), data = pval_df, size = 3, y = xmax - 0.1) +
    guides(fill = guide_legend("I values", reverse = TRUE))

ggsave("CTCLct400_cluster_moran_signif.pdf", gp)
ggsave("CTCLct400_cluster_moran_signif.png", gp)


# IMC

misty.results <- read_rds("BCexpr200.rds")

bmeta <- read_csv("data/Basel_PatientMetadata.csv")
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

sm.repr <- sm_labels(misty.results, 0.3, 0.8)
sm.repr.ext <- sm_labels(misty.results, 0.3, 0.8, freq = FALSE)

freq.sm <- sm.repr %>%
  left_join(resp, by = c("id" = "core")) %>%
  rename(target = response) %>%
  mutate(target = as.factor(make.names(target))) %>%
  select(-id)

freq_all_clusters <- sapply(sm.repr.ext, 
      function(x) x %>%
           select(-id, -x, -y) %>%
           colSums)

reliance <- model_reliance(freq.sm)
ggsave("mr_docker2.pdf", last_plot())

# Compile cluster stats

cluster_stats <- data.frame(samples = rowSums(freq_all_clusters > 0))

cluster_stats$names <- sm.repr.ext[[1]] %>%
           select(-id, -x, -y) %>%
           names

cluster_stats$windows <- rowSums(freq_all_clusters)
cluster_stats$selected <- cluster_stats$names %in% names(sm.repr)

# Number of samples, pop size per class, number of draws
expected_samples_for_x_windows <- function(x, repet = 50000){
    X = rmvhyper(repet, sapply(sm.repr.ext, nrow), x)
    X = rowSums(X > 0)
    return(quantile(x = X, probs = c(0.1,0.5,0.9)))
}

with_seed(
  42,
    expected_samples_per_windows <- sapply(min(cluster_stats$windows):max(cluster_stats$windows), 
            expected_samples_for_x_windows)
)

expected_samples_per_windows = data.frame(t(expected_samples_per_windows))
expected_samples_per_windows$windows = min(cluster_stats$windows):max(cluster_stats$windows)
expected_samples_per_windows$selected = "Expectation"

cluster_stats$names <- str_remove(string = cluster_stats$names, pattern = "...")
cluster_stats <- cluster_stats %>%
    left_join(reliance, by = join_by(names == Cluster))

# Label clusters only if they have a significant sMR
cluster_stats$label <- ifelse(abs(cluster_stats$sMR) > 1, cluster_stats$names, NA)
# and if they have an unusual number of samples compared to their number of windows
lowsamples <- (cluster_stats %>% 
    left_join(expected_samples_per_windows, by = "windows") %>%
    mutate(lowsamples = X10. > samples))$lowsamples
cluster_stats$label <- ifelse(lowsamples, cluster_stats$label, NA)

# Plot and save
gp <- ggplot(cluster_stats, aes(y = samples, x = windows, color = sMR)) +
    geom_point() +
    geom_hline(yintercept = 5, lty = 2, color = "gray50") +
    theme_classic() + 
    geom_line(data = expected_samples_per_windows,
               aes(y = X50., x = windows), color = "gray50") + 
     geom_line(data = expected_samples_per_windows,
               aes(y = X90., x = windows), color = "gray75", lty = 2) + 
     geom_line(data = expected_samples_per_windows,
               aes(y = X10., x = windows), color = "gray75", lty = 2) +
     geom_text_repel(aes(label = label)) +
     ggtitle("Frequency of window clusters") +
     xlab("Number of windows") + 
     ylab("Number of samples") + 
    scale_color_steps2(low = "#008837", mid = "white", high = "#7b3294", na.value="gray75")
ggsave("BCexpr200_cluster_frequency.pdf", gp)
ggsave("BCexpr200_cluster_frequency.png", gp)

# For all the selected clusters
# (Non-randomly spread across samples and useful for classification)
contig_scores <- cluster_stats$label %>%
    na.omit %>%
    sapply(function(x) na.omit(moran_score(sm.repr.ext, paste0("...", x), stride = 100, pval = TRUE)), simplify = FALSE) 

contig_clusters <- lapply(names(contig_scores), function(name) {
  rep(name, length(na.omit(unlist(contig_scores[[name]])))/3)
})

contiguity_df <- data.frame(clusters = rep(unlist(contig_clusters), 2))
unlist_scores <- na.omit(unlist(contig_scores))
contiguity_df$scores <- c(unlist_scores[3*1:(length(unlist_scores)/3) - 1], unlist_scores[3*1:(length(unlist_scores)/3) - 2])
contiguity_df$type <- rep(c("Expected", "Observed"), each = length(unlist(contig_clusters)))
contiguity_pval <- unlist_scores[3*1:(length(unlist_scores)/3)]

pval_df <- data.frame(pval = contiguity_pval, clusters = unlist(contig_clusters))
pval_df <- pval_df %>% group_by(clusters) %>%
    mutate(cor_p = p.adjust(pval, "bonferroni")) %>%
    summarize(mincorp = min(cor_p)) %>%
    mutate(mincorp = p.adjust(mincorp, "BH"))
pval_df$text <- paste("p =", signif(pval_df$mincorp, 3))

# Define the color mapping and the plot scales
color_mapping <- c("Expected" = "#CCCCCC", "Observed" = "#8FA1CC")
xmin <- min(contiguity_df$scores)
xmax <- max(contiguity_df$scores) + 0.2

gp <- ggplot(contiguity_df, aes(x = clusters, y = scores)) +
    geom_violin(aes(fill = type), draw_quantiles = c(0.25, 0.5, 0.75), scale = "width") +
    scale_fill_manual(values = color_mapping) +
    theme_classic() +
    xlab("Clusters") + ylab("Moran's I") +
    coord_flip() +
    ylim(xmin,xmax) +
    geom_text(aes(label = text), data = pval_df, size = 3, y = xmax - 0.1) +
    guides(fill = guide_legend("I values", reverse = TRUE))

ggsave("BCexpr200_cluster_moran_signif.pdf", gp)
ggsave("BCexpr200_cluster_moran_signif.png", gp)