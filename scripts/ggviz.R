#!/usr/bin/Rscript

# Instructions:
# Prepare input files from a GFA file using the following commands:
#
# Set the graph name (example: chr6.C4)
# name=chr6.C4
#
# Extract path coverage information:
# odgi paths -i $name.gfa -H | cut -f 1,4- | gzip > $name.coverage.tsv.gz
#
# Extract node lengths:
# zgrep '^S' $name.gfa | awk '{print("node."$2,length($3))}' OFS="\t" > $name.node-lengths.tsv
#
# Usage:
# Rscript ggviz.R <coverage_file> <node_length_file> <output_file> [plot_mode] [cluster_file]
#
# Parameters:
#   coverage_file: Path coverage data (TSV with header, columns: path.name, node.1, node.2, ...)
#   node_length_file: Node lengths (TSV, no header, columns: node_id, length)
#   output_file: Output visualization file (PDF/PNG/JPEG)
#   plot_mode: "random" (hash-based colors) or "coverage" (coverage-based colors), default: "random"
#   cluster_file: Optional clustering file (see format below)
#
# Coverage file format:
# - Tab-separated with header row
# - Column 1: path.name (e.g., "chm13#chr6:31825251-31908851")
# - Remaining columns: node.N where N is node ID, values are coverage (0, 1, 2, ...)
#
# Node length file format:
# - Tab-separated, no header
# - Column 1: node ID (e.g., "node.1", "node.2")
# - Column 2: node length in base pairs
#
# Cluster file format (optional):
# - Tab-separated file with NO header
# - Column 1: path name (must match path.name in coverage file)
# - Column 2: cluster identifier (e.g., "cluster1", "HaploGroup1", etc.)
# - Example:
#   chm13#chr6:31825251-31908851	cluster1
#   grch38#chr6:31972046-32055647	cluster1
#   HG00438#2#JAHBCA010000042.1:24398231-24449090	cluster2
# - When provided, paths will be grouped by cluster in faceted plot panels

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)
library(data.table)
library(digest)

setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE)

# Required parameters
coverage_file <- args[1]
node_length_file <- args[2]
output_file <- args[3]

# Optional parameters
plot_mode <- ifelse(length(args) >= 4 && !is.na(args[4]) && args[4] != "", args[4], "random")
cluster_file <- ifelse(length(args) >= 5 && !is.na(args[5]) && args[5] != "", args[5], NA)

#mod nodes
pad_node_ids <- function(ids) {
  node_nums <- as.integer(sub("node\\.", "", ids))
  sprintf("node.%06d", node_nums)
}

#hash-based color generation (similar to odgi viz)
get_path_color <- function(path_name) {
  hash <- digest(path_name, algo="sha256", serialize=FALSE, raw=TRUE)

  # Extract specific bytes like odgi does (indices 24, 8, 16)
  r <- as.numeric(hash[25]) / 255  # index 24 (0-based) = position 25 (1-based)
  g <- as.numeric(hash[9]) / 255   # index 8 = position 9
  b <- as.numeric(hash[17]) / 255  # index 16 = position 17

  # Normalize
  sum_rgb <- r + g + b
  r <- r / sum_rgb
  g <- g / sum_rgb
  b <- b / sum_rgb

  # Brighten the color
  f <- min(1.5, 1.0 / max(r, g, b))
  r <- round(255 * min(r * f, 1.0))
  g <- round(255 * min(g * f, 1.0))
  b <- round(255 * min(b * f, 1.0))

  return(rgb(r, g, b, maxColorValue=255))
}

#read inputs
#pangenome coverage over nodes 
coverage_data <- read.table(coverage_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "?")
coverage_long <- coverage_data %>%
    pivot_longer(cols = -path.name, names_to = "node_id", values_to = "coverage") %>%
    rename(path_name = path.name)
coverage_long$node_id <- pad_node_ids(coverage_long$node_id)

#node lengths
node_length_data <- fread(node_length_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
length_df<-data.frame(
    node_id = node_length_data$V1,
    length = node_length_data$V2,
    stringsAsFactors = FALSE
  )
length_df$node_id<-pad_node_ids(length_df$node_id)

#clustering info (optional)
has_clusters <- !is.na(cluster_file) && cluster_file != "" && file.exists(cluster_file)

if (has_clusters) {
  clustering_data <- fread(cluster_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  clustering_df <- data.frame(
      path_name = clustering_data$V1,
      cluster = as.character(clustering_data$V2),
      stringsAsFactors = FALSE
  )
}

#merge
viz_data <- merge(coverage_long, length_df, by = "node_id")
if (has_clusters) {
  viz_data <- merge(viz_data, clustering_df, by = "path_name", all.x = TRUE)
}

viz_data <- viz_data %>%
  arrange(node_id, path_name) %>%
  group_by(path_name) %>%
  mutate(
    cumulative_length = cumsum(length),
    start_pos = lag(cumulative_length, default = 0),
    end_pos = cumulative_length
  ) %>%
  ungroup()

#mod viz data based on plot mode
if (plot_mode == "coverage") {
  spectral_colors <- brewer.pal(11, "Spectral")
  max_cov_palette <- 2 + length(spectral_colors) - 1

  color_map <- c(
    "0" = "white",
    "1" = "grey60"
  )

  for (i in 2:max_cov_palette) {
    color_map[as.character(i)] <- spectral_colors[i - 1]
  }

  viz_data$coverage_clamped <- as.character(
    ifelse(viz_data$coverage >= max_cov_palette, max_cov_palette, viz_data$coverage)
  )
} else {
  # Random/hash mode: assign colors based on path name
  path_colors <- sapply(unique(viz_data$path_name), get_path_color)
  names(path_colors) <- unique(viz_data$path_name)
}

if (has_clusters) {
  # Sort clusters naturally (alphanumeric order)
  cluster_levels <- sort(unique(viz_data$cluster))
  viz_data$cluster <- factor(viz_data$cluster, levels = cluster_levels)
}

if (plot_mode == "coverage") {
  p <- ggplot(viz_data, aes(x = start_pos, xend = end_pos,
                            y = path_name,
                            color = coverage_clamped)) +
    geom_segment(linewidth = 5) +
    scale_color_manual(
      values = color_map,
      na.value = "transparent"
    ) +
    scale_y_discrete(
      expand = expansion(add=c(1,1))
    ) +
    labs(
      x = "Pangenomic position",
      y = "Haplotype",
      color = "Coverage"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position="none"
    )
} else {
  # Random/hash mode
  p <- ggplot(viz_data, aes(x = start_pos, xend = end_pos,
                            y = path_name,
                            color = path_name)) +
    geom_segment(linewidth = 5) +
    scale_color_manual(
      values = path_colors,
      na.value = "transparent"
    ) +
    scale_y_discrete(
      expand = expansion(add=c(1,1))
    ) +
    labs(
      x = "Pangenomic position",
      y = "Haplotype"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position="none"
    )
}

if (has_clusters) {
  p <- p +
    facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
}

ggsave(output_file, width=24, height=max(6, 0.17*length(unique(viz_data$path_name))), limitsize=FALSE)
