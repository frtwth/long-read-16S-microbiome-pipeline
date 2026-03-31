##### Code by S. Jang
##### Last updated: 2026-03-27
##### Code definition: Alpha diversity analysis

library(data.table)
library(dplyr)
library(tibble)
library(vegan)
library(yaml)

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
alpha_cfg <- config$alpha_diversity

path_counts <- file.path(project_root, alpha_cfg$dir_in_rel$counts)
path_metadata <- file.path(project_root, alpha_cfg$dir_in_rel$metadata)
dir_out <- file.path(project_root, alpha_cfg$dir_out_rel)

tax_col <- alpha_cfg$tax_col
sample_id_col <- alpha_cfg$sample_id_col

file_alpha_table <- alpha_cfg$output_files$alpha_table
file_alpha_rds <- alpha_cfg$output_files$alpha_rds
file_depth_table <- alpha_cfg$output_files$depth_table

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

###Step2.Load data
genus_cts <- fread(path_counts)
metadata.df <- readRDS(path_metadata)

###Step3.Preprocess genus count table
genus_cts <- genus_cts %>%
  filter(!is.na(.data[[tax_col]]), .data[[tax_col]] != "")

colnames(genus_cts) <- colnames(genus_cts) %>%
  gsub("_both$", "", .) %>%
  gsub("^BC_", "repBC", .) %>%
  gsub("_PL", "_Pl", .)

genus_cts[is.na(genus_cts)] <- 0
genus_cts <- as.data.frame(genus_cts)

###Step4.Prepare count matrix
count_df <- genus_cts
rownames(count_df) <- count_df[[tax_col]]
count_df[[tax_col]] <- NULL

count_mat <- as.matrix(count_df)

storage.mode(count_mat) <- "numeric"
count_mat <- t(count_mat)

###Step5.Calculate sample depth
depth_df <- data.frame(
  sample_id = rownames(count_mat),
  total_counts = rowSums(count_mat),
  stringsAsFactors = FALSE
)

###Step6.Calculate alpha diversity
shannon <- diversity(count_mat, index = "shannon")
simpson <- diversity(count_mat, index = "simpson")
richness <- rowSums(count_mat > 0)
evenness <- ifelse(richness > 1, shannon / log(richness), NA_real_)

alpha_df <- data.frame(
  sample_id = rownames(count_mat),
  Shannon = shannon,
  Simpson = simpson,
  Richness = richness,
  Evenness = evenness,
  stringsAsFactors = FALSE
)

###Step7.Prepare metadata
metadata <- as.data.frame(metadata.df)
rownames(metadata) <- metadata[[sample_id_col]]
metadata <- metadata[alpha_df$sample_id, , drop = FALSE]

stopifnot(all(alpha_df$sample_id == rownames(metadata)))

alpha_meta <- alpha_df %>%
  left_join(metadata, by = c("sample_id" = sample_id_col))

depth_meta <- depth_df %>%
  left_join(metadata, by = c("sample_id" = sample_id_col))

###Step8.Save outputs
write.csv(
  alpha_meta,
  file = file.path(dir_out, file_alpha_table),
  row.names = FALSE
)

write.csv(
  depth_meta,
  file = file.path(dir_out, file_depth_table),
  row.names = FALSE
)

saveRDS(
  alpha_meta,
  file = file.path(dir_out, file_alpha_rds)
)

message("Alpha diversity calculation completed.")
message("Saved outputs to: ", dir_out)
