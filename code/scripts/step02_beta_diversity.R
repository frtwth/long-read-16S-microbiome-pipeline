#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-27
##### Code definition: Beta diversity analysis

library(data.table)
library(dplyr)
library(vegan)
library(yaml)

set.seed(2026)

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
beta_cfg <- config$beta_diversity

path_abundance <- file.path(project_root, beta_cfg$dir_in_rel$abundance)
path_metadata <- file.path(project_root, beta_cfg$dir_in_rel$metadata)
dir_out <- file.path(project_root, beta_cfg$dir_out_rel)

sample_id_col <- beta_cfg$sample_id_col
tax_col <- beta_cfg$tax_col

min_abundance <- beta_cfg$abundance_filter$min_abundance
renormalize <- beta_cfg$abundance_filter$renormalize

distance_method <- beta_cfg$distance$method

permanova_formula <- as.formula(
  paste("bray_dist ~", beta_cfg$permanova$formula)
)
permanova_by <- beta_cfg$permanova$by
permanova_permutations <- beta_cfg$permanova$permutations

dispersion_group_vars <- beta_cfg$dispersion$group_vars
dispersion_permutations <- beta_cfg$dispersion$permutations

ordination_k <- beta_cfg$ordination$k
ordination_add <- beta_cfg$ordination$add

file_beta_rds <- beta_cfg$output_files$beta_rds

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

###Step2.Load data
genus_abundance <- fread(path_abundance)
metadata.df <- readRDS(path_metadata)

###Step3.Preprocess genus abundance table
genus_abundance <- genus_abundance %>%
  filter(!is.na(.data[[tax_col]]), .data[[tax_col]] != "")

colnames(genus_abundance) <- colnames(genus_abundance) %>%
  gsub("_both$", "", .) %>%
  gsub("^BC_", "repBC", .) %>%
  gsub("_PL", "_Pl", .)

genus_abundance[is.na(genus_abundance)] <- 0
genus_abundance <- as.data.frame(genus_abundance)

###Step4.Filter rare taxa
abundance_mat <- genus_abundance[, setdiff(colnames(genus_abundance), tax_col), drop = FALSE]
abundance_mat[abundance_mat < min_abundance] <- 0

genus_abundance[, setdiff(colnames(genus_abundance), tax_col)] <- abundance_mat
genus_abundance <- genus_abundance[
  rowSums(genus_abundance[, setdiff(colnames(genus_abundance), tax_col), drop = FALSE]) > 0,
]

###Step5.Re-normalize within each sample
if (isTRUE(renormalize)) {
  abundance_mat <- genus_abundance[, setdiff(colnames(genus_abundance), tax_col), drop = FALSE]
  sample_sums <- colSums(abundance_mat)

  stopifnot(all(sample_sums > 0))

  genus_abundance[, setdiff(colnames(genus_abundance), tax_col)] <- sweep(
    abundance_mat,
    2,
    sample_sums,
    "/"
  )
}

###Step6.Prepare abundance matrix for beta diversity
abundance_df <- genus_abundance
rownames(abundance_df) <- abundance_df[[tax_col]]
abundance_df[[tax_col]] <- NULL

abundance_mat <- as.matrix(abundance_df)
storage.mode(abundance_mat) <- "numeric"

###Step7.Prepare metadata
metadata <- as.data.frame(metadata.df)
rownames(metadata) <- metadata[[sample_id_col]]
metadata <- metadata[colnames(abundance_mat), , drop = FALSE]

metadata <- metadata %>%
  mutate(
    medium = factor(medium),
    time = factor(time),
    stirred = factor(stirred),
    treatment = factor(treatment)
  )

stopifnot(all(colnames(abundance_mat) == rownames(metadata)))

###Step8.Calculate beta-diversity distance
bray_dist <- vegdist(t(abundance_mat), method = distance_method)

###Step9.PERMANOVA
permanova_res <- adonis2(
  permanova_formula,
  data = metadata,
  by = permanova_by,
  permutations = permanova_permutations
)

print(permanova_res)

###Step10.Homogeneity of dispersion
dispersion_res <- list()

for (group_var in dispersion_group_vars) {
  bd <- betadisper(bray_dist, metadata[[group_var]])

  dispersion_res[[group_var]] <- permutest(
    bd,
    permutations = dispersion_permutations
  )
}

print(dispersion_res)

###Step11.PCoA coordinates
pcoa_res <- wcmdscale(
  bray_dist,
  k = ordination_k,
  eig = TRUE,
  add = ordination_add
)

pcoa_df <- data.frame(
  sample_id = rownames(pcoa_res$points),
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(
    metadata,
    by = c("sample_id" = sample_id_col)
  )

pcoa_var <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)

###Step12.Save outputs
saveRDS(
  list(
    abundance_mat = abundance_mat,
    metadata = metadata,
    bray_dist = bray_dist,
    permanova_res = permanova_res,
    dispersion_res = dispersion_res,
    pcoa_df = pcoa_df,
    pcoa_var = pcoa_var
  ),
  file = file.path(dir_out, file_beta_rds)
)

message("Beta diversity analysis completed.")
message("Saved beta diversity results to: ", file.path(dir_out, file_beta_rds))
