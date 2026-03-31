#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-27
##### Code definition: Differential abundance analysis

library(DESeq2)
library(dplyr)
library(readr)
library(yaml)

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
da_cfg <- config$da_analysis

dir_in <- file.path(project_root, da_cfg$dir_in_rel)
dir_out_tables <- file.path(project_root, da_cfg$dir_out_rel$tables)
dir_out_summary <- file.path(project_root, da_cfg$dir_out_rel$summary)

sample_id_col <- da_cfg$sample_id_col
tax_col <- da_cfg$tax_col
condition_col <- da_cfg$condition_col

min_total_count <- da_cfg$filtering$min_total_count
min_taxa_after_filter <- da_cfg$filtering$min_taxa_after_filter
round_counts <- da_cfg$filtering$round_counts

deseq2_design <- as.formula(da_cfg$deseq2$design)

file_counts <- da_cfg$output_files$counts
file_metadata <- da_cfg$output_files$metadata
file_deseq2_results <- da_cfg$output_files$deseq2_results
file_sig_padj_005 <- da_cfg$output_files$sig_padj_005
file_sig_padj_001 <- da_cfg$output_files$sig_padj_001
file_normalized_counts <- da_cfg$output_files$normalized_counts
file_size_factors <- da_cfg$output_files$size_factors
file_summary_table <- da_cfg$output_files$summary_table
file_skipped_table <- da_cfg$output_files$skipped_table

dir.create(dir_out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out_summary, recursive = TRUE, showWarnings = FALSE)

###Step2.List comparison folders
comp_dirs <- list.dirs(dir_in, recursive = FALSE, full.names = TRUE)

###Step3.Initialize summary lists
summary_list <- list()
skip_list <- list()

###Step4.Run DESeq2 for each comparison
for (dir_comp in comp_dirs) {

  comp_name <- basename(dir_comp)

  message("Running DESeq2 for: ", comp_name)

  ###Step4-1.Load data
  counts_path <- file.path(dir_comp, file_counts)
  meta_path <- file.path(dir_comp, file_metadata)

  counts_df <- read_tsv(counts_path, show_col_types = FALSE)
  meta_df <- read_tsv(meta_path, show_col_types = FALSE)

  ###Step4-2.Prepare count matrix
  counts_df <- as.data.frame(counts_df)
  rownames(counts_df) <- counts_df[[tax_col]]
  counts_df[[tax_col]] <- NULL

  counts_mat <- as.matrix(counts_df)
  storage.mode(counts_mat) <- "numeric"

  ###Step4-3.Round estimated counts if requested
  if (isTRUE(round_counts)) {
    counts_mat <- round(counts_mat)
  }

  ###Step4-4.Prepare metadata
  meta_df <- as.data.frame(meta_df)
  meta_df[[condition_col]] <- factor(meta_df[[condition_col]], levels = c("ctrl", "case"))
  rownames(meta_df) <- meta_df[[sample_id_col]]

  ###Step4-5.Match sample order
  meta_df <- meta_df[colnames(counts_mat), , drop = FALSE]
  stopifnot(all(colnames(counts_mat) == rownames(meta_df)))

  ###Step4-6.Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = meta_df,
    design = deseq2_design
  )

  ###Step4-7.Filter low-abundance taxa
  keep <- rowSums(counts(dds)) > min_total_count

  if (sum(keep) < min_taxa_after_filter) {
    message("Skipping ", comp_name, ": too few taxa after filtering.")

    skip_list[[length(skip_list) + 1]] <- data.frame(
      comparison = comp_name,
      reason = "Too few taxa after filtering",
      stringsAsFactors = FALSE
    )

    next
  }

  dds <- dds[keep, ]

  ###Step4-8.Run DESeq2
  dds <- DESeq(dds)
  res <- results(dds)

  ###Step4-9.Convert result table
  res_df <- as.data.frame(res)
  res_df[[tax_col]] <- rownames(res_df)
  res_df <- res_df %>%
    arrange(is.na(padj), padj)

  ###Step4-10.Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame()

  norm_counts[[tax_col]] <- rownames(norm_counts)

  ###Step4-11.Save size factors
  sf <- sizeFactors(dds)
  size_factors <- data.frame(
    sample = names(sf),
    size_factor = as.numeric(sf),
    stringsAsFactors = FALSE
  )

  ###Step4-12.Make output directory
  dir_comp_out <- file.path(dir_out_tables, comp_name)
  dir.create(dir_comp_out, recursive = TRUE, showWarnings = FALSE)

  ###Step4-13.Save full results
  write_csv(
    res_df,
    file.path(dir_comp_out, file_deseq2_results)
  )

  ###Step4-14.Save significant taxa (padj < 0.05)
  sig_df_005 <- res_df %>%
    filter(!is.na(padj), padj < 0.05)

  write_csv(
    sig_df_005,
    file.path(dir_comp_out, file_sig_padj_005)
  )

  ###Step4-15.Save significant taxa (padj < 0.01)
  sig_df_001 <- res_df %>%
    filter(!is.na(padj), padj < 0.01)

  write_csv(
    sig_df_001,
    file.path(dir_comp_out, file_sig_padj_001)
  )

  ###Step4-16.Save normalized counts
  write_tsv(
    norm_counts,
    file.path(dir_comp_out, file_normalized_counts)
  )

  ###Step4-17.Save size factors
  write_csv(
    size_factors,
    file.path(dir_comp_out, file_size_factors)
  )

  ###Step4-18.Store summary info
  summary_list[[comp_name]] <- data.frame(
    comparison = comp_name,
    normalization = "DESeq2_default",
    n_taxa_tested = nrow(res_df),
    sig_taxa_padj0.05 = sum(!is.na(res_df$padj) & res_df$padj < 0.05),
    sig_taxa_padj0.01 = sum(!is.na(res_df$padj) & res_df$padj < 0.01),
    up_padj0.05 = sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0),
    down_padj0.05 = sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < 0),
    min_size_factor = min(sf),
    median_size_factor = median(sf),
    max_size_factor = max(sf),
    stringsAsFactors = FALSE
  )
}

###Step5.Make summary tables
if (length(summary_list) > 0) {
  summary_table <- bind_rows(summary_list) %>%
    arrange(desc(sig_taxa_padj0.05))
} else {
  summary_table <- data.frame(
    comparison = character(),
    normalization = character(),
    n_taxa_tested = numeric(),
    sig_taxa_padj0.05 = numeric(),
    sig_taxa_padj0.01 = numeric(),
    up_padj0.05 = numeric(),
    down_padj0.05 = numeric(),
    min_size_factor = numeric(),
    median_size_factor = numeric(),
    max_size_factor = numeric(),
    stringsAsFactors = FALSE
  )
}

write_csv(
  summary_table,
  file.path(dir_out_summary, file_summary_table)
)

###Step6.Save skipped comparisons
if (length(skip_list) > 0) {
  skip_table <- bind_rows(skip_list)
} else {
  skip_table <- data.frame(
    comparison = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

write_csv(
  skip_table,
  file.path(dir_out_summary, file_skipped_table)
)

message("Differential abundance analysis completed.")
message("Saved DA result tables to: ", dir_out_tables)
message("Saved summary files to: ", dir_out_summary)
