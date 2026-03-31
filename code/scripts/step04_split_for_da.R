#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-27
##### Code definition: Table split for differential abundance analysis

library(data.table)
library(dplyr)
library(readr)
library(yaml)

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
da_cfg <- config$da_split

path_counts <- file.path(project_root, da_cfg$dir_in_rel$counts)
path_metadata <- file.path(project_root, da_cfg$dir_in_rel$metadata)

dir_out_split <- file.path(project_root, da_cfg$dir_out_rel$split)
dir_out_summary <- file.path(project_root, da_cfg$dir_out_rel$summary)

sample_id_col <- da_cfg$sample_id_col
tax_col <- da_cfg$tax_col

grouping_vars <- da_cfg$grouping_vars
metadata_cols <- da_cfg$metadata_cols

file_counts <- da_cfg$output_files$counts
file_metadata <- da_cfg$output_files$metadata
file_manifest <- da_cfg$output_files$manifest
file_skipped <- da_cfg$output_files$skipped

dir.create(dir_out_split, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out_summary, recursive = TRUE, showWarnings = FALSE)

###Step2.Load data
genus_cts <- fread(path_counts)
metadata.df <- readRDS(path_metadata)

###Step3.Preprocess genus counts table
genus_cts <- genus_cts %>%
  filter(!is.na(.data[[tax_col]]), .data[[tax_col]] != "")

colnames(genus_cts) <- colnames(genus_cts) %>%
  gsub("_both$", "", .) %>%
  gsub("^BC_", "repBC", .) %>%
  gsub("_PL", "_Pl", .)

genus_cts[is.na(genus_cts)] <- 0
genus_cts <- as.data.frame(genus_cts)

###Step4.Check sample IDs
if (!all(metadata.df[[sample_id_col]] %in% colnames(genus_cts))) {
  missing_samples <- setdiff(metadata.df[[sample_id_col]], colnames(genus_cts))
  stop(
    "Some sample IDs in metadata are missing from genus count table:\n",
    paste(missing_samples, collapse = ", ")
  )
}

###Step5.Define DA comparison pairs
comparison_list <- da_cfg$comparison_pairs

###Step6.Get unique grouping combinations
comb_df <- metadata.df %>%
  distinct(across(all_of(grouping_vars))) %>%
  arrange(across(all_of(grouping_vars)))

###Step7.Initialize summary tables
manifest_list <- list()
skip_list <- list()

###Step8.Generate DA input folders
for (i in seq_len(nrow(comb_df))) {

  current_comb <- comb_df[i, , drop = FALSE]

  message(
    "Processing: ",
    paste(unlist(current_comb), collapse = " / ")
  )

  meta_base <- metadata.df

  for (group_var in grouping_vars) {
    meta_base <- meta_base %>%
      filter(.data[[group_var]] == current_comb[[group_var]])
  }

  for (comp_name in names(comparison_list)) {

    ctrl_name <- comparison_list[[comp_name]]$control
    case_name <- comparison_list[[comp_name]]$treatment

    meta_sub <- meta_base %>%
      filter(treatment %in% c(ctrl_name, case_name))

    ###Step8-1.Skip if either group is missing
    group_counts <- table(meta_sub$treatment)

    if (!(ctrl_name %in% names(group_counts)) || !(case_name %in% names(group_counts))) {
      comparison_id <- paste(
        unlist(current_comb),
        comp_name,
        sep = "_",
        collapse = "_"
      )

      skip_list[[length(skip_list) + 1]] <- data.frame(
        comparison_id = comparison_id,
        medium = current_comb$medium,
        time = current_comb$time,
        stirred = current_comb$stirred,
        comparison = comp_name,
        reason = "One or both treatment groups missing",
        stringsAsFactors = FALSE
      )
      next
    }

    ###Step8-2.Add DESeq2 condition column
    meta_sub <- meta_sub %>%
      mutate(
        condition = case_when(
          treatment == ctrl_name ~ "ctrl",
          treatment == case_name ~ "case"
        )
      )

    ###Step8-3.Reorder samples by condition
    meta_sub <- meta_sub %>%
      arrange(condition, replicate)

    sample_vec <- meta_sub[[sample_id_col]]

    ###Step8-4.Subset count table
    counts_sub <- genus_cts %>%
      select(all_of(c(tax_col, sample_vec)))

    ###Step8-5.Make output folder
    comparison_id <- paste(
      current_comb$medium,
      paste0("t", current_comb$time),
      current_comb$stirred,
      comp_name,
      sep = "_"
    )

    dir_comp <- file.path(dir_out_split, comparison_id)
    dir.create(dir_comp, recursive = TRUE, showWarnings = FALSE)

    ###Step8-6.Save metadata
    meta_save <- meta_sub %>%
      select(all_of(metadata_cols))

    write_tsv(
      meta_save,
      file.path(dir_comp, file_metadata)
    )

    ###Step8-7.Save counts
    write_tsv(
      counts_sub,
      file.path(dir_comp, file_counts)
    )

    ###Step8-8.Save manifest info
    manifest_list[[length(manifest_list) + 1]] <- data.frame(
      comparison_id = comparison_id,
      medium = current_comb$medium,
      time = current_comb$time,
      stirred = current_comb$stirred,
      comparison = comp_name,
      control = ctrl_name,
      treatment = case_name,
      n_ctrl = sum(meta_sub$condition == "ctrl"),
      n_case = sum(meta_sub$condition == "case"),
      n_total = nrow(meta_sub),
      stringsAsFactors = FALSE
    )
  }
}

###Step9.Save summary files
manifest_df <- bind_rows(manifest_list)

write_csv(
  manifest_df,
  file.path(dir_out_summary, file_manifest)
)

if (length(skip_list) > 0) {
  skip_df <- bind_rows(skip_list)
} else {
  skip_df <- data.frame(
    comparison_id = character(),
    medium = character(),
    time = character(),
    stirred = character(),
    comparison = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

write_csv(
  skip_df,
  file.path(dir_out_summary, file_skipped)
)

message("Table split for differential abundance analysis completed.")
message("Saved split inputs to: ", dir_out_split)
message("Saved summary files to: ", dir_out_summary)
