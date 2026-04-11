#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-04-10
##### Code definition: Plot differential abundance heatmap results

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(ggh4x)
  library(yaml)
  library(grid)
})

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
fig_cfg <- config$figures$da_analysis

###Step2.Define input and output paths
dir_in <- file.path(project_root, fig_cfg$dir_in_rel$heatmap_data)
dir_out <- file.path(project_root, fig_cfg$dir_out_rel)
file_tax <- file.path(project_root, fig_cfg$dir_in_rel$taxonomy)

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

###Step3.Load custom plotting theme
path_theme <- file.path(project_root, "code/src/figure/plot_theme.R")
source(path_theme)

###Step4.Read heatmap datasets
heatmap_quercetin <- read_csv(
  file.path(dir_in, "heatmap_data_quercetin.csv"),
  show_col_types = FALSE
)

heatmap_cga <- read_csv(
  file.path(dir_in, "heatmap_data_cga.csv"),
  show_col_types = FALSE
)

###Step5.Read taxonomy table
taxonomy_df <- read_tsv(
  file_tax,
  show_col_types = FALSE
) %>%
  distinct(genus, .keep_all = TRUE)

###Step6.Define common factor levels
medium_levels <- c("GAM", "mBHI", "Schaedler")
agitation_levels <- c("Static", "Stirred")
time_levels <- c("0", "24", "48")
phylum_levels <- c(
  "Actinobacteria",
  "Bacteroidetes",
  "Firmicutes",
  "Planctomycetes",
  "Proteobacteria",
  "Synergistetes",
  "Tenericutes",
  "Verrucomicrobia"
)

###Step7.Define heatmap preparation function
prepare_heatmap_df <- function(df, taxonomy_df) {

  ###Step7-1.Join taxonomy information
  df_joined <- df %>%
    left_join(taxonomy_df, by = "genus")

  ###Step7-2.Fill missing phylum labels
  df_joined <- df_joined %>%
    mutate(
      phylum = ifelse(is.na(phylum), "Unclassified", phylum)
    )

  ###Step7-3.Set factor levels
  df_joined <- df_joined %>%
    mutate(
      medium = factor(as.character(medium), levels = medium_levels),
      agitation = factor(as.character(agitation), levels = agitation_levels),
      time = factor(as.character(time), levels = time_levels),
      phylum = factor(
        as.character(phylum),
        levels = c(phylum_levels, "Unclassified")
      )
    )

  ###Step7-4.Store genus-to-phylum mapping
  genus_phylum_df <- df_joined %>%
    distinct(genus, family, order, class, phylum)

  ###Step7-5.Define genus order
  genus_levels <- df_joined %>%
    arrange(phylum, genus) %>%
    pull(genus) %>%
    unique()

  ###Step7-6.Complete all medium/agitation/time combinations
  df_complete <- df_joined %>%
    select(
      genus,
      medium,
      agitation,
      time,
      stat,
      padj,
      log2FoldChange,
      family,
      order,
      class,
      phylum
    ) %>%
    complete(
      genus,
      medium = factor(medium_levels, levels = medium_levels),
      agitation = factor(agitation_levels, levels = agitation_levels),
      time = factor(time_levels, levels = time_levels)
    ) %>%
    left_join(genus_phylum_df, by = "genus", suffix = c("", "_joined")) %>%
    mutate(
      family = coalesce(family, family_joined),
      order = coalesce(order, order_joined),
      class = coalesce(class, class_joined),
      phylum = coalesce(phylum, phylum_joined)
    ) %>%
    select(
      genus,
      medium,
      agitation,
      time,
      stat,
      padj,
      log2FoldChange,
      family,
      order,
      class,
      phylum
    )

  ###Step7-7.Restore factor levels after completion
  df_complete <- df_complete %>%
    mutate(
      medium = factor(as.character(medium), levels = medium_levels),
      agitation = factor(as.character(agitation), levels = agitation_levels),
      time = factor(as.character(time), levels = time_levels),
      phylum = factor(
        as.character(phylum),
        levels = c(phylum_levels, "Unclassified")
      ),
      genus = factor(as.character(genus), levels = genus_levels)
    )

  ###Step7-8.Recalculate significance flags
  df_complete <- df_complete %>%
    mutate(
      sig_padj0.01 = !is.na(padj) & padj < 0.01,
      sig_padj0.05 = !is.na(padj) & padj < 0.05,
      panel_col = interaction(medium, agitation, sep = "_", lex.order = TRUE)
    )

  return(df_complete)
}

###Step8.Prepare heatmap datasets
heatmap_quercetin_plot <- prepare_heatmap_df(
  df = heatmap_quercetin,
  taxonomy_df = taxonomy_df
)

heatmap_cga_plot <- prepare_heatmap_df(
  df = heatmap_cga,
  taxonomy_df = taxonomy_df
)

###Step9.Define symmetric color scale range
max_abs_stat <- max(
  abs(c(heatmap_quercetin_plot$stat, heatmap_cga_plot$stat)),
  na.rm = TRUE
)

###Step10.Define heatmap plotting function
plot_da_heatmap <- function(df, max_abs_stat) {
  ggplot(
    df,
    aes(x = time, y = genus, fill = stat)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_point(
      data = df %>% filter(sig_padj0.01),
      aes(x = time, y = genus),
      shape = 16,
      size = 1.5,
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_nested(
      phylum ~ medium + agitation,
      scales = "free_y",
      space = "free_y",
      nest_line = element_line(linetype = 1, color = "black")
    ) +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(-max_abs_stat, max_abs_stat),
      na.value = "grey95",
      name = "Wald statistic"
    ) +
    labs(
      x = "Time (hrs)",
      y = "Genus"
    ) +
    theme_16S() +
    theme(
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size = 12, face = "italic"),
      axis.title.x = element_text(size = 15, face = "plain", margin = margin(t = 15)),
      axis.title.y = element_text(size = 15, face = "plain", margin = margin(r = 20)),
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.6),
      panel.spacing.x = unit(0.2, "lines"),
      panel.spacing.y = unit(0.45, "lines"),
      plot.title = element_text(size = 15),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
}

###Step11.Generate heatmap plots
p_quercetin <- plot_da_heatmap(
  df = heatmap_quercetin_plot,
  max_abs_stat = max_abs_stat
)

p_cga <- plot_da_heatmap(
  df = heatmap_cga_plot,
  max_abs_stat = max_abs_stat
)

###Step12.Save heatmap plots
ggsave(
  filename = file.path(dir_out, "da_heatmap_quercetin_genus.pdf"),
  plot = p_quercetin,
  width = 13,
  height = 16,
  units = "in"
)

ggsave(
  filename = file.path(dir_out, "da_heatmap_cga_genus.pdf"),
  plot = p_cga,
  width = 13,
  height = 16,
  units = "in"
)

###Step13.Print completion message
message("DA heatmap plotting with phylum grouping completed.")
message("Quercetin heatmap rows: ", nlevels(heatmap_quercetin_plot$genus))
message("CGA heatmap rows: ", nlevels(heatmap_cga_plot$genus))
message(
  "Quercetin phyla: ",
  paste(unique(na.omit(as.character(heatmap_quercetin_plot$phylum))), collapse = ", ")
)
message(
  "CGA phyla: ",
  paste(unique(na.omit(as.character(heatmap_cga_plot$phylum))), collapse = ", ")
)
