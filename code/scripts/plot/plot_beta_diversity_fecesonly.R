#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-04-09
##### Code definition: Plot beta diversity results (FecesOnly subgroup)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(yaml)
})

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
fig_cfg <- config$figures$beta_diversity

###Step2.Load custom plotting theme
path_theme <- file.path(project_root, "code/src/figure/plot_theme.R")
source(path_theme)

###Step3.Define input and output paths
path_beta <- file.path(project_root, fig_cfg$dir_in_rel$beta_rds_fecesonly)

dir_fig <- file.path(project_root, fig_cfg$dir_out_rel$figures)
dir_tab <- file.path(project_root, fig_cfg$dir_out_rel$tables)

dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

###Step4.Load beta-diversity results
beta_res <- readRDS(path_beta)

pcoa_df <- beta_res$pcoa_df
pcoa_var <- beta_res$pcoa_var
permanova_res <- beta_res$permanova_res

###Step5.Convert variables to factors
pcoa_df <- pcoa_df %>%
  mutate(
    medium = factor(medium, levels = c("mBHI", "GAM", "Schaedler")),
    time = factor(time, levels = c("0", "24", "48")),
    stirred = factor(stirred, levels = c("Static", "Stirred"))
  )

###Step6.Prepare PERMANOVA table
perma_df <- as.data.frame(permanova_res) %>%
  rownames_to_column("Factor")

p_col <- grep("Pr", colnames(perma_df), value = TRUE)

perma_table <- perma_df %>%
  mutate(
    SumOfSqs = sprintf("%.3f", SumOfSqs),
    R2 = sprintf("%.4f", R2),
    F = ifelse(is.na(F), "-", sprintf("%.2f", F)),
    p = ifelse(
      is.na(.data[[p_col]]),
      "-",
      ifelse(
        .data[[p_col]] <= 0.001,
        "<0.001",
        sprintf("%.3f", .data[[p_col]])
      )
    )
  ) %>%
  select(Factor, Df, SumOfSqs, R2, F, p)

###Step7.Define PCoA plotting function
plot_pcoa <- function(df, pcoa_var) {
  ggplot(
    df,
    aes(x = PCoA1, y = PCoA2, color = medium, shape = time)
  ) +
    geom_point(size = 2.8, alpha = 0.9) +
    stat_ellipse(
      data = df,
      mapping = aes(x = PCoA1, y = PCoA2, group = medium, color = medium),
      inherit.aes = FALSE,
      level = 0.90,
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = medium_cols,
      drop = FALSE,
      name = "Medium"
    ) +
    scale_shape_manual(
      values = time_shapes,
      drop = FALSE,
      name = "Time [hrs]"
    ) +
    labs(
      x = paste0("PCoA1 (", pcoa_var[1], "%)"),
      y = paste0("PCoA2 (", pcoa_var[2], "%)")
    ) +
    theme_16S()
}

###Step8.Define PERMANOVA bar plot function
plot_permanova_bar <- function(perma_df, p_col) {

  perma_bar_df <- perma_df %>%
    filter(Factor %in% c("medium", "time", "stirred")) %>%
    mutate(
      Factor = factor(
        Factor,
        levels = c("medium", "time", "stirred")
      ),
      p_sig = case_when(
        .data[[p_col]] <= 0.001 ~ "***",
        .data[[p_col]] <= 0.01  ~ "**",
        .data[[p_col]] <= 0.05  ~ "*",
        .data[[p_col]] <= 0.1   ~ ".",
        TRUE ~ ""
      )
    )

  y_max <- max(perma_bar_df$R2, na.rm = TRUE)

  ggplot(
    perma_bar_df,
    aes(x = Factor, y = R2, fill = Factor)
  ) +
    geom_col(width = 0.72, color = "black", linewidth = 0.4) +
    geom_text(
      aes(label = sprintf("%.3f", R2)),
      vjust = -0.35,
      size = 3.6
    ) +
    geom_text(
      aes(label = p_sig),
      vjust = -1.2,
      size = 5.5
    ) +
    scale_fill_manual(
      values = factor_cols,
      breaks = c("medium", "time", "stirred"),
      drop = FALSE,
      name = "Factor"
    ) +
    coord_cartesian(ylim = c(0, y_max * 1.18)) +
    labs(
      x = NULL,
      y = expression(R^2)
    ) +
    theme_16S() +
    theme(
      legend.position = "none"
    )
}

###Step9.Generate plots
p_pcoa_global <- plot_pcoa(
  df = pcoa_df,
  pcoa_var = pcoa_var
)

p_permanova <- plot_permanova_bar(
  perma_df = perma_df,
  p_col = p_col
)

###Step10.Save outputs
ggsave(
  filename = file.path(dir_fig, "beta_pcoa_global_genus_fecesonly.pdf"),
  plot = p_pcoa_global,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(dir_fig, "beta_permanova_barplot_genus_fecesonly.pdf"),
  plot = p_permanova,
  width = 4,
  height = 6,
  dpi = 300
)

write.csv(
  perma_table,
  file = file.path(dir_tab, "beta_permanova_table_genus_fecesonly.csv"),
  row.names = FALSE
)
