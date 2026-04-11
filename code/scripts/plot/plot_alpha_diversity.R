#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-04-11
##### Code definition: Plot alpha diversity results

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(yaml)
})

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
fig_cfg <- config$figures$alpha_diversity

###Step2.Load custom plotting theme
path_theme <- file.path(project_root, "code/src/figure/plot_theme.R")
source(path_theme)

###Step3.Define input and output paths
path_alpha <- file.path(project_root, "results/alpha/alpha_diversity_genus.rds")
dir_fig <- file.path(project_root, "results/figures/alpha")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

###Step4.Load alpha-diversity results
alpha_meta <- readRDS(path_alpha)

###Step5.Convert variables to factors
alpha_meta <- alpha_meta %>%
  mutate(
    medium = factor(medium, levels = c("mBHI", "GAM", "Schaedler")),
    time = factor(time, levels = c("0", "24", "48")),
    stirred = factor(stirred, levels = c("Static", "Stirred")),
    treatment = factor(
      treatment,
      levels = c("CGA", "DMSOControl", "FecesOnly", "Quercetin", "WaterControl")
    )
  )

###Step6.Define FaecesOnly subgroup
alpha_faeces <- alpha_meta %>%
  filter(treatment == "FecesOnly")

###Step7.Define plotting function
plot_alpha_box <- function(df, metric, ylab_text, facet_formula, use_treatment_color = TRUE) {

  p <- ggplot(df, aes(x = time, y = .data[[metric]])) +
    geom_boxplot(
      outlier.shape = NA,
      fill = "grey90",
      color = "black",
      linewidth = 0.4
    ) +
    facet_grid(facet_formula) +
    labs(
      x = "Time [hrs]",
      y = ylab_text
    ) +
    theme_16S()

  if (use_treatment_color) {
    p <- p +
      geom_jitter(
        aes(color = treatment),
        width = 0.12,
        alpha = 0.75,
        size = 1.8
      ) +
      scale_color_manual(
        values = treatment_cols,
        drop = FALSE,
        name = "Treatment"
      )
  } else {
    p <- p +
      geom_jitter(
        color = "black",
        width = 0.12,
        alpha = 0.75,
        size = 1.8
      ) +
      theme(
        legend.position = "none"
      )
  }

  return(p)
}

###Step8.Define metrics
metric_info <- list(
  shannon = list(
    column = "Shannon",
    ylab = "Shannon index"
  ),
  simpson = list(
    column = "Simpson",
    ylab = "Simpson index"
  ),
  richness = list(
    column = "Richness",
    ylab = "Observed richness"
  ),
  evenness = list(
    column = "Evenness",
    ylab = "Pielou's evenness"
  )
)

###Step9.Define plot settings
plot_sets <- list(
  all = list(
    df = alpha_meta,
    facet_formula = ~ medium,
    use_treatment_color = TRUE,
    width = 8,
    height = 4
  ),
  all_stirred = list(
    df = alpha_meta,
    facet_formula = stirred ~ medium,
    use_treatment_color = TRUE,
    width = 11,
    height = 7
  ),
  faeces = list(
    df = alpha_faeces,
    facet_formula = ~ medium,
    use_treatment_color = FALSE,
    width = 8,
    height = 4
  ),
  faeces_stirred = list(
    df = alpha_faeces,
    facet_formula = stirred ~ medium,
    use_treatment_color = FALSE,
    width = 11,
    height = 7
  )
)

###Step10.Generate plots
plot_list <- list()

for (set_name in names(plot_sets)) {
  for (metric_name in names(metric_info)) {

    plot_key <- paste(metric_name, set_name, sep = "_")

    plot_list[[plot_key]] <- plot_alpha_box(
      df = plot_sets[[set_name]]$df,
      metric = metric_info[[metric_name]]$column,
      ylab_text = metric_info[[metric_name]]$ylab,
      facet_formula = plot_sets[[set_name]]$facet_formula,
      use_treatment_color = plot_sets[[set_name]]$use_treatment_color
    )
  }
}

###Step11.Save plots
for (set_name in names(plot_sets)) {
  for (metric_name in names(metric_info)) {

    plot_key <- paste(metric_name, set_name, sep = "_")
    file_name <- paste0("alpha_", metric_name, "_", set_name, "_genus.pdf")

    ggsave(
      filename = file.path(dir_fig, file_name),
      plot = plot_list[[plot_key]],
      width = plot_sets[[set_name]]$width,
      height = plot_sets[[set_name]]$height,
      dpi = 300
    )
  }
}
