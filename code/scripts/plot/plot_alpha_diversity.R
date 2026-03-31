#!/usr/bin/Rscript

#### Code by S. Jang
#### Last updated: 2026-03-24
#### Code definition: Plot Alpha diversity results

library(dplyr)
library(ggplot2)
library(yaml)

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
fig_cfg <- config$figures$alpha_diversity

###Step2.Load custom plotting theme
path_theme <- file.path(project_root, "code/src/figure/plot_theme.R")
source(path_theme)

###Step3.Load alpha-diversity results
path_alpha <- file.path(project_root, "results/alpha/alpha_diversity_genus.rds")
alpha_meta <- readRDS(path_alpha)

###Step4.Convert variables to factors
alpha_meta <- alpha_meta %>%
  mutate(
    medium = factor(medium, levels = c("mBHI", "GAM", "Schaedler")),
    time = factor(time, levels = c("0", "24", "48")),
    stirred = factor(stirred, levels = c("Static", "Stirred")),
    treatment = factor(treatment)
  )

###Step5.Define FaecesOnly subgroup
alpha_faeces <- alpha_meta %>%
  filter(treatment == "FecesOnly")

###Step6.Define plotting function
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

###Step7.Generate plots for all samples
p_shannon_all <- plot_alpha_box(
  df = alpha_meta,
  metric = "Shannon",
  ylab_text = "Shannon index",
  facet_formula = ~ medium,
  use_treatment_color = TRUE
)

p_simpson_all <- plot_alpha_box(
  df = alpha_meta,
  metric = "Simpson",
  ylab_text = "Simpson index",
  facet_formula = ~ medium,
  use_treatment_color = TRUE
)

p_richness_all <- plot_alpha_box(
  df = alpha_meta,
  metric = "Richness",
  ylab_text = "Observed richness",
  facet_formula = ~ medium,
  use_treatment_color = TRUE
)

p_evenness_all <- plot_alpha_box(
  df = alpha_meta,
  metric = "Evenness",
  ylab_text = "Pielou's evenness",
  facet_formula = ~ medium,
  use_treatment_color = TRUE
)

###Step8.Generate plots for all samples facetted by stirred
p_shannon_all_stirred <- plot_alpha_box(
  df = alpha_meta,
  metric = "Shannon",
  ylab_text = "Shannon index",
  facet_formula = stirred ~ medium,
  use_treatment_color = TRUE
)

p_simpson_all_stirred <- plot_alpha_box(
  df = alpha_meta,
  metric = "Simpson",
  ylab_text = "Simpson index",
  facet_formula = stirred ~ medium,
  use_treatment_color = TRUE
)

p_richness_all_stirred <- plot_alpha_box(
  df = alpha_meta,
  metric = "Richness",
  ylab_text = "Observed richness",
  facet_formula = stirred ~ medium,
  use_treatment_color = TRUE
)

p_evenness_all_stirred <- plot_alpha_box(
  df = alpha_meta,
  metric = "Evenness",
  ylab_text = "Pielou's evenness",
  facet_formula = stirred ~ medium,
  use_treatment_color = TRUE
)

###Step9.Generate plots for FaecesOnly subgroup
p_shannon_faeces <- plot_alpha_box(
  df = alpha_faeces,
  metric = "Shannon",
  ylab_text = "Shannon index",
  facet_formula = ~ medium,
  use_treatment_color = FALSE
)

p_simpson_faeces <- plot_alpha_box(
  df = alpha_faeces,
  metric = "Simpson",
  ylab_text = "Simpson index",
  facet_formula = ~ medium,
  use_treatment_color = FALSE
)

p_richness_faeces <- plot_alpha_box(
  df = alpha_faeces,
  metric = "Richness",
  ylab_text = "Observed richness",
  facet_formula = ~ medium,
  use_treatment_color = FALSE
)

p_evenness_faeces <- plot_alpha_box(
  df = alpha_faeces,
  metric = "Evenness",
  ylab_text = "Pielou's evenness",
  facet_formula = ~ medium,
  use_treatment_color = FALSE
)

###Step10.Generate plots for FaecesOnly subgroup facetted by stirred
p_shannon_faeces_stirred <- plot_alpha_box(
  df = alpha_faeces,
  metric = "Shannon",
  ylab_text = "Shannon index",
  facet_formula = stirred ~ medium,
  use_treatment_color = FALSE
)

p_simpson_faeces_stirred <- plot_alpha_box(
  df = alpha_faeces,
  metric = "Simpson",
  ylab_text = "Simpson index",
  facet_formula = stirred ~ medium,
  use_treatment_color = FALSE
)

p_richness_faeces_stirred <- plot_alpha_box(
  df = alpha_faeces,
  metric = "Richness",
  ylab_text = "Observed richness",
  facet_formula = stirred ~ medium,
  use_treatment_color = FALSE
)

p_evenness_faeces_stirred <- plot_alpha_box(
  df = alpha_faeces,
  metric = "Evenness",
  ylab_text = "Pielou's evenness",
  facet_formula = stirred ~ medium,
  use_treatment_color = FALSE
)

###Step11.Save outputs
dir_fig <- file.path(project_root, "results/figures/alpha")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

## All samples
ggsave(file.path(dir_fig, "alpha_shannon_all_genus.pdf"), p_shannon_all, width = 8, height = 4, dpi = 300)
ggsave(file.path(dir_fig, "alpha_simpson_all_genus.pdf"), p_simpson_all, width = 8, height = 4, dpi = 300)
ggsave(file.path(dir_fig, "alpha_richness_all_genus.pdf"), p_richness_all, width = 8, height = 4, dpi = 300)
ggsave(file.path(dir_fig, "alpha_evenness_all_genus.pdf"), p_evenness_all, width = 8, height = 4, dpi = 300)

## All samples facetted by stirred
ggsave(file.path(dir_fig, "alpha_shannon_all_stirred_genus.pdf"), p_shannon_all_stirred, width = 11, height = 7, dpi = 300)
ggsave(file.path(dir_fig, "alpha_simpson_all_stirred_genus.pdf"), p_simpson_all_stirred, width = 11, height = 7, dpi = 300)
ggsave(file.path(dir_fig, "alpha_richness_all_stirred_genus.pdf"), p_richness_all_stirred, width = 11, height = 7, dpi = 300)
ggsave(file.path(dir_fig, "alpha_evenness_all_stirred_genus.pdf"), p_evenness_all_stirred, width = 11, height = 7, dpi = 300)

## FaecesOnly subgroup
ggsave(file.path(dir_fig, "alpha_shannon_faeces_genus.pdf"), p_shannon_faeces, width = 8, height = 4, dpi = 300)
ggsave(file.path(dir_fig, "alpha_simpson_faeces_genus.pdf"), p_simpson_faeces, width = 8, height = 4, dpi = 300)
ggsave(file.path(dir_fig, "alpha_richness_faeces_genus.pdf"), p_richness_faeces, width = 8, height = 4, dpi = 300)
ggsave(file.path(dir_fig, "alpha_evenness_faeces_genus.pdf"), p_evenness_faeces, width = 8, height = 4, dpi = 300)

## FaecesOnly subgroup facetted by stirred
ggsave(file.path(dir_fig, "alpha_shannon_faeces_stirred_genus.pdf"), p_shannon_faeces_stirred, width = 11, height = 7, dpi = 300)
ggsave(file.path(dir_fig, "alpha_simpson_faeces_stirred_genus.pdf"), p_simpson_faeces_stirred, width = 11, height = 7, dpi = 300)
ggsave(file.path(dir_fig, "alpha_richness_faeces_stirred_genus.pdf"), p_richness_faeces_stirred, width = 11, height = 7, dpi = 300)
ggsave(file.path(dir_fig, "alpha_evenness_faeces_stirred_genus.pdf"), p_evenness_faeces_stirred, width = 11, height = 7, dpi = 300)
