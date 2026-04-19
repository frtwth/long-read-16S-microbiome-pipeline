#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-04-15
##### Code definition: Plot phylum composition

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggh4x)
  library(scales)
  library(yaml)
})

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
fig_cfg <- config$figures$phylum_composition

###Step2.Load custom plotting theme
path_theme <- file.path(project_root, "code/src/figure/plot_theme.R")
source(path_theme)

###Step3.Define input and output paths
path_abundance <- file.path(project_root, fig_cfg$dir_in_rel$abundance_phylum)
path_metadata <- file.path(project_root, fig_cfg$dir_in_rel$metadata)

dir_fig <- file.path(project_root, fig_cfg$dir_out_rel)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

###Step4.Load data
phylum_abundance_df <- fread(path_abundance)
metadata_df <- readRDS(path_metadata)

###Step5.Preprocess phylum abundance table
phylum_abundance_df <- phylum_abundance_df %>%
  filter(!is.na(phylum) & phylum != "")

phylum_long_df <- phylum_abundance_df %>%
  pivot_longer(
    cols = -phylum,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  mutate(
    sample_id = gsub("_both$", "", sample_id),
    sample_id = gsub("^BC_", "repBC", sample_id),
    sample_id = gsub("_PL", "_Pl", sample_id)
  )

###Step6.Order phyla by mean relative abundance across all samples
phylum_order <- phylum_long_df %>%
  group_by(phylum) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  pull(phylum)

phylum_long_df <- phylum_long_df %>%
  mutate(
    phylum = factor(phylum, levels = phylum_order)
  )

###Step7.Merge metadata
phylum_long_df <- phylum_long_df %>%
  left_join(metadata_df, by = "sample_id")

###Step8.Calculate replicate mean
phylum_mean_df <- phylum_long_df %>%
  group_by(medium, treatment, stirred, time, condition, phylum) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_abundance = ifelse(is.na(mean_abundance), 0, mean_abundance)
  )

###Step9.Prepare plotting variables
phylum_mean_df <- phylum_mean_df %>%
  mutate(
    block = case_when(
      treatment == "FecesOnly" ~ "Baseline",
      treatment %in% c("WaterControl", "CGA") ~ "CGA",
      treatment %in% c("DMSOControl", "Quercetin") ~ "Quercetin"
    ),
    display_group = paste(treatment, stirred, sep = "_")
  )

phylum_mean_df <- phylum_mean_df %>%
  mutate(
    display_group = factor(
      display_group,
      levels = c(
        "FecesOnly_Static", "FecesOnly_Stirred",
        "WaterControl_Static", "WaterControl_Stirred",
        "CGA_Static", "CGA_Stirred",
        "DMSOControl_Static", "DMSOControl_Stirred",
        "Quercetin_Static", "Quercetin_Stirred"
      )
    ),
    medium = factor(
      medium,
      levels = c("mBHI", "GAM", "Schaedler")
    ),
    block = factor(
      block,
      levels = c("Baseline", "CGA", "Quercetin")
    ),
    time = factor(
      time,
      levels = c("0", "24", "48")
    )
  )

###Step10.Define labels and colors
x_labels <- c(
  "FecesOnly_Static" = "Base-S",
  "FecesOnly_Stirred" = "Base-T",
  "WaterControl_Static" = "Water-S",
  "WaterControl_Stirred" = "Water-T",
  "CGA_Static" = "CGA-S",
  "CGA_Stirred" = "CGA-T",
  "DMSOControl_Static" = "DMSO-S",
  "DMSOControl_Stirred" = "DMSO-T",
  "Quercetin_Static" = "Quercetin-S",
  "Quercetin_Stirred" = "Quercetin-T"
)

phylum_colors <- c(
  "Firmicutes" = "#6E94A6",
  "Bacteroidetes" = "#243B5A",
  "Actinobacteria" = "#7E8F4D",
  "Proteobacteria" = "#1F5B4A",
  "Tenericutes" = "#E08B6D",
  "Synergistetes" = "#B74C2E",
  "Planctomycetes" = "#D89A9A",
  "Verrucomicrobia" = "#8C6A4F"
)

###Step11.Define phylum composition plotting function
plot_phylum_composition <- function(df) {
  ggplot(
    df,
    aes(x = display_group, y = mean_abundance, fill = phylum)
  ) +
    geom_bar(
      stat = "identity",
      width = 0.95,
      colour = NA
    ) +
    ggh4x::facet_nested(
      rows = vars(time),
      cols = vars(medium, block),
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_manual(
      values = phylum_colors,
      drop = FALSE,
      name = "Phylum"
    ) +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(
      labels = scales::percent,
      expand = c(0.01, 0.01)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = "Relative abundance"
    ) +
    theme_16S() +
    theme(
      panel.spacing.x = unit(0.15, "lines"),
      panel.spacing.y = unit(0.45, "lines"),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text = element_text(
        face = "bold",
        size = 10,
        margin = margin(2, 2, 2, 2)
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 9.5
      ),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(face = "bold", size = 11),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.5, "cm")
    )
}

###Step12.Generate plot
p_phylum_composition <- plot_phylum_composition(
  df = phylum_mean_df
)

###Step13.Save output
ggsave(
  filename = file.path(dir_fig, "phylum_composition_barplot.pdf"),
  plot = p_phylum_composition,
  width = 11.5,
  height = 8,
  dpi = 300
)
