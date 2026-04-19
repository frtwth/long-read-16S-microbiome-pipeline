#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-04-17
##### Code definition: Plot Firmicutes-to-Bacteroidetes ratio based on phylum-level abundance

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggh4x)
  library(grid)
  library(yaml)
})

###Step1.Load config
config <- yaml::read_yaml("config/params.yaml")

project_root <- config$project$root
fig_cfg <- config$figures$fb_ratio_phylum

###Step2.Load custom plotting theme
path_theme <- file.path(project_root, "code/src/figure/plot_theme.R")
source(path_theme)

###Step3.Define input and output paths
path_abundance <- file.path(project_root, fig_cfg$dir_in_rel$abundance_phylum)
path_metadata <- file.path(project_root, fig_cfg$dir_in_rel$metadata)

dir_fig <- file.path(project_root, fig_cfg$dir_out_rel$figures)
dir_tab <- file.path(project_root, fig_cfg$dir_out_rel$tables)

dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

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

###Step6.Merge metadata
phylum_long_df <- phylum_long_df %>%
  left_join(metadata_df, by = "sample_id")

###Step7.Calculate F/B ratio for each replicate
fb_ratio_df <- phylum_long_df %>%
  filter(phylum %in% c("Firmicutes", "Bacteroidetes")) %>%
  select(sample_id, medium, treatment, stirred, time, condition, phylum, abundance) %>%
  pivot_wider(
    names_from = phylum,
    values_from = abundance,
    values_fill = 0
  ) %>%
  mutate(
    FB_ratio = ifelse(Bacteroidetes == 0, NA, Firmicutes / Bacteroidetes)
  )

###Step8.Calculate condition mean and standard error
fb_ratio_mean_df <- fb_ratio_df %>%
  group_by(medium, treatment, stirred, time, condition) %>%
  summarise(
    mean_FB_ratio = mean(FB_ratio, na.rm = TRUE),
    sd_FB_ratio = sd(FB_ratio, na.rm = TRUE),
    n = sum(!is.na(FB_ratio)),
    se_FB_ratio = sd_FB_ratio / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    mean_FB_ratio = ifelse(is.na(mean_FB_ratio), 0, mean_FB_ratio),
    sd_FB_ratio = ifelse(is.na(sd_FB_ratio), 0, sd_FB_ratio),
    se_FB_ratio = ifelse(is.na(se_FB_ratio), 0, se_FB_ratio)
  )

###Step9.Prepare plotting variables
fb_ratio_mean_df <- fb_ratio_mean_df %>%
  mutate(
    block = case_when(
      treatment == "FecesOnly" ~ "Baseline",
      treatment %in% c("WaterControl", "CGA") ~ "CGA",
      treatment %in% c("DMSOControl", "Quercetin") ~ "Quercetin"
    ),
    display_group = paste(treatment, stirred, sep = "_")
  ) %>%
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

###Step10.Define x-axis labels
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

###Step11.Define F/B ratio plotting function
plot_fb_ratio <- function(df) {
  ggplot(
    df,
    aes(x = display_group, y = mean_FB_ratio)
  ) +
    geom_point(size = 2.3) +
    geom_errorbar(
      aes(
        ymin = mean_FB_ratio - se_FB_ratio,
        ymax = mean_FB_ratio + se_FB_ratio
      ),
      width = 0.2,
      linewidth = 0.4
    ) +
    ggh4x::facet_nested(
      rows = vars(time),
      cols = vars(medium, block),
      scales = "free_x",
      space = "free_x"
    ) +
    scale_x_discrete(labels = x_labels) +
    labs(
      x = NULL,
      y = "Firmicutes / Bacteroidetes ratio"
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
      axis.title.y = element_text(face = "bold", size = 11)
    )
}

###Step12.Generate plot
p_fb_ratio <- plot_fb_ratio(
  df = fb_ratio_mean_df
)

###Step13.Save figure
ggsave(
  filename = file.path(dir_fig, "fb_ratio_phylum.pdf"),
  plot = p_fb_ratio,
  width = 11.5,
  height = 6.5,
  dpi = 300
)

###Step14.Save replicate-level F/B ratio table
fwrite(
  fb_ratio_df,
  file = file.path(dir_tab, "fb_ratio_by_replicate.tsv"),
  sep = "\t"
)

###Step15.Save condition-level summary table
fb_ratio_summary_export_df <- fb_ratio_mean_df %>%
  select(
    medium,
    treatment,
    stirred,
    time,
    condition,
    n,
    mean_FB_ratio,
    sd_FB_ratio,
    se_FB_ratio
  ) %>%
  arrange(medium, time, treatment, stirred)

fwrite(
  fb_ratio_summary_export_df,
  file = file.path(dir_tab, "fb_ratio_summary.tsv"),
  sep = "\t"
)
