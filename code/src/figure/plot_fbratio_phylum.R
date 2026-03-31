#!/usr/bin/Rscript

#### Code by S. Jang
#### Last updated: 2026-03-21
#### Code definition: Plot phylum composition

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(grid)

###Step1.Load data
phylum_abundance.df <- fread("/Users/celina/Desktop/project/16s_metagenomics_longread/data/abundance/emu-combined-abundance-phylum.tsv")
metadata.df <- readRDS("/Users/celina/Desktop/project/16s_metagenomics_longread/data/metadata/metadata.df.rds")

###Step2.Preprocess phylum abundance table
phylum_abundance.df <- phylum_abundance.df %>%
  filter(!is.na(phylum) & phylum != "")

phylum_long.df <- phylum_abundance.df %>%
  pivot_longer(
    -phylum,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  mutate(
    sample_id = gsub("_both$", "", sample_id),
    sample_id = gsub("^BC_", "repBC", sample_id),
    sample_id = gsub("_PL", "_Pl", sample_id)
  )

###Step3.Merge metadata
phylum_long.df <- phylum_long.df %>%
  left_join(
    metadata.df,
    by = "sample_id"
  )

###Step4.Calculate F/B ratio for each replicate
fb_ratio.df <- phylum_long.df %>%
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

###Step5.Calculate condition mean
fb_ratio_mean.df <- fb_ratio.df %>%
  group_by(medium, treatment, stirred, time, condition) %>%
  summarise(
    mean_FB_ratio = mean(FB_ratio, na.rm = TRUE),
    sd_FB_ratio = sd(FB_ratio, na.rm = TRUE),
    n = sum(!is.na(FB_ratio)),
    se_FB_ratio = sd_FB_ratio / sqrt(n),
    .groups = "drop"
  )

fb_ratio_mean.df$mean_FB_ratio[is.na(fb_ratio_mean.df$mean_FB_ratio)] <- 0
fb_ratio_mean.df$sd_FB_ratio[is.na(fb_ratio_mean.df$sd_FB_ratio)] <- 0
fb_ratio_mean.df$se_FB_ratio[is.na(fb_ratio_mean.df$se_FB_ratio)] <- 0

###Step6.Prepare plotting variables
fb_ratio_mean.df <- fb_ratio_mean.df %>%
  mutate(
    block = case_when(
      treatment == "FecesOnly" ~ "Baseline",
      treatment %in% c("WaterControl", "CGA") ~ "CGA",
      treatment %in% c("DMSOControl", "Quercetin") ~ "Quercetin"
    ),
    display_group = paste(treatment, stirred, sep = "_")
  )

fb_ratio_mean.df$display_group <- factor(
  fb_ratio_mean.df$display_group,
  levels = c(
    "FecesOnly_Static", "FecesOnly_Stirred",
    "WaterControl_Static", "WaterControl_Stirred",
    "CGA_Static", "CGA_Stirred",
    "DMSOControl_Static", "DMSOControl_Stirred",
    "Quercetin_Static", "Quercetin_Stirred"
  )
)

fb_ratio_mean.df$medium <- factor(
  fb_ratio_mean.df$medium,
  levels = c("mBHI", "GAM", "Schaedler")
)

fb_ratio_mean.df$block <- factor(
  fb_ratio_mean.df$block,
  levels = c("Baseline", "CGA", "Quercetin")
)

fb_ratio_mean.df$time <- factor(
  fb_ratio_mean.df$time,
  levels = c("0", "24", "48")
)

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

###Step7.Plot F/B ratio
p <- ggplot(
  fb_ratio_mean.df,
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
  theme_bw() +
  theme(
    panel.grid = element_blank(),

    # x/y panel spacing separately
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

ggsave("/Users/celina/Desktop/project/16s_metagenomics_longread/results/figures/community/fb_ratio_plot.pdf", p, width = 11.5, height = 6.5, dpi = 300)

###Step8.Save replicate-level F/B ratio table
fwrite(fb_ratio.df,"/Users/celina/Desktop/project/16s_metagenomics_longread/results/tables/community/fb_ratio_by_replicate.tsv", sep = "\t")

###Step9.Save condition-level summary table
fb_ratio_summary_export.df <- fb_ratio_mean.df %>%
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

fwrite(fb_ratio_summary_export.df, "/Users/celina/Desktop/project/16s_metagenomics_longread/results/tables/community/fb_ratio_summary.tsv", sep = "\t")
