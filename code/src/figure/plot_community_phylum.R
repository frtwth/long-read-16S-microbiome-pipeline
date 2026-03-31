#!/usr/bin/Rscript

#### Code by S. Jang
#### Last updated: 2026-03-21
#### Code definition: Plot phylum composition

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)

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

#Order phyla by mean relative abundance across all samples
phylum_order <- phylum_long.df %>%
  group_by(phylum) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  pull(phylum)

phylum_long.df$phylum <- factor(
  phylum_long.df$phylum,
  levels = phylum_order
)

###Step3.Merge metadata
phylum_long.df <- phylum_long.df %>%
  left_join(
    metadata.df,
    by = "sample_id"
  )

###Step4.Calculate replicate mean
phylum_mean.df <- phylum_long.df %>%
  group_by(medium, treatment, stirred, time, condition, phylum) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = "drop"
  )

phylum_mean.df$mean_abundance[is.na(phylum_mean.df$mean_abundance)] <- 0

#Check sum
#phylum_mean.df %>%
#  group_by(condition) %>%
#  summarise(total = sum(mean_abundance))

###Step5.Bar plot
phylum_mean.df <- phylum_mean.df %>%
  mutate(
    block = case_when(
      treatment == "FecesOnly" ~ "Baseline",
      treatment %in% c("WaterControl", "CGA") ~ "CGA",
      treatment %in% c("DMSOControl", "Quercetin") ~ "Quercetin"
    ),
    display_group = paste(treatment, stirred, sep = "_")
  )

phylum_mean.df$display_group <- factor(
  phylum_mean.df$display_group,
  levels = c(
    "FecesOnly_Static", "FecesOnly_Stirred",
    "WaterControl_Static", "WaterControl_Stirred",
    "CGA_Static", "CGA_Stirred",
    "DMSOControl_Static", "DMSOControl_Stirred",
    "Quercetin_Static", "Quercetin_Stirred"
  )
)

phylum_mean.df$medium <- factor(
  phylum_mean.df$medium,
  levels = c("mBHI", "GAM", "Schaedler")
)

phylum_mean.df$block <- factor(
  phylum_mean.df$block,
  levels = c("Baseline", "CGA", "Quercetin")
)

phylum_mean.df$time <- factor(
  phylum_mean.df$time,
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

p <- ggplot(
  phylum_mean.df,
  aes(x = display_group, y = mean_abundance, fill = phylum)
) +
  geom_bar(stat = "identity", width = 0.95, colour = NA) +
  ggh4x::facet_nested(
    rows = vars(time),
    cols = vars(medium, block),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(values = phylum_colors, drop = FALSE) +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(
    labels = percent,
    expand = c(0.01, 0.01)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Phylum"
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
    axis.title.y = element_text(face = "bold", size = 11),

    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")
  )

ggsave("/Users/celina/Desktop/project/16s_metagenomics_longread/results/figures/community/community_barplot_phylum.pdf", p, width = 11.5, height = 8, dpi = 300)
