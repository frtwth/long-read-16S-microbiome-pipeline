#!/usr/bin/env Rscript

##### Code by S. Jang
##### Last updated: 2026-03-08
##### Code definition: ggplot2 theme setting.

library(ggplot2)

### Colors
medium_cols <- c(
  "mBHI" = "#1b9e77",
  "GAM" = "#d95f02",
  "Schaedler" = "#7570b3"
)

factor_cols <- c(
  medium    = "#8ECae6",  # soft sky blue
  time      = "#F4A261",  # pastel orange
  treatment = "#A8DADC",  # pastel mint
  stirred   = "#CED4DA"   # light gray
)

time_shapes <- c(
  "0" = 15,
  "24" = 16,
  "48" = 17
)

treatment_cols <- c(
  CGA          = "#F8766D",  # red
  DMSOControl  = "#7CAE00",  # yellow-green
  FecesOnly    = "#00BFC4",  # green-cyan
  Quercetin    = "#619CFF",  # blue 
  WaterControl = "#C77CFF"   # purple
)

### Common theme
theme_16S <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", size = 12),
 
      text = element_text(size = 11),
      strip.text = element_text(size = 11),

      legend.position = "right",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10)
    )
}
