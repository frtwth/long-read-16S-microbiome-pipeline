[README.md](https://github.com/user-attachments/files/26371891/README.md)
# Long-read 16S rRNA microbiome analysis pipeline (EMU-based)

This repository contains analysis scripts for long-read 16S rRNA microbiome data based on EMU-derived taxonomic profiles.  
The pipeline includes diversity analysis, differential abundance analysis, and visualization.

---

## Project Structure

The analysis is performed using the following directory layout:

```
project_root/
|-- code/
|-- data/
`-- results/
```

### Input data structure

```
data/
|-- abundance/   # Relative abundance tables
|-- counts/      # Estimated count tables
|-- metadata/    # Sample metadata
`-- taxonomy/    # Taxonomic annotation tables
```

- The abundance, count, and taxonomy tables follow a predefined format.  
- The metadata file (`metadata.rds`) is required for the analysis.  
- The `results/` directory is generated automatically during execution.  

---

## How to Run

All commands should be executed from the `code/` directory.

Example:

```
Rscript scripts/step01_alpha_diversity.R
```

All parameters, file paths, and settings are defined in the configuration file.

---

## Analysis Workflow

### Core analyses

1. Alpha diversity
```
Rscript scripts/step01_alpha_diversity.R
```

2. Beta diversity
```
Rscript scripts/step02_beta_diversity.R
Rscript scripts/step03_beta_diversity_fecesonly.R
```

3. Differential abundance
```
Rscript scripts/step04_split_for_da.R
Rscript scripts/step05_da_analysis.R
Rscript scripts/step06_prepare_heatmap_data.R
```

---

### Additional analyses

4. Phylum-level community composition
```
Rscript scripts/plot/plot_phylum_composition.R
```

5. Firmicutes-to-Bacteroidetes (F/B) ratio
```
Rscript scripts/plot/plot_fb_ratio.R
```

---

### Visualization scripts

```
Rscript scripts/plot/plot_alpha_diversity.R
Rscript scripts/plot/plot_beta_diversity.R
Rscript scripts/plot/plot_beta_diversity_fecesonly.R
Rscript scripts/plot/plot_da_analysis.R
```

---

## Notes

- The pipeline is tailored to a specific EMU-based input format used in our analysis workflow.  
- Direct application to external datasets may require modification of input formats and file paths.  
- Raw sequencing data are not included in this repository.  
- Designed for structured microbiome data analysis workflows.
