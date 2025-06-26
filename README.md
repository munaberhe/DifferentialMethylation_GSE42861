# Differential Methylation Analysis Pipeline (GSE42861)

This repository provides a complete pipeline for analyzing differential methylation using Illumina 450k array data from GEO accession [GSE42861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42861). The workflow is implemented in R and covers data acquisition, preprocessing, quality control, normalization, differential analysis using the `limma` package, annotation, and visualization.

---

## 1. Project Structure

```
DifferentialMethylationPipeline/
├── data/                      # Raw IDAT files (downloaded from GEO)
├── results/                   # Output results: CSVs, plots
├── scripts/                   # R scripts for each analysis step
├── plots/                     # PCA, Volcano, and MA plots (PDF/PNG)
├── pheno.csv                 # Sample metadata (group assignments)
├── README.md
```

---

## 2. Getting Started

### Prerequisites

Make sure R and the following packages are installed:

```r
install.packages(c("tidyverse", "ggplot2", "ggfortify"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("minfi", "limma", "GEOquery",
                       "IlluminaHumanMethylation450kmanifest",
                       "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
```

---

## 3. Pipeline Overview

### Step 1: Data Download & Preparation

Download and extract the IDAT files from GEO:

```bash
# Inside scripts/download_data.sh
mkdir -p data
cd data
# (Manually download and extract GSE42861 supplementary tarball)
tar -xvf GSE42861_RAW.tar
```

### Step 2: Load IDATs

```r
library(minfi)
setwd("data/")
rgSet <- read.metharray.exp(base = getwd())
```

### Step 3: Quality Control & Normalization

```r
detP <- detectionP(rgSet)
sample_qc <- colMeans(detP)
rgSet_filtered <- rgSet[, sample_qc < 0.01]
mSet <- preprocessQuantile(rgSet_filtered)

keep_probes <- rowSums(detP < 0.01) == ncol(rgSet_filtered)
mSet <- mSet[keep_probes, ]
```

### Step 4: Create `pheno.csv`

Ensure `pheno.csv` matches your sample IDs and contains a `group` column (e.g., Case, Control):

```
SampleID,group
GSM1052203,Case
GSM1052204,Control
...
```

Then:

```r
pheno <- read.csv("pheno.csv", row.names = 1)
mSet <- mSet[, rownames(pheno)]
```

---

## 4. Differential Methylation Analysis

### Step 5: Run Limma

```r
m_values <- getM(mSet)
design <- model.matrix(~ pheno$group)
colnames(design) <- c("Intercept", "GroupCase")

library(limma)
fit <- lmFit(m_values, design)
fit <- eBayes(fit)
top_dmps <- topTable(fit, coef = "GroupCase", number = Inf, adjust.method = "BH")
write.csv(top_dmps, "results/DMP_results.csv")
```

---

## 5. Annotation & Visualization

### Volcano Plot

```r
library(ggplot2)
top_dmps$significant <- top_dmps$adj.P.Val < 0.05 & abs(top_dmps$logFC) > 1

ggplot(top_dmps, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Methylated Positions")
```

### Annotation

```r
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
top_dmps_annotated <- merge(top_dmps, ann450k, by.x = "row.names", by.y = "Name")
colnames(top_dmps_annotated)[1] <- "CpG"
write.csv(top_dmps_annotated, "results/DMP_results_annotated.csv", row.names = FALSE)
```

### Beta Distribution of Top DMPs

```r
beta <- getBeta(mSet)
sig_cpgs <- head(top_dmps_annotated$CpG, 5)

for (cpg in sig_cpgs) {
  df <- data.frame(Beta = beta[cpg, ], Group = pheno$group)
  p <- ggplot(df, aes(x = Group, y = Beta, fill = Group)) +
    geom_boxplot() +
    labs(title = paste("Beta Values for", cpg), y = "Beta Value") +
    theme_minimal()
  print(p)
}
```

### PCA Plot

```r
library(ggfortify)
top_var_probes <- head(order(apply(beta, 1, var), decreasing = TRUE), 5000)
pca_res <- prcomp(t(beta[top_var_probes, ]), scale. = TRUE)

autoplot(pca_res, data = pheno, colour = "group") +
  labs(title = "PCA of Top 5000 Variable Probes") +
  theme_minimal()
```

---

## 6. Outputs

- `results/DMP_results.csv`: Raw DMP table
- `results/DMP_results_annotated.csv`: Annotated CpG table
- `plots/`: PCA, volcano, and boxplots for significant CpGs

---

## 7. Best Practices

- Run all scripts in a clean R environment
- Ensure your pheno.csv matches your sample IDs exactly
- Visualize and validate QC steps carefully
- Adjust significance thresholds as needed


---

## Author

Muna Berhe  
