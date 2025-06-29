#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(minfi)
  library(limma)
  library(ggplot2)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

# Set working directory
setwd("~/Desktop/DifferentialMethylationPipeline/data")

# Load normalized methylation data (you should already have mSet saved in memory)
# If not, load RGChannelSet and preprocess again

# Extract M-values and Beta values
m_values <- getM(mSet)
beta <- getBeta(mSet)

# Load phenotype metadata
pheno_full <- read.csv("pheno.csv", row.names = 1, stringsAsFactors = FALSE)

# Match and clean phenotype data
sample_ids <- sapply(strsplit(colnames(mSet), "_"), `[`, 1)
rownames(pheno_full) <- pheno_full$geo_accession
pheno <- pheno_full[sample_ids, ]
rownames(pheno) <- sample_ids

# Extract and clean group variable
pheno$group <- sub("disease state: ", "", pheno$characteristics_ch1.1)
pheno$group <- factor(pheno$group)

# Ensure two levels exist
if (length(levels(pheno$group)) != 2) {
  stop("Group variable must have exactly 2 levels.")
}

# Create design matrix
design <- model.matrix(~ pheno$group)
colnames(design) <- c("Intercept", paste0("Group", levels(pheno$group)[2]))

# Fit linear model
fit <- lmFit(m_values, design)
fit <- eBayes(fit)

# Get all DMPs
top_dmps <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

# Save results
dir.create("../results", showWarnings = FALSE)
write.csv(top_dmps, "../results/DMP_results.csv")

# Volcano plot
top_dmps$significant <- top_dmps$adj.P.Val < 0.05 & abs(top_dmps$logFC) > 1

volcano <- ggplot(top_dmps, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Methylated Positions")

ggsave("../plots/volcano_plot.png", plot = volcano, width = 8, height = 6)

# Annotate top CpGs
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
top_dmps_annotated <- merge(top_dmps, ann450k, by = "row.names")
colnames(top_dmps_annotated)[1] <- "CpG"

write.csv(top_dmps_annotated, "../results/DMP_results_annotated.csv", row.names = FALSE)

# Optional: Plot top 5 CpGs
top_cpgs <- head(top_dmps_annotated$CpG, 5)
for (cpg in top_cpgs) {
  df <- data.frame(
    Beta = beta[cpg, ],
    Group = pheno$group
  )
  p <- ggplot(df, aes(x = Group, y = Beta, fill = Group)) +
    geom_boxplot() +
    labs(title = paste("Beta Value Distribution -", cpg)) +
    theme_minimal()
  ggsave(paste0("../plots/", cpg, "_boxplot.png"), plot = p, width = 6, height = 4)
}

