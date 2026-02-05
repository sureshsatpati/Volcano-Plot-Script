# Volcano-Plot-Script
Script to plot volcano plot for DEG analysis

# Load required libraries
library(ggplot2)
library(dplyr)

# Example: DESeq2 results dataframe
# res <- results(dds)
# res_df <- as.data.frame(res)
# res_df$gene <- rownames(res_df)

# If already have a DEG table, load it
# res_df <- read.csv("DEG_results.csv")

# Remove NA values
res_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

# Define significance thresholds
logFC_cutoff <- 1
padj_cutoff  <- 0.05

# Add DEG classification
res_df <- res_df %>%
  mutate(
    DEG_status = case_when(
      padj < padj_cutoff & log2FoldChange >= logFC_cutoff  ~ "Upregulated",
      padj < padj_cutoff & log2FoldChange <= -logFC_cutoff ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = DEG_status), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot of Differential Gene Expression",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Status"
  ) +
  theme_classic()
