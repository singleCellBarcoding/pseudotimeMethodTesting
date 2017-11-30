library(monocle)

data_path <- "../data/MelanomaDrugResistance/SingleCellmRNAFISH/"

# Read data matrix and metadata
data <- readRDS(paste0(data_path, "WM9_data_500.rds"))
metadata <- readRDS(paste0(data_path, "WM9_metadata_500.rds"))

# Remove lowest 5% of the data
#keep <- rowSums(data) > quantile(rowSums(data), 0.05)
#data <- data[keep,]
#metadata <- metadata[keep,]

# Make CellDataSet object from data
expr <- as.matrix(t(data))
pd <- as.data.frame(metadata)
rownames(pd) <- rownames(data)
pd <- new("AnnotatedDataFrame", data = pd)
fd <- as.data.frame(cbind(colnames(data), colnames(data)))
colnames(fd) <- c("gene_id", "gene_short_name")
rownames(fd) <- fd$gene_id
fd <- new("AnnotatedDataFrame", data = fd)
D <- newCellDataSet(expr, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size(), lowerDetectionLimit = 1)

# Estimate size factors and dispersions
D <- estimateSizeFactors(D)
D <- estimateDispersions(D)

# Reduce dataset dimensionality
D <- reduceDimension(D, max_components = 2, method = "DDRTree")

# Order cells in pseudotime
D <- orderCells(D)

# Plot cell trajectory in the reduced dimensional space
pdf("WM9_500_monocle2_plots.pdf")
plot_cell_trajectory(D, color_by = "timepoint")
plot_cell_trajectory(D, color_by = "Pseudotime")
graphics.off()

colours <- c("#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c")

df_class <- data.frame(D@phenoData@data$Pseudotime, factor(metadata$timepoint, levels = c("noDrug", "48hrs", "1week", "4weeks")))
colnames(df_class) <- c("pseudotime", "timepoint")

saveRDS(df_class, "WM9_500_monocle2_pseudotime.rds")

pdf("WM9_500_monocle2_pseudotime.pdf")
ggplot(df_class, aes(x = timepoint, y = pseudotime, fill = timepoint)) +
  geom_violin(show.legend = F) +
  xlab("Cell timepoint") +
  ylab("Pseudotime") +
  ggtitle("Monocle2") +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = colours)
graphics.off()
