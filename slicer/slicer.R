library(SLICER)
library(lle)

data_path <- "../data/MelanomaDrugResistance/SingleCellmRNAFISH/"

# Read data matrix and metadata
data <- readRDS(paste0(data_path, "WM9_data_500.rds"))
metadata <- readRDS(paste0(data_path, "WM9_metadata_500.rds"))

# Remove lowest 5% of the data
#keep <- rowSums(data) > quantile(rowSums(data), 0.05)
#data <- data[keep,]
#metadata <- metadata[keep,]

# Find "optimal" k for LLE
k <- select_k(log2(data+1))

data.lle <- lle(log2(data+1), m = 2, k)$Y
data.graph <- conn_knn_graph(data.lle, 5)
ends <- find_extreme_cells(data.graph, data.lle) # 604, 190
start = 190 # WM9_noDrug_20150810_1194
data.order <- cell_order(data.graph, 190)

colours <- c("#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c")
tmp <- factor(metadata[,2], levels = c("noDrug", "48hrs", "1week", "4weeks"))
pdf("WM9_500_slicer_model_LLE_map.pdf")
plot(data.lle, xlab = "LLE component 1", ylab = "LLE component 2", pch = 16, col = colours[tmp], main = "SLICER LLE")
graphics.off()
plot(data.order, tmp, xlab = "Pseudotime", ylab = "Timepoint", pch = 16, col = colours[tmp], main = "SLICER")

df_class <- data.frame(data.order, factor(metadata$timepoint, levels = c("noDrug", "48hrs", "1week", "4weeks")))
colnames(df_class) <- c("pseudotime", "timepoint")

pdf("WM9_500_slicer_pseudotime.pdf")
ggplot(df_class, aes(x = timepoint, y = pseudotime, fill = timepoint)) +
  geom_violin(show.legend = F) +
  xlab("Cell timepoint") +
  ylab("Pseudotime") +
  ggtitle("SLICER") +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = colours)
graphics.off()
