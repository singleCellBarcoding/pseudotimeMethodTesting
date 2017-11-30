library(TSCAN)

data_path <- "../data/MelanomaDrugResistance/SingleCellmRNAFISH/"

# Read data matrix and metadata
data <- readRDS(paste0(data_path, "WM9_data_500.rds"))
metadata <- readRDS(paste0(data_path, "WM9_metadata_500.rds"))

# Remove lowest 5% of the data
#keep <- rowSums(data) > quantile(rowSums(data), 0.05)
#data <- data[keep,]
#metadata <- metadata[keep,]

# Run TSCAN
mclustdata <- exprmclust(log2(t(data)+1), clusternum = 2, reduce = T) # clusternum = 2 seems the best (with reduce); reduce = F doesn't work.

tscan_order <- TSCANorder(mclustdata, orderonly = F, flip = F)
tscan_order <- tscan_order[rownames(metadata),]

pdf("WM9_500_tscan_model_plots.pdf")
plotmclust(mclustdata, show_cell_names = F)
colours <- c("#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c")
tmp <- factor(metadata[,2], levels = c("noDrug", "48hrs", "1week", "4weeks"))
plot(mclustdata$pcareduceres[,1:2], xlab = "PC1", ylab = "PC2", pch = 16, col = colours[tmp], main = "TSCAN PCA")
#plot(tscan_order$Pseudotime, tmp, xlab = "Pseudotime", ylab = "Timepoint", pch = 16, col = colours[tmp], main = "TSCAN")
graphics.off()

df_class <- data.frame(tscan_order$Pseudotime, tmp)
colnames(df_class) <- c("pseudotime", "timepoint")

pdf("WM9_500_tscan_pseudotime.pdf")
ggplot(df_class, aes(x = timepoint, y = pseudotime, fill = timepoint)) +
  geom_violin(show.legend = F) +
  xlab("Cell timepoint") +
  ylab("Pseudotime") +
  ggtitle("TSCAN") +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = colours)
graphics.off()
