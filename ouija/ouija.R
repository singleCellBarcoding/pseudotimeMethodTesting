library(ouija)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = 6)

data_path <- "../data/MelanomaDrugResistance/SingleCellmRNAFISH/"

# Read data matrix and metadata
data <- readRDS(paste0(data_path, "WM9_data_500.rds"))
metadata <- readRDS(paste0(data_path, "WM9_metadata_500.rds"))

# Remove lowest 5% of the data
# keep <- rowSums(data) > quantile(rowSums(data), 0.05)
# data <- data[keep,]
# metadata <- metadata[keep,]

oui <- ouija(log2(as.matrix(data)+1), iter = 1000)
saveRDS(oui, "WM9_500_ouija_model.rds")

pdf("WM9_500_ouija_model_plots.pdf")
plot_diagnostics(oui)
plot_expression(oui)
plot_switch_times(oui)
#plot_peak_times(oui) # needs transient genes
graphics.off()

colours <- c("#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c")
#tmp <- factor(metadata[,2], levels = c("noDrug", "48hrs", "1week", "4weeks"))
#plot(map_pseudotime(oui), tmp, xlab = "Pseudotime", ylab = "Timepoint", pch = 16, col = colours[tmp], main = "Ouija")

df_class <- data.frame(map_pseudotime(oui), factor(metadata$timepoint, levels = c("noDrug", "48hrs", "1week", "4weeks")))
colnames(df_class) <- c("pseudotime", "timepoint")

pdf("WM9_500_ouija_pseudotime.pdf")
ggplot(df_class, aes(x = timepoint, y = pseudotime, fill = timepoint)) +
  geom_violin(show.legend = F) +
  xlab("Cell timepoint") +
  ylab("MAP pseudotime") +
  ggtitle("Ouija") +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = colours)
graphics.off()
