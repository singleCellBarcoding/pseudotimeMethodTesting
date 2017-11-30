library(destiny)

data_path <- "../data/MelanomaDrugResistance/SingleCellmRNAFISH/"

# Read data matrix and metadata
data <- readRDS(paste0(data_path, "WM9_data_500.rds"))
metadata <- readRDS(paste0(data_path, "WM9_metadata_500.rds"))

# Remove lowest 5% of the data
#keep <- rowSums(data) > quantile(rowSums(data), 0.05)
#data <- data[keep,]
#metadata <- metadata[keep,]

# Create diffusion map
dm <- DiffusionMap(log2(data+1))
saveRDS(dm, "WM9_500_destiny_model.rds")

colours <- c("#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c")
tmp <- factor(metadata[,2], levels = c("noDrug", "48hrs", "1week", "4weeks"))
pdf("WM9_500_destiny_model_map.rds")
plot(eigenvectors(dm)[,1:2], xlab = "Diffusion component 1", ylab = "Diffusion component 2", pch = 16, col = colours[tmp], main = "Destiny diffusion map")
graphics.off()

dpt <- DPT(dm)

# Get pseudotime from dpt object...
###
branch_divide <- function(dpt, divide = integer(0L)) {

  if (length(divide) == 0L) return(dpt)
  
  for (b in divide) {
    super_rows <- dpt@branch[, 1] == b & !is.na(dpt@branch[, 1])
    if (!any(super_rows)) {
      available <- na.omit(unique(dpt@branch[, 1]))
      stop('invalid branch to divide ', b, ' not in ', available)
    }
    
    # shift sub branches/tips to the left
    dpt@branch[super_rows, ] <- cbind(dpt@branch[super_rows, -1], NA)
    dpt@tips  [super_rows, ] <- cbind(dpt@tips  [super_rows, -1], NA)
    
    # TODO: maybe also modify DPT?
  }
  
  vacant_levels <- apply(dpt@branch, 2L, function(col) all(is.na(col)))
  dpt@branch <- dpt@branch[, !vacant_levels]
  dpt@tips   <- dpt@tips  [, !vacant_levels]
  
  dpt
}

dpt_for_branch <- function(dpt, branch_id) {
  branch_idx <- dpt@branch[, 1L] == branch_id
  stopifnot(any(branch_idx))
  tip_cells <- which(branch_idx & dpt@tips[, 1L])
  if (length(tip_cells) == 0L) tip_cells <- which(branch_idx)
  dpt[tip_cells[[1L]], ]
}

dpt_flat <- branch_divide(dpt)
root <- min(dpt_flat@branch[, 1], na.rm = TRUE)
pt_vec <- dpt_for_branch(dpt_flat, root)
######

df_class <- data.frame(pt_vec, factor(metadata$timepoint, levels = c("noDrug", "48hrs", "1week", "4weeks")))
colnames(df_class) <- c("pseudotime", "timepoint")

pdf("WM9_500_destiny_pseudotime.pdf")
ggplot(df_class, aes(x = timepoint, y = pseudotime, fill = timepoint)) +
  geom_violin(show.legend = F) +
  xlab("Cell timepoint") +
  ylab("Pseudotime") +
  ggtitle("destiny") +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = colours)
graphics.off()
