options(scipen = 999, stringsAsFactors = F)

setwd("data/MelanomaDrugResistance/SingleCellmRNAFISH/")
seed <- 20170708

files <- c("WM9_noDrug_20150618.txt",
           "WM9_noDrug_20150810.txt",
           "WM9_48hrs_20150624.txt",
           "WM9_48hrs_20151023.txt",
           "WM9_48hrs_20160302.txt",
           "WM9_1week_20150715.txt",
           "WM9_1week_20150821.txt",
           "WM9_1week_20160225.txt",
           "WM9_4weeks_20150701.txt",
           "WM9_4weeks_20150831.txt")

data <- do.call("rbind", lapply(files, function(x) {
  y <- read.csv(x)
  n <- gsub("[.]txt", "", x)
  rownames(y) <- paste(n, y[,1], sep = "_")
  y <- y[,-(1:3)]
  return(y)
}))

metadata <- as.data.frame(stringr::str_split_fixed(rownames(data), "_", 4))
colnames(metadata) <- c("cell_line", "timepoint", "date", "id")
rownames(metadata) <- rownames(data)

# Remove all-zero cells
metadata <- metadata[rowSums(data) > 0,]
data <- data[rowSums(data) > 0,]

saveRDS(data, "WM9_data.rds")
saveRDS(metadata, "WM9_metadata.rds")

# Downsampled sets
set.seed(seed)
x2000 <- do.call(c, lapply(unique(metadata$timepoint), function(t) sample(rownames(metadata[metadata$timepoint==t,]), 2000)))
x1000 <- do.call(c, lapply(unique(metadata$timepoint), function(t) sample(rownames(metadata[metadata$timepoint==t,]), 1000)))
x500 <- do.call(c, lapply(unique(metadata$timepoint), function(t) sample(rownames(metadata[metadata$timepoint==t,]), 500)))
x100 <- do.call(c, lapply(unique(metadata$timepoint), function(t) sample(rownames(metadata[metadata$timepoint==t,]), 100)))

saveRDS(data[x2000,], "WM9_data_2000.rds")
saveRDS(metadata[x2000,], "WM9_metadata_2000.rds")

saveRDS(data[x1000,], "WM9_data_1000.rds")
saveRDS(metadata[x1000,], "WM9_metadata_1000.rds")

saveRDS(data[x500,], "WM9_data_500.rds")
saveRDS(metadata[x500,], "WM9_metadata_500.rds")

saveRDS(data[x100,], "WM9_data_100.rds")
saveRDS(metadata[x100,], "WM9_metadata_100.rds")
