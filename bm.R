library(Matrix)
install.packages('Seurat')
library(Seurat)
library(dplyr)
library(patchwork)

#Load the data
#  Based on the files we have, we'll need to load multiple samples and combine them into a single Seurat object.
sample_files <- c("niche-col23.txt", "niche-lepr.txt", "niche-vecad.txt", 
                  "5FU-cntrl.txt", "5FU-cntrl-col23.txt", 
                  "5FU-treat-col23.txt", "5FU-treat-lepr.txt", "5FU-treat-vecad.txt")

#read each file and create a Seurat object for each:
seurat_objects <- lapply(X, FUN, ...)apply(sample_files, function(file) {
  data <- read.table(file, header = TRUE, row.names = 1, sep = "\t")
  sample_name <- tools::file_path_sans_ext(file)
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
  seurat_obj$sample <- sample_name
  return(seurat_obj)
})

#Merge all seurat objects into one:
combined <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)], 
                  add.cell.ids = names(seurat_objects), project = "BM_Niche")

#Add metadata:
metadata <- read.csv("metadata.csv", row.names = 1)
combined <- AddMetaData(combined, metadata)

#Load gene annotations:
genes <- read.table("genes.tsv", header = TRUE, sep = "\t")
rownames(combined) <- genes$gene_symbol[match(rownames(combined), genes$gene_id)]

##ANALYSIS

# 1. Quality Control Violin Plots
pdf("qc_violin_plots.pdf", width = 12, height = 6)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1)
dev.off()

# 2. PCA Plot
pdf("pca_plot.pdf", width = 10, height = 8)
DimPlot(combined, reduction = "pca")
dev.off()

# 3. Elbow Plot (to help determine significant PCs)
pdf("elbow_plot.pdf", width = 8, height = 6)
ElbowPlot(combined, ndims = 50)
dev.off()

# 4. UMAP Plot
pdf("umap_plot.pdf", width = 10, height = 8)
DimPlot(combined, reduction = "umap", group.by = "seurat_clusters")
dev.off()

# 5. UMAP Plot colored by sample
pdf("umap_plot_by_sample.pdf", width = 10, height = 8)
DimPlot(combined, reduction = "umap", group.by = "sample")
dev.off()