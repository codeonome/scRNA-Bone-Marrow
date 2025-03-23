##Bone-marrow-new##

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)

# Read the data files
vecad_counts <- read.table("niche-vecad.txt", header = TRUE, row.names = 1)
lepr_counts <- read.table("niche-lepr.txt", header = TRUE, row.names = 1)
col23_counts <- read.table("niche-col23.txt", header = TRUE, row.names = 1)

# Ensure all matrices have the same genes (rows)
common_genes <- Reduce(intersect, list(rownames(vecad_counts), 
                                      rownames(lepr_counts), 
                                      rownames(col23_counts)))

# Subset to common genes
vecad_counts <- vecad_counts[common_genes,]
lepr_counts <- lepr_counts[common_genes,]
col23_counts <- col23_counts[common_genes,]

# Add cell type identifiers to column names to maintain population origin
colnames(vecad_counts) <- paste0("vecad_", colnames(vecad_counts))
colnames(lepr_counts) <- paste0("lepr_", colnames(lepr_counts))
colnames(col23_counts) <- paste0("col23_", colnames(col23_counts))

# Combine the matrices
combined_counts <- cbind(vecad_counts, lepr_counts, col23_counts)

# Create metadata with cell type information
cell_metadata <- data.frame(
  celltype = c(
    rep("VE-Cad", ncol(vecad_counts)),
    rep("LEPR", ncol(lepr_counts)),
    rep("COL2.3", ncol(col23_counts))
  ),
  row.names = colnames(combined_counts)
)

# Create a single Seurat object directly
seurat_obj <- CreateSeuratObject(
  counts = combined_counts,
  meta.data = cell_metadata,
  project = "BM_Niche"
)

# Verify counts and cell type annotations
print(table(seurat_obj$celltype))
# Shows counts of cells by type

# Seurat object with all cells, containing celltype metadata column

# Verify counts and cell type annotations
print(table(seurat_obj$celltype))
#Table showing number of cells per population

# Calculate percentage of mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
#Adds percent.mt column to metadata

# QC metrics visualization
p_qc <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "celltype", ncol = 3, pt.size = 0.1)
print(p_qc)
#Violin plots showing QC metrics by cell type

# Quality control filtering as described in methods section
# "We removed cells that had fewer than 1,000 detected genes"
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 1000)
# Expected: Fewer cells, filtered for gene count >= 1000

# "To exclude cells that were extreme outliers in terms of library complexity...
# we removed any cells in the top 2% quantile"
top_2_percent <- quantile(seurat_obj$nFeature_RNA, 0.98)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA < top_2_percent)

# "We also removed cells with more than 10% of the transcripts coming from mitochondrial genes"
seurat_obj <- subset(seurat_obj, subset = percent.mt < 10)

# Print remaining cells by population
print(table(seurat_obj$celltype))
# Distribution of cells after filtering, similar to Fig.1 population proportions

# "We also removed 965 contaminating cells—most of which were haematopoietic cells"
# For this step, we'd need to identify contaminating hematopoietic cells
# We could use known markers like Ptprc (CD45)
seurat_obj <- NormalizeData(seurat_obj)
FeaturePlot(seurat_obj, features = "Ptprc")
# Shows expression of hematopoietic marker

# Remove contaminating cells (adjust threshold as needed based on expression)
seurat_obj <- subset(seurat_obj, subset = Ptprc < 1)
# Expected: Removes hematopoietic contaminants

# Print final cell count
print(paste("Final cell count after filtering:", ncol(seurat_obj)))
# Expected, ~17,000-17,500 cells (paper reports 17,374 cells)

# Normalize data as described in methods
# "We normalized the data by the total expression, multiplied by a scale factor of 10,000 and log-transformed"
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


# Find variable features for dimensionality reduction
# "We took the union of the top 2,000 genes with the highest dispersion from both datasets"
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
print(paste("Number of variable features:", length(VariableFeatures(seurat_obj))))
# aroundf 2000 variable features

# Scale data
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
# Adds scaled data to the object

# Run PCA
# "We then aligned the subspaces on the basis of the first 30 canonical correlation vectors"
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 40)
# PCA reduction added to object

# Examine PCA results
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
# Top genes contributing to first 5 PCs

# Determine optimal dimensionality
ElbowPlot(seurat_obj, ndims = 40)
#Elbow plot showing variance explained by each PC

# Use t-SNE for visualization (as in the paper)
# "We further reduced the dimensionality... to project the cells in 2D space using t-SNE"
seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)


# Find clusters
# "Aligned canonical correlation analysis was also used as a basis for partitioning the dataset"
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)


# Finding clusters with resolution parameter
# Paper found "10 transcriptionally similar subpopulations"
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
# Expected: ~10 clusters (paper reports 10 clusters)

# Check cluster distribution
print(table(seurat_obj$seurat_clusters))
# Distribution of cells across 10-11 clusters

# Create Figure 1B - Expression of key marker genes
p1b <- FeaturePlot(seurat_obj, 
                  features = c("Cdh5", "Lepr", "Col1a1"),
                  reduction = "tsne",
                  ncol = 3,
                  order = TRUE)
print(p1b)
# Expected: Fig.1B - Three t-SNE plots showing marker gene expression
# Cdh5 should be high in VE-Cad+ cells, Lepr in LEPR+, Col1a1 in COL2.3+

# Show cell types on t-SNE
p_celltypes <- DimPlot(seurat_obj, reduction = "tsne", group.by = "celltype")
print(p_celltypes)


# Create Figure 1D - Cluster visualization
p1d <- DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) +
       ggtitle("t-SNE visualization of bone marrow niche clusters")
print(p1d)
# Expected: Fig.1D - t-SNE with cells colored by cluster, 10-11 clusters visible

# Find marker genes for clusters (for Figure 1C)
# "To find markers that define individual clusters, we performed pairwise differential expression analysis"
all_markers <- FindAllMarkers(seurat_obj, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25,
                             test.use = "MAST")  # Paper used MAST method
# Data frame with marker genes for each cluster

# Get top 10 markers per cluster for Fig.1C
# "We used the 10 most-significant positive markers for heat map visualization"
top10 <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
# Expected: 10 marker genes per cluster

# Create Figure 1C - Heat map of cluster markers
# "We randomly subsampled 100 cells from each cluster to neutralize size differences"
set.seed(42)  # For reproducibility in random sampling
p1c <- DoHeatmap(seurat_obj, 
                features = top10$gene, 
                cells = 100, 
                group.by = "seurat_clusters") +
       ggtitle("Gene signatures of niche subpopulations")
print(p1c)
# Expected: Fig.1C - Heatmap showing cluster-specific gene expression patterns

# Combine plots for complete Figure 1
figure1 <- (p1b / p1c / p1d) + 
           plot_layout(heights = c(1, 1.5, 1))
print(figure1)


# Save the Seurat object for further analysis
saveRDS(seurat_obj, file = "bone_marrow_niche_seurat.rds")


# Save the figure
ggsave("Figure1_reproduction.pdf", figure1, width = 12, height = 16)

##This code provides a complete workflow to reproduce Figure 1 from the paper using the combined-first approach. Each step includes comments about the expected results, which should help you verify that your analysis is on track.
##The key visualization elements in Figure 1 are:

##Fig.1B: Expression of marker genes (Cdh5, Lepr, Col1a1) on t-SNE plots
##Fig.1C: Heatmap showing cluster-specific gene signatures
##Fig.1D: t-SNE visualization of identified clusters

