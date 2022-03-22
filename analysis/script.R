library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(harmony)
library(ggplot2)

Convert("/Users/charlotte/Desktop/Jessica 254/254_original.h5ad", dest = "h5seurat", overwrite = FALSE)

immune_all <- LoadH5Seurat("254_original.h5seurat")
immune_all

immune_all$nCount_RNA = colSums(x = immune_all, slot = "counts")  # nCount_RNA
immune_all$nFeature_RNA = colSums(x = GetAssayData(object = immune_all, slot = "counts") > 0)  # nFeatureRNA

#remove mitochondrial
immune_all[["percent.mt"]] <- PercentageFeatureSet(immune_all, pattern = "^MT-")
VlnPlot(immune_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(immune_all, feature1 = "nCount_RNA", feature2 = "percent.mt")
?geom_vline
plot2 <- FeatureScatter(immune_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") & geom_vline(xintercept = 200, linetype = "dashed") & geom_vline(xintercept = 5000, linetype = "dashed") & geom_hline(yintercept = 200, linetype = "dashed")
plot2
# removing the cells with low or extreme high gene counts
immune_all <- subset(immune_all, subset = nCount_RNA > 200 & nCount_RNA < 5000 & nFeature_RNA > 200)
plot2 <- FeatureScatter(immune_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") & geom_hline(yintercept = 200) & geom_hline(yintercept = 5000)
plot2
dim(immune_all)
# start normalization
immune_all = NormalizeData(immune_all, verbose = FALSE)
immune_all = FindVariableFeatures(immune_all, selection.method = "vst", nfeatures = 2000)
?FindVariableFeatures
immune_all = ScaleData(immune_all, verbose = FALSE)
# run PCA to denoise
immune_all = RunPCA(immune_all, verbose = TRUE) # default is 50
# determine the dimensionality of the dataset
ElbowPlot(immune_all)
# run harmony to remove donor batch effect
unique(immune_all$batch)
immune_all = RunHarmony(immune_all, group.by.vars = "batch", dims.use = 1:30)
ElbowPlot(immune_all)

immune_all = RunUMAP(immune_all, reduction = "harmony", dims = 1:10)
immune_all = FindNeighbors(immune_all, reduction = "harmony", dims = 1:10)
immune_all = FindClusters(immune_all, resolution = 0.6)
immune_all = identity(immune_all)

# save the sample
saveRDS(immune_all, file = "processed_immune_final.rds")
DimPlot(immune_all, reduction = "umap", label=TRUE)
DimPlot(immune_all, reduction = "umap", group.by="species", label=TRUE, shuffle = TRUE)
DimPlot(immune_all, reduction = "umap", group.by="final_annotation", shuffle = TRUE)

# Find all markers
immune_all.markers <- FindAllMarkers(immune_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune_all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
write.csv(immune_all.markers, file = "immune_species.csv")

#immune_all = readRDS("processed_immune.rds")

# plot stacked barplot for annotating the clusters
immune_all$orig.ident = immune_all$species
orig = immune_all$orig.ident
orig
u_orig = sort(unique(orig))
u_orig
cluster = as.character(Idents(immune_all))
uct = sort(unique(cluster))
z = matrix(nrow=length(uct), ncol=length(u_orig))
rownames(z) = uct
colnames(z) = u_orig
z
for (i in 1:nrow(z)){
  for (j in 1:ncol(z)){
    z[i,j] = sum(orig[cluster==uct[i]] == u_orig[j])
  }
}
# Make a plot about num cells for each cluster by orig
z = as.matrix(z)
z1 = data.frame(ncells=as.vector(z),
                type=rep(as.integer(rownames(z)), ncol(z)),
                orig=rep(colnames(z), each=nrow(z)))
ggplot(z1, aes(x=type, y=ncells, fill=orig, group=orig)) +
  geom_bar(stat="identity", width=0.7) +
  theme_classic() +
  labs(x="Cluster Number", y="Numer of Cells", fill=NULL)+
  scale_x_continuous(breaks=0:max(z1$type))

?top_n

