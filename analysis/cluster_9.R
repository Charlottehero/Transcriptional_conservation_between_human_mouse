immune_all = readRDS("processed_immune_final.rds")
unique(Idents(immune_all))

subsample = subset(immune_all, idents = c(9))
unique(subsample$final_annotation)

subsample$orig.ident = subsample$final_annotation
orig = subsample$orig.ident
orig
uni_orig = sort(unique(orig))
uni_orig
z = matrix(ncol= length(uni_orig), nrow = 1)
colnames(z) = uni_orig
for (j in 1:ncol(z)){
  z[1,j] = sum(orig == uni_orig[j]) / length(orig)
}
z
z = as.matrix(z)
lbls <- uni_orig
pie(z, labels = lbls, main="Pie Chart of Celltypes")
?pie

Idents(immune_all) = immune_all$species
Idents(immune_all)
human = subset(immune_all, idents = c("Human"))
Idents(human) = human$seurat_clusters
human = subset(human, idents = c(9))
human$orig.ident = human$final_annotation
orig = human$orig.ident
orig
uni_orig = sort(unique(orig))
uni_orig
z = matrix(ncol= length(uni_orig), nrow = 1)
colnames(z) = uni_orig
for (j in 1:ncol(z)){
  z[1,j] = sum(orig == uni_orig[j]) / length(orig)
}
z
z = as.matrix(z)

immune_all$orig.ident = immune_all$seurat_clusters
immune_all$final_annotation
DimPlot(immune_all, group.by = "final_annotation", shuffle = TRUE)

DimPlot(subsample, group.by = "species", shuffle = TRUE)

Idents(subsample) = subsample$species
# Find all markers
subsample.markers <- FindAllMarkers(subsample, only.pos = TRUE, min.pct = 0.25)
subsample.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(subsample.markers, file = "immune_species_c10.csv")

human.markers = read.csv("human_cluster_markers.csv", row.names = 1)
mouse.markers = read.csv("mouse_cluster_markers.csv", row.names = 1)

length(rownames(mouse.markers)[mouse.markers$cluster==9])
length(rownames(human.markers)[human.markers$cluster==9])
common = intersect(rownames(mouse.markers)[mouse.markers$cluster==9][1:20], rownames(human.markers)[human.markers$cluster==9][1:20])
length(common)

dim(immune_all)
