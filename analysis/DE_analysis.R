# Read in the Seurat Object
immune_all = readRDS("processed_immune_final.rds")
Idents(immune_all)
Idents(immune_all) = immune_all$species
human = subset(immune_all, idents = c("Human"))
mouse = subset(immune_all, idents = c("Mouse"))
Idents(immune_all) = immune_all$seurat_clusters
Idents(human) = human$seurat_clusters
Idents(mouse) = mouse$seurat_clusters
# Find all markers
human.markers <- FindAllMarkers(human, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
human.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(human.markers, file = "human_cluster_markers.csv")

# Find all markers
mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mouse.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(mouse.markers, file = "mouse_cluster_markers.csv")

# checking cluster 9 
Idents(mouse) = mouse$seurat_clusters
sub_human = subset(human, idents = c(9))
sub_mouse = subset(mouse, idents = c(9))
subsample = subset(immune_all, idents = c(9))
Idents(subsample) = subsample$final_annotation
B_cells_9 = subset(subsample, idents = c("B cells"))

# start counting the percentage of species for each cell type
orig = B_cells_9$species
orig
uni_orig = sort(unique(orig))
uni_orig
z = matrix(ncol= length(uni_orig), nrow = 1)
colnames(z) = uni_orig
for (j in 1:ncol(z)){
  z[1,j] = sum(orig == uni_orig[j]) / length(orig)
}
z

# checking cluster 0
Idents(human) = human$seurat_clusters
Idents(mouse) = mouse$seurat_clusters
sub_human = subset(human, idents = c(0))
sub_mouse = subset(mouse, idents = c(0))
subsample = subset(immune_all, idents = c(0))
# pring out the number of samples in cluster 0
length(colnames(subsample))
Idents(subsample) = subsample$final_annotation

monocytes = subset(subsample, idents = c("Monocytes"))
neutrophils = subset(subsample, idents = c("Neutrophils"))

orig = monocytes$species
orig
uni_orig = sort(unique(orig))
uni_orig
z = matrix(ncol= length(uni_orig), nrow = 1)
colnames(z) = uni_orig
for (j in 1:ncol(z)){
  z[1,j] = sum(orig == uni_orig[j]) / length(orig)
}
z

orig = neutrophils$species
orig
uni_orig = sort(unique(orig))
uni_orig
z = matrix(ncol= length(uni_orig), nrow = 1)
colnames(z) = uni_orig
for (j in 1:ncol(z)){
  z[1,j] = sum(orig == uni_orig[j]) / length(orig)
}
z

# print out number of cells in 16
subsample = subset(immune_all, idents = c(16))
length(colnames(subsample))
