# Read in gene markers by species
human.markers = read.csv("human_cluster_markers.csv", row.names = 1)
mouse.markers = read.csv("mouse_cluster_markers.csv", row.names = 1)

# Find intersection
common = intersect(rownames(mouse.markers)[mouse.markers$cluster==9][1:100], rownames(human.markers)[human.markers$cluster==9][1:100])
length(common)
length(row.names())

# Start kuncheva analysis
kuncheva = function(k) {
  common = length(intersect(rownames(mouse.markers)[mouse.markers$cluster==9][1:k], rownames(human.markers)[human.markers$cluster==9][1:k]))
  result = (common - (k**2)/8135)/(k-(k**2)/8135)
  return(result)
}
results = c()
for (k in 1:50) {
  results = c(results, kuncheva(k))
}
plot(1:50, results, type = "b", frame = TRUE, pch = 19, 
     col = "red", xlab = "K", ylab = "Kuncheva Index", ylim = c(0, 1))
