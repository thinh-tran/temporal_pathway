list.of.packages <- c('reshape','reshape2','pheatmap','tsne','igraph','ggplot2','mclust','grid','Rtsne','cccd', "Vennerable", "dplyr", "GO.db", "org.Hs.eg.db")
lapply(list.of.packages, require, character.only = TRUE)


source("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/cellrouter-master/CellRouter_Class.R")
libdir <- "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/cellrouter-master/CellRouter/"
maindir <- "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/cell_router_d20/"
load(paste0(maindir, "cellrouter_object.RData"))

load("/Volumes/Thinh/Bader_Lab/heartdev/Data/Seurat_objects/Ventricular_D20_allcells.RData")
exprMatrix <- as(eb1S@data, "matrix")
exprMatrix <- as.data.frame(exprMatrix)

cellrouter <- CellRouter(exprMatrix, min.cells=0, min.genes=0)

cellrouter <- addInfo(cellrouter, metadata = eb1S@meta.data, colname = "seurat_cluster", metadata.column = 'res.0.6')
cellrouter <- scaleData(cellrouter)
cellrouter <- computePCA(cellrouter, num.pcs = 50, seed=42)
plot(cellrouter@pca$sdev, xlab='PC', ylab='Standard deviation of PC')
cellrouter <- computeTSNE(cellrouter, num.pcs = 20, seed=42, max_iter = 500)
cellrouter <- customSpace(cellrouter, cellrouter@tsne$cell.embeddings)
plotReducedDimension(cellrouter, reduction.type = 'tsne', dims.use = c(1,2), annotation = "seurat_cluster", annotation.color = 'seurat_cluster_color', showlabels = TRUE, 4.5, 3.5, filename="tnse.pdf")

cellrouter <- buildKNN(cellrouter, k = 10, column.ann = 'seurat_cluster', num.pcs = 20, sim.type = 'jaccard')
filename <- paste0(maindir, "cell_edge_weighted_network.txt")
write.table(cellrouter@graph$edges, file=filename, sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE)

sources <- c('4')
targets <- setdiff(as.vector(cellrouter@sampTab$seurat_cluster), sources)
methods <- c("euclidean", "maximum", "manhattan","canberra","binary", 'graph')
cellrouter <- findPaths(cellrouter, column='seurat_cluster', libdir, "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/cell_router_d20", method="euclidean")

library(igraph)
#Preprocess trajectories
ranks <- c('path_cost', 'path_flow', 'rank', 'length')
cellrouter <- processTrajectories(cellrouter, rownames(cellrouter@ndata), path.rank=ranks[3], num.cells = 3, neighs=3, column.ann = 'seurat_cluster', column.color = 'seurat_cluster_color')

names <- unique(names(cellrouter@pathsinfo$distr))
##Identify genes regulated along each trajectory (up or down-regulated)
cellrouter <- correlationPseudotime(cellrouter, type='spearman')
cellrouter <- topGenes(cellrouter, 0.8, 0.1)

cellrouter <- smoothDynamics(cellrouter, names(cellrouter@pathsinfo$distr))
cellrouter <- clusterGenesPseudotime(cellrouter, 5)

plotReducedDimension(cellrouter, reduction.type = 'tsne', annotation="seurat_cluster", annotation.color = 'seurat_cluster_color',showlabels = TRUE, width = 4.5, height = 3.5, filename=paste0(maindir, 'tSNE_graphClustering_clusters.pdf'))


grn.data <- buildGRN(cellrouter, species = 'Hs', genes.use = rownames(cellrouter@ndata), zscore = 5, filename = paste0(maindir, 'GRN.R'))
plotClusterHeatmap(cellrouter, names(cellrouter@pathsinfo$distr), 10, 10, 2, paste0(maindir, 'dynamics.pdf'))


transitions <- sort(names(cellrouter@pathsinfo$distr)) #transitions to be analyzed
grn.scores <- grnscores(cellrouter, ggrn = grn.data$GRN, tfs = grn.data$tfs, transitions = transitions[9], direction='both', flip=TRUE, q.up=0.80, q.down=0.0,
                        dir.targets='up', columns=3, width=5, height=7, filename=paste0(maindir, 'dynamics_transitions_', paste0(transitions[9], collapse=", "), '.pdf'))

library('geomnet')
p <- '4.9'
plottr(cellrouter, p, grn.scores[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste0(maindir, p, 'Regulators_transition.pdf',sep=''))


#save(cellrouter, file=paste0(maindir, "cellrouter_object.RData"))



