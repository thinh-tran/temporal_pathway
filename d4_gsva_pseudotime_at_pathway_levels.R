library(monocle)
load("/Volumes/Thinh/Bader_Lab/heartdev/Data/d4_gsva.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/Data/Seurat_objects/Ventricular_D4.RData")

d4_gsva <- as(d4_gsva, "matrix")

#### create CellDataSet object
pd <- colnames(d4_gsva)
pd <- unlist(lapply(pd, function(x) paste(strsplit(x, '_')[[1]][1], strsplit(x, '_')[[1]][2], sep = '_')))
clusters <- as.character(unlist(eb1S@meta.data$res.0.8))
pd_df <- data.frame(type=pd,seurat_clust=clusters)
rownames(pd_df) <- colnames(d4_gsva)
pda <- new('AnnotatedDataFrame', data = pd_df)

mes_d4 <- newCellDataSet(d4_gsva,
                         phenoData = pda,
                         expressionFamily=negbinomial.size())

mes_d4 <- estimateSizeFactors(mes_d4)
mes_d4 <- estimateDispersions(mes_d4)

fData(mes_d4)$gene_short_name <- sub("%.*", "", rownames(fData(mes_d4)))
# 

###### Clustering

#filter for genes with expression above 0.1
mes_d4_clust <- detectGenes(mes_d4, min_expr = 0.1) 

#only take genes that are expressed in at least 5% of the cells to order
fData(mes_d4_clust)$use_for_ordering <- fData(mes_d4_clust)$num_cells_expressed > 0.05 * ncol(mes_d4_clust)

#pca
plot_pc_variance_explained(mes_d4_clust, return_all = F) 

#tSNE
mes_d4_clust <- reduceDimension(mes_d4_clust,
                                max_components = 2,
                                norm_method = 'log',
                                num_dim = 12,
                                reduction_method = 'tSNE',
                                verbose = T) 

mes_d4_clust <- clusterCells(mes_d4_clust,
                             rho_threshold = 15,
                             delta_threshold = 16,
                             skip_rho_sigma = T,
                             verbose = F)

plot_cell_clusters(mes_d4_clust, color_by = 'as.factor(Cluster)')
plot_cell_clusters(mes_d4_clust, color_by = 'seurat_clust')

#### Pseudotime
#take genes thar are expressed in more than 10 cells
mes_expgene <-  row.names(subset(fData(mes_d4_clust),
                                 num_cells_expressed >= 10))

#differential expression by Monocle/seurat cluster
clustering_DEG_genes <- differentialGeneTest(mes_d4_clust[mes_expgene,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 1)

#take the top most 1000 DE genes to build a pseudotime trajectory
ordering_gene <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

mes_d4_clust <-
  setOrderingFilter(mes_d4_clust,
                    ordering_genes = ordering_gene)

#build DDRTree 
mes_d4_clust <-
  reduceDimension(mes_d4_clust, method = 'DDRTree')

#build pseudotime trajectory
mes_d4_clust <-
  orderCells(mes_d4_clust)

#plot trajectory
plot_cell_trajectory(mes_d4_clust, color_by = 'as.factor(Cluster)')
plot_cell_trajectory(mes_d4_clust, color_by = 'seurat_clust')
plot_cell_trajectory(mes_d4_clust, color_by = 'State')

### Pseudotime DEA

pseudotime_de <- differentialGeneTest(mes_d4_clust,
                                      fullModelFormulaStr = "~Pseudotime)",
                                      cores = 1)

pseudotime_de <- pseudotime_de[order(pseudotime_de$qval),]
pseudotime_gene <- rownames(pseudotime_de)[1:10]

plot_genes_in_pseudotime(mes_d4_clust[pseudotime_gene,])

de_state <- differentialGeneTest(mes_d4_clust,
                                      fullModelFormulaStr = "~State",
                                      cores = 1)
save(mes_d4, mes_d4_clust, pseudotime_de, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D4_GSVA_Monocle.RData")

#BEAM



