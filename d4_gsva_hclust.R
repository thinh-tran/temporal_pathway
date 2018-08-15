rm(list=ls())
load("/Volumes/Thinh/Bader_Lab/heartdev/Data/d4_gsva_pseudotime.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D4_GSVA_Monocle.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/Data/Seurat_objects/Ventricular_D4.RData")
d4_hallmark <- d4_gsva_pseudotime[which(grepl("^HALLMARK", rownames(d4_gsva_pseudotime))==TRUE),]
rownames(d4_hallmark) <- sub("%.*", "", rownames(d4_hallmark))
rownames(d4_hallmark) <- sub("^HALLMARK_", "", rownames(d4_hallmark))
d4_hallmark <- as.matrix(d4_hallmark)

cell_order <- seq(1:ncol(d4_gsva_pseudotime))
colnames(d4_hallmark) <-cell_order


d4_hm_dist <- dist(t(d4_hallmark))
d4_hm_clust <- hclust(d4_hm_dist)
hm_subgroup <- cutree(d4_hm_clust, k=5)
plot(d4_hm_clust, xlab="", ylab="", xaxt='n', cex=0.55, main="D4 Ventricular Cells")
rect.hclust(d4_hm_clust, k=5)

d4_hallmark_m <- matrix(as.numeric(unlist(d4_hallmark)), nrow=nrow(d4_hallmark), 
                        dimnames=list(rownames(d4_hallmark), colnames(d4_hallmark)))

heatmap(d4_hallmark_m)

d4_hclust_group <- data.frame(cell=colnames(d4_gsva_pseudotime), order=names(hm_subgroup), hclust_cluster=hm_subgroup)

d4_seurat_clust <- as.data.frame(eb1S@meta.data[, which(colnames(eb1S@meta.data) == 'res.0.8')])
d4_seurat_clust$cell <- rownames(eb1S@meta.data)
d4_hclust_group <- merge(d4_hclust_group, d4_seurat_clust, by="cell")
colnames(d4_hclust_group)[[4]] <- "seurat_cluster"

state <- as.data.frame(mes_d4_clust@phenoData@data$State)
state$cell <- rownames(mes_d4_clust@phenoData@data)
d4_hclust_group <- merge(d4_hclust_group, state, by="cell")
colnames(d4_hclust_group)[[4]] <- "state"

d4_hallmark_m_state <- rbind(d4_hallmark_m, hm_subgroup)
pearson_test <- apply(d4_hallmark_m, 1, function(x) cor.test(x, as.numeric(unlist(d4_hclust_group$state)), method="p"))
pearson_pvals <- as.data.frame(unlist(lapply(pearson_test, function(x) x$p.value)))
