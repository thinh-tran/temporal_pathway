load("/Volumes/Thinh/Bader_Lab/heartdev/Data/Seurat_objects/Ventricular_D4.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/Data/d4_gsva.RData")

exprMatrix <- eb1S@data
rm(pDat, uniqDE, minDEgenes, neighbDE, maxPCt)

d4_pdata <- d4_pdata[order(d4_pdata$Pseudotime), ]
pseudotime <- subset(d4_pdata, select = "Pseudotime")
pseudotime <- t(pseudotime)
pseudotime <- as.data.frame(pseudotime)


d4_pseudotime <- rbind(exprMatrix, pseudotime)
d4_pseudotime_ordered <- d4_pseudotime[, order(d4_pseudotime[18510,])]
colnames(d4_pseudotime_ordered) <- seq(1:ncol(d4_pseudotime_ordered))

d4_clust <- hclust(dist(d4_pseudotime_ordered[-c(18510), ]))

d4_gsva_cor <- cor(d4_gsva)
d4_gene_cor <- cor(as(exprMatrix, "matrix"))

library(reshape2)

d4_gsva_cor_single <- melt(d4_gsva_cor)
d4_gene_cor_single <- melt(d4_gene_cor)

cor_test <- cortest.jennrich(d4_gsva_cor, d4_gene_cor, n1=1862, n2=1862)

plot(x=d4_gene_cor_single$value, y=d4_gsva_cor_single$value, xlab="D4_Gene_Exp", ylab= "D4_GSVA", main="GSVA and gene expression correlation", xlim=c(0, 1), ylim=c(-.5, 1))

