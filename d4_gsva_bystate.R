library(GSVA)
library(tm)
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D20_GSVA_Monocle_pathways.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/Data/Seurat_objects/Ventricular_D20_allcells.RData")
rm(mes_d4, pseudotime_de, minDEgenes, neighbDE, pDat, uniqDE, maxPCt)

exprMatrix <- eb1S@data
state <- as.data.frame(mes_d4_clust@phenoData@data$State)
rownames(state) <- rownames(mes_d4_clust@phenoData@data)
state$exprMatrix <- colnames(exprMatrix)

exprMatrix_bystate <- list()
for (i in levels(state$`mes_d4_clust@phenoData@data$State`)){
  exprMatrix_bystate[[i]] <- as.matrix(exprMatrix[, which(state$`mes_d4_clust@phenoData@data$State` == i)])
}

if (length(which(lapply(exprMatrix_bystate, ncol) == 1)) > 0){
  exprMatrix_bystate_f <- exprMatrix_bystate[-which(lapply(exprMatrix_bystate, ncol) == 1)]
} else { exprMatrix_bystate_f <- exprMatrix_bystate }


geneexp_bystate <- lapply(exprMatrix_bystate_f, rowSums)
geneexp_state_m <- Reduce(cbind, geneexp_bystate)
colnames(geneexp_state_m) <- names(geneexp_bystate)
geneexp_state_m <- cbind(geneexp_state_m, unlist(exprMatrix_bystate[which(lapply(exprMatrix_bystate, ncol) == 1)]))
#colnames(geneexp_state_m)[11] <- names(exprMatrix_bystate)[which(lapply(exprMatrix_bystate, ncol) == 1)]

#d4_exprMatrix_bystate <- geneexp_state_m[, c(1, 11, seq(2, 10))]

d20_exprMatrix_bystate <- geneexp_state_m

save(d20_exprMatrix_bystate, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D20_exprMatrix_bystate.RData")


library(GSVA)
library(GSEABase)
library(wordcloud)
library(tm)

load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D20_exprMatrix_bystate.RData")
pathways <- getGmt("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Human_GO_AllPathways_no_GO_iea_July_01_2018_symbol.gmt")

d20_gsva_bystate <- gsva(d20_exprMatrix_bystate, pathways, method="gsva", min.sz=5, max.sz=350, parallel.sz=detectCores()-2)

save(d20_gsva_bystate, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D20_gsva_bystate.RData")


load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D4_gsva_bystate.RData")
rownames(d20_gsva_bystate) <- sub("%.*", "", rownames(d20_gsva_bystate))

upreg <- apply(d20_gsva_bystate, 2, function(x) x[which(x > 0)])
downreg <- apply(d20_gsva_bystate, 2, function(x) x[which(x < 0)])

upreg_signaling <- lapply(upreg, function(x) x[which(grepl("SIGNALING", names(x))==TRUE)])
downreg_signaling <- lapply(downreg, function(x) x[which(grepl("SIGNALING", names(x))==TRUE)])

source("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/gsva_annotation_functions.R")

upreg_bystate <- FindRepPathways(lapply(upreg, function(x) names(x)), plot=T, col="Dark2", min_freq = 4, max_words = 10)
