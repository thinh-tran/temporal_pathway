
### GSVA by cluster ###

library(GSEABase)
library(GSVA)
library(iprior)
library(shape)
library(reshape2)
library(scales)


load("/Volumes/Thinh/Bader_Lab/heartdev/Data/Seurat_objects/Ventricular_D20_CTNT.RData")
exprMatrix <- as(eb1S@data, "matrix")
tSNE <- eb1S@dr$tsne@cell.embeddings

rm(minDEgenes, neighbDE, pDat, uniqDE, maxPCt)
idents <- eb1S@meta.data$res.0.6
names(idents) <- rownames(eb1S@meta.data)

exprMatrix_bycluster <- list()
for (i in levels(idents)){
  exprMatrix_bycluster[[i]] <- rowMeans(exprMatrix[, which(colnames(exprMatrix) %in% names(idents)[which(idents == i)])])
}

exprMatrix_bycluster <- do.call(cbind, exprMatrix_bycluster)
#save(exprMatrix_bycluster, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/d4_exprMatrix_bycluster.RData")

pathways <- getGmt("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Human_GO_AllPathways_no_GO_iea_July_01_2018_symbol.gmt")
gsva_bycluster <- gsva(exprMatrix_bycluster, pathways, method="gsva", min.sz=5, max.sz=350, parallel.sz=detectCores()-2)
#gsva_bycluster <- as.data.frame(gsva_bycluster)


### Find paths between clusters ###
#pathway_var <- apply(gsva_bycluster, 1, var)
#hivar_pathways <- which(pathway_var > median(pathway_var)+2*mad(pathway_var, center = median(pathway_var)))

temp_cor <- list()
for ( i in 1:ncol(gsva_bycluster)){
  temp_matrix <- gsva_bycluster[, setdiff(seq(1:ncol(gsva_bycluster)), i)]
  temp_cor[[i]] <- apply(temp_matrix, 2, function(x) cor(gsva_bycluster[,i], x))
}

# temp_cor_gene <- list()
# for ( i in 1:ncol(exprMatrix_bycluster)){
#   temp_matrix <- exprMatrix_bycluster[hivar_pathways, setdiff(seq(1:ncol(exprMatrix_bycluster)), i)]
#   temp_cor_gene[[i]] <- apply(temp_matrix, 2, function(x) cor(exprMatrix_bycluster[hivar_pathways,i], x))
# }

scores <- trajectories <- list()
for (i in 1:length(temp_cor)){
  temp <- temp_cor
  avail <- seq(1:length(temp))
  total_score <- 0
  startpt <- i
  trajectory <- c(i)
  while (length(avail) > 0) {
    #print(startpt)
    target <- names(temp[[startpt]])[which(temp[[startpt]]==max(temp[[startpt]]))]
    total_score <- total_score + max(temp_cor[[startpt]])
    trajectory <- c(trajectory, target)
    avail <- avail[-which(avail %in% trajectory)]
    #print(paste0("available: ", paste(avail, collapse=", ")))
    temp <- lapply(temp_cor, function(x) x[-which(names(x) %in% trajectory)])
    startpt <- as.numeric(target)
  }
  scores[[i]] <- total_score
  trajectories[[i]] <- trajectory
}

names(scores) <- lapply(trajectories, function(x) paste(x, collapse="-"))

#chosen_trajectory <- trajectories[which(unlist(scores)==max(unlist(scores)))][[1]]

### finding strong edges ###
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/d20_gsva_cell_cell_correlation.RData")

# load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D20_gsva.RData")
# d20_gsva_cor <- cor(gsva_d20)
# d20_gsva_cor_norep <- replace(d20_gsva_cor, lower.tri(d20_gsva_cor, TRUE), NA)
# 
# for (i in 1:nrow(d20_gsva_cor_norep)){
#   if (length(d20_gsva_cor_norep[i,][which(!is.na(d20_gsva_cor_norep[i,]))]) > 1) {
#     cutoff <- max(d20_gsva_cor_norep[i,], na.rm=T)-2*mad(d20_gsva_cor_norep[i,], na.rm=T)
#     d20_gsva_cor_norep[i,][which(d20_gsva_cor_norep[i,] < cutoff)] <- NA
#   }
# }
# 
# d20_gsva_cor_strong <- melt(d20_gsva_cor_norep, na.rm=T)
# 
# d20_gsva_cor_strong_Cell1_clust <- sapply(d20_gsva_cor_strong$Var1, function(x) as.numeric(idents[which(names(idents)==x)]))
# d20_gsva_cor_strong_Cell2_clust <- sapply(d20_gsva_cor_strong$Var2, function(x) as.numeric(idents[which(names(idents)==x)]))
# 
# d20_gsva_cor_strong$Cell1_clust <- d20_gsva_cor_strong_Cell1_clust
# d20_gsva_cor_strong$Cell2_clust <- d20_gsva_cor_strong_Cell2_clust 
# 
# d20_gsva_cor <- melt(replace(d20_gsva_cor, lower.tri(d20_gsva_cor, TRUE), NA), na.rm = TRUE)
# 
# d20_gsva_cor_Cell1_clust <- sapply(d20_gsva_cor$Var1, function(x) as.numeric(idents[which(names(idents)==x)]))
# d20_gsva_cor_Cell2_clust <- sapply(d20_gsva_cor$Var2, function(x) as.numeric(idents[which(names(idents)==x)]))
# # d20_gsva_cor$Cell1_clust <- d20_gsva_cor_Cell1_clust
# # d20_gsva_cor$Cell2_clust <- d20_gsva_cor_Cell2_clust

combination <- combn(as.numeric(levels(idents)), 2)
combination <- lapply(seq_len(ncol(combination)), function(i) combination[,i])

edge_counts <- all_edge_counts <- strong_edge_percent <- list()
for (i in 1:length(combination)){
  edge_counts[[i]] <- length(which(d20_gsva_cor_strong_Cell1_clust == combination[[i]][1] & d20_gsva_cor_strong_Cell2_clust == combination[[i]][2]))
  all_edge_counts[[i]] <- length(which(d20_gsva_cor_Cell1_clust == combination[[i]][1] & d20_gsva_cor_Cell2_clust == combination[[i]][2]))
  strong_edge_percent[[i]] <- edge_counts[[i]]/all_edge_counts[[i]]
}

names(strong_edge_percent) <- lapply(combination, function(x) paste0(x, collapse="-"))

strong_edges_df <- data.frame(Edge=names(strong_edge_percent), Strength=unlist(strong_edge_percent))

duplicated_edges <- list()
for (i in 1:nrow(strong_edges_df)){
  duplicated_edges[[i]] <- strong_edges_df$Strength[[i]]
  names(duplicated_edges)[[i]] <- paste(rev(unlist(strsplit(as.character(strong_edges_df$Edge[[i]]),split="-"))),collapse="-")
}

duplicated_edges_df <- data.frame(Edge=names(duplicated_edges), Strength=unlist(duplicated_edges))
strong_edges_df <- rbind(strong_edges_df, duplicated_edges_df)

#save(d20_gsva_cor_Cell1_clust, d20_gsva_cor_Cell2_clust, d20_gsva_cor_strong_Cell1_clust, d20_gsva_cor_strong_Cell2_clust, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/d20_gsva_cell_cel_correlation.RData")


### pruning ###
edge_pairs <- lapply(trajectories, function(x) apply(cbind(x[-length(x)], x[-1]), 1, function(y) paste0(y, collapse="-")))
edge_strength_score <- lapply(edge_pairs, function(x) sum(strong_edges_df$Strength[which(rownames(strong_edges_df) %in% x)]))

# chosen_trajectory <- unlist(trajectories[which(unlist(edge_strength_score)==max(unlist(edge_strength_score[which(unlist(scores)==max(unlist(scores)))])))])
# chosen_trajectory_index <- which(unlist(edge_strength_score)==max(unlist(edge_strength_score[which(unlist(scores)==max(unlist(scores)))])))
chosen_trajectory <- trajectories[[1]]
chosen_trajectory_index <- 1
chosen_trajectory_edge_strength <- strong_edges_df$Strength[which(rownames(strong_edges_df) %in% edge_pairs[[chosen_trajectory_index]])]
chosen_trajectory_edge_type <- ifelse(chosen_trajectory_edge_strength < quantile(strong_edges_df$Strength, 0.6), 2, 1)
chosen_trajectory_edge_strength <- rescale(chosen_trajectory_edge_strength, to=c(2,5))
chosen_trajectory_edge_pairs <- edge_pairs[[chosen_trajectory_index]]

### Plot ###
cluster_coord <- list()
for (i in levels(idents)){
  cluster_coord[[i]] <- apply(tSNE[which(rownames(tSNE) %in% names(idents)[which(idents == i)]), ], 2, mean)
}
cluster_coord <- Reduce(rbind, cluster_coord)
rownames(cluster_coord) <- levels(idents)

from <- to <- list()
for (i in 1:(length(chosen_trajectory)-1)){
  from[[i]] <- cluster_coord[as.numeric(chosen_trajectory[[i]]), ] 
  to[[i]] <- cluster_coord[as.numeric(chosen_trajectory[[i+1]]), ] 
}
from <- Reduce(rbind, from)
to <- Reduce(rbind, to)

if (length(levels(eb1S@ident)) <= 8) {
  clustCols <- brewer.pal(length(levels(eb1S@ident)),"Dark2")
} else {
  clustCols <- gg_colour_hue(length(levels(eb1S@ident)))
}

plot(eb1S@dr$tsne@cell.embeddings,pch=20,xlab="tSNE_1",ylab="tSNE_2",
     col=alpha(clustCols[eb1S@ident],.7),
     bg=alpha(clustCols[eb1S@ident],.3), main="D20 ventricular cardiomyocytes")

arrow_col <- rgb(0,0,0, max=255, alpha=125)

for (i in 1:(nrow(cluster_coord)-1)){
  par(new=T)
  Arrows(from[i, 1], from[i, 2], to[i, 1], to[i, 2], lwd=chosen_trajectory_edge_strength[[i]], lty = chosen_trajectory_edge_type[[i]], col=arrow_col)
}

for (i in 1:nrow(cluster_coord)){
  par(new=T)
  text(cluster_coord[i, 1], cluster_coord[i, 2]+2, labels=as.character(i), col="black", cex=1.2, font=2)
}


### Find differentially expressed pathways in transitions ###

wilcox_pvals <- difference <- list()
for (i in 1:length(chosen_trajectory_edge_pairs)){
  edges <- as.numeric(strsplit(chosen_trajectory_edge_pairs[[i]], split="-")[[1]])
  temp_gsva_1 <- gsva_d20[, which(colnames(gsva_d20) %in% names(idents)[which(idents == edges[[1]])])]
  temp_gsva_2 <- gsva_d20[, which(colnames(gsva_d20) %in% names(idents)[which(idents == edges[[2]])])]
  temp_wilcox_pvals <-  list()
  for (j in 1:nrow(gsva_d20)){
    temp_wilcox_pvals[[j]] <- wilcox.test(temp_gsva_1[j, ], temp_gsva_2[j, ])
  }
  wilcox_pvals[[i]] <- unlist(lapply(temp_wilcox_pvals, function(x) x$p.value))
  difference[[i]] <- rowMeans(temp_gsva_1)-rowMeans(temp_gsva_2)
  names(wilcox_pvals[[i]]) <- names(difference[[i]]) <- rownames(gsva_d20)
}

upreg_first <- upreg_later <- list()
for (i in 1:length(wilcox_pvals)){
  result_df <- data.frame(Pathway=names(wilcox_pvals[[i]]), P_value=wilcox_pvals[[i]], Difference=difference[[i]])
  upreg_first[[i]] <- result_df[which(result_df$P_value < 0.05 & result_df$Difference > 0),]
  upreg_first[[i]] <- upreg_first[[i]][order(upreg_first[[i]]$Difference, decreasing = T), ]
  upreg_later[[i]] <- result_df[which(result_df$P_value < 0.05 & result_df$Difference < 0),]
  upreg_later[[i]] <- upreg_later[[i]][order(upreg_later[[i]]$Difference, decreasing = F), ]
}

names(upreg_first) <- names(upreg_later) <- chosen_trajectory_edge_pairs

save(upreg_first, upreg_later, file="/Volumes/Thinh/Bader_Lab/heartdev/Data/Seurat_objects/d20_seurat_trajectory_de_pathways.RData")

upreg_first_viz <- lapply(upreg_first, function(x) inner_join(x,  pathway_df, by=c("Pathway_name"="Pathway")))
upreg_later_viz <- lapply(upreg_later, function(x) inner_join(x, pathway_df, by=c("Pathway"="Pathway_name")))

upreg_first_viz <- lapply(upreg_first_viz, function(x) x[, c(1, 4, 2)])
upreg_later_viz <- lapply(upreg_later_viz, function(x) x[, c(1, 4, 2)])

filepath <- "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/d20_seurat_trajectory/"

# for (i in 1:length(upreg_first_viz)){
#   write.table(upreg_first_viz[[i]], file=paste0(filepath, i, "_", names(upreg_first_viz)[[i]], "_upreg_first.txt"), quote=F, row.names = F, sep="\t")
#   write.table(upreg_later_viz[[i]], file=paste0(filepath, i, "_", names(upreg_later_viz)[[i]], "_upreg_later.txt"), quote=F, row.names = F, sep="\t")
# }

load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/d20_seurat_trajectory_de_pathways.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/d20_seurat_trajectory/em_themes.RData")

source("/Volumes/Thinh/Bader_Lab/summarize_em_themes.R")

library(RCy3)
library(httr)
library(RJSONIO)
library(purrr)

port.number = 1234
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")

pvalue_threshold <- "0.05"

similarity_threshold <- "0.25"
similarity_metric = "COMBINED"

generic_gmt_file <- "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Human_GO_AllPathways_no_GO_iea_July_01_2018_symbol.gmt"
# 
# network_suid_first <- summary_first <- list()
# for (i in 3:length(upreg_first)){
#   print(i)
  results_filename <- paste0(filepath, i, "_", names(upreg_first)[[i]], "_upreg_first.txt")
  
  em_command = paste('enrichmentmap build analysisType="generic" gmtFile=',generic_gmt_file,
                     'pvalue=',pvalue_threshold,
                     'similaritycutoff=',similarity_threshold,'enrichmentsDataset1=',results_filename,
                     sep=" ")
  
  #current_network_name <- paste0("PC_", i)
  network_suid_first[[i]] <- commandsRun(em_command)
  summary_first[[i]] <- summarize_em_themes(as.numeric(network_suid_first[[i]]))
  summary_first[[i]] <- summary_first[[i]][order(summary_first[[i]]$number_of_nodes, decreasing=T), ]
} 

network_suid_later <- summary_later <- list()
for (i in 1:length(upreg_later)){
  print(i)
  results_filename <- paste0(filepath, i, "_", names(upreg_later)[[i]], "_upreg_later.txt")
  
  em_command = paste('enrichmentmap build analysisType="generic" gmtFile=',generic_gmt_file,
                     'pvalue=',pvalue_threshold,
                     'similaritycutoff=',similarity_threshold,'enrichmentsDataset1=',results_filename,
                     sep=" ")
  
  #current_network_name <- paste0("PC_", i)
  network_suid_later <- commandsRun(em_command)
  summary_later[[i]] <- summarize_em_themes(as.numeric(network_suid_later))
  summary_later[[i]] <- summary_later[[i]][order(summary_later[[i]]$number_of_nodes, decreasing=T), ]
} 

names(summary_first) <- names(summary_later) <- names(upreg_first)
save(summary_first, summary_later, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/d20_seurat_trajectory/em_themes.RData")



### Find changing themes at each point ###

load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/d20_seurat_trajectory/em_themes.RData")

gsva_cluster_signaling <- gsva_bycluster[which(grepl("signaling", rownames(gsva_bycluster), ignore.case = T)==T), ]
signaling_cor <- cor(t(gsva_cluster_signaling))
signaling_cor[!lower.tri(signaling_cor)] <- 0
signaling_cor_f <- signaling_cor[, apply(signaling_cor,2,function(x) all(x<=0.85))]

pathway_cor <- cor(t(gsva_bycluster))
pathway_cor[!lower.tri(pathway_cor)] <- 0
pathway_cor_f <- pathway_cor[, apply(pathway_cor,2,function(x) all(x<=0.85))]

gsva_bycluster_f <- gsva_bycluster[which(rownames(gsva_bycluster) %in% colnames(signaling_cor_f)), ]
gsva_var_f <- sort(apply(gsva_bycluster_f, 1, var), decreasing=T)[1:50]
gsva_tbg <- gsva_bycluster[which(rownames(gsva_bycluster) %in% names(gsva_var_f)), as.numeric(chosen_trajectory)]
rownames(gsva_tbg) <- sub("%.*", "", rownames(gsva_tbg))

heatmap.2(gsva_tbg, dendrogram = "row",Colv = F, main="Signaling pathways in D20 ventricular data", margins=c(5, 25), cexRow=0.75)

characteristic_pathway <- list()
for (i in 1:(length(upreg_first)-1)){
  characteristic_pathway[[i]] <- Reduce(intersect, list(upreg_first[[i+1]]$Pathway, upreg_later[[i]]$Pathway, rownames(gsva_bycluster_f)))
}

for (i in 1:length(chosen_trajectory_edge_pairs)){
  edges <- as.numeric(strsplit(chosen_trajectory_edge_pairs[[i]], split="-")[[1]])
  gf <- union(intersect(upreg_first[[i]]$Pathway, rownames(gsva_bycluster_f)), intersect(upreg_later[[i]]$Pathway, rownames(gsva_bycluster_f)))
  temp_gsva_bycluster <- gsva_bycluster[which(rownames(gsva_bycluster) %in% gf), edges]
  temp_gsva_var_f <- sort(apply(temp_gsva_bycluster, 1, var), decreasing=T)[1:50]
  temp_gsva_tbg <- temp_gsva_bycluster[which(rownames(temp_gsva_bycluster) %in% names(temp_gsva_var_f)), ]
  rownames(temp_gsva_tbg) <- sub("%.*", "", rownames(temp_gsva_tbg))
  heatmap.2(temp_gsva_tbg, dendrogram = "row", main=paste0("Changing signaling pathways between clusters ", chosen_trajectory_edge_pairs[[i]]), margins=c(5, 20), cexRow=0.8)
}


### Find growth factors that are differentially expressed

growth_factor <- read.table("/Volumes/Thinh/Bader_Lab/growth_factor_list.txt")
gf_index <- which(rownames(exprMatrix) %in% growth_factor$V1)

wilcox_pvals_gf <- difference_gf <- list()
for (i in 1:length(chosen_trajectory_edge_pairs)){
  edges <- as.numeric(strsplit(chosen_trajectory_edge_pairs[[i]], split="-")[[1]])
  temp_exprMatrix_1 <- exprMatrix[gf_index, which(colnames(exprMatrix) %in% names(idents)[which(idents == edges[[1]])])]
  temp_exprMatrix_2 <- exprMatrix[gf_index, which(colnames(exprMatrix) %in% names(idents)[which(idents == edges[[2]])])]
  temp_wilcox_pvals <-  list()
  for (j in 1:nrow(exprMatrix[gf_index, ])){
    temp_wilcox_pvals[[j]] <- wilcox.test(temp_exprMatrix_1[j, ], temp_exprMatrix_2[j, ])
  }
  wilcox_pvals_gf[[i]] <- unlist(lapply(temp_wilcox_pvals, function(x) x$p.value))
  difference_gf[[i]] <- rowMeans(temp_exprMatrix_1)-rowMeans(temp_exprMatrix_2)
  names(wilcox_pvals_gf[[i]]) <- names(difference_gf[[i]]) <- rownames(exprMatrix[gf_index, ])
}

upreg_first_gf <- upreg_later_gf <- list()
for (i in 1:length(wilcox_pvals_gf)){
  result_df <- data.frame(Pathway=names(wilcox_pvals_gf[[i]]), P_value=wilcox_pvals_gf[[i]], Difference=difference_gf[[i]])
  upreg_first_gf[[i]] <- result_df[which(result_df$P_value < 0.05 & result_df$Difference > 0),]
  upreg_first_gf[[i]] <- upreg_first_gf[[i]][order(upreg_first_gf[[i]]$Difference, decreasing = T), ]
  upreg_later_gf[[i]] <- result_df[which(result_df$P_value < 0.05 & result_df$Difference < 0),]
  upreg_later_gf[[i]] <- upreg_later_gf[[i]][order(upreg_later_gf[[i]]$Difference, decreasing = F), ]
}

names(upreg_first_gf) <- names(upreg_later_gf) <- chosen_trajectory_edge_pairs

characteristic_gf <- list()
for (i in 1:(length(upreg_first_gf)-1)){
  characteristic_gf[[i]] <- intersect(upreg_first_gf[[i+1]]$Pathway, upreg_later_gf[[i]]$Pathway)
}

exprMatrix_bycluster_df <- exprMatrix_bycluster[gf_index, ]
gf_var <- apply(exprMatrix_bycluster_df, 1, var)




library("gplots")

heatmap.2(exprMatrix_bycluster_df[which(gf_var > quantile(gf_var, 0.7)), as.numeric(chosen_trajectory)], dendrogram = "row",Colv = F, main="Growth factors in D20 ventricular data")
par(mfrow=c(2,4))
for (i in 1:length(chosen_trajectory_edge_pairs)){
  edges <- as.numeric(strsplit(chosen_trajectory_edge_pairs[[i]], split="-")[[1]])
  gf <- union(upreg_first_gf[[i]]$Pathway, upreg_later_gf[[i]]$Pathway)
  temp_exprMatrix_bycluster <- exprMatrix_bycluster[which(rownames(exprMatrix_bycluster) %in% gf), edges]
  heatmap.2(temp_exprMatrix_bycluster, dendrogram = "row", main=paste0("Changing growth factors between clusters ", chosen_trajectory_edge_pairs[[i]]))
}







