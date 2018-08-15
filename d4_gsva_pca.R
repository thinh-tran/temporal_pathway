rm(list=ls())
library(GSEABase)
library(httr)
library(RCy3)
library(RJSONIO)
library(scales)
library(dplyr)

##PCA
load("/Volumes/Thinh/Bader_Lab/heartdev/Data/d4_gsva.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/Data/Seurat_objects/Ventricular_D4.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D4_GSVA_Monocle.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D4_gsva_bystate.RData")
state <- subset(mes_d4_clust@phenoData@data, select = "Pseudotime")
state$cell_names <- rownames(state)
state <- state[order(state$Pseudotime, decreasing = F),]
state$cell_order <- seq(1:nrow(state))
d4_gsva_m <- matrix(as.numeric(unlist(d4_gsva)), nrow=nrow(d4_gsva), 
                    dimnames = list(rownames(d4_gsva), colnames(d4_gsva)))

#rownames(d4_gsva_m) <- sub("%.*", "", rownames(d4_gsva_m))
d4_gsva_pca <- prcomp(t(d4_gsva_m), center=T, scale=T)
screeplot(d4_gsva_pca, npcs=100, type="lines", main="PCA on D4 GSVA data")

pc_var_explained <- (((d4_gsva_pca$sdev)^2)/sum((d4_gsva_pca$sdev)^2))*100
#plot(cumsum(pc_var_explained[1:100]))
pc_used <- which(cumsum(pc_var_explained) >= 80)[1]
pc_used <- 50


###Try to remove the first 3 PCs
d4_gsva_denoised <- d4_gsva_pca$x[,3:pc_used] %*% t(d4_gsva_pca$rotation[,3:pc_used])
if(any(d4_gsva_pca$scale != 0) == TRUE){
  d4_gsva_denoised <- scale(d4_gsva_denoised, center = FALSE , scale=1/d4_gsva_pca$scale)
}
if(any(d4_gsva_pca$center != 0) == TRUE){
  d4_gsva_denoised <- scale(d4_gsva_denoised, center = -1 * d4_gsva_pca$center, scale=FALSE)
}





###Annotate significant pathways

significant_pathways <- pc_cor <- significant_pvals <- list()
for (i in 1:pc_used){
  genes_scaled <- scale(d4_gsva_pca$rotation[,i])
  significant_pathways[[i]] <- names(which(genes_scaled[,1] > 1.5 | genes_scaled[,1] < -1.5))
  significant_pvals[[i]] <- rescale(d4_gsva_pca$rotation[,i][which(genes_scaled[,1] > 1.5 | genes_scaled[,1] < -1.5)], to=c(0, 0.05))
  #pc_cor[[i]] <- cor(t(d4_gsva_m[which(rownames(d4_gsva_m) %in% significant_pathways[[i]]), ]))
}

#significant_pathways <- unique(unlist(significant_pathways))

pathways_db <- getGmt("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Human_GO_AllPathways_with_GO_iea_July_01_2018_symbol.gmt")

pathway_names <- pathway_description <- list()
for (i in 1:length(pathways_db)){
  pathway_names[[i]] <- pathways_db[[i]]@setName
  pathway_description[[i]] <- pathways_db[[i]]@shortDescription
}

pathway_names <- as.character(Reduce(rbind, pathway_names))
pathway_description <- as.character(Reduce(rbind, pathway_description))

pathways_db_table <- data.frame(Pathway_name = pathway_names, Description = pathway_description)

filepath <- "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Significant_signaling_pathways_by_PCs_enrichment/"

significant_pvals_signaling <- lapply(significant_pvals, function(x) x[which(grepl("SIGNALING", names(x))==TRUE)])
significant_pathways_signaling <- lapply(significant_pathways, function(x) grep("SIGNALING", x, value=TRUE))


for (i in 1:length(significant_pathways_signaling)){
  temp_enrichment_file <-  pathways_db_table[which(pathways_db_table$Pathway_name %in% significant_pathways_signaling[[i]]), ]
  p_values <- data.frame(Pathway_name = names(significant_pvals_signaling[[i]]), p_value=significant_pvals_signaling[[i]])
  temp_enrichment_file <- inner_join(temp_enrichment_file, p_values)  
  write.table(temp_enrichment_file, file=paste0(filepath, "PC_", i, ".txt"), quote=F, sep="\t", row.names=F)
}


####Cytoscape

source("/Volumes/Thinh/Bader_Lab/summarize_em_themes.R")

port.number = 1234
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")

pvalue_threshold <- "0.05"

similarity_threshold <- "0.25"
similarity_metric = "COMBINED"

generic_gmt_file <- "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Human_GO_AllPathways_with_GO_iea_June_01_2018_symbol.gmt"

network_suid <- summary <- themes <- list()
for (i in 1:pc_used){
  print(i)
  results_filename <- paste0("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Significant_pathways_by_PCs_enrichment/PC_", i, ".txt")
  
  em_command = paste('enrichmentmap build analysisType="generic" gmtFile=',generic_gmt_file,
                     'pvalue=',pvalue_threshold,
                     'similaritycutoff=',similarity_threshold,'enrichmentsDataset1=',results_filename,
                     sep=" ")
  
  #current_network_name <- paste0("PC_", i)
  network_suid[[i]] <- commandRun(em_command)
  #response <- renameNetwork(current_network_name, network = current_network_suid)
  
  summary[[i]] <- summarize_em_themes(as.numeric(network_suid[[i]]))
  summary[[i]] <- summary[[i]][order(summary[[i]]$number_of_nodes, decreasing=T), ]
  themes[[i]] <- as.character(unlist(summary[[i]]$label))[1:10]
}


# #cur_model_name <- "d4_gsva_pca_1" 
# results_filename <- "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Significant_pathways_by_PCs_enrichment/PC_1.txt"
# 
# em_command = paste('enrichmentmap build analysisType="generic" gmtFile=',generic_gmt_file,
#                    'pvalue=',pvalue_threshold,
#                    'similaritycutoff=',similarity_threshold,'enrichmentsDataset1=',results_filename,
#                    sep=" ")
# 
# current_network_name <- "PC_1"
# current_network_suid <- commandsRun(em_command)
# response <- renameNetwork(current_network_name, network = current_network_suid)
# 
# summary <- summarize_em_themes(as.numeric(current_network_suid))
# summary <- summary[order(summary$number_of_nodes, decreasing=T), ]
# themes <- as.character(unlist(summary$label))[1:10]
# save(summary, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/EM_theme_50PCs.RData")
# save(theme, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/EM_themes_only_50PCs.RData")


##plot trees and save

d4_gsva_significant_genes <- d4_gsva_m[which(rownames(d4_gsva_m) %in% important_genes), ]

matching_states <- data.frame(cell_names = colnames(d4_gsva_significant_genes))
matching_states <- merge(state, matching_states)

colnames(d4_gsva_significant_genes) <- matching_states$cell_order
d4_hclust <- hclust(dist(t(d4_gsva_significant_genes)))
d4_subgroup <- cutree(d4_hclust, k=10)

svg(file="/Volumes/Thinh/Bader_Lab/GSVA_Temporal_tools/d4_pca_heatmap.svg", width=60, type="cairo")
heatmap.2(d4_gsva_significant_genes, ColSideColors = as.character(d4_subgroup), trace="none", cexRow = 0.1, cexCol = 0.1)
dev.off()
early <- list()
for (i in seq(1,10)){
  early[[i]] <- length(which(as.numeric(names(d4_subgroup)[which(d4_subgroup == i)]) < 950))/length(which(d4_subgroup == i))
}
early_df <- Reduce(rbind, early)
rownames(early_df) <- seq(1:10)


###Annotate PCs
for (i in 1:length(summary)){
  summary[[i]]$PC <- rep(i, nrow(summary[[i]]))
}
sig_pathways <- Reduce(rbind, lapply(summary, function(x) x[1:10,]))
sig_pathways <- sig_pathways[order(sig_pathways$number_of_nodes, decreasing=T), ]
sig_pathways <- sig_pathways[!duplicated(as.character(unlist(sig_pathways$label))), ]

corpus <- Corpus(VectorSource(as.character(unlist(sig_pathways$label))))
toSpace <- content_transformer(function (x , pattern) gsub(pattern, " ", x))


corpus <- tm_map(corpus, removeWords, c("pathway", "signaling", "receptor", "factor", "members",
                                        "downstream", "signalling", "cell", "type", "network",
                                        "mediated", "dependent"))

dtm <-TermDocumentMatrix(corpus)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)



assocs <- list()
for (i in 1:nrow(d[which(d$freq >= 4),])){
  assocs[[i]] <- names(findAssocs(dtm, as.character(d[which(d$freq >= 4),]$word[[i]]), 0.3)[[1]])
}
assocs1 <- unlist(assocs)

d_subset <- d[which(d$freq > 2), ]
d_subset_index <- which(d_subset$word %in% assocs1 & d_subset$freq < 4)
d_subset <- d_subset[-d_subset_index, ]

pop_list <- c()
for (i in 1:(length(nodes_in_themes)-1)){
  if (length(intersect(nodes_in_themes[[i]], nodes_in_themes[[i+1]]))/length(nodes_in_themes[[i]]) > 0.8){
    pop_list <- c(pop_list, i)
  }
}

d_subset <- d_subset[-pop_list, ]


pathways_in_themes <- nodes_in_themes <- list()
for (i in 1:nrow(d_subset)){
  pathways_in_themes[[i]] <- grep(as.character(unlist(sig_pathways$label)), pattern = d_subset$word[[i]], value=TRUE)
  nodes_in_themes[[i]] <- Reduce(unlist, lapply(as.character(unlist(sig_pathways$node_labels[grep(as.character(unlist(sig_pathways$label)), pattern = d_subset$word[[i]])])), 
                                                function(x) strsplit(x, split=";")))
}
names(pathways_in_themes) <- names(nodes_in_themes) <- d_subset$word

theme_expression <- list()
for (i in 1:length(nodes_in_themes)){
  theme_expression[[i]] <- colSums(d4_gsva_bystate[which(rownames(d4_gsva_bystate) %in% nodes_in_themes[[i]]), ])
}
theme_expression <- Reduce(rbind, theme_expression)
rownames(theme_expression) <- names(nodes_in_themes)

correlation_theme_exp <- cor(t(theme_expression))

##remove redundancy
redundancy <- apply(correlation_theme_exp, 2, function(x) rownames(correlation_theme_exp)[which(x > 0.95)])

redundancy_filtered <- redundancy[which(names(redundancy) %in% rownames(theme_exp_unique))]

##list out pathways that are perfectly overlapping
redundant_pathways <- list()
for (i in 1:length(high_info)){
  redundant_pathways[[i]] <- names(redundancy)[which(sapply(redundancy, function(x) isTRUE(all.equal(x, high_info[[i]]))))]
}
names(redundant_pathways) <- high_info

keep_list <- c("notch", "respiratory", "wnt", "integrin", "senescence", "degradation", "metabolism",
               "smad2", "glucose", "insulin", "nucleus", "oxidation", "translation", "braf", "nectin")

keep_list <- c(keep_list, as.character(redundant_pathways[which(lengths(redundant_pathways) == 1)]))

theme_exp_unique <- theme_expression[which(rownames(theme_expression) %in% keep_list), ]
theme_exp_unique <- theme_exp_unique[-which(rownames(theme_exp_unique) %in% setdiff(unique(unlist(redundancy)), keep_list)), ]

redundancy_filtered <- redundancy[which(names(redundancy) %in% rownames(theme_exp_unique))]

for (i in 1:length(redundancy_filtered)){
  if (length(intersect(rownames(theme_exp_unique), redundancy_filtered[[i]])) > 1){
    temp_intersect <- intersect(rownames(theme_exp_unique), redundancy_filtered[[i]])
    print(temp_intersect)
    important_word <- readline(prompt="Identify a key word: ")
    theme_exp_unique <- theme_exp_unique[-(which(rownames(theme_exp_unique) 
                                                 %in% temp_intersect[-which(temp_intersect == as.character(important_word))])), ]
  }
}

lowinfo <- c("pid", "defective", "complex", "nucleus", "diseases", "fatty", "rho", "collagen", "glucagon", "naba", "a6b1")
theme_exp_unique <- theme_exp_unique[-which(rownames(theme_exp_unique) %in% lowinfo), ]

pathways_in_themes_subset <- pathways_in_themes[which(names(pathways_in_themes) %in% rownames(theme_exp_unique))]
pathway_names <- c()
for (i in 1:length(pathways_in_themes_subset)){
  print(pathways_in_themes_subset[[i]])
  temp_label <- readline(prompt="Identify main pathway name: ")
  pathway_names <- c(pathway_names, temp_label)
}


##Plot
colors <- colorRampPalette(brewer.pal(8, "Spectral"))
colors_count <- nrow(theme_expression)


par(new=FALSE)
plot(theme_expression[1,], type='l', xaxt="n", bty="n", yaxt="n", xlab = "State", ylab="Enrichment", ylim=c(-120, 80), col=colors(colors_count)[i])
for (i in 2:nrow(theme_expression)){
  par(new=TRUE)
  plot(theme_expression[i,], type='l', xaxt="n", bty="n", yaxt="n", xlab = "", ylab="",  ylim=c(-120, 80), col=colors(colors_count)[i])
}
axis(side=1, at=seq(1,11), pos=c(0,0), outer=TRUE, lwd.ticks = 0)
axis(side=2, at=seq(-60, 80, 10))
title(main="Theme enrichment over pseudotime state")
legend("bottom", legend=rownames(theme_expression), fill=NULL, lty=rep(1, nrow(theme_expression)), col=colors(colors_count), bty="n", cex=0.4, ncol=5, xjust=0)


###compare with seurat cluster
nodes_in_themes_subset <- nodes_in_themes[which(names(nodes_in_themes) %in% rownames(theme_exp_unique))]

clusters <- subset(eb1S@meta.data, select="res.0.8")
expr_bycluster <- list()
for (i in levels(clusters$res.0.8)){
  expr_bycluster[[i]] <- rowMeans(d4_gsva_m[, which(colnames(d4_gsva_m) %in% rownames(clusters)[which(clusters$res.0.8 == i)])])
}
expr_bycluster <- do.call(cbind, expr_bycluster)

theme_expr_bycluster <- lapply(nodes_in_themes_subset, function(x) colSums(expr_bycluster[which(rownames(expr_bycluster) %in% x), ]))
theme_expr_bycluster <- do.call(rbind, theme_expr_bycluster)

par(new=FALSE)
plot(theme_expr_bycluster[1,], type='l', xaxt="n", bty="n", yaxt="n", xlab = "Cluster", ylab="Enrichment", ylim=c(-25, 10), col=colors(colors_count)[i])
for (i in 2:nrow(theme_exp_unique)){
  par(new=TRUE)
  plot(theme_expr_bycluster[i,], type='l', xaxt="n", bty="n", yaxt="n", xlab = "", ylab="",  ylim=c(-25, 10), col=colors(colors_count)[i])
}
axis(side=1, at=seq(1,6), pos=c(0,0), outer=TRUE, lwd.ticks = 0)
axis(side=2, at=seq(-15, 10, 5))
title(main="Theme enrichment across Seurat clusters")
legend("bottom", legend=pathway_names, fill=NULL, lty=rep(1, nrow(theme_exp_unique)), col=colors(colors_count), bty="n", cex=0.6, ncol=2, xjust=0)
''


########Most variable pathways##########
d4_gsva_m <- as.matrix(d4_gsva_pseudotime)
rownames(d4_gsva_m) <- sub("%.*", "", rownames(d4_gsva_m))

#rownames(d4_gsva_m) <- sub("%.*", "", rownames(d4_gsva_m))
pathway_variance <- sort(apply(d4_gsva_m[-c(nrow(d4_gsva_m)-1, nrow(d4_gsva_m)), ], 1,
                               var), decreasing=T)
pathway_variance <- sort(pathway_variance, decreasing = T)
hist(pathway_variance)
significant_genes <- names(pathway_variance)[which(scale(pathway_variance) > 1.5)]
signficant_genes_matrix <- d4_gsva_m[which(rownames(d4_gsva_m) %in% significant_genes), ]
colnames(signficant_genes_matrix) <- seq(1:ncol(signficant_genes_matrix))
d4_hclust_1 <- hclust(dist(t(signficant_genes_matrix)))
d4_subgroup_1 <- cutree(d4_hclust, k=5)

early_1 <- list()
for (i in seq(1,5)){
  early_1[[i]] <- length(which(as.numeric(names(d4_subgroup_1)[which(d4_subgroup_1 == i)]) < 950))/length(which(d4_subgroup_1 == i))
}
early_df_1 <- Reduce(rbind, early_1)
rownames(early_df_1) <- seq(1:5)

significant_gene_matrix <- apply(signficant_genes_matrix, 1, as.numeric)
heatmap.2(t(significant_gene_matrix), trace="none", ColSideColors=as.character(d4_subgroup), cexRow = 0.1, cexCol = 0.1)

hallmark_matrix <- d4_gsva_m[which(grepl("^HALLMARK", rownames(d4_gsva_m))== TRUE), ]
hallmark_matrix_num <- apply(hallmark_matrix, 2, as.numeric)
rownames(hallmark_matrix_num) <- sub("^HALLMARK_", "", rownames(hallmark_matrix))
colnames(hallmark_matrix_num) <- as.character(seq(1:ncol(hallmark_matrix_num)))
d4_hm_clust <- hclust(dist(t(hallmark_matrix_num)))
d4_hm_subgroup <- cutree(d4_hm_clust, k=5)
heatmap.2(hallmark_matrix_num, trace="none", ColSideColors=as.character(d4_hm_subgroup), cexRow = 0.3, cexCol = 0.1)


#line plot for each hallmark pathway 
library(RColorBrewer)
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D4_gsva_bystate.RData")
colors <- colorRampPalette(brewer.pal(11, "Spectral"))
colors_count <- nrow(gsva_bystate_hallmark)

gsva_bystate_hallmark <- d4_gsva_bystate[which(grepl("^HALLMARK", rownames(d4_gsva_bystate))== TRUE), ]
rownames(gsva_bystate_hallmark) <- sub("^HALLMARK_", "", rownames(gsva_bystate_hallmark))
rownames(gsva_bystate_hallmark) <- sub("%.*", "", rownames(gsva_bystate_hallmark))

for (i in 1:nrow(gsva_bystate_hallmark)){
  par(new=TRUE)
  plot(gsva_bystate_hallmark[i, c(1, 3:ncol(gsva_bystate_hallmark))], type='b', xlab="", ylab="", xaxt="n", yaxt="n")
}

gsva_hm_without2 <- gsva_bystate_hallmark[, c(1, 3:ncol(gsva_bystate_hallmark))]

plot(gsva_hm_without2[1,], type='l', xaxt="n", bty="n", yaxt="n", xlab = "Pseudotime state", ylab="Enrichment", ylim=c(-1, 0.8), col=colors(colors_count)[1], 
     main="Hallmark pathway enrichment over pseudotime")
for (i in 2:nrow(gsva_bystate_hallmark)){
  par(new=TRUE)
  plot(gsva_hm_without2[i,], type='l', xaxt="n", bty="n", yaxt="n", xlab = "", ylab="", ylim=c(-1, 0.8), col=colors(colors_count)[i])
}
axis(side=1, at=seq(1,10), labels=c(1, 3:ncol(gsva_bystate_hallmark)), pos=c(0,0), outer=TRUE, lwd.ticks = 0)
axis(side=2, at=seq(-1, 1.2, 0.2))
legend("bottom", legend=rownames(gsva_hm_without2), fill=NULL, lty=rep(1, 50), col=colors(colors_count), bty="n", cex=0.3, ncol=3, xjust=0)

save(pathway_df, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Pathways_table.RData")
