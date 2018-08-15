library(dplyr)
library(scales)
library(wordcloud)
library(RCy3)
library(httr)
library(RJSONIO)
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D20_gsva_bystate.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Pathways_table.RData")
filepath <- "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/enrichment_bystate_d20/"
pathway_df$Pathway_name <- as.character(unlist(pathway_df$Pathway_name))

upreg_pvals <- apply(d20_gsva_bystate, 2, function(x) rescale(x[which(x > 0)], to=c(0, 0.05)))
upreg_bystate <- lapply(upreg_pvals, function(x) data.frame(Pathway_name = as.character(unlist(names(x))), pvalues = as.numeric(unlist(x))))

upreg_bystate <- lapply(upreg_bystate, function(x) inner_join(x, pathway_df, by="Pathway_name"))

downreg_pvals <- apply(d20_gsva_bystate, 2, function(x) rescale(x[which(x < 0)], to=c(0, 0.05)))
downreg_bystate <- lapply(upreg_pvals, function(x) data.frame(Pathway_name = names(x), pvalues = x))
upreg_bystate <- lapply(upreg_bystate, function(x) inner_join(x, pathway_df, by="Pathway_name"))


for (i in 1:length(upreg_bystate)){
  #upreg_bystate[[i]]$Description <- as.character(unlist(upreg_bystate[[i]]$Description))
  upreg_bystate[[i]] <- upreg_bystate[[i]][, c(1,3,2)]
  write.table(upreg_bystate[[i]], file=paste0(filepath, "state_", i, ".txt"), quote=F, sep="\t", row.names=F)
}

source("/Volumes/Thinh/Bader_Lab/summarize_em_themes.R")

port.number = 1234
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")

pvalue_threshold <- "0.05"

similarity_threshold <- "0.25"
similarity_metric = "COMBINED"

generic_gmt_file <- "/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/Human_GO_AllPathways_no_GO_iea_July_01_2018_symbol.gmt"

network_suid_state <- summary_state <- themes_state <- list()
for (i in 1:11){
  print(i)
  results_filename <- paste0("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/enrichment_bystate_d20/state_", i, ".txt")
  
  em_command = paste('enrichmentmap build analysisType="generic" gmtFile=',generic_gmt_file,
                     'pvalue=',pvalue_threshold,
                     'similaritycutoff=',similarity_threshold,'enrichmentsDataset1=',results_filename,
                     sep=" ")
  
  #current_network_name <- paste0("PC_", i)
  network_suid_state[[i]] <- commandsRun(em_command)
  #response <- renameNetwork(current_network_name, network = current_network_suid)
  
  summary_state[[i]] <- summarize_em_themes(as.numeric(network_suid_state[[i]]))
  summary_state[[i]] <- summary_state[[i]][order(summary_state[[i]]$number_of_nodes, decreasing=T), ]
  themes_state[[i]] <- as.character(unlist(summary_state[[i]]$label))[1:10]
}


# results_filename <- paste0("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/enrichment_bystate/state_", i, ".txt")
# 
# em_command = paste('enrichmentmap build analysisType="generic" gmtFile=',generic_gmt_file,
#                    'pvalue=',pvalue_threshold,
#                    'similaritycutoff=',similarity_threshold,'enrichmentsDataset1=',results_filename,
#                    sep=" ")
# 
# #current_network_name <- paste0("PC_", i)
# test1 <- commandsRun(em_command)
# #response <- renameNetwork(current_network_name, network = current_network_suid)
# 
# summary_test <- summarize_em_themes(as.numeric(test1))
# summary_test <- summary_test[order(summary_test$number_of_nodes, decreasing=T), ]
# themes_state[[i]] <- as.character(unlist(summary_state[[i]]$label))[1:10]
# 
# 
# edgetable_colnames <- getTableColumnNames(table="edge",  network = as.numeric(test1))
# 
# #get the correct attribute names
# similarity_attrib <- edgetable_colnames[grep(edgetable_colnames, pattern = "similarity_coefficient")]
# 
# #get the column from the nodetable and node table
# nodetable_colnames <- getTableColumnNames(table="node",  network = as.numeric(test1))
# 
# descr_attrib <- nodetable_colnames[grep(nodetable_colnames, pattern = "_GS_DESCR")]
# 
# aa_label_url <- paste("autoannotate label-clusterBoosted labelColumn=", descr_attrib," maxWords=3 nodeList=SUID:192515", sep="")
# current_name <-commandsGET(aa_label_url)
# 

for (i in 1:length(summary_state)){
  summary_state[[i]] <- summary_state[[i]][order(summary_state[[i]]$number_of_nodes, decreasing=T),]
  min_freq <- summary_state[[i]]$number_of_nodes[5]
  filepath <- "/Users/thinhtran/Desktop/Bader_Lab/Temp_WC/state/"
  png(file=paste0(filepath, i, ".png"), width=200, height=200, res=100, bg='transparent')
  par(mar = rep(0, 4))
  wordcloud(summary_state[[i]]$label, summary_state[[i]]$number_of_nodes, min.freq=min_freq, max.words=20, scale=c(1.1, .3),colors = brewer.pal(8, "Dark2"), rot.per=0)
  dev.off()
}

save(summary_state, file="/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/EM_themes_bystate.RData")
load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/EM_themes_bystate.RData")

node_list <- lapply(summary_state, function(x) strsplit(as.character(unlist(x$node_labels)), split=";"))


node_gsva_score <- tbg <- list()
for (i in 1:length(node_list)){
  node_gsva_score[[i]] <- lapply(node_list[[i]], function(t) 
    median(d4_gsva_bystate[which(rownames(d4_gsva_bystate) %in% t), i]))
  summary_state[[i]]$node_score <- unlist(node_gsva_score[[i]])
  tbg[[i]] <- summary_state[[i]][which(summary_state[[i]]$number_of_nodes > 5 
                                       & summary_state[[i]]$label != "No annotation returned"), c(2, 5)]
  tbg[[i]] <- tbg[[i]][order(tbg[[i]]$node_score, decreasing = T), ]
  tbg[[i]] <- tbg[[i]][c(1:5), ]
  tbg[[i]]$freq <- rescale(tbg[[i]]$node_score, to=c(2, 5))
}


for (i in 1:length(tbg)){
  #summary_state[[i]] <- summary_state[[i]][order(summary_state[[i]]$number_of_nodes, decreasing=T),]
  #min_freq <- summary_state[[i]]$number_of_nodes[5]
  filepath <- "/Volumes/Thinh/Bader_Lab/Temp_WC/state_1/"
  png(file=paste0(filepath, i, ".png"), width=400, height=200, res=90, bg='transparent')
  par(mar = rep(0, 4))
  wordcloud(tbg[[i]]$label, tbg[[i]]$freq, min.freq=2, max.words=10, scale=c(1, .3),colors = brewer.pal(8, "Dark2"), rot.per=0)
  dev.off()
}



