load("/Volumes/Thinh/Bader_Lab/heartdev/GSVA_Temporal_tools/D4_GSVA_Monocle.RData")
library(Matrix)
library(magrittr)
library(dplyr)
library(monocle)
library(tibble)
library(ggimage)
library(grid)
library(png)
reduced_dim_coords <- mes_d4_clust@reducedDimK
lib_info_with_pseudo <- pData(mes_d4_clust)

ica_space_df <- Matrix::t(reduced_dim_coords) %>%
  as.data.frame() %>%
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))

dp_mst <- mes_d4_clust@minSpanningTree

edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  select_(source = "from", target = "to") %>%
  left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
  left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")

sample_state <- mes_d4_clust@phenoData$State
sample_name <- NA
  
data_df <- t(mes_d4_clust@reducedDimS) %>%
  as.data.frame() %>%
  select_(data_dim_1 = 1, data_dim_2 = 2) %>%
  rownames_to_column("sample_name") %>%
  mutate(sample_state) %>%
  left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")


mst_branch_nodes <- mes_d4_clust@auxOrderingData[[mes_d4_clust@dim_reduce_type]]$branch_points

states <- list()
for (i in levels(data_df$State)){
  states[[i]] <- apply(data_df[which(data_df$State == as.numeric(i)), c(2,3)], 2, mean)
}

states <- Reduce(rbind, states)
rownames(states) <- unlist(levels(data_df$State))
# plot(t(mes_d4_clust@reducedDimS) )
# points(states, pch=15, col='red')



g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) 
g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=0.75, linetype="solid", na.rm=TRUE, data=edge_df) + theme(axis.line=element_blank(),
                                                                                                                                                                                                                           axis.text.x=element_blank(),
                                                                                                                                                                                                                           axis.text.y=element_blank(),
                                                                                                                                                                                                                           axis.ticks=element_blank(),
                                                                                                                                                                                                                           axis.title.x=element_blank(),
                                                                                                                                                                                                                           axis.title.y=element_blank(),
                                                                                                                                                                                                                           legend.position="none",
                                                                                                                                                                                                                           panel.background=element_blank(),
                                                                                                                                                                                                                           panel.border=element_blank(),
                                                                                                                                                                                                                           panel.grid.major=element_blank(),
                                                                                                                                                                                                                           panel.grid.minor=element_blank(),
                                                                                                                                                                                                                           plot.background=element_blank(),
                                                                                                                                                                                                                           plot.margin = unit(c(1, 1, 1, 1), "cm")) 


branch_point_df <- ica_space_df %>%
  slice(match(mst_branch_nodes, sample_name)) %>%
  mutate(branch_point_idx = seq_len(n()))

g <- g +
  geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
             size=5, na.rm=TRUE, branch_point_df) +
  geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2", label="branch_point_idx"),
            size=4, color="white", na.rm=TRUE, branch_point_df)
print(g)

filepath <- "/Volumes/Thinh/Bader_Lab/Temp_WC/state_1/"
wordcloud_files <- list.files(path=filepath, all.files=FALSE)
tmp_down <- list()
state_num <- length(wordcloud_files)
y_lim <- layer_scales(g)$y$range$range
x_lim <- layer_scales(g)$x$range$range
for (i in 1:state_num){
  tmp_down[[i]] <- rasterGrob(readPNG(paste0(filepath, wordcloud_files[i])), interpolate=TRUE)
  #tmp_up[[i]] <- rasterGrob(readPNG(paste0(filepath, wordcloud_files[(state_num+i)])), interpolate=TRUE)
  width <- 6
  if (states[i,2] < y_lim[1]) {
    g <- g + geom_subview(subview=tmp_down[[i]], x=states[i,1]+(0.1*width), y=states[i,2]+(0.1*width), width=width*1.15, heigh=width) + theme_transparent()
    #g <- g + geom_subview(subview=tmp_up[[i]], x=states[i,1]-(0.2*width), y=states[i,2]+(0.25*width), width=width, heigh=width) + theme_transparent()
  } else if (states[i,2] > y_lim[2]) {
    g <- g + geom_subview(subview=tmp_down[[i]], x=states[i,1]-(0.03*width), y=states[i,2]-(0.02*width), width=width*1.15, heigh=width) + theme_transparent()
    #g <- g + geom_subview(subview=tmp_up[[i]], x=states[i,1]-(0.2*width), y=states[i,2]-(0.15*width), width=width, heigh=width) + theme_transparent()
  } else {
    g <- g + geom_subview(subview=tmp_down[[i]], x=states[i,1]+(0.015*width), y=states[i,2]+(0.05*width), width=width*1.1, heigh=width) + theme_transparent()
    #g <- g + geom_subview(subview=tmp_up[[i]], x=states[i,1]-(0.2*width), y=states[i,2]+(0.1*width), width=width, heigh=width) + theme_transparent()
  }
}

print(g)



# tempfile_down <- readPNG(paste0(filepath, wordcloud_files[1]))
# tmp_grob <- rasterGrob(tempfile_down, interpolate=TRUE)
# tempfile_up <- readPNG(paste0(filepath, wordcloud_files[12]))
# tmp_grob_up <- rasterGrob(tempfile_up, interpolate=TRUE)
# width <- 2.5
# g <- g + geom_subview(subview=tmp_grob, x=states[1,1]+(0.4*width), y=states[1,2]+(0.2*width), width=width, heigh=width) + theme_transparent()
# g <- g + geom_subview(subview=tmp_grob_up, x=states[1,1]-(0.2*width), y=states[1,2], width=width, heigh=width) + theme_transparent()


