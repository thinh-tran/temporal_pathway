#take a list of gsva matrix and a monocle output with cell state 
#do a wilcoxon rank-sum test to find pathways that are DE between one state and the rest
#of the cells
#split the result into pathways that are upregulated (higher average expression in one
#specific state) and downregulated (lower average expression in one specific state)
#returns two lists: upregulated and downregulated pathways in each state
findDEPathways <- function(exprMatrix, monocle){
  upreg <- downreg <- list()
  idents <- levels(monocle@phenoData@data$State)
  for (i in idents){
    print(i)
    ioi <- i
    ioi_exprMatrix <- exprMatrix[, which(monocle@phenoData@data$State == ioi)]
    other <- setdiff(idents, ioi)
    other_exprMatrix <- exprMatrix[, which(monocle@phenoData@data$State %in% other)]
    
    wilcox_stats <- pvals <- list()
    if (is.matrix(ioi_exprMatrix) == FALSE) {
      for (j in 1:length(ioi_exprMatrix)) {
        wilcox_stats[[j]] <- wilcox.test(as.numeric(ioi_exprMatrix[j]), as.numeric(other_exprMatrix[j,]), p.adjust.methods="fdr")
        pvals[[j]] <- wilcox_stats[[j]]$p.value #output p-value in a list
      }
    } else {  
      for (j in 1:nrow(ioi_exprMatrix)){
        wilcox_stats[[j]] <- wilcox.test(as.numeric(ioi_exprMatrix[j,]), as.numeric(other_exprMatrix[j,]), p.adjust.methods="fdr")
        pvals[[j]] <- wilcox_stats[[j]]$p.value #output p-value in a list
      }
    }     
    
    names(wilcox_stats) <- names(pvals) <- sub("%.*", "", rownames(exprMatrix)) 
    result_df <- t(as.data.frame(pvals))
    colnames(result_df) <- c("Wilcox_test_adj_pvalue")
    result_df <- as.data.frame(result_df)
    
    mean_diff <- mapply(function(x,y) x-y, rowMeans(as.matrix(ioi_exprMatrix)), rowMeans(as.matrix(other_exprMatrix)))
    result_df$Mean_difference_of_interest_vs_other_myocyte <- mean_diff
    
    upreg[[i]] <- result_df[which(result_df$Wilcox_test_adj_pvalue < 0.05 & result_df$Mean_difference_of_interest_vs_other_myocyte > 0), ]
    downreg[[i]] <- result_df[which(result_df$Wilcox_test_adj_pvalue < 0.05 & result_df$Mean_difference_of_interest_vs_other_myocyte < 0), ]
  }
  tbr <- list(upreg, downreg)
  names(tbr) <- c("Upregulated pathways", "Downregulated pathways")
  return(tbr)
}



FindRepPathways <- function(pathwaylist, plot=FALSE, col=NULL, min_freq=NULL, max_words=NULL){
  corpus <- lapply(pathwaylist, function(x) Corpus(VectorSource(x)))
  toSpace <- content_transformer(function (x , pattern) gsub(pattern, " ", x))
  corpus <- lapply(corpus, function(x) tm_map(x, toSpace, "\\_"))
  corpus <- lapply(corpus, function(x) tm_map(x, toSpace, "\\."))
  
  corpus <- lapply(corpus, function(x) tm_map(x, removeWords, c("SIGNALING", "PATHWAY", "IN", "AND", "MEDIATED", "EVENTS", "RECEPTOR","MUTANTS", "BIOCARTA","FAMILY", "DOWNSTREAM", "HALLMARK", "THE", "TYPE", "CLASS", "NONCANONICAL", "CELL", "INTERLEUKIN", "REGULATION")))
  
  dtm <- lapply(corpus, function(x) TermDocumentMatrix(x))
  m <- lapply(dtm, as.matrix)
  v <- lapply(m, function(x) sort(rowSums(x),decreasing=TRUE))
  d <- lapply(v, function(x) data.frame(word = names(x),freq=x))
  
  if (plot == TRUE){
    for (i in 1:length(d)){
      if (nrow(d[[i]]) == 0){
        plot(0,type='n',axes=FALSE,ann=FALSE)
      } else {
        wordcloud(d[[i]]$word, d[[i]]$freq, min.freq=min_freq, max.words=max_words, scale=c(1.5, .3),colors = brewer.pal(8, col), rot.per=0)
      }
    }
  }
  
  freq_upreg_terms <- lapply(dtm, function(x) findFreqTerms(x, lowfreq = 4))
  common_freq <- Reduce(intersect, freq_upreg_terms)
  freq_upreg_terms_f <- lapply(freq_upreg_terms, function(x) setdiff(x, common_freq))
  return(freq_upreg_terms_f)
}

