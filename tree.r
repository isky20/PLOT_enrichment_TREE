library(ape)
library(dendextend)
library(tidyverse)
library(ggtree)
library(gridExtra)
library(ggplot2)
library(ggtreeExtra)

ref <- read.csv("ref_enrichment_endometriosis_norm.csv")
data <- read.csv("aim.csv")

RCTM <- merge(data,ref, by = c("term","category","description"), all.x = TRUE)
RCTM_rows <- RCTM[RCTM$category == "RCTM", ]
RCTM_rows <- RCTM_rows[!is.na(RCTM_rows$MergedGenes), ]

KEGG <- merge(data,ref, by = c("term","category","description"), all.x = TRUE)
KEGG_rows <- KEGG[KEGG$category == "KEGG", ]
KEGG_rows <- KEGG_rows[!is.na(KEGG_rows$MergedGenes), ] # remove empty row in gene list




for (df in list(RCTM_rows, KEGG_rows)){
  # Get unique genes across all terms
  title <- paste(unique(df$category), collapse = ",")
  
  all_genes <- unique(unlist(strsplit(gsub(" ", "", df$MergedGenes), ",")))
  # Create a binary matrix indicating gene presence for each term
  
  binary_matrix <- matrix(0, nrow = nrow(df), ncol = length(all_genes), dimnames = list(df$description, all_genes)) #df$term
  for (i in seq_along(df$MergedGenes)) {
    gene_list <- unlist(strsplit(gsub(" ", "", df$MergedGenes[i]), ","))
    binary_matrix[df$description[i], intersect(all_genes, gene_list)] <- 1 #df$term
  }
  binary_position_matrix <- (binary_matrix == 1)
  
  # Compute Jaccard similarity matrix
  matrix <- dist(binary_position_matrix,  method = "binary")
  hclust_object <- hclust(as.dist(matrix), method = "ward.D2" )
  # Create a dendrogram object
  dendrogram_object <- as.dendrogram(hclust_object)
  newick <- as.phylo(dendrogram_object)
  
  p1 <- ggtree(newick, branch.length="none", layout="circular", open.angle=15) + 
    geom_tiplab(align=TRUE, linesize = 0, size = 3.5, offset = 7) + ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+ 
    xlim(0, 20)
  selected_data <- data[, c(1, grep("^count_g", names(data)))] # select columns starting with "count_g"
  
  normalized_df <- apply(selected_data[, 2:ncol(selected_data)], 1, function(row) {
    normalized_row <- (row / max(row)) * 100  # Normalize each value in the row
    return(normalized_row)
  })
  data_normalized <- as.data.frame(t(normalized_df))
  
  data_normalized$description <- selected_data$description
  p2 <- p1 + 
    geom_fruit(
      data = data_normalized,
      geom = geom_col,
      mapping = aes(y=description, x=count_g1, fill = "C")
    ) + scale_size_continuous(range = c(1, 20))
  
  p3 <- p2 + 
    geom_fruit(
      data = data_normalized,
      geom = geom_col,
      mapping = aes(y=description, x=count_g2, fill = "EII")
    ) + scale_size_continuous(range = c(1, 20))
  
  p4 <- p3 + 
    geom_fruit(
      data = data_normalized,
      geom = geom_col,
      mapping = aes(y=description, x=count_g3, fill = "EIV")
    ) + scale_size_continuous(range = c(1, 20)) + scale_fill_manual(values=c("red","blue","green"),
                                                                    guide=guide_legend(title = "Case",keywidth=1, keyheight=1)) +
    theme(legend.position=c(1, 0.51),
          legend.title=element_text(size=7),
          legend.text=element_text(size=7),
          legend.spacing.y = unit(0.1, "cm")) 
  
  ggsave(paste(gsub(" ", "", title),".svg",collapse = "_"), p4, device = "svg",
         width = 18, height = 16, units = "in", dpi = 500)
}

