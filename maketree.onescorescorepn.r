library(ape)
library(dendextend)
library(tidyverse)
library(ggtree)
library(gridExtra)
library(ggplot2)
library(ggtreeExtra)
library(ggplot2)
library(ggstar)


df <- read.csv("/home/cnr/Desktop/THEWORK2803/10.Amyloidosis(heart)_H.sapiens/data_for_tree/tabs_first5_coppi_chord/TTR_control.csv")
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

p1 <- ggtree(newick, branch.length="none", layout="circular", open.angle=10)+
  geom_tiplab(align=TRUE, linesize = 0, size = 4, offset = 4) + 
  theme(plot.title = element_text(hjust = 0.5)) 


 
p2 <- p1+ 
  geom_fruit(
    data = df,
    geom = geom_col,
    mapping = aes(y=description, x=score, fill = correlated.to),
    axis.params=list(
      axis = "x"),
    grid.params=list()) + scale_fill_manual(values=c("#49FFFF","#F46D43")) 


p3 <- p2 +
  geom_fruit(
    data = df,
    geom = geom_star,
    mapping = aes(y=description, starshape=category),
    size = 2)+ scale_starshape_manual(values=c(1,  12))

p3 = p3 +
  theme(
    plot.background = element_rect(fill = "white", color = NA, size = 3)  # Adjust size as needed
  )

# Save the plot with a white background covering the whole plot area
ggsave("TTR_control.svg", p3, device = "svg",
       width = 8, height = 6, units = "in", dpi = 1200)

