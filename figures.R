library(synapser)
synLogin()
library(data.table)
library(fastcluster)
library(ggdendro)
library(dendextend)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# mRNA dendrogram
# ----------------------------------------------------------------------------------
gene_types <- synGet('syn9846154')
dfg <- read.table(gene_types$path,stringsAsFactors = TRUE, 
                  sep=',', header=TRUE)
coding_genes = as.character(dfg[dfg$gene_biotype == 'protein_coding', ]$hgnc_symbol)
# Get RNASeq file
file <- synGet('syn7080753')
dfr <- read.table(file$path, stringsAsFactor = TRUE, 
                 header=TRUE, row.names='id')
df <- dfr[rownames(dfr) %in% coding_genes, ]
dft <- transpose(df)
rownames(dft) <- colnames(df)
colnames(dft) <- rownames(df)

cell_line_order <- rownames(dft)#[hc$order]
# Get cell line categories
dc <- as.data.frame(read_csv('cell_lines_categories.csv'))
dc2 <- dc[dc$Sample %in% cell_line_order, ]
# Sort dc2 by cell line order
dc3 <- dc2[match(cell_line_order, dc2$Sample), ]
colors = add.alpha(as.character(dc3$color_map), 0.7)
dc4 <-  as.data.frame(dc3$Receptor_Status)
colnames(dc4) <- "Receptor_status"
rownames(dc4) <- dc3$Sample
rs_colors <- add.alpha(c("red", "#1F78B4","orange2", '#33A02C'), 0.9)
#rs_colors <- add.alpha(c("gray", "#1F78B4","#FB9A99", '#33A02C'), 0.9)
#rs_colors <- c('blue', 'black', 'red', 'green')
names(rs_colors) <- c('NM', 'HER2amp', 'TNBC', 'HR+')
anno_colors <- list(Receptor_status=rs_colors)

gg1 <- pheatmap(dft,
               color=0, cellwidth=0, legend=FALSE,
               cluster_cols=FALSE, show_colnames=FALSE,
               annotation_row=dc4, annotation_colors=anno_colors,
               annotation_names_row=FALSE, annotation_legend = FALSE,
               main='mRNA')

# protein (iBAQ) dendrogram
# ----------------------------------------------------------------------------------
file <- synGet('syn7437201')
dfr <- read.table(file$path, stringsAsFactor = TRUE, 
                  header=TRUE, row.names='Uniprot_Id', sep=',')
df <- dfr[, 2:37]
dft <- t(df)
rownames(dft) <- colnames(df)
colnames(dft) <- rownames(df)

cell_line_order <- rownames(dft)#[hc$order]
# Get cell line categories
# dc <- as.data.frame(read_csv('cell_lines_categories.csv'))
dc2 <- dc[dc$Sample %in% cell_line_order, ]
# Sort dc2 by cell line order
dc3 <- dc2[match(cell_line_order, dc2$Sample), ]
colors = add.alpha(as.character(dc3$color_map), 0.7)
dc4 <-  as.data.frame(dc3$Receptor_Status)
colnames(dc4) <- "Receptor_status"
rownames(dc4) <- dc3$Sample
rs_colors <- add.alpha(c("red", "#1F78B4","orange2", '#33A02C'), 0.9)
#rs_colors <- add.alpha(c("gray", "#1F78B4","#FB9A99", '#33A02C'), 0.9)
#rs_colors <- c('blue', 'black', 'red', 'green')
names(rs_colors) <- c('NM', 'HER2amp', 'TNBC', 'HR+')
anno_colors <- list(Receptor_status=rs_colors)

gg2 <- pheatmap(dft,
                color=0, cellwidth=0, legend=FALSE,
                cluster_cols=FALSE, show_colnames=FALSE,
                annotation_row=dc4, annotation_colors=anno_colors,
                annotation_names_row=FALSE, annotation_legend=FALSE,
                main='protein')

# phosphoprotein dendrogram
# ----------------------------------------------------------------------------------
file <- synGet('syn7437201')
dfr <- read.table(file$path, stringsAsFactor = TRUE, 
                  header=TRUE, row.names='Uniprot_Id', sep=',')
df <- dfr[, 2:37]
dft <- t(df)
rownames(dft) <- colnames(df)
colnames(dft) <- rownames(df)

cell_line_order <- rownames(dft)#[hc$order]
# Get cell line categories
# dc <- as.data.frame(read_csv('cell_lines_categories.csv'))
dc2 <- dc[dc$Sample %in% cell_line_order, ]
# Sort dc2 by cell line order
dc3 <- dc2[match(cell_line_order, dc2$Sample), ]
colors = add.alpha(as.character(dc3$color_map), 0.7)
dc4 <-  as.data.frame(dc3$Receptor_Status)
colnames(dc4) <- "Receptor_status"
rownames(dc4) <- dc3$Sample
rs_colors <- add.alpha(c("red", "#1F78B4","orange2", '#33A02C'), 0.9)
#rs_colors <- add.alpha(c("gray", "#1F78B4","#FB9A99", '#33A02C'), 0.9)
#rs_colors <- c('blue', 'black', 'red', 'green')
names(rs_colors) <- c('NM', 'HER2amp', 'TNBC', 'HR+')
anno_colors <- list(Receptor_status=rs_colors)

gg3 <- pheatmap(dft,
                color=0, cellwidth=0, legend=FALSE,
                cluster_cols=FALSE, show_colnames=FALSE,
                annotation_row=dc4, annotation_colors=anno_colors,
                annotation_names_row=FALSE, annotation_legend=TRUE,
                main='phosphoprotein')
