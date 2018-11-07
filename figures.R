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
dd <- dist(dft, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
hcd <- as.dendrogram(hc)
par(cex=0.9) # Set dendrogram label size
plot(hcd, main='RNASeq')

cell_line_order <- rownames(dft)[hc$order]
# Get cell line categories
dc <- read_csv('cell_lines_categories.cgp1sv')
dc2 <- dc[dc$Sample %in% cell_line_order, ]
# Sort dc2 by cell line order
dc3 <- dc2[match(cell_line_order, dc2$Sample), ]
colors = add.alpha(as.character(dc3$color_map), 0.7)
colored_bars(colors=colors, dend=hcd, rowLabels="")

#------------------
p1 <- ggplot(dend_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3) +ylim(-100, 450) + theme_void() + scale_x_continuous(expand=c(0.02,0.02))

p2 <- ggplot(dc3, aes(x=Sample, y=0.2, fill=Receptor_Status)) +geom_tile() + 
  theme(legend.position="none") + theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none")

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

grid.arrange(gp1, gp2, ncol=1,heights=c(9/10,1/10))
