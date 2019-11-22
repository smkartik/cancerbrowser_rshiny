library(ggpubr)
library(ggrepel)
#library(tcltk)


ms_colors <- c("#0073C2B2", "#EFC000B2", "#868686B2", "#CD534CB2", "#7AA6DCB2")
names(ms_colors) <- c('Basal A', 'Basal B', 'Bridge', 'Luminal', 'Non malignant, Basal')
#---Figure 1a------------------------------------------------------------------
dfr <- read.csv(file='rnaseq_pca.csv', stringsAsFactors = TRUE,header=TRUE)

f1a <- ggplot(dfr, aes(x=pc_1, y=pc_2, color=Molecular_subtype)) + 
  geom_point(size=3, alpha=0.7) + theme_bw() + 
  scale_color_manual(values=ms_colors) +
  #scale_color_jco(name='transcriptional subtype') + 
  xlab('PC 1 (20.43%)') + ylab('PC 2 (9.98%)') +
  ggtitle("transcriptome") +
  geom_label_repel(aes(label = sample),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme( axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

#---Figure 1b------------------------------------------------------------------
dfm <- read.csv(file='ms_pca.csv', stringsAsFactors = TRUE,header=TRUE)

f1b <- ggplot(dfm, aes(x=pc_1, y=pc_2, color=Molecular_subtype)) + 
  geom_point(size=3, alpha=0.7) + theme_bw() +
  scale_color_manual(values=ms_colors) +
  #scale_color_jco(name='transcriptional subtype') +
  xlab('PC 1 (18.61%)') + ylab('PC 2 (10.32%)') +
  ggtitle("proteome") +
  theme( axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position='none') 

#---Figure 1c------------------------------------------------------------------
dfp <- read.csv(file='pMS_pca.csv', stringsAsFactors = TRUE,header=TRUE)

f1c <- ggplot(dfp, aes(x=pc_1, y=pc_2, color=Molecular_subtype)) + 
  geom_point(size=3, alpha=0.7) + theme_bw() +
  scale_color_manual(values=ms_colors) +
  #scale_color_jco(name='transcriptional subtype') +
  xlab('PC 1 (14.57%)') + ylab('PC 2 (9.27%)') +
  ggtitle("phosphoproteome") +
  theme( axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text( size = 10, face = "bold" ),
    axis.title.y = element_text( size = 10, face = "bold" ),
    plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
    panel.spacing = unit(c(1,1,1,1),"cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position='none') 

#---Figure 1d------------------------------------------------------------------
dfg <- read.csv(file='batch1_GR_AOC.csv', stringsAsFactors = TRUE,header=TRUE)

f1d <- ggplot(dfg, aes(x=GR_AOC, y=DrugName, color=Molecular_subtype)) + 
    geom_point(alpha=0.7) + theme_bw() + 
    scale_color_manual(values=ms_colors, name='') +
    #scale_color_jco(name='transcriptional subtype') +
    ylab('') + xlab('GR AOC') + 
    ggtitle("Drug response across 35 breast cancer cell lines") +
    theme( #axis.text.x = element_blank(),
          axis.text.y = element_text( size = 10, face = "bold" ),
          #axis.ticks.x = element_blank(),
          axis.title.x = element_text( size = 12, face = "bold" ),
          axis.title.y = element_text( size = 12, face = "bold" ),
          plot.title = element_text( hjust = 0.5, face = "bold" ),
          panel.spacing = unit(c(1,1,1,1),"cm"),
          legend.position='bottom')  

#----Figure legend --------------------------
text <- paste("Figure 1: Datasets collected",
              "(A) Principal component analyses of baseline expression data across 35 cell lines.",
              "20,840 transcripts measured by RNASeq, ",
              "11,710 proteins measured by TMT-based shotgun proteomics, and",
              "13, 295 phosphopeptides measured in common across all cell lines by TMT-based shotgun phoshoproteomics.",
              "(B) Drug response for each drug across 35 cell lines measured as area over the GR curve (GR_AOC)." , sep = " ")
text.p <- ggparagraph(text = text, size = 11, color = "black")
#----Combine subplots to single figure 1--------
p1 <- ggarrange(f1a, f1b, f1c,
                ncol = 1, nrow = 3)
p1 <-annotate_figure(p1, top = text_grob("baseline expression profiles", face = "bold", size = 12))

f <- ggarrange(p1, f1d, ncol=2, nrow=1, widths=c(0.35, 1)) #, labels=c('A', 'B')) 
ggarrange(f, text.p, ncol=1, nrow=2, heights=c(1, 0.15))
#print(f)
#ggsave("Figure1.pdf")

#x11()
#print(f)
#prompt  <- "hit spacebar to close plots"
#extra   <- "some extra comment"
#capture <- tk_messageBox(message = prompt, detail = extra)
