library( Matrix )
library( tidyverse )
library( reshape2 )

## Platform ordering
plord <- c( "RNAseq", "MS", "pMS", "Kinase_activity" )

## Drug groupings
cdk12 <- c( "R0-3306", "BMS-265246")
cdk46 <- c("Abemaciclib", "Palbociclib", "LEE011", "BSJ-03-123", "BSJ-03-124", "Flavopiridol")
cdk79 <- c("YKL-5-124", "SY-1365", "ZZ1-33B")
cdk_high <- c("THZ1", "THZ531", "BSJ-01-175", "MFH-2-90", "E17", "FMF-04-107-2")
#grp2 <- c( "Alpelisib", "Pictilisib", "Taselisib", "Buparlisib" )

## Retrieve the table
X <- read_csv( "auc_cdks.csv" ) %>%
  rename( Drug = drug ) %>%
  arrange( desc(Drug) ) %>%
  mutate( Drug = factor( Drug, Drug ) ) %>%
  gather( Platform, AUC, -Drug ) %>%
  mutate( AUC = pmax( AUC, 0.5 ) ) %>%
  mutate( Platform = factor( Platform, levels = plord ) ) %>%
  mutate( Group = ifelse( Drug %in% cdk12, "CDK1/2i", "Other" ) ) %>%
  mutate( Group = ifelse( Drug %in% cdk46, "CDK4/6i", Group ) ) %>%
  mutate( Group = ifelse( Drug %in% cdk79, "CDK7/9i", Group ) ) %>%
  mutate( Group = ifelse( Drug %in% cdk_high, "CDK12/13/14i", Group ) ) %>%
  mutate( Group = factor( Group, levels = c("CDK1/2i", "CDK4/6i", "CDK7/9i", "CDK12/13/14i", "Other") ) )

#Plot only YKL
X_ykl <- X[X$Drug == 'YKL-5-124',]
## Clip the headings to avoid overlapping
levels(X$Platform)[4] <- c( "Kinase_Act." )

## Define plot elements
cl1 <- "#203c52"
cl2 <- "#a5d8ff"
fpal <- colorRampPalette( c(cl1,cl2) )
sg <- guide_legend( override.aes = list( color = fpal(6) ) )
etxt <- function(s) element_text( face="bold", size=s )

## Make a dot plot
ykl_aoc <- ggplot( X_ykl, aes( x = Platform, y = Drug, size = AUC, color = AUC ) ) +
  geom_point() + theme_bw() + ylab( "Drug" ) +
  facet_grid( Group ~ ., scales="free", space="free" ) +
  scale_radius( range = c(2,8), breaks=seq(0.5,1,0.1), limits=c(0.5,1) ) +
  scale_color_gradient( limits=c(0.5,1), low=cl1, high=cl2 ) +
  #+ ggtitle("Predictive perforamce of each dataset") +
  theme( axis.text = etxt(11), axis.title = etxt(12),
         panel.grid.major = element_line( linetype = "dashed" ),
         strip.text.y = element_blank() ) +
  guides( size = sg, color=FALSE )


## Importances
dfi <- read.csv('ykl_kinase.csv', stringsAsFactors = TRUE,header=TRUE)


f3a <- ggbarplot(dfi, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab('Kinase activity feature importance') + xlab('') +
  ggtitle("YKL-5-124") +
  annotate(geom='text', label='AUC = 0.77', x=5, y=max(dfi$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

dfi <- read.csv('ykl_rna_imp.csv', stringsAsFactors = TRUE,header=TRUE)


f3b <- ggbarplot(dfi, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman_rho), size=2) + scale_colour_gradient2() +
  ylab('RNA feature importance') + xlab('') +
  ggtitle("YKL-5-124") +
  annotate(geom='text', label='AUC = 0.79', x=5, y=max(dfi$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

f3ab <- ggarrange(f3b, f3a, nrow = 1, ncol=2)
f3 <- ggarrange(ykl_aoc, f3ab, nrow=2, ncol=1,heights = c(0.3, 1))
print(f3)


####
dfi <- read.csv('azd1775_kinase.csv', stringsAsFactors = TRUE,header=TRUE)


f3c <- ggbarplot(dfi, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab('Kinase activity feature importance') + xlab('') +
  ggtitle("AZD1775") +
  annotate(geom='text', label='AUC = 0.58', x=5, y=max(dfi$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

dfi <- read.csv('azd1775_rna_imp.csv', stringsAsFactors = TRUE,header=TRUE)


f3d <- ggbarplot(dfi, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman_rho), size=2) + scale_colour_gradient2() +
  ylab('RNA feature importance') + xlab('') +
  ggtitle("AZD1775") +
  annotate(geom='text', label='AUC = 0.66', x=5, y=max(dfi$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 





