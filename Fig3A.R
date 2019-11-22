library( Matrix )
library( tidyverse )
library( reshape2 )

## Platform ordering
plord <- c( "RNAseq", "MS", "pMS", "TF_activity", "Kinase_activity" )

## Drug groupings
grp1 <- c( "Abemaciclib", "Palbociclib", "Ribociclib" )
grp2 <- c( "Alpelisib", "Pictilisib", "Taselisib", "Buparlisib" )

## Retrieve the table
X <- read_csv( "auc_alldrugs.csv" ) %>%
  rename( Drug = X1 ) %>%
  arrange( desc(Drug) ) %>%
  mutate( Drug = factor( Drug, Drug ) ) %>%
  gather( Platform, AUC, -Drug ) %>%
  mutate( AUC = pmax( AUC, 0.5 ) ) %>%
  mutate( Platform = factor( Platform, levels = plord ) ) %>%
  mutate( Group = ifelse( Drug %in% grp1, "CDK4/6i", "Other" ) ) %>%
  mutate( Group = ifelse( Drug %in% grp2, "PI3Ki", Group ) ) %>%
  mutate( Group = factor( Group, levels = c("CDK4/6i", "PI3Ki", "Other") ) )

## Clip the headings to avoid overlapping
levels(X$Platform)[4:5] <- c( "TF_Act.", "Kinase_Act." )

## Define plot elements
cl1 <- "#203c52"
cl2 <- "#a5d8ff"
fpal <- colorRampPalette( c(cl1,cl2) )
sg <- guide_legend( override.aes = list( color = fpal(6) ) )
etxt <- function(s) element_text( face="bold", size=s )

## Make a dot plot
gg <- ggplot( X, aes( x = Platform, y = Drug, size = AUC, color = AUC ) ) +
  geom_point() + theme_bw() + ylab( "Drug" ) +
  facet_grid( Group ~ ., scales="free", space="free" ) +
  scale_radius( range = c(2,8), breaks=seq(0.5,1,0.1), limits=c(0.5,1) ) +
  scale_color_gradient( limits=c(0.5,1), low=cl1, high=cl2 ) +
  #+ ggtitle("Predictive perforamce of each dataset") +
  theme( axis.text = etxt(11), axis.title = etxt(12),
         panel.grid.major = element_line( linetype = "dashed" ),
         strip.text.y = element_blank() ) +
  guides( size = sg, color=FALSE )

###--------------importances------------------------
dfi <- read.csv('phosphoms_cdk46_importance.csv', stringsAsFactors = TRUE,header=TRUE)

dfi1 <- dfi[dfi$drug == 'Abemaciclib/LY2835219', ]
dfi2 <- dfi[dfi$drug == 'Palbociclib/PD0332991', ]
dfi3 <- dfi[dfi$drug == 'LEE011/Ribociclib', ]

f3a <- ggbarplot(dfi1, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab('pMS feature importance') + xlab('') +
  #ggtitle(" ") +
  annotate(geom='text', label='AUC = 0.88', x=5, y=max(dfi1$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

f3b <- ggbarplot(dfi2, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab(' ') + xlab('') +
  #ggtitle(" ") +
  annotate(geom='text', label='AUC = 0.92', x=5, y=max(dfi1$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

f3c <- ggbarplot(dfi3, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab(' ') + xlab('') +
  #ggtitle(" ") +
  annotate(geom='text', label='AUC = 0.91', x=5, y=max(dfi1$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

pms <- ggarrange(f3a, f3b, f3c, nrow=1, ncol=3)

# -------feature importance protein -------------
dfi <- read.csv('ms_cdk46_importance.csv', stringsAsFactors = TRUE,header=TRUE)

dfi1 <- dfi[dfi$drug == 'Abemaciclib/LY2835219', ]
dfi2 <- dfi[dfi$drug == 'Palbociclib/PD0332991', ]
dfi3 <- dfi[dfi$drug == 'LEE011/Ribociclib', ]

f3a <- ggbarplot(dfi1, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab('MS feature importance') + xlab('') +
  ggtitle("Abemaciclib/LY2835219") +
  annotate(geom='text', label='AUC = 0.84', x=5, y=max(dfi1$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

f3b <- ggbarplot(dfi2, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab(' ') + xlab('') +
  ggtitle("Palbociclib/PD0332991") +
  annotate(geom='text', label='AUC = 0.87', x=5, y=max(dfi1$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

f3c <- ggbarplot(dfi3, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab(' ') + xlab('') +
  ggtitle("LEE011/Ribociclib") +
  annotate(geom='text', label='AUC = 0.87', x=5, y=max(dfi1$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=10, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' )

ms <- ggarrange(f3a, f3b, f3c, nrow=1, ncol=3)


#---------viper kinase feature importace---------------------------- 
dfi <- read.csv('viper_kinase_cdk46_importance.csv', stringsAsFactors = TRUE,header=TRUE)

dfi1 <- dfi[dfi$drug == 'Abemaciclib/LY2835219', ]
dfi2 <- dfi[dfi$drug == 'Palbociclib/PD0332991', ]
dfi3 <- dfi[dfi$drug == 'LEE011/Ribociclib', ]

f3a <- ggbarplot(dfi1, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab('Kinase activity feature importance') + xlab('') +
  #ggtitle("Abemaciclib/LY2835219 (0.81)") +
  annotate(geom='text', label='AUC = 0.81', x=5, y=max(dfi1$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

f3b <- ggbarplot(dfi2, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab(' ') + xlab('') +
  #ggtitle("Palbociclib/PD0332991 (0.81)") +
  annotate(geom='text', label='AUC = 0.93', x=5, y=max(dfi1$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

f3c <- ggbarplot(dfi3, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab(' ') + xlab('') +
  #ggtitle("LEE011/Ribociclib") +
  annotate(geom='text', label='AUC = 0.94', x=5, y=max(dfi1$importance), vjust=1, fontface='bold') +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' )

kinact <- ggarrange(f3a, f3b, f3c, nrow=1, ncol=3)

# ------- group stuff -------------
importance <- ggarrange(ms, pms, kinact, nrow=3, ncol=1)
importance <-annotate_figure(importance, top = text_grob("Predictors of response to CDK4/6 inhibitors", face = "bold", size = 12))

gg <-annotate_figure(gg, top = text_grob("Predictive performace of each dataset", face = "bold", size = 12))
fig <- ggarrange(gg, importance, nrow=1, ncol=2, widths=c(0.7, 1))
