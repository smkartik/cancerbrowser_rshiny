dfi <- read.csv('pi3k_feaure_importances.csv', stringsAsFactors = TRUE,header=TRUE)

dfi1 <- dfi[dfi$drug == 'Alpelisib/BYL719', ]
dfi2 <- dfi[dfi$drug == 'Pictilisib/GDC0941', ]
dfi3 <- dfi[dfi$drug == 'Taselisib/GDC0032', ]
dfi4 <- dfi[dfi$drug == 'Buparlisib/NVP-BKM120', ]

label_colors <- ifelse(dfi1$validated == 'yes', "red", "black")
f3a <- ggbarplot(dfi1, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab('feature importance') + xlab('') +
  ggtitle("Alpelisib/BYL719 (0.81)") +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.text.x = element_text(colour = label_colors),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

label_colors <- ifelse(dfi2$validated == 'yes', "red", "black")
f3b <- ggbarplot(dfi2, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab(' ') + xlab('') +
  ggtitle("Pictilisib/GDC0941 (0.77)") +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.text.x = element_text(colour = label_colors),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

label_colors <- ifelse(dfi3$validated == 'yes', "red", "black")
f3c <- ggbarplot(dfi3, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab('feature importance') + xlab('') +
  ggtitle("Taselisib/GDC0032 (0.71)") +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.text.x = element_text(colour = label_colors),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 

label_colors <- ifelse(dfi4$validated == 'yes', "red", "black")
f3d <- ggbarplot(dfi4, x='features', y='importance', 
                 x.text.angle=90, sort.val='desc') +
  geom_point(aes(y=imp2, color=spearman), size=2) + scale_colour_gradient2() +
  ylab(' ') + xlab('') +
  ggtitle("Buparlisib/NVP-BKM120 (0.88)") +
  theme( axis.title.x = element_text( size = 10, face = "bold" ),
         axis.text.x = element_text(colour = label_colors),
         axis.title.y = element_text( size = 10, face = "bold" ),
         plot.title = element_text( size=12, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 
r1 <- ggarrange(f3a, f3b, nrow=1, ncol=2)
r2 <- ggarrange(f3c, f3d, nrow=1, ncol=2)
f3 <-ggarrange(r1, r2, nrow=2, ncol=1)



  