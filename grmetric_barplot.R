library(shiny)
library(ggplot2)
library(scatterpie)

# have on OWL dataset behind this drug response
# ----------------------------------------------

#gr_metric <- read.csv(file = "grmetrics_phase1.csv",
#                      stringsAsFactors = TRUE, header = TRUE)
#gr_values <- read.csv(file = "grvalues_phase1.csv",
#                      stringsAsFactors = TRUE, header = TRUE)

gr_metric <- read.csv(file = "allBatches_grmetrics.csv",  #grmetrics_pooled.csv",
                      stringsAsFactors = TRUE, header = TRUE)
gr_values <- read.csv(file = "allBatches_GRvalues_mean.csv", #"grvalues_pooled.csv",
                      stringsAsFactors = TRUE, header = TRUE)
#gr_cellcycle <- read.csv(file = "alldata.csv",
#                         stringsAsFactors = TRUE, header = TRUE)


#---------------------------------------------------------------
# define function that adds transperancy to colors
# Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
rs_colors <- add.alpha(c('#c85f61', '#4dab9b', '#ca7832', '#7da245'), 0.7)
names(rs_colors) <- c('NM', 'HER2amp', 'TNBC', 'HR+')
pie_colors <- add.alpha(c('#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69'), 0.7)
names(pie_colors) <- c('beyondG2', 'S_dropout', 'subG1', 'G1', 'G2', 'M', 'S') 


# -------------------------------------------------
# Define UI for application that draws dose repsone and GR metric histogram
ui <- fluidPage(
  
  # Give the page a title
  titlePanel("GR dose response & metric across BRCA cell lines"),
  
  # Generate a row with a sidebar
  sidebarLayout(
    
    # Define the sidebar with three inputs: drug, GR metric, and receptor status
    sidebarPanel(
      selectInput("drug", "Select drug treatment:",
                  choices=as.character(unique(gr_metric$agent)),
                  selected='Flavopiridol'),
      selectInput("cell_line", "Select cell line:",
                  choices=as.character(unique(gr_values$cell_line)),
                  selected='MCF7'),
      selectInput("metric", "Choice of GR metric:",
                  choices=c("GR_AOC", "GRmax", "GR50"), #, "GEC50", "GRinf"),
                  selected='GR_AOC'),
      checkboxGroupInput("receptor", "Select Receptor Status:", 
                         choices = list("TNBC" = "TNBC", 
                                        "Non-Malignant (NM)" = "NM", 
                                        "HR+" = "HR+",
                                        "HER2amp" = "HER2amp"),
                         selected = c("TNBC", "NM", "HR+", "HER2amp")),
      hr(),
      helpText("Summary GR metric for the selected drug across 39 breast cancer cell lines")
    ),
    
    # Create a spot for the barplot
    mainPanel(
      #plotOutput("dosePlot"),
      tabPanel("GR metric", plotOutput("phonePlot")),
      tabPanel("GR curves and cell cycle distribution", plotOutput("ccPlot"))
    )
  )
)

# ----------------------------------------------
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Fill in the spot we created for a plot
  output$phonePlot <- renderPlot({
    
    # Subset dataframe base on selected receptor status
    gr_metric <- gr_metric[gr_metric$receptor_status %in% c(input$receptor), ]
    
    # Subset dataframe based on selected drug/agent 
    gr_metric <- gr_metric[gr_metric$agent == input$drug, ]
    rownames <- gr_metric$cell_line
    gr_metric <- gr_metric[order(gr_metric[input$metric]),]
    
    dimnames(gr_metric) <- list(as.character(gr_metric$cell_line),
                                colnames(gr_metric))
    # Add transperancy to colors
    alpha_colors = add.alpha(as.character(gr_metric$color_map), 0.7)
    
    # Render a barplot
    barplot(gr_metric[,input$metric], 
            main=input$metric,
            ylab="",
            xlab="",
            names.arg=gr_metric$cell_line, las=2,
            col=alpha_colors)
    #horiz=TRUE)
    legend("topleft", inset=.02, title="Receptor status",
           c("TNBC", "NM", "HR+", "HER2amp"), 
           fill=add.alpha(c('#ca7832', '#c85f61', '#7da245', '#4dab9b'), 0.7), 
           horiz=TRUE, cex=0.8)
  })
  
  output$ccPlot <- renderPlot({
    grc <- gr_values[gr_values$agent == input$drug &
                          gr_values$cell_line == input$cell_line, ]
    #grc <- aggregate(. ~ concentration, grc[-2], mean)
    grc$dummyconc <- 0.4*c(1:length(grc$concentration))
    as.integer(grc$dummyconc)
    grc$log10c <- log10(grc$concentration)
    labels <- formatC(grc$log10c)
    #labels <- formatC(grc$concentration, format='e', digits=1)
    #sp <- barplot(grc[, 'GRvalue'],
    #              main='GRvalue',
    #              ylab="",
    #              xlab="",
    #              names.arg=grc$concentration, las=2)
    # Sigmoidal fit function
    sig1 <- smooth.spline(x=grc$dummyconc, y=grc$GRvalue, df=3)
    sig2 <- smooth.spline(x=grc$dummyconc, y=grc$GR_static, df=3)
    sig3 <- smooth.spline(x=grc$dummyconc, y=grc$GR_toxic, df=3)
    cmin <- min(grc$dummyconc)
    cmax <- max(grc$dummyconc)
    xi <- seq(cmin, cmax, by=0.05)
    sig_predict1 <- data.frame(predict(sig1, xi))
    sig_predict2 <- data.frame(predict(sig2, xi))
    sig_predict3 <- data.frame(predict(sig3, xi))
    ggplot() + 
      geom_line(data=sig_predict1, aes(x=x, y=y, group=1), size=0.7, alpha=0.9, color='blue') +
      geom_line(data=sig_predict2, aes(x=x, y=y, group=1), size=0.7, alpha=0.7, color='green') +
      geom_line(data=sig_predict3, aes(x=x, y=y, group=1), size=0.7, alpha=0.7, color='red') +
      geom_scatterpie(data=grc, aes(x=dummyconc, y=GRvalue, r=0.15),
                      cols=c('subG1', 'G1', 'S', 'S_dropout',
                             'G2', 'beyondG2', 'M')) + ylim(-1, 1.5) +
      scale_x_continuous(breaks=grc$dummyconc,labels=labels)+
      #guides( fill=guide_legend(label.position ="top")) +
      scale_fill_manual(values=pie_colors, guide=guide_legend(nrow=1))+
      xlab("concentration log10(uM)") + ylab("GR value") +
      theme(
        axis.title=element_text(size=16, face="bold"),
        legend.title=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(),
        axis.ticks=element_line(),
        axis.text=element_text(size=14, face='bold'),
        #legend.direction="horizontal",
        legend.position="top",
        legend.text = element_text(colour="black", size=10), 
        #                           face="bold"),
        #legend.text.align = 0,
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.spacing.y = unit(0.2, "cm") 
      ) +
      coord_equal() +
      geom_hline(yintercept=1, alpha=0.5) + geom_hline(yintercept=0, alpha=0.5)
    #scale_color_manual('cell cycle fractions', values=pie_colors) +
    
  }, width=400, height=400)
  
  ##  #sp <- ggplot(grc, aes(x=dummyconc, y=GRvalue, 
  ##  #                       group=cell_line))#, color=receptor_status))
  ##  sp <- ggplot() + geom_scatterpie(aes(x=dummyconc, y=GRvalue, group=well, r=0.15),
  ##                            data=grc,
  ##                            cols=c('subG1', 'G1', 'S', 'S_dropout',
  ##                                   'G2', 'beyondG2', 'M'))
  ##  sp + #geom_line() + geom_point() #+ coord_fixed() +
  ##      #scale_color_manual(values=rs_colors) +
  ##      scale_x_continuous(breaks=c(1:9), labels=grc$concentration) +
  ##      ylim(-0.6, 1.2) +
  ##      theme(panel.grid.major = element_blank(), 
  ##           panel.grid.minor = element_blank(),
  ##            panel.background = element_blank(), 
  ##            axis.line = element_line(colour = "black"),
  ##            #aspect.ratio = 0.7,
  ##           text = element_text(size=12, face='bold')) +
  ##  ##    #facet_grid(. ~ title) +
  ##        geom_hline(yintercept=1, alpha=0.5) + geom_hline(yintercept=0, alpha=0.5)
  
  ## }, height = 400)
  
  
  ##   sp <- ggplot() + geom_scatterpie(aes(x=log10c, y=GRvalue, group=well, r=0.15),
  ##                                    data=grc,
  ##                                    cols=c('subG1', 'G1', 'S', 'S_dropout', 'G2', 'beyondG2', 'M'))
  ##   sp + geom_line() + geom_point() + coord_fixed() + #scale_x_continuous(trans='log10') +
  ##     theme(panel.grid.major = element_blank(), 
  ##           panel.grid.minor = element_blank(),
  ##           panel.background = element_blank(), 
  ##           axis.line = element_line(colour = "black"),
  ##           #aspect.ratio = 0.5,
  ##           text = element_text(size=12, face='bold')) +
  ##     #facet_grid(. ~ title) +
  ##     geom_hline(yintercept=1, alpha=0.5) + geom_hline(yintercept=0, alpha=0.5)# + 
  ##     #coord_equal()
  ## }, height = 500)#, width=400)
  
  # All curves----------------------   
  
  output$dosePlot <- renderPlot({
    
    # Subset GR values dataframe based on chosen drug and receptor status
    grv <- gr_values[gr_values$agent == input$drug &
                       gr_values$receptor_status %in% c(input$receptor), ]
    
    grv$title <- input$drug
    # Render a dose response plot
    colormap <- unique(grv[, c("receptor_status", "color_map")])
    sp <- ggplot(grv, aes(x=concentration, y=GRvalue, 
                          group=cell_line, color=receptor_status))
    sp + geom_line() + geom_point() + 
      scale_color_manual(values=rs_colors) +
      scale_x_continuous(trans='log10') +
      coord_cartesian(ylim=c(-0.6, 1.2)) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            aspect.ratio = 0.7,
            text = element_text(size=12, face='bold')) +
      facet_grid(. ~ title) +
      geom_hline(yintercept=1, alpha=0.5) + geom_hline(yintercept=0, alpha=0.5)
    
  }, height = 800, width = 400)
}

# Run the application 
shinyApp(ui = ui, server = server)
