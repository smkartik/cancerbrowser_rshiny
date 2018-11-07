library(shiny)
library(ggplot2)

#gr_metric <- read.csv(file = "grmetrics_phase1.csv",
#                      stringsAsFactors = TRUE, header = TRUE)
#gr_values <- read.csv(file = "grvalues_phase1.csv",
#                      stringsAsFactors = TRUE, header = TRUE)

gr_metric <- read.csv(file = "grmetrics_pooled.csv",
                      stringsAsFactors = TRUE, header = TRUE)
gr_values <- read.csv(file = "grvalues_pooled.csv",
                      stringsAsFactors = TRUE, header = TRUE)


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
                  choices=as.character(unique(gr_metric$agent))),
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
      plotOutput("dosePlot"),
      plotOutput("phonePlot")
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
    
  }, height = 400, width = 400)
}

# Run the application 
shinyApp(ui = ui, server = server)