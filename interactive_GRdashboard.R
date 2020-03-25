library(shinydashboard)
library(shiny)
library(ggplot2)
library(dplyr)
library(ggrepel)


# Load data files
# ---------------
gr_metric <- read.csv(file = "allBatches_grmetrics.csv",  #grmetrics_pooled.csv",
                      stringsAsFactors = TRUE, header = TRUE)
gr_values <- read.csv(file = "allBatches_GRvalues_mean.csv", #"grvalues_pooled.csv",
                      stringsAsFactors = TRUE, header = TRUE)


# Define color palette
# ---------------------
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

body <- dashboardBody(
  fluidRow(
  box(title='Flavopiridol GR metrics', width=8,
    plotOutput("gr_metrics", click = "plot_click", 
               hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce")),
    uiOutput("hover_info")
  )
  )
)

server <- function(input, output) {
  output$gr_metrics <- renderPlot({
    ds$cell_line <- factor(ds$cell_line, levels = ds$cell_line[order(-ds$GR_AOC)])
    ggplot(data=ds) + 
      geom_point(aes_string(x='cell_line', y='GR_AOC', 
                            fill='receptor_status', color='receptor_status'),
                 size=5, alpha=0.7) +
      # geom_label_repel(data=ds, aes_string(x='rank', y='GR_AOC', label='cell_line')) +
      scale_color_manual(values = rs_colors) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle=90, size=10, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=16, face='bold'),
        axis.title.y =element_text(size=16, face='bold'),
        plot.title =element_text(size=14, face='bold', hjust=0.5)
      ) +
      xlab('') + ylab('GR AOC')# +
      #ggtitle(sprintf('%s GR metrics', unique(ds$agent)))
    
    #plot(ds$GR_AOC, ds$mpg)
  })
  
  output$hover_info <- renderUI({
    ds$cell_line <- factor(ds$cell_line, levels = ds$cell_line[order(-ds$GR_AOC)])
    hover <- input$plot_hover
    point <- nearPoints(ds, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    
    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    right_pct <- (hover$x - hover$domain$right) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    right_px <- hover$range$right + right_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "right:", 0.5 * right_pct, "px; top:", top_pct, "px;")
    
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> cell line: </b>", point$cell_line, "<br/>",
                    "<b> GR AOC: </b>",sprintf("%.02f",  point$GR_AOC), "<br/>",
                    "<b> Receptor_status: </b>", point$receptor_status, "<br/>")))
    )
  })
}

shinyApp(
  ui = dashboardPage(
    dashboardHeader(title = 'LINCS Cancer Browser', titleWidth = 400),
    dashboardSidebar(disable = TRUE),
    body
  ),
  server = server
)