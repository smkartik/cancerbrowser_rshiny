library(shiny)
library(ggplot2)
library(tidyverse)

dfp <- read.csv(file = "tnbc_molecules_cell.csv",
                stringsAsFactors = TRUE, header = TRUE)
dfp <- dfp[, colnames(dfp) != 'MCF10Arep2']
dc <- as.data.frame(read_csv('cell_lines_categories.csv'))


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
  titlePanel("Protein abundance across BRCA cell lines"),
  
  # Generate a row with a sidebar
  sidebarLayout(
    
    # Define the sidebar with three inputs: drug, GR metric, and receptor status
    sidebarPanel(
      textInput("protein", h3("HGNC symbol"), 
                value = "ESR1"),
      #selectInput("protein", "Select protein:",
      #            choices=as.character(unique(dfp$GENE_NAME))),
      hr(),
      helpText("Protein abundance (molecules/cell) across 35 breast cancer cell lines")
    ),
    
    # Create a spot for the barplot
    mainPanel(
      plotOutput("phonePlot")
    )
  )
)

# ----------------------------------------------
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Fill in the spot we created for a plot
  output$phonePlot <- renderPlot({
  
    # Subset dataframe based on selected drug/agent 
    dfp <- dfp[dfp$GENE_NAME == input$protein, ]
    pseries <- as.data.frame(t(dfp)[-1, ])
    colnames(pseries) <- input$protein
    rownames(dc) <- dc$Sample
    dc2 <- dc[rownames(dc) %in% rownames(pseries), ]
    pseries <- merge(pseries, dc2, by=0, all=T)
    # pseries[, input$protein] <- as.numeric(as.vector(pseries[input$protein]))
    #pseries <- pseries[order(pseries[input$protein]), ]
    rownames(pseries) <- pseries$Sample

    # Add transperancy to colors
    alpha_colors = add.alpha(as.character(pseries$color_map), 0.7)
    
    # Render a barplot
    barplot(as.numeric(as.vector(pseries[, input$protein])), 
            main=input$protein,
            ylab="",
            xlab="",
            names.arg=rownames(pseries), las=2,
            col=alpha_colors)
    #horiz=TRUE)
    legend("topleft", inset=.02, title="Receptor status",
           c("TNBC", "NM", "HR+", "HER2amp"), 
           fill=add.alpha(c('#ca7832', '#c85f61', '#7da245', '#4dab9b'), 0.7), 
           horiz=TRUE, cex=0.8)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
