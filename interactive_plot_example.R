library(shiny)
library(ggplot2)
library(ggrepel)


ui <- basicPage(
  plotOutput("plot1", click = "plot_click"),
  plotOutput("plot2")
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    ggplot(data=mtcars) + 
      geom_point(aes(x=wt, y=mpg)) +
      geom_label_repel(data=mtcars, aes(x=wt, y=mpg, label=model)) +
      theme_bw() + 
      xlab('WT') + ylab('MPG')
  })
  
  output$plot2 <- renderPlot({
    # With base graphics, need to tell it what the x and y variables are.
    ds <- nearPoints(mtcars, input$plot_click, xvar = "wt", yvar = "mpg") 
    if (dim(ds)[1] >= 1){
      plot(ds$wt, ds$mpg)
    }
    
    # nearPoints() also works with hover and dblclick events
  })
}

shinyApp(ui, server)
