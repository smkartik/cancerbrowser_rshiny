library(shinydashboard)
library(shiny)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringr)
library(tidyr)


# Load data files
# ---------------
#gr_metrics <- read.csv('data/gr_metrics_all.csv', row.names = 'X')
#gr_values <- read.csv('data/grvalues_cellcycle_DDD.csv')
#metadata <- read.csv('data/metadata.csv')

#dfm = metadata[metadata$molecular_subtype != 'Bridge' & metadata$molecular_subtype != 'ovarian', ]
#dfm$cell_line <- sapply(dfm$cell_line,  function(x){unlist(str_split(x, '[ (]'))[1]})
#dfm$cell_line <- sapply(dfm$cell_line, function(x){str_replace(x, '-', '')})
#dfm$cell_line <- sapply(dfm$cell_line, function(x){str_replace(x, '-', '')})
#dfm$cell_line <- sapply(dfm$cell_line, function(x){str_replace(x, 'Hs578T', 'HS578T')})
#dfm$cell_line <- sapply(dfm$cell_line, function(x){str_replace(x, 'hME1', 'HTERT_HME1')})
#dfm <- dfm[c('cell_line', 'molecular_subtype', 'receptor_status')]
#dfm <- rbind(dfm, data.frame(cell_line='MCF10A (GM)', molecular_subtype='Non malignant, Basal',
#                             receptor_status='NM'))
#gr_metric <- merge(gr_metrics, dfm, by.x='cell_line', by.y='cell_line')
#gr_metric <- gr_metric %>% distinct()
#gr_metric$cell_line <- as.character(gr_metric$cell_line)

#gr_values <- merge(gr_values, dfm, by.x='cell_line', by.y='cell_line')
#gr_values$cell_line <- as.character(gr_values$cell_line)

load('data.RData')


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
names(rs_colors) <- c('NM', 'Her2amp', 'TNBC', 'HR+')
pie_colors <- add.alpha(c('#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69'), 0.85)
names(pie_colors) <- c('beyondG2', 'S_dropout', 'subG1', 'G1', 'G2', 'M', 'S')

# Genetic functions
# -----------------
sigmoidal_fit <- function(drug, cell_line){    
  grc <- gr_values[gr_values$agent == drug & 
                     gr_values$cell_line == cell_line, ]
  grc$log10_conc <- log10(grc$concentration)
  
  grm <- gr_metric[gr_metric$agent == drug & 
                     gr_metric$cell_line == cell_line, ]

  cc <- grc$log10_conc
  param_df = tibble::tribble(
    ~parameter,                 ~lower,     ~prior,    ~upper,
    "GRinf",                       -1,        0.1,         1, 
    "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
    "h_GR",                         0.1,          2,         5
  )
  
  opfct_sig = function(x, p) {
    p[1] + (1 - p[1])/(1 + (10^x / (10^p[2])) ^ p[3])
  }
  
  sum_square_sig = function(x, y, p) {sum((y - opfct_sig(x, p))^2)}
  
  startVec = param_df$prior
  psVec <- abs(startVec)
  psVec[psVec < 1e-4] <- 1
  
  xx <- grc$concentration
  yy <- grc$GRvalue
    
  fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),                    
                                   function(p, x, y) sum_square_sig(x = log10(xx), y = yy, p = p),
                                   hessian = TRUE, method = "L-BFGS-B",
                                   lower = param_df$lower, upper = param_df$upper)))
  
  minc <- min(grc$concentration)
  maxc <- max(grc$concentration)
  len = (log10(maxc) - log10(minc))*20
  concentration = 10^(seq(log10(minc) - 1, log10(maxc) + 1, length.out = len))
  
  sig_fit <- unlist(sapply(concentration,
                          function(x){fit$par[1] + (1-fit$par[1]) / (1 + (x / (10** fit$par[2])) ** fit$par[3])}))
                           
  dr <- data.frame(cc=concentration, yfit=sig_fit, log10_conc=log10(concentration))
  return(dr)
}

biphasic_fit <- function(drug, cell_line){    
  grc <- gr_values[gr_values$agent == drug & 
                     gr_values$cell_line == cell_line, ]
  grc$log10_conc <- log10(grc$concentration)
  
  grm <- gr_metric[gr_metric$agent == drug & 
                     gr_metric$cell_line == cell_line, ]
  
  conc <- grc$concentration
  cc <- grc$log10_conc
  
  ec50_low = log10(max(c(min(unlist(conc)) * 1e-4, 1e-7))) 
  ec50_high = log10(min(c(max(unlist(conc)) * 1e2, 1e2))) 
  
  param_df = tibble::tribble(
    ~parameter,                 ~lower,     ~prior,             ~upper,
    "GRinf_1",                     -.05,        0.1,                  1, 
    "log10_GEC50_1",          ec50_low, median(unlist(cc)),   log10(1),
    "h_GR_1",                       0.025,          2,                  5,
    "GRinf_2",                       -1,       -0.1,                0.5, 
    "log10_GEC50_2",        log10(0.3),   log10(1),          ec50_high,
    "h_GR_2",                       0.025,          2,                  5
)
  
  opfct_biphasic = function(x, p) {
    term1 = 1 + (p[1] + (1 - p[1])/(1 + (10^x / (10^p[2])) ^ p[3]))
    term2 = 1 + (p[4] + (1 - p[4])/(1 + (10^x / (10^p[5])) ^ p[6]))
    2^( 0.5*( log2(term1) + log2(term2) ) ) - 1
  }
  
  sum_square_biphasic = function(x, y, p) {sum((y - opfct_biphasic(x, p))^2)}
  
  startVec = param_df$prior
  psVec <- abs(startVec)
  psVec[psVec < 1e-4] <- 1
  
  xx <- grc$concentration
  yy <- grc$GRvalue
  
  fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),                    
                                   function(p, x, y) sum_square_biphasic(x = log10(xx), y = yy, p = p),
                                   hessian = TRUE, method = "L-BFGS-B",
                                   lower = param_df$lower, upper = param_df$upper)))
  
  minc <- min(grc$concentration)
  maxc <- max(grc$concentration)
  len = (log10(maxc) - log10(minc))*20
  concentration = 10^(seq(log10(minc) - 1, log10(maxc) + 1, length.out = len))
  
  create_biphasic_data <- function(conc){
    term1 = 1 + (fit$par[1] + (1 - fit$par[1])/(1 + (conc / (10^fit$par[2])) ^ fit$par[3]))
    term2 = 1 + (fit$par[4] + (1 - fit$par[4])/(1 + (conc / (10^fit$par[5])) ^ fit$par[6]))
    biphasic_fit_data <- 2^( 0.5*( log2(term1) + log2(term2) ) ) - 1
    return(biphasic_fit_data)
  }
  
  biphasic_fit <- unlist(sapply(concentration, create_biphasic_data))
  
  dr <- data.frame(cc=concentration, yfit=biphasic_fit, log10_conc=log10(concentration))
  return(dr)
}


make_cc_long_table <- function(agent, cell_line){
  grc <- gr_values[gr_values$agent == agent & gr_values$cell_line == cell_line, ]
  cc_cols <- c("G1", "S", "G2", "M", "S_dropout","subG1", "beyondG2")
  long <- gather(grc, phase, fraction, cc_cols, factor_key=T)
  long$log10_conc <- log10(long$concentration)
  long$phase_count <- mapply(function(x, y){x * y}, long$cell_count, long$fraction)
  return(long)
}
 

body <- dashboardBody(
  # Select Drug
  fluidRow(
    column(3,
           box(status='primary', width=12,
               selectInput("agent", "Select drug of interest:",
                           choices=as.character(unique(gr_metric$agent)),
                           selected='Flavopiridol')
               )
           ),
    column(5,
           box(status='info', width = 12,
               title='Project information',
               verbatimTextOutput("project_info")
               #tags$head(tafs$style("project_info{font-size:16px,}")),
               #uiOutput("project_info")
           )),
    column(4, img(src='LSP_Horiz-Logo.pdf'))
  ),
  
  fluidRow(
  box(title=textOutput('boxtitle'), width=8,
    plotOutput("gr_metrics", click = "plot_click", 
               hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce")),
    uiOutput("hover_info")
  )
  ),
  
  fluidRow(
    box(width=4,
        fluidPage(
          fluidRow(
            column(12, selectInput('fit_option', '',
                                   choices=c('sigmoidal', 'biphasic'),
                                   selected='sigmoidal')
            )
          ),
          fluidRow(
            column(12, plotOutput("GRdose_response"))
          )
    )),
    box(width=4,
        fluidPage(
          fluidRow(
            column(12, selectInput('cc_option', '',
                                   choices=c('fraction', 'cell count'),
                                   selected='fraction')
                   )
          ),
          fluidRow(
            column(12, plotOutput("cellcycle_response"))
          )
        )
    )
  )
)

server <- function(input, output) {
  
  output$boxtitle <- renderText({
    paste(input$agent, ' GR metrics across breast cancer cell lines')
  })
  
  output$project_info <- renderText({
    paste('The datasets visualized in this app were collected and analyzed as part of the "The Library of Integrated Network-Based Cellular Signatures" (LINCS) programe and funded by the NIH Common Fund program (U54 grant HL127365).\nThe datasets are currently available under a Creative Commons License CC BY 4.0')
  })
  
  output$gr_metrics <- renderPlot({
    ds <- gr_metric[gr_metric$agent == input$agent, ]
    ds$cell_line <- factor(ds$cell_line, levels = ds$cell_line[order(-ds$GR_AOC)])
    p <- ggplot(data=ds) + 
      geom_point(aes_string(x='cell_line', y='GR_AOC', 
                            color='receptor_status'),
                 size=5, alpha=0.7) +
      # geom_label_repel(data=ds, aes_string(x='rank', y='GR_AOC', label='cell_line')) +
      scale_color_manual(name='Receptor status', values = rs_colors) +#, name='Receptor status') +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle=90, size=10, face='bold', hjust=0.95),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=16, face='bold'),
        axis.title.y =element_text(size=16, face='bold'),
        #plot.title =element_text(size=14, face='bold', hjust=.05),
        legend.title = element_text(colour="black", size=14, face='bold'),
        legend.text = element_text(colour="black", size=12)
      ) +
      xlab('') + ylab('GR AOC')
    return(p)
    
    #plot(ds$GR_AOC, ds$mpg)
  })
  
  output$hover_info <- renderUI({
    ds <- gr_metric[gr_metric$agent == input$agent, ]
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
  
  
  output$GRdose_response <- renderPlot({
    ds <- gr_metric[gr_metric$agent == input$agent, ]
    ds$cell_line <- factor(ds$cell_line, levels = ds$cell_line[order(-ds$GR_AOC)])
    
    point <- nearPoints(ds, input$plot_click, xvar='cell_line', yvar='GR_AOC')
    if (dim(point)[1] >=1){
      cell_line <- as.character(point$cell_line)
      gr50 <- point$GR50
      gr_max <- point$GRmax
      
      grc <- gr_values[gr_values$agent == input$agent & 
                         gr_values$cell_line == cell_line, ]
      grc$log10_conc <- log10(grc$concentration)
      
      if (input$fit_option == 'sigmoidal'){
        dr <- sigmoidal_fit(input$agent, cell_line)
      } else if (input$fit_option == 'biphasic'){
        dr <- biphasic_fit(input$agent, cell_line)
      }
      
      p <- ggplot() + geom_point(data=grc, aes_string(x='log10_conc', y='GRvalue'), color='blue', size=2, alpha=0.5) +
        geom_line(data=dr, aes_string(x='log10_conc', y='yfit'), color='black', size=2, alpha=0.5) +
        xlab(expression(paste('log10 ', mu, 'M'))) + ylab('GR value') +
        ggtitle(cell_line) +
        ylim(-1, 1.1) +
        #theme_bw() + 
        theme(
          panel.background=element_blank(),
          axis.title=element_text(size=16, face="bold"),
          #legend.title=element_blank(),
          axis.line=element_line(),
          axis.ticks=element_line(),
          axis.text.x = element_text(size=12, face='bold'),
          axis.text.y = element_text(size=12, face='bold'),
          axis.text=element_text(size=14, face='bold')
        ) +
        #coord_equal() +
        geom_hline(yintercept=1, alpha=0.5) + 
        geom_hline(yintercept=0, alpha=0.5) +
        geom_hline(yintercept=-1, alpha=0.5) +
        #geom_hline(yintercept=gr_max, alpha=0.5) +
        geom_segment(aes(x = log10(gr50), y = -1, xend = log10(gr50), yend = -0.8), 
                     color='purple', alpha=0.7, size=2)
      return(p)
    }
  },width=300, height=300)
  
  
  output$cellcycle_response <- renderPlot({
    ds <- gr_metric[gr_metric$agent == input$agent, ]
    ds$cell_line <- factor(ds$cell_line, levels = ds$cell_line[order(-ds$GR_AOC)])
   
    point <- nearPoints(ds, input$plot_click, xvar='cell_line', yvar='GR_AOC')
    if (dim(point)[1] >=1){
      cell_line <- as.character(point$cell_line)
      cl <- make_cc_long_table(input$agent, cell_line)
      if (input$cc_option == 'fraction'){
        pos='fill'
        ylabel='cell cycle fraction'} else if (
        input$cc_option == 'cell count'){
          pos='stack'
          ylabel='# cells'}
      p <- ggplot(cl, aes(fill=phase, y=phase_count, x=log10_conc)) + 
        geom_bar(position=pos, stat='identity') +
        scale_fill_manual(name='cell cycle phase', values=pie_colors) +
        guides( fill=guide_legend(title.position ="top")) +
        xlab('') + ylab(ylabel) +
        theme(
          panel.background=element_blank(),
          axis.title=element_text(size=16, face="bold"),
        #legend.title=element_blank(),
          axis.line=element_line(),
          axis.ticks=element_line(),
          axis.text.x = element_text(size=12, face='bold'),
          axis.text.y = element_text(size=12, face='bold'),
          axis.text=element_text(size=14, face='bold'),
          legend.title = element_text(colour="black", size=14, face='bold'),
          legend.text = element_text(colour="black", size=12),
          legend.position='bottom'
      )
      return(p)
    }
  },width=350, height=350)
}


shinyApp(
  ui = dashboardPage(
    dashboardHeader(title = "LINCS Cancer Browser", titleWidth = 400),
    dashboardSidebar(disable = TRUE),
    body
  ),
  server = server
)