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
#control_drugs <- c('Actinomycin D', 'DMSO', 'GSK2126458', 
#                   'Paclitaxel_pos', 'Storausporin', 'Vincristin' )
#gr_metric <- gr_metric[!gr_metric$agent %in% control_drugs, ]

#cm <- read.csv('compound_metadata.csv')



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
sigmoidal_fit <- function(drug, cell_line, measure='GRvalue', low_high=NULL){    
  grc <- gr_values[gr_values$agent == drug & 
                     gr_values$cell_line == cell_line, ]
  grc$log10_conc <- log10(grc$concentration)
  
  grm <- gr_metric[gr_metric$agent == drug & 
                     gr_metric$cell_line == cell_line, ]
  
  conc <- grc$concentration
  cc <- grc$log10_conc
  
  ec50_low = log10(max(c(min(unlist(conc)) * 1e-4, 1e-7))) 
  ec50_high = log10(min(c(max(unlist(conc)) * 1e2, 1e2))) 
  
  p_sig_GR = tibble::tribble(
    ~parameter,                 ~lower,     ~prior,    ~upper,
    "GRinf",                       -1,        0.1,         1, 
    "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
    "h_GR",                         0.1,          2,         5
  )
  
  p_sig_GR_static =  tibble::tribble(
    ~parameter,                 ~lower,     ~prior,    ~upper,
    "GRinf",                       0,        0.5,         1, 
    "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
    "h_GR",                         0.1,          2,         5
  )
  
  p_sig_GR_toxic = tibble::tribble(
    ~parameter,                 ~lower,     ~prior,    ~upper,
    "GRinf",                       -1,        -0.5,         0, 
    "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
    "h_GR",                         0.1,          2,         5
  )
  
  p_sig_GR_low = tibble::tribble(
    ~parameter,                 ~lower,     ~prior,    ~upper,
    "GRinf",                     -.05,        0.1,            1, 
    "log10_GEC50",          ec50_low, median(unlist(cc)),   log10(1),
    "h_GR",                       0.025,          2,          5
  )
  
  ### priors and upper/lower bounds for high/late sigmoid fit (GR)
  p_sig_GR_high = tibble::tribble(
    ~parameter,                 ~lower,     ~prior,    ~upper,
    "GRinf",                       -1,       -0.1,        0, 
    "log10_GEC50",            log10(0.3),   log10(1),   ec50_high,
    "h_GR",                       0.025,          2,      5
  )
  
  
  param_list <- list(p_sig_GR, p_sig_GR_static, p_sig_GR_toxic, p_sig_GR_low, p_sig_GR_toxic)
  names(param_list) <- c('GRvalue', 'GR_static', 'GR_toxic', 'GRvalue_low', 'GRvalue_high')
  
  if (!is.null(low_high)){
    param_choice = sprintf("%s_%s", measure, low_high)
  } else {param_choice=measure}
  param_df <- param_list[param_choice][[1]]
  
  
  opfct_sig = function(x, p, bound=1) {
    p[1] + (bound - p[1])/(1 + (10^x / (10^p[2])) ^ p[3])
  }
  
  sum_square_sig = function(x, y, p, bound=1) {sum((y - opfct_sig(x, p, bound))^2)}
  
  startVec = param_df$prior
  psVec <- abs(startVec)
  psVec[psVec < 1e-4] <- 1
  
  xx <- grc$concentration
  yy <- grc[measure]
  
  if (measure == 'GR_toxic'){bound=0} else {bound=1}  
  fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),                    
                                   function(p, x, y) sum_square_sig(x = log10(xx), y = yy, p = p, bound=bound),
                                   hessian = TRUE, method = "L-BFGS-B",
                                   lower = param_df$lower, upper = param_df$upper)))
  
  minc <- min(grc$concentration)
  maxc <- max(grc$concentration)
  len = (log10(maxc) - log10(minc))*20
  concentration = 10^(seq(log10(minc) - 1, log10(maxc) + 1, length.out = len))
  
  sig_fit <- unlist(sapply(concentration,
                          function(x){fit$par[1] + (bound-fit$par[1]) / (1 + (x / (10** fit$par[2])) ** fit$par[3])}))
                           
  dr <- data.frame(cc=concentration, yfit=sig_fit, log10_conc=log10(concentration))
  return(dr)
}

biphasic_fit <- function(drug, cell_line, measure='GRvalue'){    
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
  yy <- grc[measure]
  
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

# Taken from http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


make_cc_long_table <- function(agent, cell_line){
  grc <- gr_values[gr_values$agent == agent & gr_values$cell_line == cell_line, ]
  cc_cols <- c("G1", "S", "G2", "M", "S_dropout","subG1", "beyondG2")
  long <- gather(grc, phase, fraction, cc_cols, factor_key=T)
  long$log10_conc <- log10(long$concentration)
  long$phase_count <- mapply(function(x, y){x * y}, long$cell_count, long$fraction)
  return(long)
}
 
header = dashboardHeader(
  # from https://stackoverflow.com/questions/48978648/align-header-elements-in-shiny-dashboard
  tags$li(
    class = "dropdown",
    tags$style(
      ".main-header {max-height: 120px;
      font-size:36px; 
      font-weight:bold; 
      line-height:100px;}"),
    tags$style(
      ".main-header .logo {height: 120px;
      font-size:36px; 
      font-weight:bold; 
      line-height:100px;align}"
    )
    ),
  title = HTML(
    "<div style = 'background-color:white; vertical-align:middle'>
    <img src = 'LSP_Horiz-Logo.jpg' align = 'left' height = '100px'>
    HMS LINCS Breast Cancer Browser
    </div>"),
  titleWidth = "92%"
)
#<img src = 'http://www.clipartbest.com/cliparts/nTX/8nj/nTX8njyEc.jpeg' align = 'left' height = '100px'>
#{{imageOutput('lsp_logo')}}
#<img src = '', align = 'left' height = '100px'>

body <- dashboardBody(
  fluidRow(
    box(status='info', width = 12,
        #title=Project information',
        #verbatimTextOutput("project_info"),
        uiOutput('project_info'),
        tags$head(tags$style("#project_info{font-size:12px; font-style:italic; max-height:100px; background-color: solidwhite; text-align:center;}"))
        #uiOutput("project_info")
    )
  ),
  # Select Drug
  fluidRow(
    column(4,
           box(status='primary', width=12,
               selectInput("agent", "Select drug of interest:",
                           choices=as.character(unique(gr_metric$agent)),
                           selected='Flavopiridol')
               )
           ),
    column(4,
           box(title='Nominal targets', width=12,
               verbatimTextOutput('nominal_targets'))
),
    column(4,
           box(title='HMS LINCS small molecules DB', width=12,
           uiOutput('hmsl_id'),
           tags$head(tags$style("#hmsl_id{font-size:20px;}")))) 
  #column(4, imageOutput('lsp_logo'))
  ),
  
  fluidRow(
  box(title=textOutput('boxtitle'), width=12,
    plotOutput("gr_metrics", click = "plot_click", 
               hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce")),
    #uiOutput("hover_info")
    #uiOutput("mouse_pointer"),
    tags$head(tags$style("#gr_metrics{cursor:pointer;}"))
    #tags$head(HTML('<style>{{uiOutput("mouse_pointer")}}</style>'))
    
    #gr_metrics{cursor:{{
  
  )
  ),
  
  fluidRow(
    box(width=4, title='GR dose response',
        fluidPage(
          fluidRow(
            column(6, selectInput('fit_option', 'Curve fit:',
                                   choices=c('sigmoidal', 'biphasic'),
                                   selected='sigmoidal')
            ),
            column(6, 
                   conditionalPanel(
                     condition = "input.fit_option == 'sigmoidal'",
                     checkboxInput("st_tox", "Show GR static & toxic response", FALSE))
          )),
          fluidRow(
            column(12, plotOutput("GRdose_response"))
          )
    )),
    box(width=4, title='Cell viability',
        plotOutput("incfrac_dead")
    ),
    box(width=4, title='Cell cycle distribution',
        fluidPage(
          fluidRow(
            column(12, selectInput('cc_option', 'Display as:',
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

server <- function(input, output, session) {
  
  output$boxtitle <- renderText({
    drug_label = str_replace(input$agent, '_', '/')
    paste(drug_label, ' GR metrics across breast cancer cell lines')
  })
  
  output$project_info <- renderUI({
    url <- a("HMS LINCS Terms of Use", href="http://lincs.hms.harvard.edu/terms/", target="_blank")
    tagList(HTML(paste('<b>This open-access Shiny app is funded by NIH U54 grant HL127365. Please see the ', url, 'regarding use and citation of the published and unpublished data presented here.</br>',
            '<b>© 2020 Sorger Lab, Harvard Medical School</br>')))
  })
  
  output$lsp_logo <- renderImage({
    filename <- normalizePath('./www/LSP_Horiz-Logo.png')
    list(src=filename,
         width = 400,
         height = 100)
  })
  
  output$nominal_targets <- renderText({
    drug_label = str_replace(input$agent, '_', '/')
    targets <- cm[cm$agent == drug_label, 'nominal_target']
    paste(targets)
  })
  
  output$hmsl_id <- renderUI({
    drug_label = str_replace(input$agent, '_', '/')
    lincs_id <- cm[cm$agent == drug_label, 'hmsl_id']
    url <- a(lincs_id, href=sprintf('http://lincs.hms.harvard.edu/db/sm/%s/', lincs_id),
             target='_blank')
    tagList(HTML(paste(url)))
  })
  
  output$gr_metrics <- renderPlot({
    ds <- gr_metric[gr_metric$agent == input$agent, ]
    ds$cell_line <- factor(ds$cell_line, levels = ds$cell_line[order(-ds$GR_AOC)])
    p <- ggplot(data=ds, aes_string(x='cell_line', y='GR_AOC', fill='receptor_status')) + 
      geom_bar( stat='identity')+#,
       #          size=5, alpha=0.7) +
      # geom_label_repel(data=ds, aes_string(x='rank', y='GR_AOC', label='cell_line')) +
      scale_fill_manual(name='Receptor status', values = rs_colors) +#, name='Receptor status') +
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
  
  # https://gitlab.com/snippets/16220
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
    style <- paste0("cursor: pointer; position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "right:", 0.5 * right_px, "px; top:", top_px, "px;")
    
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> cell line: </b>", point$cell_line, "<br/>",
                    "<b> GR AOC: </b>",sprintf("%.02f",  point$GR_AOC), "<br/>",
                    "<b> Receptor_status: </b>", point$receptor_status, "<br/>")))
    )
  })
  
  output$mouse_pointer <- renderUI({
    hover <- input$plot_hover
    ds <- gr_metric[gr_metric$agent == input$agent, ]
    tagList(
      HTML(paste('#gr_metrics{cursor:pointer;}')))
    #if (round(hover$x) == as.numeric(ds$cell_line)){
    #  keeprows <- round(input$plot_click$x) == as.numeric(ds$cell_line)
    #  point <- ds[keeprows, ]
    #  gr_aoc <- point$GR_AOC
    #}
    #if ((round(hover$x) == as.numeric(ds$cell_line)) &
    #    ((between(hover$y, 0, gr_aoc)) | (between(hover$y, gr_aoc, 0)))){
    #  tags$style( '#{ cursor: pointer; }')
    #}
  })
  
  output$GRdose_response <- renderPlot({
    ds <- gr_metric[gr_metric$agent == input$agent, ]
    ds$cell_line <- factor(ds$cell_line, levels = ds$cell_line[order(-ds$GR_AOC)])
    drug_label = str_replace(input$agent, '_', '/')
    
    if (is.null(input$plot_click$x)) return()
    else{ 
      keeprows <- round(input$plot_click$x) == as.numeric(ds$cell_line)
      point <- ds[keeprows, ]
      cell_line <- as.character(point$cell_line)
      gr50 <- point$GR50
      gr_max <- point$GRmax
      gr_aoc <- point$GR_AOC
      if ((between(input$plot_click$y, 0, gr_aoc)) | (between(input$plot_click$y, gr_aoc, 0))){
        grc <- gr_values[gr_values$agent == input$agent & 
                           gr_values$cell_line == cell_line, ]
        grc$log10_conc <- log10(grc$concentration)
        tgc <- summarySE(grc, measurevar="GRvalue", groupvars=c("log10_conc"))
        sc_conc <- unique(sapply(grc$log10_conc, function(x){formatC(10**x, format='e', digit=0)}))
        
        if (input$fit_option == 'sigmoidal'){
          dr <- sigmoidal_fit(input$agent, cell_line)
        } else if (input$fit_option == 'biphasic'){
          dr <- biphasic_fit(input$agent, cell_line)
          #dr_low <- sigmoidal_fit(input$agent, cell_line, low_high='low')
          #dr_high <- sigmoidal_fit(input$agent, cell_line, low_high='high')#, measure='GR_toxic')
        }
        
        p <- ggplot() + 
          geom_line(data=dr, aes_string(x='log10_conc', y='yfit'), color='black', size=2, alpha=0.5) +
          geom_errorbar(data=tgc, aes(ymin=GRvalue-se, 
                            ymax=GRvalue+se, x=log10_conc), width=.1) +
          geom_point(data=tgc, aes_string(x='log10_conc', y='GRvalue'), fill='white', size=3, shape=21) +
          
          xlab(paste(drug_label, ' (µM)')) + ylab('GR value') +
          ggtitle(cell_line) +
          ylim(-1, 1.1) +
          #theme_bw() + 
          theme(
            panel.background=element_blank(),
            axis.title=element_text(size=16, face="bold"),
            #legend.title=element_blank(),
            axis.line=element_line(),
            axis.ticks=element_line(),
            axis.text.x = element_text(size=10),
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
        cl <- unique(grc$log10_conc)
        breaks = sort(cl[cl == round(cl)])
        breaks = breaks[seq(1, length(breaks), 2)]
        xlabels = sapply(breaks, function(x){formatC(10**x, format='e', digit=0)})
        #p <- p + xlim(min(breaks)-0.25, max(breaks)+0.25) 
        p <- p + scale_x_continuous(breaks=breaks, 
                                    labels=xlabels)
        if (input$st_tox & input$fit_option=='sigmoidal'){
          dr_static <- sigmoidal_fit(input$agent, cell_line, measure='GR_static')
          dr_toxic <- sigmoidal_fit(input$agent, cell_line, measure='GR_toxic')
          p <- p + geom_line(data=dr_static, aes_string(x='log10_conc', y='yfit'), color='green', size=2, alpha=0.5) +
            geom_line(data=dr_toxic, aes_string(x='log10_conc', y='yfit'), color='red', size=2, alpha=0.5)
        }
        
        #if (input$fit_option=='biphasic'){
        #  p <- p + geom_line(data=dr_low, aes_string(x='log10_conc', y='yfit'), color='black', linetype='dotted', size=0.5, alpha=0.5) +
        #    geom_line(data=dr_high, aes_string(x='log10_conc', y='yfit'), color='black', linetype='dashed', size=0.5, alpha=0.5)
        #}
        return(p)
      } else return()
    }
    #point <- nearPoints(ds, input$plot_click, xvar='cell_line', yvar='GR_AOC')
    #if (dim(point)[1] >=1){
  },width=300, height=300)
  
  output$incfrac_dead <- renderPlot({
    ds <- gr_metric[gr_metric$agent == input$agent, ]
    ds$cell_line <- factor(ds$cell_line, levels = ds$cell_line[order(-ds$GR_AOC)])
    drug_label = str_replace(input$agent, '_', '/')
    
    if (is.null(input$plot_click$x)) return()
    else{ 
      keeprows <- round(input$plot_click$x) == as.numeric(ds$cell_line)
      point <- ds[keeprows, ]
      cell_line <- as.character(point$cell_line)
      gr50 <- point$GR50
      gr_max <- point$GRmax
      gr_aoc <- point$GR_AOC
      if ((between(input$plot_click$y, 0, gr_aoc)) | (between(input$plot_click$y, gr_aoc, 0))){
        grc <- gr_values[gr_values$agent == input$agent & 
                           gr_values$cell_line == cell_line, ]
        grc$log10_conc <- log10(grc$concentration)
        tgc <- summarySE(grc, measurevar="increase_fraction_dead", groupvars=c("log10_conc"))
        
        p <- ggplot(tgc, aes_string(x='log10_conc', y='increase_fraction_dead')) + 
          geom_errorbar(aes(ymin=increase_fraction_dead-se, 
                            ymax=increase_fraction_dead+se), width=.1) +
          geom_line() +
          geom_point(size=3, shape=21, fill="white") +
          xlab(paste(drug_label, ' (µM)')) + ylab('increase fraction dead') +
          ggtitle(cell_line) +
          ylim(-0.05, 1) +
          #theme_bw() + 
          theme(
            panel.background=element_blank(),
            axis.title=element_text(size=16, face="bold"),
            #legend.title=element_blank(),
            axis.line=element_line(),
            axis.ticks=element_line(),
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=12, face='bold'),
            axis.text=element_text(size=14, face='bold')
          ) + 
          geom_hline(yintercept=1, alpha=0.5) +
          geom_hline(yintercept=0, alpha=0.5)
        cl <- unique(grc$log10_conc)
        breaks = sort(cl[cl == round(cl)])
        xlabels = sapply(breaks, function(x){formatC(10**x, format='e', digit=0)})
        p <- p + scale_x_continuous(breaks=breaks, 
                                    labels=xlabels)
        return(p)
      } else return()
    }
  },width=300, height=300)
  
  
  output$cellcycle_response <- renderPlot({
    ds <- gr_metric[gr_metric$agent == input$agent, ]
    ds$cell_line <- factor(ds$cell_line, levels = ds$cell_line[order(-ds$GR_AOC)])
    drug_label = str_replace(input$agent, '_', '/')
    
    if (is.null(input$plot_click$x)) return()
    else{ 
      keeprows <- round(input$plot_click$x) == as.numeric(ds$cell_line)
      point <- ds[keeprows, ]
      cell_line <- as.character(point$cell_line)
      gr50 <- point$GR50
      gr_max <- point$GRmax
      gr_aoc <- point$GR_AOC
      if ((between(input$plot_click$y, 0, gr_aoc)) | (between(input$plot_click$y, gr_aoc, 0))){
        cl <- make_cc_long_table(input$agent, cell_line)
        cl$phase <- factor(cl$phase, levels = c('M', 'beyondG2', 'G2', 'S_dropout', 'S', 'G1', 'subG1'))
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
          xlab(paste(drug_label, ' (µM)')) + ylab(ylabel) +
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
        if (!is.na(cl$phase_count)){
          lc <- unique(cl$log10_conc)
          breaks = sort(lc[lc == round(lc)])
          xlabels = sapply(breaks, function(x){formatC(10**x, format='e', digit=0)})
          p <- p + scale_x_continuous(breaks=breaks, 
                                      labels=xlabels)
        }
        return(p)
      } else return()
    }
  },width=350, height=350)
}


shinyApp(
  ui = dashboardPage(title='HMS LINCS BRCA Browser',
    header,
    #dashboardHeader(title = "HMS LINCS Breast Cancer Browser", titleWidth = 600),
    dashboardSidebar(disable = TRUE),
    body,
    skin='black'
  ),
  server = server
)