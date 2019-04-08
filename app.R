library(shiny)
library(plotly)
library(dplyr)
library(shinyjs)
library(shinydashboard)
library(lavaan)
library(rms)
library(rsvg)
library(shinyBS)
library(data.table)
library(DT)
library(heatmaply)

load("Shiny Workspace.RData")
options(shiny.error = T, shiny.trace=F)

cols.GISTIC <- c("blue", "lightblue", "gray", "pink", "red")
names(cols.GISTIC) <- c("Deep loss", "Shallow loss", "No change", "Gain", "Amp")
cols.factor <- c("Gray", "Orange", "Green", "Blue", "Red")

getTable <- function(graph){
  nodes <- sapply(graph$nodes, function(x) x[1])
  tab <- sapply(graph$edges, function(x){
    strsplit(x[1], split = ":", fixed = T)[[1]]
  })
  types <- sapply(graph$edgeTypes, function(x) sapply(x, function(y) y[1]))
  rbind(tab, types)
}

coefVar <- function(x) {
  out <- tryCatch( (x - mean(x, na.rm = T))/sd(x, na.rm = T),
                   error=function(cond) {
                     return(NA)
                   })    
  return(out)
  
}

make.formula <- function(mat, ftab){
  
  factors <- mat[mat[,3] == "Factor",,drop=F]
  regress <- mat[mat[,3] == "Causal",,drop=F]
  
  ##Parsing factor loadings
  mat <- factors
  innodes <- unique(mat[,2])
  
  validNodes = union(innodes, colnames(ftab))
  
  model <- ""
  for(i in innodes){
    onodes <- intersect(mat[mat[,2] == i,1], validNodes)
    if(length(onodes) == 0) next
    out <- paste(onodes, collapse = " + ")
    expr <- paste(i, "=~", out, collapse = " ")
    model <- paste(model, "\n", expr, collapse = "")
  }
  
  ##Parsing regressions
  mat <- regress
  innodes <- intersect(mat[,2], validNodes)
  
  for(i in innodes){
    onodes <- intersect(mat[mat[,2] == i,1], validNodes)
    if(length(onodes) == 0) next
    out <- paste(onodes, collapse = " + ")
    expr <- paste(i, "~", out, collapse = " ")
    model <- paste(model, "\n", expr, collapse = "")
  }
  
  return(model)
}

makeFTableFun <- function(nodes, dataset){
  types <- sapply(colnames(dataset), function(x) strsplit(x, split = ".", fixed = T)[[1]][2])
  
  res <- dataset[, !(types %in% c("Factor", "Surv"))]
  
  res <- droplevels(res)
  facs <- intersect(colnames(res), names(factors2levels))
  
  if(length(facs) > 0){
    for(i in facs) res[[i]] <- factor(res[[i]], levels = names(sort(table(res[[i]]), decreasing = T)))
    opt <- options("na.action")
    options(na.action = "na.pass")
    mat <- model.matrix(~ ., data = res[facs])
    options(na.action = opt$na.action)
    m <- match(facs, colnames(res))
    res <- res[,-m, drop=F]
    res <- cbind(res, mat)
  }
  
  for(i in 1:ncol(res)){
    res[[i]] <- coefVar(as.numeric(res[[i]]))
  }
  
  dr <- apply(res, 2, function(y) all(is.na(y)))
  res <- subset(res, select = colnames(res)[!dr])
  
  colnames(res) <- sapply(colnames(res), function(x) gsub(" ", "_", x, fixed = T))
  
  return(res)
}

calculateSurvival <- function(mat, ftab, clin){
  nodes <- mat[,2]
  types <- sapply(nodes, function(x) strsplit(x, split = ".", fixed = T)[[1]][2])
  m <- mat[types == "Surv", ,drop = F]
  if(length(types[types == "Surv"]) == 0) return(NULL)
  if(is.vector(m)) m <- matrix(m, nrow = 1, dimnames = list(1, names(m)))
  vars <- unique(m[,1])
  var <- paste(vars, collapse = "+");
  f <- paste("Surv(T,event) ~", var, collapse = "+")
  f <- as.formula(f)
  a <- intersect(rownames(ftab), rownames(clin))
  
  tryCatch({
    tab <- data.frame(ftab[a,vars], T = clin[a,"OS_time"], event = clin[a,"OS"])
    colnames(tab)[1:length(vars)] <- vars
    cox <- coxph(f, data = tab)
    s <- as.data.frame(summary(cox)$coef)
    s <- s[, c(1,3:5)]
    colnames(s) <- c("std.all", "se", "z", "pvalue")
    s$lhs <- rep("Overall_Survival.Surv", nrow(s))
    s$rhs <- vars
    s$op <- " ~"
    s <- s[,c("lhs", "op", "rhs", "z", "pvalue", "std.all")]
    return(list(aic = round(extractAIC(cox)[2], digits = 0), params = s))
  }, error = function(e) return(NULL)
  )
}

collapseParams <- function(mat, params, factorlist){
  m <- apply(mat, 1, function(x){
    if(x[1] %in% names(factorlist)) a <- factorlist[[x[1]]] else a <- x[1]
    if(x[2] %in% names(factorlist)) b <- factorlist[[x[2]]] else b <- x[2]
    pm <- params[params$op == "~",]
    pm <- params[params[,1] %in% b & params[,3] %in% a,]
    z <- which.max(abs(pm[,"z"]))
    if(length(z) == 0) z <- NA else z <- pm[z, "z"]
    std.all <- which.max(abs(pm[,"std.all"]))
    if(length(std.all) == 0) std.all <- NA else std.all <- pm[std.all, "std.all"]
    pvalue <- min(pm[,"pvalue"], na.rm = T)
    return(c(lhs = x[2], op = "~", rhs = x[1], z = round(z, digits = 2), pvalue = pvalue, std.all = round(std.all, digits = 2)))
  })
  colnames(m) <- 1:ncol(m)
  m <- as.data.frame(t(m))
  m[m == Inf | m == -Inf] <- NA
  return(m)
}

factorizeMat <- function(mat, col = 1, factorlist, ftab){
  i <- mat[,col] %in% names(factorlist)
  if(length(i[i]) == 0) return(mat)
  m <- mat[i,,drop=F]
  
  mat <- mat[!(mat[,col] %in% names(factorlist)), , drop = F]
  m2 <- apply(m, 1, function(x){
    f <- intersect(colnames(ftab), factorlist[[x[col]]])
    s <- sapply(f, function(y) 
      ifelse(col == 1, return(c(y, x[2:3])), return(c(x[1], y, x[3]))))
    return(list(s))
  })
  
  mm <- t(m2[[1]][[1]])
  if(length(m2) > 1){
    for(i in m2[2:length(m2)]) mm <- rbind(mm, t(i[[1]]))
  }
  
  if(nrow(mat) > 0) return(rbind(mm,mat))
  else return(mm)
}


downloadGraph <- function(file, what){
  writeLines(what, con = file)
}

showCorMatrixDialog <- modalDialog(
  fluidRow(div(id = "corHeatmapDialogDiv", title = "Select dataset for the correlation heatmap",
      selectInput("selectDatasetHeatmap", label = "Select dataset", choices = as.character(disease.map), selected = "BLCA"))),
  fluidRow(div(id = "corHeatmapDiv", plotlyOutput("corHeatmap")))
  , title = "Correlation matrix for the model variables", easyClose = T)

showCorsDialog <- modalDialog({
  wellPanel(
    style = "background-color: #FFFFFF;",
    fluidRow(
      p("View how your selected variables in the model correlate with all other variables in the selected datasets below"),
      div(id = "corsDialogDiv1", title = "Select data type to view correlations for",
          selectInput("corsDialogDataTypeSelect", label = "Select data type", choices = names(dataNames)[-8], selected = "RPPA")),
      div(id = "corsDialogDiv2", title = "Select dataset to view the correlations for",
          selectInput("corsDialogDatasetSelect", label = "Select dataset", choices = as.character(disease.map), selected = "BLCA"))
    ),
    hr(),
    h4("Table of genome-wide correlations*"),
    dataTableOutput("corsTable"),
    p("*: Table reports Cox t-statistic for correlations with survival, and Pearson's r values for others"),
    hr(),
    actionButton("importVars", label = strong("Import selected variables to model"))
  )}, title = "Global pair-wise correlations of selected variables", easyClose = T)
  

tableStatsDialog <- bsModal(
  id = "statsModalDialog",
  title = "Fit statistics for evaluated models with SEM",
  trigger = "showStatsBtn",
  selectInput(
    "modelSelect",
    label = "Dataset to display the fit statistics for",
        choices = as.character(disease.map),
        selected = as.character(disease.map)[1],
        width = "180px"
    ),
  wellPanel(
    id = "modelTablePanel",
    shiny::dataTableOutput("modelTable"),
    style = "padding: 0px; background-color: #ffffff; height: 370px;"
  ),
  size = "large"
)

graphBox <- box(
  id = "graphBox",
  width = 12,
  height = 676,
  
  wellPanel(
    id = "graphWellPanel",
    style = "background-color: #fffaf4;",
    fluidRow(
      
      div(id = "panelElements",
          div(id= "left",
              div(id = "edgeSelectDiv",
                  selectInput(
                    "edgeSelect",
                    width = "80px",
                    label = "Link type",
                    choices = c("Causal", "Factor"),
                    selected = "Causal"
                  )
              )
              ,
              div(id = "layoutDiv",
                  #disabled(
                  selectInput("layoutSelect",
                              width = "110px",
                              label = "Layout type",
                              choices = c("Custom", 
                                          "Hierarchical",
                                          "Organic",
                                          "Orthogonal",
                                          "Tree", "Radial"
                              ),
                              selected = "Hierarchical")
                  
                  #)
                  
              )
          ),
          div(id = "center",
              div(id = "edgeLabelDiv", title = "Select parameters to mark links",
                  disabled(
                    selectInput(
                      "edgeLabel",
                      width = "90px",
                      label = "Mark links by",
                      choices = c("Nothing", "z", "p-value", "std.all"),
                      selected = "Nothing"
                    )
                    
                  )
              ),
              div(id = "datasetChooseDiv", title = "Select dataset results to display",
                  disabled(selectizeInput("datasetSelect", 
                                          width = "90px",
                                          label = "Dataset",
                                          choices = c("Average", as.character(disease.map)),
                                          selected = "Average"
                  ))
              ),
              div(id = "modelChooseDiv", title = "Select model graph to display",
                  disabled(selectInput("modelSelectMenu",
                                          width = "80px",
                                          label = "Model", choices = NULL)))
          ),
          div(id = "right",
              div(id = "submitButtonDiv",
                  actionButton("submit", "Evaluate")
              ),
              div(id = "exportBtnDiv", title = "Export graph",
                  actionButton("exportBtn", label = "", icon = icon("upload"), class = "greyButton")
              ),
              div(id = "printGraphDiv", title = "Print graph",
                  actionButton("btnPrintGraph", label = "", icon = icon("print"), class = "greyButton")
              ),
              div(id = "binarizeDiv", title = "Categorize the variable (e.g. high and low)",
                  disabled(
                    actionButton("binarizeBtn", label = "Bin", class = "greyButton")
                  )
              ),
              div(id = "combineButtonDiv", title = "Combine variables into one",
                  disabled(
                    actionButton("combineBtn", label = "", icon = icon("object-group"), class = "greyButton")
                  )
              ),
              div(id = "showStatsButtonDiv", title = "Display the fit statistics for models",
                  disabled(
                    actionButton("showStatsBtn", label = "Stats", class = "greyButton")
                  ),
                  tableStatsDialog
              ),
              div(id = "showCorsDialogBtnDiv", title = "View global correlations of selected variables",
                  actionButton("showCorsDialogBtn", "", icon = icon("globe"), class = "greyButton")
                  ),
              div(id = "showCorMatrixDialogBtnDiv", title = "View the covariance matrix of data in the model",
                  disabled(actionButton("showCorMatrixDialogBtn", "", icon = icon("th-large"), class = "greyButton"))
              )
          )
          
          
      )
      ,
      hidden(
        div(class = "loader", id = "loading2")  
      )
      
    )
    
    
  ),
  fluidRow(
    div(id = "graphComponent", style = "height: 513px;"),
    style = "padding-left:0px; padding-right:0px; padding-top: 0px; border: 2px solid lightgray"
    
  )
  
  ,
  # tags$hr(),
  
  
  style = "padding-left: 20px; padding-right: 20px;"
  
)

addNodesBox <- box(
  width = 30,
  strong("Add nodes to model:", style = "font-size:120%"),
  fluidRow(div(id = "shapeDiv",
      selectInput(
        inputId = "shape",
        label = "",
        choices = c(
          "Choose node type" = "",
          "mRNA (RNAseq)",
          "miRNA (miRNA-seq)",
          "Protein (RPPA)",
          "Somatic mutation (WXS)",
          "Germline variation (WXS)",
          "Copy Number Variation (SNP6)",
          "Copy Number Variation (GISTIC)",
          "Tumor features (PanCan iAtlas)",
          "Clinical",
          "Factor"
        )
      )
  )),
  fluidRow(
    div(id = "nameDiv",
      #Use selectizeinput unless shape is not Selected Factor
      conditionalPanel(
        condition = "input.shape != 'Factor'",
        
        selectizeInput(
          "dynamic",
          label = NULL,
          choices = NULL,
          # width = 305, 
          options = list(
            placeholder = "Start typing...",
            maxItems = NULL,
            maxOptions = 20,
            closeAfterSelect = T,
            addPrecedence = T
          )
        )
      ),
      conditionalPanel(
        condition = "input.shape == 'Factor'",
        textInput(
          "factor",
          #width = 305,
          placeholder = "Enter factor name...",
          value = "",
          label = NULL
        )
      )
  ),
  actionButton("nodeButton", "Add", class = "greyButton")
  )
)

plotOptionsBox <- box(id = "plotOptionsBox",
  width = "60px",
  selectInput("datasets", label = "Choose datasets to plot", 
              choices = datasetNames, selected = datasetNames[1], multiple = TRUE, size = length(datasetNames), selectize = F)
  #  downloadButton("downloadPlot", "Download as PDF")
)

exportDialog <- modalDialog({
  fluidRow(column(
    width = 6,
    radioButtons(
      "exportRadio",
      label = "Options",
      choices = c("png", "pdf", "svg", "raw data"),
      inline = T
    ),
    downloadButton("exportDialogBtn", "export")
  ))
},
title = "Export graph", size = "m", easyClose = T)

binarizeDialog <- modalDialog({
  fluidRow(
    p('Choose a cutoff to categorize the selected variable into two. Label the categories (e.g. "Low" and "High",
      or "Neutral" and "Amplified" for the two categories)'),
    p("Analyze the distribution of the data in the plot window to choose a proper data cutoff"),
    wellPanel(
      fluidRow(
        column(width = 5, textInput("leftCategory", "Label 1st Category", value = "Low")),
        column(width = 2, numericInput("cutoff", "Cutoff", 0)),
        column(width = 5, textInput("rightCategory", "Label 2nd Category", "High"))
      )
    ),
    style = "padding-left:20px; padding-right:20px; padding-top: 10px;"
  )
}, title = "Categorize variable", footer = wellPanel(actionButton("okModalBtn", "OK"),modalButton("Cancel")))


logx <- div(id = "logxCheckDiv",
            checkboxInput(
              "xaxislog",
              label = "log-scale on X"
            )
)
logy <- div(id = "logyCheckDiv",
            checkboxInput(
              "yaxislog",
              label = "log-scale on Y"
            )
)

logy2 <- div(id = "logyCheckDiv2",
            checkboxInput(
              "yaxislog2",
              label = "log-scale on Y"
            )
)

exportPDFDiv <- div(id = "pdfdiv", title = "Export Plot",
                    downloadButton("exportPDF", "PDF"),
                    downloadButton("exportPNG", "PNG")
                 )

heatmapPanel <- div(
  id = "heatmapPanelDiv",
  div(id = "heatmapDatasetSelectDiv",
      selectInput("heatmapDatasetSelect", label = "Select dataset", choices = as.character(disease.map), selected = "BLCA"))
)

scatterPanel <- div(id = "scatterPanelDiv",
                logx,
                logy2,
                div(id = "Charttype",
                    selectInput(
                      "chartType",
                      width = "90px",
                      label = "Chart type",
                      choices = c("Scatter", "Histogram 2D"),
                      selected = "Scatter"
                    )
                )
              )
boxPanel <- div(id = "boxPanelDiv",
                logy,
                div(id = "showDataPoints",
                    checkboxInput(
                      "showPointsCheck",
                      label = "Show all data points"
                    )
                )
              )
factorPanel <- div(id = "proportions",
                selectInput(
                  "factorType",
                  width = "130px",
                  label = "Y-axis options",
                  choices = c("Total numbers", "Proportion"),
                  selected = "Total numbers"
                )
            )
survPanel <- disabled(div(id = "SurvOptions",
                                       div(id = "groupNumberDiv",
                                           selectInput("numberGroupsSurv", "Number of groups", choices = c(2, 3),
                                                       width = "130px", selected = 2)
                                           ),
                                       div(id = "cutoffSliderDiv",
                                           conditionalPanel(
                                             condition = "input.numberGroupsSurv == 2",
                                             sliderInput(
                                               "cutoffSurv",
                                               label = "Percentile cutoff for groups",
                                               min = 1, max = 99, value = 50)
                                           ),
                                           conditionalPanel(
                                             condition = "input.numberGroupsSurv == 3",
                                             sliderInput(
                                               "cutoffSurv2",
                                               label = "Percentile cutoff for groups",
                                               min = 1, max = 99, value = c(25, 75))
                                           )
                                           
                                       ))
                         )

ui <- fluidPage(
                 tags$head(
                   HTML('<link rel="icon", href="sema icon 5.png", 
                                   type="image/png" />'),
                   HTML('<title>sema-cancer</title>'),
                   tags$style(includeHTML("Styles.css")),
                   tags$style(type="text/css",
                              ".shiny-output-error { visibility: hidden; }",
                              ".shiny-output-error:before { visibility: hidden; }"
                   ),
                   div(class = "loader", id = "loading"),
                   useShinyjs(),
                   tags$link(rel = "stylesheet", href = "yFiles2/lib/yfiles.css"),
                   tags$script(src = "yWorks.yFilesHTML.License.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/lang.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/core-lib.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/graph-core.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/graph-binding.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/graph-input.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/graph-styles-core.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/graph-styles-default.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/graph-styles-other.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/graph-styles-template.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/layout-core.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/algorithms.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/layout-hierarchic.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/layout-radial.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/layout-orthogonal.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/layout-organic.js"),
                   tags$script(src = "yFiles2/lib/yfiles/impl/graph-layout-bridge.js"),
                   tags$script(src = "yFiles2/ide-support/yfiles-typeinfo.js"),
                   tags$script(
                     'Shiny.addCustomMessageHandler("refocus",
                        function(NULL) {
                         document.getElementById("factor").focus();
                        });'
                   )
                 ),
                 fluidPage(
                   fluidRow(column(width = 1, div(img(src = "sema icon 3.png", width = 150))), column(width = 1, div(img(src = "CCHMC-2016-Logo.png", width = 200)))),
                   hr(),
                            div(
                              id = "editPanel",
                              fluidRow(
                                column(
                                  width = 8,
                                  graphBox
                                ),
                                column(
                                  width = 4,
                                  fluidRow(
                                    column(width = 6, addNodesBox)
                                  ),
                                  hr(),
                                  fluidRow(
                                    disabled(actionButton("plotButton", "Plot selection", icon = icon("bar-chart-o"))),
                                    disabled(actionButton("heatmapButton", "Heatmap", icon = icon("th")))),
                                  fluidRow(
                                    plotOptionsBox
                                  )
                                )
                                
                              ),
                              fluidRow(
                                column(
                                  width = 10,
                                  column(
                                    width = 9,
                                    box(
                                        id = "largePlotBox",
                                        width = 12,
                                        wellPanel(
                                          id = "plotOptionsPanel",
                                          style = "background-color: white;",
                                          fluidRow(
                                          div(id = "edgePlotDiv",
                                            hidden(
                                              selectInput("edgePlot", label = "Plot for the link:", width = '170px',
                                                      choices = c("Raw data", 
                                                                  "Beta values", "Z-scores", "P-values (-log10)"),
                                                      selected = "Raw data"))
                                          ),
                                          div(id = "plotOptionsDiv",
                                          conditionalPanel(
                                            condition = "output.pType == 'Box'",
                                            boxPanel
                                          ),
                                          conditionalPanel(
                                            condition = "output.pType == 'Scatter'",
                                            scatterPanel
                                          ),
                                          conditionalPanel(
                                            condition = "output.pType == 'Factor'",
                                            factorPanel
                                          ),
                                          conditionalPanel(
                                            condition = "output.pType == 'Surv'",
                                            survPanel
                                          ),
                                          conditionalPanel(
                                            condition = "output.pType == 'Heatmap'",
                                            heatmapPanel
                                          ),
                                          hidden(exportPDFDiv))
                                        )),
                                        plotlyOutput("largePlot")
                                    )
                                  )
                                )
                              ),
                              fluidRow(
                                style = "padding-left:0px; padding-right:10px; padding-top: 110px")
                            )),
                 tags$script(src = "Script.js"),
                 extendShinyjs(script = "shapeFun.js"),
                 tags$script(
                   "window.onload = function(){ 
                        document.getElementById('loading').style.display = 'none';
                     }"
                 )
)

combineNumericDialog <- modalDialog(
  title = strong("Combine variables"),
  textInput("combineNumericText", "Name the group variable:"),
  hr(),
  p("The new variable will be defined based on the average or the sum of the selected variables"),
  hr(),
  wellPanel(selectInput("combineNumericSelect", "Define function to combine the variables with",
              choices = c("Average", "Sum"), selected = "Average")),
  footer = wellPanel(actionButton("combineNumericBtn", "OK"),modalButton("Cancel"))
)

combineFactorDialog <- modalDialog(
  title = strong("Combine variables"),
  p("The new variable will be a categorical (binary) variable. It will be assigned to the second category if \n
    all or any of the conditions below are satisfied, and to the first category otherwise."),
  hr(),
  textInput("combineFactorText", "Name the group variable:"),
  wellPanel(
    id = "factorWellPanel",
    style = "background-color: #fffaf4;",
    fluidRow(
      div(id = "factor1", textInput("firstCat", "Name first category")),
      div(id = "factor2", textInput("secondCat", "Name second category"))
    )
  ),
  hr(),
  p("Select the rules for the new group variable:"),
  checkboxInput("allorany", "Satisfy all (checked) or any (unchecked) of the conditions below"),
  uiOutput("combineCategoricalTable"),
  footer = wellPanel(actionButton("combineFactorBtn", "OK"),modalButton("Cancel"))
)

server <- function(input, output, session) {
  #Keeps the name of the chosen dataset
  submitted <- reactiveVal(FALSE)
  chosenDataset <- reactiveVal(value = NULL)
  layoutType <- reactiveVal(value = NULL)
  edgeLabelType <- reactiveVal(value = "Nothing");
  arrowType <- reactiveVal(value = "Causal");
  selectedNode <- reactiveVal(value = NULL)
  plotType <- reactiveVal(value = NULL);
  modelCount <- reactiveVal(value = 0)
  svg_img <- reactiveVal()
  currentPlot <- reactiveVal()
  
  survSlider1 <- reactive(input$cutoffSurv) %>% debounce(1000)
  survSlider2 <- reactive(input$cutoffSurv2) %>% debounce(1000)
  
  output$pType <- reactive({
    plotType()
  })
  
  outputOptions(output, "pType", suspendWhenHidden = FALSE)
  
  #initialize tempData
  tempData <- reactiveValues();
  tempData$DF <- Model.DF
  
  #initialize modelsList
  modelsList <- reactiveValues();
  
  getNumericData <- function(dtype){
    require(plyr)
    
    m <- tempData$DF[tempData$DF$Type == dtype,4:ncol(tempData$DF)]
    types <- sapply(colnames(m), function(x) strsplit(x, split = ".", fixed = T)[[1]][2])
    m <- m[, !(types %in% c("Clin", "Surv"))]
    m <- apply(m, 2, function(x) 
      as.numeric(mapvalues(x, c("ALT", "REF", "WT", "Mut", "Deep loss", "Shallow loss", "No change", "Gain", "Amp"), 
                           c(2, 0, 0, 2, -3, -1.5, 0, 1.5, 3))))
    rownames(m) <- rownames(tempData$DF[tempData$DF$Type == dtype,])
    return(m)
  }
  
  #After shape is selected the dataset is assigned to be used in selectizeInput
  
  output$combineCategoricalTable <- renderUI({
    name <- paste0(input$combineNumericText, ".Group")
    selection <- strsplit(input$multiSelection, split = ",")[[1]]
    lapply(selection, function(x)
      selectInput(x, label = paste(x, "should be:"), choices = levels(tempData$DF[[x]])))
  })
  
  output$exportPDF <- downloadHandler(filename = function() {
    paste("myPlot", Sys.time(), "pdf", sep = ".")}, content = function(file) {
      pp <- isolate(currentPlot())
      export(pp, file = file, cliprect = "viewport")
  })
  
  output$exportSVG <- downloadHandler(filename = function() {
    paste("myPlot", Sys.time(), "svg", sep = ".")}, content = function(file) {
      pp <- isolate(currentPlot())
      export(pp, file = file)
    })
  
  output$exportPNG <- downloadHandler(filename = function() {
    paste("myPlot", Sys.time(), "png", sep = ".")}, content = function(file) {
      pp <- isolate(currentPlot())
      export(pp, file = file)
    })
  
  observeEvent(input$combineFactorBtn, {
    name <- gsub(".", "", input$combineFactorText, fixed = T)
    name <- gsub(" ", "", name, fixed = T)
    name <- gsub(",", "", name, fixed = T)
    name <- gsub("/", "", name, fixed = T)
    name <- gsub("-", "", name, fixed = T)
    name <- paste0(name, ".", "Group")
    if(name %in% colnames(tempData$DF)){
      showNotification("Name already exists in the model")
      return(NULL)
    }
    if(input$combineFactorText == ""){
      showNotification("Name cannot be empty")
      return(NULL)
    }
    if(input$firstCat == "" | input$secondCat == ""){
      showNotification("Category labels cannot be empty")
      return(NULL)
    }
    selection <- strsplit(input$multiSelection, split = ",")[[1]]
    options <- sapply(selection, function(x) input[[x]])
    if(input$allorany) 
      s <- apply(tempData$DF[,selection], 1, function(x) all(x == options))
    else s <- apply(tempData$DF[,selection], 1, function(x) any(x == options))
    dat <- factor(s)
    levels(dat) <- c(input$firstCat, input$secondCat)
    tempData$DF[[name]] <- dat
    js$addNode(type = "Group",
               label = name,
               index = 0)
    removeModal()
    js$layout(layoutType())
  })
  
  observeEvent(input$showCorsDialogBtn, {
    showModal(showCorsDialog)
  })
  
  observeEvent(input$showCorMatrixDialogBtn, {
    showModal(showCorMatrixDialog)
  })
  
  heatmapDataset <- reactive(input$selectDatasetHeatmap)
  
  output$corHeatmap <- renderPlotly({
    require(plyr)
    m <- getNumericData(heatmapDataset())
    
    m <- cor(m, use = "pairwise.complete.obs")
    m[is.na(m)] <- 0
    m[m < -0.75] <- -0.75
    m[m > 0.75] <- 0.75
    heatmaply(m, colors = c("blue", "white", "red"), limits = c(-0.75, 0.75)) %>%
      config(displayModeBar = F)
  })
  
  output$modelTable <- shiny::renderDataTable({
    if(modelCount() == 0) return(NULL)
    tbl <- sapply(1:modelCount(), function(x) {
      nam <- paste("Model", x)
      m <- modelsList[[nam]][["measures"]][,input$modelSelect]
      return(m)
    })
    tbl1 <- data.frame(t(tbl))
    return(tbl1)
  },
  options = list(
    searching = F,
    autowidth = T,
    lengthChange = F,
    ordering = T,
    paging = F,
    processing = T,
    scrollY = "340px",
    scrollCollapse = T,
    info = F,
    nowrap = T

  )
  )
  
  observeEvent(input$importVars,{
    rows <- input$corsTable_rows_selected
    type <- input$corsDialogDataTypeSelect
    names <- dataNames[[type]][rows]
    switch(type,
           "RNA" = shape <- "mRNA (RNAseq)",
           "RPPA" = shape <- "Protein (RPPA)",
           "CNV" = shape <- "Copy Number Variation (GISTIC)",
           "Mut" = shape <- "Somatic mutation (WXS)",
           "Germ" = shape <- "Germline variation (WXS)",
           "miRNA" = shape <- "miRNA (miRNA-seq)",
           "iAtlas" = shape <- "Tumor features (PanCan iAtlas)"
    )
    for(i in names) {
      addNode(shape = shape, name = i)
      js$addNode(type = shape, label = i, index = 0)
    }
    removeModal()
    js$layout(layoutType())
  })
  
  output$corsTable <- DT::renderDataTable({
    selection <- NA
    s <- input$multiSelection
    if(!is.null(s)) selection <- strsplit(s, split = ",")[[1]]
    if(length(selection) == 0) selection <- input$clicked
    
    dataset <- input$corsDialogDatasetSelect
    datatype <- input$corsDialogDataTypeSelect
    
    dd <- data.frame(Names = dataNames[[datatype]])
    tt <- paste0("Ave.in.", dataset)
    dd[[tt]] <- round(Means[[datatype]][dataNames[[datatype]], dataset], digits = 3)
    
    if(!is.na(selection)){
      for(i in selection){
        s <- strsplit(i, split = ".", fixed = T)[[1]]
        name <- s[1]
        type <- s[2]
        
        if(type == "GISTIC") type <- "CNV"
        
        p <- chmatch(name, dataNames[[type]])
        
        n <- paste0("f", (p-1)%/%1000)
        
        file <- paste0(type, "_to_", datatype)
        p <- (p-1) %% 1000 + 1
        
        if(type == "Surv") {
          file <- paste0(datatype, "_to_", type)
          f <- paste("Data/Correlations", dataset, file, sep = "/")
          if(file.exists(f)) dd[[i]] <- round(as.matrix(fread(f, nrows = 1, sep = ",", data.table = F, header = F))[1,],
                                              digits = 3)
        }
        else {
          f <- paste("Data/Correlations", dataset, file, n, sep = "/")
          if(file.exists(f)) dd[[i]] <- round(as.matrix(fread(f, nrows = 1, sep = ",", skip = p, data.table = F, header = F))[1,],
                                              digits = 3)
        }
      }
    }
    
    dd <- datatable(dd, rownames = F, filter = "bottom", extensions = "Buttons",
                    options = list(
                      dom = "Bltrp", 
                      pagingType = "simple",
                      pageLength = 25,
                      processing = T,
                      scrollX = T,
                      scrollY = "250px",
                      buttons = list("csv"),
                      lengthMenu = list(c(25, 50, 100, -1), list("25", "50", "100", "All"))
                    ))
      
    
    return(dd)
  })
  
  output$exportDialogBtn <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste("myModel", Sys.time(),
            
            switch(
              input$exportRadio,
              "png" = "png",
              "pdf" = "pdf",
              "svg" = "svg",
              "raw data" = "csv"
            ),
            sep = ".")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      switch(
        input$exportRadio,
        "png" = rsvg_png(charToRaw(svg_img()), file),
        "pdf" = rsvg_pdf(charToRaw(svg_img()), file),
        "svg" = writeLines(svg_img(), file),
        "raw data" = {
          types <- sapply(colnames(tempData$DF)[4:ncol(tempData$DF)], function(x) strsplit(x, split = ".", fixed = T)[[1]][2])
          df <- tempData$DF[, c(1:3, (4:ncol(tempData$DF))[!(types %in% c("Germ", "Surv"))])]
          fwrite(df, file = file, row.names = T)
        }
      )
      
    },
    contentType = switch(
      input$exportRadio,
      "png" = "image/png",
      "pdf" = "application/pdf",
      "svg" = "image/svg"
    )
  )
  
  observeEvent(input$combineNumericBtn,{
    if(input$combineNumericSelect == "Average") func <- mean else func <- sum
    varname <- gsub(".", "", input$combineNumericText, fixed = T)
    varname <- gsub(" ", "", varname, fixed = T)
    varname <- gsub(",", "", varname, fixed = T)
    varname <- gsub("/", "", varname, fixed = T)
    varname <- gsub("-", "", varname, fixed = T)
    varname <- paste0(varname, ".", "Group")
    
    if(varname %in% colnames(tempData$DF)){
      showNotification("Name already exists in the model")
      return(NULL)
    }
    if(input$combineNumericText == ""){
      showNotification("Name cannot be empty")
      return(NULL)
    }
    
    selection <- strsplit(input$multiSelection, split = ",")[[1]]
    tempData$DF[[varname]] <- apply(tempData$DF[,selection], 1, func)
    js$addNode(type = "Group",
               label = varname,
               index = 0)
    removeModal()
    js$layout(layoutType())
  })
  
  observeEvent(input$exportBtn, {
    showModal(exportDialog)
  })
  
 
  observeEvent(input$datasetSelect, {
    chosenDataset(input$datasetSelect)
    
    p <- tab()
    
    switch(edgeLabelType(),
           "z" = values <- p$zmat[,chosenDataset()],
           "std.all" = values <- p$stdmat[,chosenDataset()],
           "p-value" = values <- p$pmat[,chosenDataset()],
           "Nothing" = values <- rep("", length(p$lhs))
    )
    
    df <- data.frame(lhs = p$lhs, rhs = p$rhs, zmat = p$zmat[, chosenDataset()], 
                     std = p$stdmat[,chosenDataset()], val = values);
    
    tryCatch({
      js$results(
        df
      )
    }, error = function(e) {})
    
  })
  observeEvent(input$shape, {
    my_data <- NULL
    
    switch(
      input$shape,
      "mRNA (RNAseq)" = my_data <- dataNames[["RNA"]],
      "Protein (RPPA)" = my_data <- dataNames[["RPPA"]],
      "Copy Number Variation (SNP6)" = my_data <- dataNames[["CNV"]],
      "Copy Number Variation (GISTIC)" = my_data <- dataNames[["GISTIC"]],
      "Somatic mutation (WXS)" = my_data <- dataNames[["Mut"]],
      "Germline variation (WXS)" = my_data <- dataNames[["Germ"]],
      "miRNA (miRNA-seq)" = my_data <- dataNames[["miRNA"]],
      "Tumor features (PanCan iAtlas)" = {
        my_data <- dataNames[["iAtlas"]]
        updateSelectizeInput(session = session, "dynamic", choices = my_data, 
                             server = T, options = list(openOnFocus = TRUE,selectOnTab = FALSE,
                                                        maxOptions = 30,
                                                        render = I(
                                                          '{
                                                          function(item, escape) {
                                                          return "<div>" + escape(item.name) + "</div>"
                                                          }
                                  }'
                                                        )))
        return(NULL);
      },
      "Clinical" = {
        my_data <- c("Overall_Survival", clinical_parameters[1:7])
        updateSelectizeInput(session = session, "dynamic", selected = my_data[1], choices = my_data,
                             server = T, options = list(openOnFocus = TRUE,selectOnTab = FALSE,
                                                        maxOptions = 30,
                                                        render = I(
                                                          '{
                                    function(item, escape) {
                                        return "<div>" + escape(item.name) + "</div>"
                                    }
                                  }'
                                                        )))
        return(NULL);
      }
    )
    
    
    
    #Renders the dataset and returns the options while the user typing
    updateSelectizeInput(
      session = session,
      "dynamic",
      choices = my_data,
      server = TRUE,
      options = list(
        openOnFocus = FALSE,
        selectOnTab = FALSE,
        render = I(
          '{
          function(item, escape) {
          return "<div>" + escape(item.name) + "</div>"
          }
  }'
        )
        
        
      )
      
    )
    
  })
  
 
  #clears the factor when user clicks on it
  onclick("factor", {
    updateTextInput(session, "factor", value = "")
  })
  
  
  observeEvent(input$btnPrintGraph, {
    js$printGraph();
  })
  
  observeEvent(input$jsImage, {
    svg_img(input$jsImage)
  })
  
  onclick("exportBtn", {
    js$toImage()
  })
  
  observeEvent(input$layout, {
    js$layout(layoutType())
  })
  
  observeEvent(input$layoutSelect, {
    layoutType(input$layoutSelect)
    if(input$layoutSelect != "Custom")
      js$layout(layoutType())
  })
  
  observeEvent(input$submit, {
    js$submit(chosenDataset())
    enable("showStatsBtn")
  })
  
  
  observeEvent(input$edgeSelect, {
    #js$edgeChange(input$edgeSelect)
    txt <- input$edgeSelect
    if (txt != arrowType()){
      arrowType(txt);
      js$edgeChange(txt)
    }
  })
  
  observeEvent(input$edgeLabel, {
    edgeLabelType(input$edgeLabel)
    res <- tab()
    if(is.null(res)) return(NULL);
    switch(edgeLabelType(),
           "z" = values <- data.frame(lhs=res$lhs, rhs=res$rhs, val = res$zmat[,chosenDataset()]),
           "std.all" = values <- data.frame(lhs=res$lhs, rhs=res$rhs, val=res$stdmat[,chosenDataset()]),
           "p-value" = values <- data.frame(lhs=res$lhs, rhs=res$rhs, val=res$pmat[,chosenDataset()]),
           "Nothing" = values <- data.frame(lhs=res$lhs, rhs=res$rhs, val=rep("", length(res$lhs)))
    )
    
    js$renderLabel(values)
  })
  
  
  observeEvent(input$nodeButton, {
    showNotification("Node Added")
    print(submitted())
    if(submitted()){
      js$results(TRUE)
    }
    
    if (input$shape == "Factor") {
      lab <- gsub(" ", "_", input$factor)
      lab <- gsub("\\.", "_", lab)
      lab <- gsub("-", "_", lab)
      js$addNode(type = input$shape,
                 label = lab,
                 index = 0)
      updateTextInput(session, "factor",
                      value = "")
      session$sendCustomMessage(type = "refocus", message = list(NULL))
    } else{
      for (name in input$dynamic) {
        addNode(input$shape, name)
        js$addNode(type = input$shape,
                   label = name,
                   index = 0)
      }
      # if(length(input$dynamic) > 1){
      #
      # }
      updateSelectInput(session, "dynamic",
                        selected = "")
    }
    enable("layoutSelect")
    enable("layout")
    js$layout(layoutType())
    js$focus()
  })
  
  addNode <- function(shape, name, type){
    if(name == "Overall_Survival") return(NULL);
    if(missing(type)){
      switch(
        shape,
        "mRNA (RNAseq)" = type <- "RNA",
        "Protein (RPPA)" = type <- "RPPA",
        "Copy Number Variation (SNP6)" = type <- "CNV",
        "Copy Number Variation (GISTIC)" = type <- "GISTIC",
        "Somatic mutation (WXS)" = type <- "Mut",
        "Germline variation (WXS)" = type <- "Germ",
        "miRNA (miRNA-seq)" = type <- "miRNA",
        "Tumor features (PanCan iAtlas)" = type <- "iAtlas",
        "Clinical" = type <- "Clin",
        "Factor" = type <- "Factor"
      )
    }
    
    node.name <- paste(name, type, sep = ".")
    
    if(type == "Factor"){
      return(NULL)
    }
    
    if(type == "Clin"){
      dd <- Clinical[,name]
      names(dd) <- rownames(Clinical)
    } else {
      p <- chmatch(name, dataNames[[type]])
      if(type == "miRNA"){
        n <- "f0"
      } else {
        n <- paste0("f", (p-1)%/%1000)
      }
      f <- paste("Data/TCGA", type, n, sep = "/")
      p <- (p-1) %% 1000 + 1
      dd <- fread(f, nrows = 1, sep = ",", skip = p, data.table = F, header = F)
      dd <- as.matrix(dd)[1,]
      names(dd) <- sampleNames[[type]]
      
      if(type %in% names(sampleNames.TARGET)){
        f <- paste("Data/TARGET", type, n, sep = "/")
        p <- (p-1) %% 1000 + 1
        dd2 <- fread(f, nrows = 1, sep = ",", skip = p, data.table = F, header = F)
        dd2 <- as.matrix(dd2)[1,]
        names(dd2) <- sampleNames.TARGET[[type]]
        dd <- c(dd, dd2)
      }
    }
    
    if(type == "Mut") dd <- factor(dd, labels = c("WT", "Mut"))
    if(type == "Germ") dd <- factor(dd, labels = c("REF", "ALT"))
    if(type == "GISTIC") dd <- factor(dd, levels = 
                                        intersect(c("Deep loss", "Shallow loss", "No change", "Gain", "Amp"),
                                                  unique(dd)))
    
    tempData$DF[[node.name]] <- dd[rownames(tempData$DF)]
    if(ncol(tempData$DF) > 4) {
      enable("showCorMatrixDialogBtn")
      enable("heatmapButton")
    }
  }
  
  observeEvent(input$binarizeBtn, {
    showModal(binarizeDialog)
  })
  
  observeEvent(input$combineBtn, {
    selection <- strsplit(input$multiSelection, split = ",")[[1]]
    if("Overall_Survival.Surv" %in% selection){
      showNotification("Survival variable cannot be combined with other variales")
      return(NULL);
    }
    
    selection <- intersect(selection, colnames(tempData$DF))
    if(length(selection) < 2) return(NULL)
    types <- unlist(sapply(tempData$DF[selection], function(x) class(x)))
    if(!all(types == types[1])){
      showNotification("Variables to be combined must be all either numeric or categorical")
      return(NULL);
    }
    
    if(types[1] == "factor"){
      showModal(combineFactorDialog)
    } else {
      showModal(combineNumericDialog)
    }
    
  })
  
  observeEvent(input$modelSelectMenu, {
    
    model = input$modelSelectMenu
    modelN = strsplit(model, split = " ", fixed = T)[[1]][2]
    
    graph <- isolate(modelsList[[model]][["graph"]])
    pars <- isolate(modelsList[[model]][["pars"]])
    
    js$loadGraph(modelN);

    updateGraph(pars)
  }, ignoreInit = TRUE)
  
  observeEvent(input$okModalBtn, {
    left <- input$leftCategory
    right <- input$rightCategory
    thresh <- input$cutoff
    removeModal()
    temp <- unlist(strsplit(selectedNode(), "\\."))
    if(left == "" | right == "") {
      showNotification("Category labels cannot be blank")
      return(NULL);
    }
    tryCatch(
      {
        nd <- selectedNode()
        dat <- tempData$DF[[nd]]
        dat <- cut(dat, breaks = c(-Inf, thresh, Inf), labels = c(left, right))
        tempData$DF[[nd]] <- dat
      }, error = function(e){
        showNotification("Wrong inputs. Check your inputs (categories and cutoff) and try again.")
        return(NULL);
      }
    )
  })

  observeEvent(input$multiSelection, {
    selection <- strsplit(input$multiSelection, split = ",")[[1]]
    if(length(selection) > 1) {
      enable("combineBtn")
      enable("showCorsDialogBtn")
    }else {
      disable("combineBtn");
    }
  })
  
  observeEvent(input$heatmapButton, {
    show("pdfdiv")
    output$largePlot <- renderPlotly({
      plotType("Heatmap")
      dtype <- input$heatmapDatasetSelect
      m <- isolate(getNumericData(dtype))
      m <- m[, !apply(m, 2, function(x) all(is.na(x)))]
      m <- na.omit(m)
      cc <- isolate(sapply(tempData$DF[,colnames(m)], is.numeric))
      types <- sapply(colnames(m), function(x) strsplit(x, split = ".", fixed = T)[[1]][2])
      cc <- colnames(m)[cc & !(types %in% c("RPPA", "CNV"))]
      if(length(cc) > 0) m2 <- apply(m[,cc, drop = F], 2, function(x) log2((x+1)/mean(x+1, na.rm = T)))
      if("CNV" %in% types) {
        cnv <- m[, types == "CNV", drop = F]
        m[,types == "CNV"] <- cnv
      }
      if("RPPA" %in% types) {
        rppa <- m[, types == "RPPA", drop = F]
        rppa <- scale(rppa)
        m[,types == "RPPA"] <- rppa
      }
      if(length(cc) == 1) m[,cc] <- m2
      if(length(cc) > 1) m[,cc] <- m2[,cc]
      m[m > 3] <- 3
      m[m < -3] <- -3
      p <- heatmaply(t(m), colors = c("blue", "white", "red"),
                     plot_method = "plotly", limits = c(-3,3)) %>% config(displayModeBar = F)
      currentPlot(p)
      return(p)
    })
  })
  
  observeEvent(input$nodeDeleted, {
    selection <- strsplit(input$multiSelection, split = ",")[[1]]
    if(length(selection) == 0) selection <- input$clicked
    for(node in selection){
      if(node %in% colnames(tempData$DF)) tempData$DF[[node]] <- NULL
    }
    if(ncol(tempData$DF) < 5) {
      disable("showCorMatrixDialogBtn")
      disable("heatmapButton")
    }
  })
  
  observeEvent(input$clicked, {
    enable("plotButton")
    boo <- grepl(pattern = ' : ', input$clicked)
    if(boo){
      disable("binarizeBtn")
      disable("showCorsDialogBtn")
    } else {
      enable("showCorsDialogBtn")
      temp <- unlist(strsplit(input$clicked, "\\."))
      if(temp[2] == "Surv"){
        disable("binarizeBtn")
      } else {
        if(class(tempData$DF[[input$clicked]]) != "factor"){
          enable("binarizeBtn")
        } else {
          disable("binarizeBtn")
        }
      }
    }
  })
  
  observeEvent(input$plotButton, {
    clickedNode <- isolate(input$clicked)
    boo <- isolate(grepl(pattern = ' : ', clickedNode))
    show("pdfdiv")
    index <- NA
    # if boo is true the item is an edge, otherwise it's a node
    if (boo) {
      temp1 <- unlist(strsplit(clickedNode, " : "))
      source <- temp1[1]
      target <- temp1[2]
      sourcePart <- unlist(strsplit(source, "\\."))
      targetPart <- unlist(strsplit(target, "\\."))
      if(submitted()) {
        res <- tab()
        for(i in 1:length(res$lhs)) 
          if(res$lhs[i] == target && res$rhs[i] == source) index <- i
          if(!is.na(index)) show("edgePlot")
      }
      
      output$largePlot <- renderPlotly({
        sel.datasets <- disease.map[input$datasets]
        
        if(input$edgePlot != "Raw data" & !is.na(index)){
          
          if(input$edgePlot == "Beta values") dat <- res$stdmat[index,]
          if(input$edgePlot == "Z-scores") dat <- res$zmat[index,]
          if(input$edgePlot == "P-values (-log10)") dat <- -log10(res$pmat[index,])
          plotType("Bar")
          p <- plot_ly(type = "bar", x = sel.datasets, y = dat[sel.datasets]) %>% 
            layout(title = input$edgePlot, yaxis = list(title = paste(source, "->", target, sep = ""))) %>% config(displayModeBar = F)
          currentPlot(p)
          return(p)
        }
        
        #        if(class(tempData[["BRCA"]][[sourcePart[2]]][,sourcePart[1]]) == "factor"){
        #          xx <- factor(xx, labels = levels(tempData[["BLCA"]][[sourcePart[2]]][,sourcePart[1]]))
        #        }
        
        if (targetPart[2] != "Surv"){
          df <- data.frame(x = tempData$DF[tempData$DF$Type %in% sel.datasets, source], 
            y = tempData$DF[tempData$DF$Type %in% sel.datasets, target], 
            z = tempData$DF[tempData$DF$Type %in% sel.datasets, "Type"])
        } else {
          df <- data.frame(OS_time = tempData$DF[tempData$DF$Type %in% sel.datasets, "OS_time"],
                           OS = tempData$DF[tempData$DF$Type %in% sel.datasets, "OS"],
                           x = tempData$DF[tempData$DF$Type %in% sel.datasets, source], 
                           z = tempData$DF[tempData$DF$Type %in% sel.datasets, "Type"])
          df <- na.omit(df)
        }
        
        df$z <- factor(df$z)
        
        if(class(df$x) == "factor" && class(df$y) == "factor"){
          plotType("Factor")
          cc <- dplyr::count(df, x, y)
          cc2 <- left_join(cc, dplyr::count(cc, x, wt = n))
          cc2 <- mutate(cc2, prop = n/nn)
          if(input$factorType == "Total numbers"){
            cc2 <- mutate(cc2, res = n)
            tit <- "Number of samples"
          } else {
            cc2 <- mutate(cc2, res = prop)
            tit <- "Proportion of"
          }
          
          if(targetPart[2] == "GISTIC") cols <- cols.GISTIC else cols <- cols.factor
          
          p <- cc2  %>% mutate(prop = n/nn) %>% plot_ly(x = ~ x, y = ~res, color = ~y, colors = cols) %>% 
            add_bars() %>% layout(barmode = "stack", xaxis = list(title = source), 
                                  yaxis = list(title = paste(tit, target))) %>% config(displayModeBar = F)
          currentPlot(p)
          return(p)
        } else if ((class(df$x) == "factor" | class(df$y) == "factor") & targetPart[2] != "Surv") {
          plotType("Box")
          bx <- ifelse(input$showPointsCheck, "all", "none")
          opts <- ifelse(input$yaxislog, "log", "linear")
          if(sourcePart[2] == "GISTIC" | targetPart[2] == "GISTIC")
            cols <- cols.GISTIC else cols <- cols.factor
          if(class(df$y) == "factor"){
            p <- plot_ly(df, x = ~z, y = ~x, color = ~y, type = "box", boxpoints = bx, jitter = 0.3, pointpos = 0, colors = cols) %>%
                layout(boxmode = "group", xaxis = list(title = target), yaxis = list(title = source, type = opts)) %>% 
                config(displayModeBar = F)
            currentPlot(p)
            return(p)
          } else {
            p <- plot_ly(df, x = ~z, y = ~y, color = ~x, type = "box", boxpoints = bx, jitter = 0.3, pointpos = 0, colors = cols) %>%
                layout(boxmode = "group", xaxis = list(title = source), yaxis = list(title = target, type = opts)) %>% config(displayModeBar = F)
            currentPlot(p)
            return(p)
          }
        } else if (targetPart[2] == "Surv"){
          plotType("Surv")
          if(sourcePart[2] == "GISTIC") cols <- cols.GISTIC[levels(df$x)] else cols <- cols.factor
          if(class(df$x) != "factor"){
            enable("SurvOptions")
            x <- as.numeric(df$x)
            
            if(input$numberGroupsSurv == 2){
              probs <- survSlider1()
              labels <- c("Low", "High")
            } else {
              probs = survSlider2()
              labels <- c("Low", "Med", "High")
            }
            
            x <- cut(x, breaks = c(-Inf, quantile(x, probs = probs/100, na.rm = T), Inf), include.lowest = T, 
                     labels = labels)
            df$x <- x
          } else disable("SurvOptions")
          
          fit <- npsurv(Surv(OS_time, OS) ~ x, data = df)
          p <- survplotp(fit, col = cols) %>% config(displayModeBar = F)
          currentPlot(p)
          return(p)
        } else {
          plotType("Scatter")
          xopts <- ifelse(input$xaxislog, "log", "linear")
          yopts <- ifelse(input$yaxislog2, "log", "linear")
          if(input$chartType == "Scatter"){
            p <- plot_ly(df, x = ~x, y = ~y, color = ~z, type = "scatter") %>% 
              layout(xaxis = list(title = source, type = xopts), yaxis = list(title = target, type = yopts)) %>% config(displayModeBar = F)
          } else {
            xbins <- cut(df$x, breaks = 50)
            ybins <- cut(df$y, breaks = 50)
            counts <- with(df, table(xbins, ybins))
            zmin <- min(counts)
            zmax <- mean(counts) + 2*sd(counts)
            p <- plot_ly(df, x = ~x, y = ~y, z = ~counts) %>% 
              add_histogram2d(zmin = zmin, zmax = zmax) %>%
              layout(xaxis = list(title = source, type = xopts), yaxis = list(title = target, type = yopts)) %>% config(displayModeBar = F)
          }
          currentPlot(p)
          return(p)
        }
      })
      
    }
    # Plot for a node
    else{
      hide("edgePlot")
      output$largePlot <- renderPlotly({
        sel.datasets <- disease.map[input$datasets]
        temp2 <- unlist(strsplit(clickedNode, "\\."))
        nd <- clickedNode
        
        if(temp2[2] == "Surv"){
          plotType("Surv")
          disable("SurvOptions")
          fit <- npsurv(Surv(OS_time,OS) ~ Type, data = tempData$DF[tempData$DF$Type %in% sel.datasets,])
          p <- survplotp(fit) %>% config(displayModeBar = F)
          currentPlot(p)
          return(p)
        } else {
          df <- data.frame(x = tempData$DF[tempData$DF$Type %in% sel.datasets, nd], 
                           z = tempData$DF[tempData$DF$Type %in% sel.datasets, "Type"])
          df$z <- factor(df$z)
          if (class(df$x) == "factor") {
            plotType("Factor")
            
            cc <- dplyr::count(df, z, x)
            cc <- na.omit(cc)
            cc2 <- left_join(cc, dplyr::count(cc, z, wt = n))
            cc2 <- dplyr::mutate(cc2, prop = n/nn)
            
            if(input$factorType == "Total numbers"){
              tit = "Number of samples"
              cc2$res <- cc2$n
            } else {
              tit = "Proportion of samples"
              cc2$res <- cc2$prop
            }
            
            if(temp2[2] == "GISTIC") {
              cols <- cols.GISTIC[levels(df$x)]
            } else cols <- c("Gray", "Orange", "Green", "Blue", "Red")
            p <- plot_ly(cc2, x = ~z, y = ~res, name = ~x, color = ~x, colors = cols) %>% add_bars() %>%
              layout(yaxis = list(title = tit), barmode = "stack") %>% config(displayModeBar = F)
            currentPlot(p)
            return(p)
          }else {
            plotType("Box")
            bx <- ifelse(input$showPointsCheck, "all", "none")
            opts <- ifelse(input$yaxislog, "log", "linear")
            selectedNode(input$clicked)
            p <- plot_ly(df, y = ~x, color =~ z, type = "box", boxpoints = bx, jitter = 0.3, pointpos = 0) %>%
                 layout(yaxis = list(title = nd, type = opts)) %>% config(displayModeBar = F)
            currentPlot(p)
            return(p)
          }
        }
      })
    }
  })
  
  tab <- eventReactive(input$graph, {
    graph <- input$graph
    mat <- getTable(graph)
    
    mat <- t(mat)
    nodes <- union(mat[, 1], mat[, 2])
    a <- setdiff(nodes, colnames(tempData$DF))
    if(length(a) > 0) tempData$DF[a] <- NA
    toggle("loading2")
    parameters <- lapply(levels(tempData$DF$Type), function(x){
      tryCatch({
        ftab <- makeFTableFun(nodes, tempData$DF[tempData$DF$Type == x,])
        
        mat2 <- factorizeMat(mat, 1, factors2levels, ftab)
        mat2 <- factorizeMat(mat2, 2, factors2levels, ftab)
        
        model <- make.formula(mat2, ftab)
        measures <- rep(NA, 7)
        names(measures) <- c("df", "cfi", "tli", "aic.sem", "srmr", "ntotal", "aic.cox")
        if(model != ""){
          fit <- sem(model, data = ftab, control = list(iter.max = 250))
          params <- parameterEstimates(fit, standardized = T, ci = F)
          measures[1:6] <- tryCatch(round(fitMeasures(fit, 
                                                      fit.measures = c("df", "cfi", "tli", "aic", "srmr", "ntotal")),
                                          digits = 3), error = function(e) rep(NA, 6))
          measures["aic.sem"] <- round(measures["aic.sem"], digits = 0)
          
          ftab2 <- na.omit(ftab)
          pred <- predict(fit, newdata= ftab2)
          rownames(pred) <- rownames(ftab2)
          if(ncol(pred) > 0) {
            ftab[colnames(pred)] <- NA
            ftab[rownames(pred), colnames(pred)] <- pred
            tempData$DF[rownames(pred), colnames(pred)] <- pred
          }
        } else params <- NA
        
        surv <- calculateSurvival(mat2, ftab, tempData$DF)
        
        if(is.na(params) & is.null(surv)) return(NA)
        
        if(!is.na(params) & !is.null(surv)) {
          params <- rbind(params[,colnames(surv$params)], surv$params)
          measures["aic.cox"] <- surv$aic
        } else if(!is.null(surv)) {
          params <- surv$params
          measures["aic.cox"] <- surv$aic
        }
        
        params <- collapseParams(mat, params, factors2levels)
        
        return(list(params = params, measures = measures))
      },
      error = function(e) NA)
    })
    toggle("loading2")
    
    tryCatch({
      names(parameters) <- disease.map
      
      z.mat <- sapply(parameters, function(x) 
        if(!is.na(x)) {
          x <- x$params
          m <- x[x[,2] != "~~", ] 
          v <- apply(mat, 1, function(y) {
            num <- as.numeric(as.vector(m[m[,1] == y[2] & m[,3] == y[1], "z"]))
            if(length(num) == 0) num <- NA
            return(num)
          })
          return(v)
        } 
        else rep(NA, nrow(mat)))
      p.mat <- sapply(parameters, function(x) 
        if(!is.na(x)) {
          x <- x$params
          m <- x[x[,2] != "~~", ] 
          v <- apply(mat, 1, function(y) {
            num <- as.numeric(as.vector(m[m[,1] == y[2] & m[,3] == y[1], "pvalue"]))
            if(length(num) == 0) num <- NA
            return(num)
          })
          return(v)
        } 
        else rep(NA, nrow(mat))
      )
      std.mat<-sapply(parameters, function(x) 
        if(!is.na(x)) {
          x <- x$params
          m <- x[x[,2] != "~~", ] 
          v <- apply(mat, 1, function(y) {
            num <- as.numeric(as.vector(m[m[,1] == y[2] & m[,3] == y[1], "std.all"]))
            if(length(num) == 0) num <- NA
            return(num)
          })
          return(v)
        } 
        else rep(NA, nrow(mat))
      )
      
      
      if(is.vector(z.mat)) z.mat <- matrix(z.mat, nrow = 1, dimnames = list(1, names(z.mat)))
      if(is.vector(p.mat)) p.mat <- matrix(p.mat, nrow = 1, dimnames = list(1, names(p.mat)))
      if(is.vector(std.mat)) std.mat <- matrix(std.mat, nrow = 1, dimnames = list(1, names(std.mat)))
      
      colnames(z.mat) <- names(parameters)
      colnames(p.mat) <- names(parameters)
      colnames(std.mat)<-names(parameters)
      
      avg.z <- apply(z.mat, 1, function(x) mean(x, na.rm = T))
      avg.p <- apply(p.mat, 1, function(x) mean(x, na.rm = T))
      avg.std <-apply(std.mat, 1, function(x) mean(x, na.rm = T))
      
      z.mat <- cbind(z.mat, Average = avg.z)
      p.mat <- cbind(p.mat, Average = avg.p)
      std.mat <- cbind(std.mat, Average = avg.std)
      
      count <- 0
      for(i in parameters) if(!is.na(i$params)) break else count <- count+1;
      if(count == length(parameters)) return(NULL);
      lhs = mat[,2]
      rhs = mat[,1]
      pars <- list(lhs = lhs, rhs = rhs,
                   zmat = z.mat, pmat = p.mat, stdmat = std.mat) 
      isolate(modelCount(modelCount()+1))
      n <- paste("Model", modelCount())
      fits <- sapply(parameters, function(x) 
        tryCatch(c(Model = modelCount(), x$measures), error = function(e) rep(NA, 8)))
      
      rownames(fits) <- toupper(rownames(fits))
      modelsList[[n]] <- list()
      modelsList[[n]][["pars"]] <- pars
      modelsList[[n]][["measures"]] <- fits
      modelsList[[n]][["tempData"]] <- tempData
      js$copyGraph();
      if(!submitted()){
        isolate({
          enable("edgeLabel")
          enable("datasetSelect")
          enable("modelSelectMenu")
        })
        submitted(TRUE)
      }
      return(pars)
    }, error = function(e) return(NULL))
    
  })
  
  observeEvent(tab(), {
    
    resTable <- tab()
    p <- resTable
    if(is.null(p)) return(NULL);

    isolate({
##      updateSelectInput(session, "edgeLabel", selected = "z")
##      updateSelectInput(session, "datasetSelect", selected = "Average")
      updateSelectInput(session, "modelSelectMenu", choices = names(modelsList), 
                        selected = paste("Model", modelCount()))
    })
    
    updateGraph(p)
  })
  
  #prevent entering node name until type of the node selected
  observeEvent(input$shape, {
    if (input$shape == "") {
      disable("dynamic")
      disable("nodeButton")
    }
    else{
      enable("dynamic")
      enable("nodeButton")
    }
    
  })
  
  updateGraph <- function(p){
    switch(edgeLabelType(),
           "z" = values <- p[["zmat"]][,chosenDataset()],
           "std.all" = values <- p[["stdmat"]][,chosenDataset()],
           "p-value" = values <- p[["pmat"]][, chosenDataset()],
           "Nothing" = values <- rep("", length(p[["lhs"]]))
    )
    
    df <- data.frame(lhs = p$lhs, rhs = p$rhs, zmat = p$zmat[, chosenDataset()], 
                     std = p$stdmat[,chosenDataset()], pmat = p$pmat[, chosenDataset()], val = values);
    
    js$results(
      df
    )
  }
  
  # output$header <- renderUI({
  #   if (is.null(chosenDataset()))
  #     h1("Please choose a dataset to begin...")
  #   else
  #     h1(paste("The chosen dataset: ", chosenDataset()))
  # })
  
  
}


shinyApp(ui = ui, server = server, options = list(port = 80, host = "10.200.42.161"))