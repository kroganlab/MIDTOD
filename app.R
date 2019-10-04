# load required packages
suppressPackageStartupMessages(library(package = "shiny"))
suppressPackageStartupMessages(library(package = "shinyjs"))
suppressPackageStartupMessages(library(package = "shinyFiles"))
# load midtod internal functions
source(file = "midtod.R")

# expand size of file loaded using the shiny interface
options(shiny.maxRequestSize = 30*1024^2)

# user interface
ui <- fluidPage(
  tags$div(
    h1("MIDTOD")
  ),

  fileInput(inputId = "resultFile",
	    label   = h4("Result file")),
  fileInput(inputId = "evidenceFile",
	    label   = h4("Evidence file"),),
  selectInput(inputId  = "species",
	      label    = h4("Species"), 
	      choices  = list("Human" = "human", "Mice" = "mouse"),
	      selected = "human"),
  h4("Output directory"),
  shinyDirButton(id    = "outputDir",
		 label = "Output directory",
		 title = "Select"),
  br(),
  br(),
  useShinyjs(),
  actionButton("action", label = "Run"),
  hr(),
  textOutput("text")
)

server <- function(input, output){
  shinyDirChoose(
    input = input,
    id    = "outputDir",
    roots = c(wd = getwd()))

  global <- reactiveValues(datapath = getwd())
  
  dir <- reactive(input$dir)
  
  output$dir <- renderText({
    global$datapath
  })

  observeEvent(ignoreNULL  = TRUE,
               eventExpr   = {
                 input$dir
               },
               handlerExpr = {
                 if (!"path" %in% names(dir()))
		   return()
                 home <- normalizePath("~")
                 global$datapath <-
                   file.path(home,
			     paste(unlist(dir()$path[-1]),
				   collapse = .Platform$file.sep))
               })
  
  observeEvent(input$action, {
    withCallingHandlers({
      html(id = "text", "")
      midtod(input$resultFile$datapath,
	     input$evidenceFile$datapath,
	     input$species,
	     input$outputDir$path[[2]])
    },
    message = function(m) {
      html(id = "text", html = paste0(m$message, "<br>"), add = TRUE)
    })
  })
}
shinyApp(ui = ui, server = server)
