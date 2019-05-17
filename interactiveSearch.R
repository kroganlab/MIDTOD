
## Load to_search from "metabolomics_pipeline.R"
head(to_search)


hmdb <- read.delim(hmdb_file, sep='\t', stringsAsFactors=F, header=T)
flu <- read.delim(flu_file, sep='\t', stringsAsFactors=F, header=T)

# INITIALIZE CONSTRAINTS
# any log2 fold change above this value is considered significant (also applies to the negative value in the oposite way)
log2FC = 1
# any p-value below this is considered significant
pvalue = 0.05
THRESH = 0.05  # this is the amount we are willing to let the masses be off for identification +/-
max_weight = 500 # The maximum weight (in Daltons) of metabolites in HMDB to be included in the search

# log2FC condition (column name) in results to check for significant fold change in
log2FC_condition = 'LungH1.LungPBS_log2FC'
# adjusted p-value condition (column name) in results file to check for significant p-values in
adj_pval_condition = 'LungH1.LungPBS_adj.pvalue'


#####################################################################################


server <- function(input, output) {
#   output$distPlot <- renderPlot({
#     hist(input$pval, col = 'darkgray', border = 'white')
#   }),
  output$text1 <- renderDataTable({ 
    # idx <- intersect( grep(paste(input$check_proteins,collapse="|"), x.clust$Proteins), which(x.clust$cluster %in% input$check_clusters) )
    # paste("You have chosen a Percentile range per bait of\n", paste(unique(x.clust$Proteins[idx]),collapse=","))
    
    tmp <- flu[which( (flu$q_value<input$pval) & (abs(flu$log2fc)>=input$log2fc) & (flu$condition_2 %in% input$condition_2) ),]
    # tmp <- head(flu)
    tmp
  })
}

ui <- basicPage(
  sidebarLayout(
    sidebarPanel(
      sliderInput("pval", 
                  label="Select p-value threshold:", min = 0, max = 1, value = 0.05),
      sliderInput("log2fc", 
                  label="Select Log2 Fold Change threshold:", min = 0, max = 10, value = 1),
      sliderInput("err_thresh", 
                  label="Select amount of error allowed in matching masses:", min = 0, max = 1, value = .05),
      sliderInput("maxweight", 
                  label="Select maximum weight to search against:", min = 0, max = 1000, value = 500),
      selectizeInput( 'condition_2',
                      label="Conditions to Search Against:", choices = unique(flu$condition_2), options=list(create=TRUE), multiple=TRUE )
    ),
    
    mainPanel( 
      tabsetPanel(
        tabPanel("Main output",
                 dataTableOutput("text1")
        )
      )
    )
  )
)


shinyApp(ui = ui, server = server)





