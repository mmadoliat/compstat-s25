library(shiny);

options(shiny.sanitize.errors = FALSE)

shinyUI(fluidPage(
  titlePanel("EM-Algorithm (Gaussian Mixture Model)"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file', 'Choose CSV File', accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), ','),
      numericInput('column','Pick Column:', value = 1, min = 1, max = 1),
      
      tags$hr(), numericInput('modes','# of Modes:', value = 1, min = 1, max = 10),
      tags$hr(), numericInput('step','EM-Step:', value = 0, min = 0),
      
      tags$hr(), checkboxInput('fixms', 'Fix Means and SDs (Nonparam. Dens. Est. for large # of modes)', FALSE),
      checkboxGroupInput('show', 'Show:', c("Initial Est." = "init", "Final Est." = "final", "Legend" = "legend"), inline = TRUE),
      sliderInput("bins", "Number of bins:", min = 1, max = 50, value = 30)),
    mainPanel(
      tabsetPanel(type = "tabs", 
                  tabPanel("Data", tableOutput('data')),
                  tabPanel("Model Selection", plotOutput("ic")),
                  tabPanel("Plot", plotOutput("plot")), 
                  tabPanel("Summary", tableOutput("summary"))))
    #tabPanel("Summary", verbatimTextOutput("summary"))))
  )
))
