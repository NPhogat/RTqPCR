library(shiny)
shinyUI(fluidPage(
  titlePanel("Analysis of RT-qPCR data"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1","Choose .txt tab separated file of fluo data to upload"),          
      tags$hr(),
      checkboxInput("header","Header", TRUE),
      radioButtons("txt.type","Type of txt file",c(dec= ".", Sep = "\t")),
      
      fileInput("file2","Select .txt tab separated file of sample info data to upload"),          
      tags$hr(),
      checkboxInput("header","Header", TRUE),
      radioButtons("txt.type","Type of txt file",c(dec= ".", Sep= "\t")),                                             
      
      selectInput("pcrtype","Select the PCR type", choices = c("LC480","Mx3005P"), selectize = FALSE),
      br(),
      selectInput("efficiencymethod","Select efficiency computing method", choices = c("sigfit","expfit",
                                                                                "sliwin","LRE")),                                                                        
      br(),
      selectInput("replaceabovecutoff","Select to apply the cut off", choices = c("None","Replace")),
      br(),
      numericInput("Cqcutoff", "Cq cut off",""),
      numericInput("effscutoff", "efficiency cut off",""),
      br(),
      selectInput("replacevalue", "Select to replace the Cq and Eff. values", choices = c("None","Replace")),
      br(),
      numericInput("cqrow","The row, where Cq value is to be replaced",""),
      numericInput("newcq", "New Cq Value", ""),
      numericInput("effsrow", "The row, where efficiency value is to be replaced",""),
      numericInput("neweffs", "New Efficiency value",""),
      br(),
      selectInput("nondetects", "Select method to remove nondetects", choices = c("None","Mean","Median")),
      br(),
      selectInput("replacenas", "Select to replace all NAs by user defined value", choices = c("None","Replace")),
      br(),
      numericInput("cqna", "Provide the Cq value to replace NAs",""),
      numericInput("effsna","Provide the efficiency value to replace NAs", ""),
      br(),
      selectInput("ctechreps", "Select method to combine Technical Replicates", choices = c("Mean",
                                                                          "Median", "Geom")),
      br(),
      textInput("Reference", "Write the name of reference gene under column 'Target Name'", value = "",
                width = NULL, placeholder = NULL),
      br(),
      actionButton("compute","COMPUTE!"),
      p("Click on COMPUTE! button to compute the results"),
      br(),
      actionButton("run.example", "Run example"),
      br(),
      textInput("directory","Write the name of the directory to download the report file","NULL"),
      downloadButton("download.Cq", "Download Results")
      ),
    mainPanel(
      uiOutput("tabset")
      )
    )
  ))