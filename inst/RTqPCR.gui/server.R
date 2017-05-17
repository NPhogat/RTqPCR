library(shiny)
library(RTqPCR)
shinyServer(function(input,output){
  pcr <- reactive({
    input$pcrtype
  })
  effmethod <- reactive({
    input$efficiencymethod
  })
  n.det <- reactive({
    input$nondetects
  })
  r.cutoff <-  reactive({
    input$replaceabovecutoff
  })
  cq.cutoff <- reactive({
    input$Cqcutoff
  })
  effs.cutoff <- reactive({
    input$effscutoff
  })
  replace.value <- reactive({
    input$replacevalue
  })
  replace.nas <- reactive({
    input$replacenas
  })
  c.techreps <- reactive({
    input$ctechreps
  })
  ref.target <- reactive({
    input$Reference
  })
  result <- reactive({
    input$selectresult
  })
  null.input <- reactive({
    is.null(input[["file1"]]) && input[["run.example"]] == 0
  })
  
  fluo <- reactive({
   if(is.null(input$file1)){
      infile1 <- NULL
    } else
      infile1 <- read.RTqPCR(input$file1$datapath, PCRtype = pcr())
      x.exprs <- as.data.frame(exprs(infile1))
      colnames(x.exprs) <- colnames(exprs(infile1))
      x.featureData <- as(featureData(infile1),"data.frame")
      x.data.frame <- as.data.frame(cbind(x.exprs, x.featureData))
      names(x.data.frame) <- c(names(x.exprs),names(x.featureData))
      x.data.frame
  })
  
  samInfo <- reactive({
    if(is.null(input$file2)){
      infile2 <- NULL
    } else{
      infile2 <- read.RTqPCRSampleInfo(input$file2$datapath, PCRtype = pcr())
      infile2
    }
  }) 
  
  Cq <- eventReactive(input$compute,{
     runif(input$file2)
    isolate({
    if(!(is.null(samInfo)))
      {
    fluofile <- fluo()
    fluofile1 <- RTqPCR.dataframe(fluofile, fun = "read.RTqPCR", pcrtype = pcr())
    M <- merge(fluofile1, samInfo())
    res <- CqEffs(M, PCRtype = pcr(),Effmethod = effmethod(), baseline = "none")
    x1 <- exprs(res)
    x2 <- effs(res)
    x3 <- fData(res)
    x4 <- as.data.frame(cbind(x1,x2,x3))
    names(x4) <- c("exprs","effs",names(x3))
    x4
      }
    else{
      message = ("Please check fluorescence data and sample information data files are uploaded!")
    }
    })
  })
  
  Cq.cut <- reactive({
    
    x <- Cq()
    res <- RTqPCR.dataframe(x, fun = "CqValues", pcrtype = pcr())
    if (r.cutoff()=="None")
    {
      res1 <- res
    }
    else
    {
    res1 <- ReplaceAboveCutOff(res, NewVal = NA, Cqcutoff = cq.cutoff(),effscutoff = effs.cutoff())
    }
    x1 <- exprs(res1)
    x2 <- effs(res1)
    x3 <- fData(res1)
    x4 <- as.data.frame(cbind(x1,x2,x3))
    names(x4) <- c("exprs","effs",names(x3))
    x4
  })
  
  Cq.re <- reactive({
    x <- Cq.cut()
    res <- RTqPCR.dataframe(x, fun = "CqValues", pcrtype = pcr())
    if (replace.value()== "None")
    {
      res1 <- res
    }
    else {
      res1 <- ReplaceValue(res, NewCq = input$newcq, Neweffs = input$neweffs,
                           Cqrow = input$cqrow, effsrow = input$effsrow)
    }
    x1 <- exprs(res1)
    x2 <- effs(res1)
    x3 <- fData(res1)
    x4 <- as.data.frame(cbind(x1,x2,x3))
    names(x4) <- c("exprs","effs",names(x3))
    x4
  })
  
  Cq.nd <- reactive({
    x <- Cq.re()
    res <- RTqPCR.dataframe(x, fun = "CqValues", pcrtype = pcr())
    if (n.det()== "None")
    {
      res1 <- res
    }
    else {
      res1 <- NonDetects(res, Calc = n.det())
    }
    x1 <- exprs(res1)
    x2 <- effs(res1)
    x3 <- fData(res1)
    x4 <- as.data.frame(cbind(x1,x2,x3))
    names(x4) <- c("exprs","effs",names(x3))
    x4
  })
  
  Cq.rena <- reactive({
      x <- Cq.nd()
      res <- RTqPCR.dataframe(x, fun = "CqValues", pcrtype = pcr())
      if (replace.nas()== "None")
      {
        res1 <- res
      }
      else {
        res1 <- ReplaceNAs(res, NewCqNA = input$cqna, NeweffsNA = input$effsna)
      }
      x1 <- exprs(res1)
      x2 <- effs(res1)
      x3 <- fData(res1)
      x4 <- as.data.frame(cbind(x1,x2,x3))
      names(x4) <- c("exprs","effs",names(x3))
      x4
  })
  
  CqReps <- reactive({
    x <- Cq.rena()
    res <- RTqPCR.dataframe(x, fun = "CqValues", pcrtype = pcr())
    res1 <- CombineTechReps(res, calc = c.techreps())
    res2 <- as.data.frame(fData(res1))
  })
  
  effsReps <- reactive({
    x <- Cq.rena()
    res <- RTqPCR.dataframe(x, fun = "CqValues", pcrtype = pcr())
    res1 <- CombineTechReps(res, calc = c.techreps(), cRepCq = FALSE)
    res2 <- as.data.frame(fData(res1))
  })
  
  deltaCq <- reactive({
    x <- CqReps()
    res <- RTqPCR.dataframe(x, fun = "CombineTechReps", pcrtype = pcr())
    res1 <- DeltaCq(res, Ref = ref.target())
    res2 <- as.data.frame(fData(res1))
    res2
  })
  
  deltadeltaCq <- reactive({
    x <- deltaCq()
    res <- RTqPCR.dataframe(x, fun = "CombineTechReps", pcrtype = pcr())
    res1 <- DeltaDeltaCqAll(res)
    res2 <- as.data.frame(fData(res1))
    res2
  })
  
  NRQeffs <- reactive({
    x <- effsReps()
    res <- RTqPCR.dataframe(x, fun = "CombineTechReps", pcrtype = pcr())
   res1 <- NRQeffsAll(res, y = ref.target())
   res2 <- as.data.frame(fData(res1))
   res2
  })
  
  output$tabset <- renderUI({
    if(null.input()) {
      tabPanel("No input detected")       
    } 
    else {
      tabsetPanel(
        tabPanel("Cq, efficiency", tableOutput("Cq")),
        tabPanel("Cq.cutoff, effs.cutoff", tableOutput("Cq.cut")),
        tabPanel("Cq.replacevalue, effs.replacevalue", tableOutput("Cq.re")),
        tabPanel("Cq.nondetects, effs.nondetects", tableOutput("Cq.nd")),
        tabPanel("Cq.replaceNAs, effs.replaceNAs", tableOutput("Cq.rena")),
        tabPanel("Cq Rep", tableOutput("CqReps")),
        tabPanel("Effs Rep", tableOutput("effsReps")),
        tabPanel("deltaCq", tableOutput("deltaCq")),
        tabPanel("deltadeltaCq", tableOutput("deltadeltaCq")),
        tabPanel("log(Relative Expression)", tableOutput("NRQeffs"))
      )
    }
  })
  
  output$Cq <- renderTable({
    x <- Cq()
    x.result <- x[,c("exprs","effs")]
    names(x.result) <- c("Cq","Efficiency")
    x.result
  })
  
  output$Cq.cut <- renderTable({
    x <- Cq.cut()
    x.result <- x[,c("exprs","effs")]
    names(x.result) <- c("Cq","Efficiency")
    x.result
  })
  
  output$Cq.re <- renderTable({
    x <- Cq.re()
    x.result <- x[,c("exprs","effs")]
    names(x.result) <- c("Cq","Efficiency")
    x.result
  })
  
  output$Cq.rena <- renderTable({
    x <- Cq.rena()
    x.result <- x[,c("exprs","effs")]
    names(x.result) <- c("Cq","Efficiency")
    x.result
  })
  
  output$Cq.nd <- renderTable({
    x <- Cq.nd()
    x.result <- x[,c("exprs","effs")]
    names(x.result) <- c("Cq","Efficiency")
    x.result
  })
  
  output$CqReps <- renderTable({
    x <- CqReps()
    x.result <- x[,c("ID","Cq","sd.Cq")]
    x.result
  })
  
  output$effsReps <- renderTable({
    x <- effsReps()
    x.result <- x[,c("ID","Cq","sd.Cq","effs","sd.effs")]
    x.result
  })
  
  output$deltaCq <- renderTable({
    x <- deltaCq()
    x.result <- x[,c("ID","deltaCq","sd.deltaCq")]
  })
  
  output$deltadeltaCq <- renderTable({
    deltadeltaCq()
  })
  
  output$NRQeffs <- renderTable({
    NRQeffs()
  })
  
  text1 <- reactive({
    "Select option other than 'None' and fill up necessary boxes to implement the method!"
  }) 
  
  text2 <- reactive({
    "Select option other than 'None' to implement the method!"
  })
  
  output$download.Cq <- downloadHandler(
    filename  = "result_report.html",
    content <- function(file) {
      knitr:::knit(input = "result_report.Rmd", 
                   output = "result_report.md", quiet = TRUE)
      markdown:::markdownToHTML("result_report.md", "result_report.html")
    }
  )
  
  
})
