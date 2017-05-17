#'@name read.RTqPCR
#'@aliases read.RTqPCR
#'@title  To read in the raw fluorescent data of LC480 light cycler and Mx3005P RT-qPCR experiment.
#'@description Function \code{read.RTqPCR} reads in the raw fluorescent data of LC480 light cycler and Mx3005P
#'RT-qPCR and use the data to populate an object of class \code{"RTqPCRBatch"}. 
#'@param file the name of the file to read in.
#'@param PCRtype the type of RT-qPCRs ("LC480" or "Mx3005P") whose file is to read in.
#'@param \dots Other parameters to be passed to downstream methods.
#'@return \code{"RTqPCRBatch"} object with exprs slots.
#'@details Function \code{read.RTqPCR} reads in the raw fluorescent data of LC480 light cycler
#'and Mx3005P RTqPCR and is based on functions \code{\link[ReadqPCR]{read.LC480}} and \code{\link{read.Mx3005P}}.
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@examples
#'##To read in the fluorescent data of LC480 light cycler
#'
#'LC480.example <- file.path(path, "LC480_Example.txt") 
#'rtData.LC480 <- read.RTqPCR(LC480.example, PCRtype = "LC480")
#'rtData.LC480 ## to visualise the overview of data
#'head(exprs(rtData.LC480)) ## to visualise all of the fluorescence data
#'head(exprs(rtData.LC480[,1:5])) ## to visualise the first five columns of fluorescence data
#'head(pData(rtData.LC480)) ## to express the phenoData 
#'head(fData(rtdata.LC480)) ## to express featureData
#'
#'##To read in the fluorescent data of Mx3005P RT-qPCR
#'
#'Mx3005P.example <- file.path(path, "Mx3005P_Example.txt") 
#'rtData.Mx <- read.RTqPCR(Mx3005P.example, PCRtype = "Mx3005P")
#'rtData.Mx ## to visualise the overview of data
#'head(exprs(rtData.Mx)) ## to visualise all of the fluorescence data
#'head(exprs(rtData.Mx[,1:5])) ## to express the first five columns of fluorescence data
#'head(pData(rtData.Mx)) ## to express the phenoData
#'head(fData(rtData.Mx)) ## to express the featureData
#'
#'@export 
read.RTqPCR <- function (file, PCRtype, ...)
{
  if (PCRtype == "LC480")
  {
    x <- read.LC480(file = file)
    x.exprs <- exprs(x)
    x.featureData <- featureData(x)
    x.phenoData <- phenoData(x)
    cycData <-  new("RTqPCRBatch", exprs = x.exprs, featureData = x.featureData, phenoData = x.phenoData)
    }
    
  else
  {
    if (PCRtype == "Mx3005P")
    {
      cycData <- read.Mx3005P(file = file)
    }
    else
    {
      stop("Define appropriate PCRtype out of LC480 or Mx3005P!")
    }
  }
  return(cycData)
}