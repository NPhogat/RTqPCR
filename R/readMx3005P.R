#'@name read.Mx3005P
#'@aliases read.Mx3005P
#'@title To read in the raw fluorescent data of Mx3005P (Agilent) RT-qPCR experiments
#'@description Function \code{read.Mx3005P} reads in txt files of raw fluorescent data of Mx3005P 
#'RT -qPCR and uses the data to populate an object of class \code{"RTqPCRBatch"}.
#'@param file the name of the file to read in.
#'@param colNames a character vector of names to be assumed for the columns.
#'@param cycleThreshold maximum number of cycles which will be read in.
#'@param fileType the type of the file.
#'@param skip an integer: the number of lines of the data file to skip before beginning to read data.
#'@param header a logical value indicating whether the file contains the names of the variables 
#'as its first line. If missing, the value is determined from the file format: header is set to TRUE
#'if and only if the first row contains one fewer field than the number of columns.
#'@param sep the field separator character. Values on each line of the file are separated by this character.
#'See \code{\link[utils]{read.table}}.
#'@param quote the set of quoting characters. To disable quoting altogether, use \code{quote = ""}.
#'See \code{\link[base]{scan}} for the behaviour on quotes embedded in quotes. Quoting is only considered 
#'for columns read as character, which is all of them unless colClasses is specified.
#'@param dec the character used in the file for decimal points.
#'@param fill logical. If TRUE then in case the rows have unequal length, blank fields are implicitly added.
#'See \code{\link[utils]{read.table}}.
#'@param comment.char character: a character vector of length one containing a single character or an empty string. 
#'Use \code{""} to turn off the interpretation of comments altogether.
#'@param \dots other parameters to be passed to downstream methods.
#'@return object of Class \code{"RTqPCRBatch"} with exprs slots.
#'@details Allows the user to read in qPCR raw fluorescence data from Mx3005P RT-qPCR which has been exported 
#'to a txt-file, alongside phenotypic data. More details on how to read in raw fluorescence data from RT-qPCR
#'or Light Cycler in the \code{RTqPCR} package vignette.
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@keywords CyclesSet
#'@examples
#'path <- system.file("exData", package = "RTqPCR")
#'Mx3005P.example <- file.path(path, "Mx3005P_Example.txt")
#'cycData <- read.Mx3005P(Mx3005P.example)
#'cycData[1:5] ##to visualise the first 5 cycdata and for all data remove [1:5] from the code
#'featureData(cycData)  ## to visualise the overview of the featureData (feature information of data) 
#'fData(cycData)[1:5] ## To visualise the content of the first five featureData
#'## fData(cycData)   To visualise content of all featureData
#'phenoData(cycData) ## To visualise the overview of phenoData (phenotypic information of data)
#'pData(cycData)[1:5] ## To visualise the content of the first five phenoData
#'## pData(cycData) To visualise the content of all phenoData
#'@export 
####################################################################
## Read Experiment Text File of Mx3005P into RTqPCRBatch
####################################################################
 
                                                                   
                       ## read data from Mx3005P
                       read.Mx3005P <- function(file, colNames = c("mxpSegment","Ramp/Plateaue", "Ramp/Plateaue#",
                                                                  "Well", "Dye",
                                                                  "Cycle#", "Fluorescence",
                                                                  "Temperature"),
                                          cycleThreshold = 45, fileType = "txt", skip = 1,
                                          header = TRUE, sep = "\t", quote = "\"", dec = ".",
                                          fill = TRUE, comment.char = ""){
                               ## txt file
                               if(fileType == "txt")
                               res <- readMx3005Ptxt(file = file, colNames = colNames,
                                          cycleThreshold = cycleThreshold, skip = skip,
                                          header = header, sep = sep, quote = quote, dec = dec,
                                          fill = fill, comment.char = comment.char)
                              else
                              res <- readMx3005Pxml(file = file, colNames = colNames, dec = dec,
                                         cycleThreshold = cycleThreshold)
                              res
                                }
                         
                              ## not exported function for reading xml files from LC480
                               readMx3005Pxml <- function(file, colNames, dec, cycleThreshold){
                                          stop ("This function is not yet implemented!")
                                        }
                                       
                            ## not exported function for reading text files from LC480 
                               readMx3005Ptxt <- function(file, colNames, cycleThreshold, header, sep, quote,
                                                   dec, fill, comment.char, skip){
                               allData <- read.table(file = file, header = header, sep = sep, quote = quote,
                                         dec = dec, fill = fill, comment.char = comment.char,
                                         skip = skip, stringsAsFactors = FALSE)
                            ## remove empty columns
                               allData <- allData[,colSums(is.na(allData)) != nrow(allData)]
                            ## rename of columns
                            if(ncol(allData) == length(colNames))
                              colnames(allData) <- colNames
                            else
                            stop("Number of non-empty columns not equal to length of 'colNames'!")
                            ## split data by Sample position
                            fac <- factor(allData[,"Well"],
                            levels = unique(allData[,"Well"]))
                            allData.split <- split(x = allData, f = fac)
                            ## extract Fluorescence data
                            fluoData <- sapply(allData.split, function(x) x[seq_len(cycleThreshold),"Fluorescence"])
                           rownames(fluoData) <- as.character(seq_len(cycleThreshold))
                          ## feature Data (here: cycles)
                           acTemp <- allData.split[[1]][seq_len(cycleThreshold),"Temperature"]
                           fData <- data.frame("Cycle#" = seq_len(cycleThreshold),
                                   "Temperature" = acTemp,
                                   row.names = seq_len(cycleThreshold),
                                   check.names = FALSE,
                                  stringsAsFactors = FALSE)
                                  fMetaData <- data.frame(labelDescription = c("Cycle#",
                                                      "Temperature"),
                                  row.names = names(fData), check.names = FALSE,
                                  stringsAsFactors = FALSE)
                                  featureData <- AnnotatedDataFrame(data = fData, varMetadata = fMetaData)
                           ## phenotypic data
                            samPos <- sapply(allData.split, function(x) x[1,"Well"])
                            samPos <- factor(samPos, levels = unique(samPos))
                            samNam <- sapply(allData.split, function(x) x[1,"Ramp/Plateaue"])
                            prgNum <- sapply(allData.split, function(x) x[1,"Ramp/Plateaue#"])
                            segNum <- sapply(allData.split, function(x) x[1,"Dye"])
                            mspseg <- sapply(allData.split, function(x) x[1,"mxpSegment"])   
                            pData <- data.frame("Well" = samPos,
                                                "Ramp/Plateaue" = samNam,
                                                "Ramp/Plateaue#" = prgNum,
                                                "Dye" = segNum,
                                                "mspSegment" = mspseg, 
                                                 row.names = samPos, check.names = FALSE,
                                                 stringsAsFactors = FALSE)
                         pMetaData <- data.frame(labelDescription = c("Well",
                                                                      "Ramp/Plateaue",
                                                                      "Ramp/Plateaue#",
                                                                      "Dye", "mspSegment"),
                                                  row.names = names(pData), check.names = FALSE,
                                                  stringsAsFactors = FALSE)
                      phenoData <- AnnotatedDataFrame(data = pData, varMetadata = pMetaData)
                      new("RTqPCRBatch", exprs = fluoData, featureData = featureData,
                            phenoData = phenoData)
                          }
                  