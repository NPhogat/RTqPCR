#'@name RTqPCR.dataframe
#'@aliases RTqPCR.dataframe
#'@title Convert the object of Class \code{"data.frame"} into an object of Class \code{"RTqPCRBatch"}.
#'@description Convert the object of Class \code{"data.frame"} into an object of Class \code{"RTqPCRBatch"}.
#'@param x An object of Class \code{"data.frame"}. 
#'@param fun Name of function whose output is to be converted into RTqPCRBatch.
#'@param \dots other parameters to be passed to the downstream methods.
#'@return An object of Class \code{"RTqPCRBatch"} with new slots.
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@keywords data.frame RTqPCRBatch
#'@export 
setMethod("RTqPCR.dataframe", signature = "data.frame", definition = 
            function(x, fun, pcrtype, ...)
            {
              if ((fun == "read.RTqPCR")|| (fun == "merge"))
              {
                if (pcrtype == "LC480"){
                  fData <- as.data.frame(cbind(x[,c("Cycle number","Acquisition time", 
                                                    "Acquisition temperature")]))
                  drops <- c("Cycle number","Acquisition time","Acquisition temperature")
                  
                }
                else if (pcrtype == "Mx3005P"){
                  fData <- as.data.frame(cbind(x[,c("Cycle#","Temperature")]))
                  drops <- c("Cycle#", "Temperature")
                }
                else
                {
                  stop ("please check the pcrtype!")
                }
                fMetaData <- data.frame(labelDescription = names(fData),
                                        row.names = names(fData), check.names = FALSE,
                                        stringsAsFactors = FALSE)
                featureData <- AnnotatedDataFrame(data = fData, varMetadata = fMetaData)
                
                fluoData <- x[ , !(names(x) %in% drops)]
                fluoData1 <- as.matrix(fluoData)
                names(fluoData1) <- names(fluoData)
                
                #### phenoData
                Sampos <- colnames(fluoData1)
                samPos <- factor(Sampos, levels = unique(Sampos))
                if (pcrtype == "LC480")
                {
                  pData <- data.frame("Sample position" = samPos, row.names = samPos, check.names = FALSE,
                                      stringsAsFactors = FALSE)
                }
                else if (pcrtype == "Mx3005P")
                {
                  pData <- data.frame("Well" = samPos, row.names = samPos, check.names = FALSE,
                                      stringsAsFactors = FALSE)
                }
                
                else
                {
                  stop ("please check the pcrtype!")
                }
                pMetaData <- data.frame(labelDescription = names(pData),
                                        row.names = names(pData), check.names = FALSE,
                                        stringsAsFactors = FALSE)
                phenoData = AnnotatedDataFrame(data = pData, varMetadata = pMetaData)
                
                raw.data <- new("RTqPCRBatch", exprs = fluoData1, featureData = featureData, phenoData = phenoData)
                
              }
                
              else if ((fun == "CqValues")|| (fun == "NonDetects") || (fun == "ReplaceNAs") || (fun == "ReplaceValue") 
              || (fun == "ReplaceAboveCutOff"))
              {
              x.exprs <- as.matrix(x[,"exprs"])
              x.effs <- as.matrix(x[,"effs"])
              x1 <- x[-1] 
              x2 <- x1[-1]
              if (pcrtype == "LC480")
              {  
              x2[,"Sample position"] <- factor(x2[,"Sample position"],levels = x2[,"Sample position"])
              }
              else if (pcrtype == "Mx3005P")
              {
                x2[,"Well"] <- factor(x2[,"Well"], levels = x2[,"Well"])
              }
              else
              {
                stop ("Please check the pcrtype!")
              }
              metaData <- data.frame(labelDescription = names(x2),
                                     row.names = names(x2), check.names = FALSE,
                                     stringsAsFactors = FALSE)
              fData <- AnnotatedDataFrame(data = x2, varMetadata = metaData)
              effs <- assayDataNew("environment", effs = x.effs)
              res <- new("RTqPCRBatch", exprs = x.exprs, featureData = fData,
                         effs = x.effs)
            }
            else
              {
              if (fun == "CombineTechReps")
              {
                x[,"ID"] <- factor(x[,"ID"],levels = x[,"ID"])
                metaData <- data.frame(labelDescription = names(x),
                                       row.names = names(x), check.names = FALSE,
                                       stringsAsFactors = FALSE)
                fData <- AnnotatedDataFrame(data = x, varMetadata = metaData)
                rep.effs <- new("RTqPCRBatch", featureData = fData) 
              }  
            }
            
            }
)
