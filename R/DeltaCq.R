#'@name DeltaCq
#'@aliases DeltaCq
#'@title Compute delta Cq value and standard deviation of delta Cq
#'@description Compute delta Cq and standard deviation of delta Cq of object of \code{"RTqPCRBatch"}.
#'@param RTqPCRBatch An object of class \code{"RTqPCRBatch"} produced as an output of CombineTechReps function.
#'@param Ref Name of reference gene under column Target name.
#'@param \dots other parameters to be passed to the downstream methods.
#'@return An object of \code{"RTqPCRBatch"} with new slots.
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@examples 
#'
#'## Read in the raw fluorescent data
#'
#'## Read in the sample information data
#'
#'## Merge the fluorescence and sample information data through merge function
#'
#'## Compute the Cq values and amplification efficiencies
#'
#'##  Implement the functions ReplaceNAs, ReplaceAboveCutOff, ReplaceValue and NonDetects as per requirement.
#'
#'## Combine technical replicates (based on Cq values) : see also CombineTechReps function and vignettes
#'
#'Cqreps <- CombineTechReps(Cqeffs.LC480)
#'Cqreps   ## to visualise the overview of the resulting data
#'
#'## To compute the delta Cq value
#'deltaCq <- DeltaCq(Cqreps, Ref = "hk")
#'
#'## For more details see the vignettes
#'
#'@export 
setMethod("DeltaCq", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, Ref, ...)
            { 
                  x <- as.data.frame(fData(RTqPCRBatch))
                  if (is.null(Ref))
                  {
                    warning("Please provide Reference gene. Otherwise, all reference values 
                            will be considered")
                    ref <- as.matrix(x[x$"Combined sample and target type" == "Ref Unknown", ])
                  }
                  
                  else
                  {
                    ref <- as.data.frame(x[x$"Target name" == Ref, ])
                    ref <- as.matrix(ref[ref$"Combined sample and target type" == "Ref Unknown", ])                   
                  }
                  
                  ref.Cq <- mean((as.numeric(as.vector(ref[,"Cq"]))), na.rm = TRUE)
                  ref.sd.Cq <- mean((as.numeric(as.vector(ref[,"sd.Cq"]))), na.rm = TRUE)
                  target.Cq <- as.numeric(as.matrix(x[,"Cq"]))
                  target.sd.Cq <- as.numeric(as.matrix(x[,"sd.Cq"]))
                  
                  deltaCq <- (target.Cq) - (ref.Cq)
                  is.na(deltaCq) <- NA
                  deltaCq <- as.matrix(format(deltaCq, digit = 3))
                  
                  sd.deltaCq <- sqrt((target.sd.Cq)^2 + (ref.sd.Cq)^2)      
                  is.na(sd.deltaCq) <- NA
                  sd.deltaCq <- as.matrix(format(sd.deltaCq, digit = 3))
                   
                  ID <- as.matrix(x$ID)
                  Target.Name <- as.matrix(x[,"Target name"])
                  CSTT <- as.matrix(x[,"Combined sample and target type"]) 
                  
                  deltaCq.table <- as.data.frame(cbind(ID, deltaCq, sd.deltaCq,Target.Name, CSTT))
                  names(deltaCq.table) <- c("ID", "deltaCq", "sd.deltaCq", "Target name",
                                            "Combined sample and target type")
                  deltaCq.table1 <- deltaCq.table[deltaCq.table$"Combined sample and target type" 
                                                                == "Target Unknown", ]
                  deltaCq.table2 <- deltaCq.table[deltaCq.table$"Combined sample and target type"
                                                  == "Target PosCalibrator", ]
                  deltaCq.table3 <- as.data.frame(rbind(deltaCq.table1, deltaCq.table2))
                  deltaCq.table3[,"ID"] <- factor(deltaCq.table3[,"ID"],levels = deltaCq.table3[,"ID"])
                  metaData <- data.frame(labelDescription = names(deltaCq.table3),
                                         row.names = names(deltaCq.table3), check.names = FALSE,
                                         stringsAsFactors = FALSE)
                  fData <- AnnotatedDataFrame(data = deltaCq.table3, varMetadata = metaData)
                  deltaCq.table4 <- new("RTqPCRBatch", featureData = fData)
                
                
                    return (deltaCq.table4)
                }
                    
                (RTqPCRBatch, ...)
            )

