#'@name DeltaDeltaCqAll
#'@aliases DeltaDeltaCqAll
#'@title Compute delta delta Cq, standard deviation of delta delta Cq and fold concentration of all genes 
#'at one time.
#'@description Compute delta delta Cq, standard deviation of delta delta Cq and fold concentration of an 
#'object of \code{"RTqPCRBatch"}.
#'@param x An object of class \code{"RTqPCRBatch"} produced as an output of DeltaCq function.
#'@param \dots other parameters to be passed to the downstream methods.
#'@return An object of \code{"RTqPCRBatch"} with new slots.
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'#'@examples 
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
#'deltaCq
#'
#'## For more details see the vignettes
#'deltadeltaCqall <- DeltaDeltaCqAll(deltaCq)
#'deltadeltaCqall
#'@export 
setMethod("DeltaDeltaCqAll", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, ...)
            {
              x <- as.data.frame(fData(RTqPCRBatch))
              Target <- x[,"Target name"]
              fac <- unique(Target)
              deltadeltaCq <- c()
              for (i in 1:length(fac)){
                deltadeltaCq <- append(DeltaDeltaCq(RTqPCRBatch, Target = fac[i]),i)
              }
              x2 <- as.data.frame(deltadeltaCq$ID)
              x3 <- as.data.frame(deltadeltaCq$deltadeltaCq)
              x4 <- as.data.frame(deltadeltaCq$"Fold Conc.")
              x5 <- as.data.frame(deltadeltaCq$sd.deltadeltaCq)
              ddCq <- cbind(x2,x3,x4,x5)
              row.names(ddCq) <- 1:length(row.names(ddCq))
              names(ddCq) <- c("ID","deltadeltaCq", "Fold Conc.[log]", "sd.deltadeltaCq")
              ddCq[,"ID"] <- factor(ddCq[,"ID"],levels = ddCq[,"ID"])
              metaData <- data.frame(labelDescription = names(ddCq),
                                     row.names = names(ddCq), check.names = FALSE,
                                     stringsAsFactors = FALSE)
              fData <- AnnotatedDataFrame(data = ddCq, varMetadata = metaData)
              ddCq.table <- new("RTqPCRBatch", featureData = fData)
              return(ddCq.table)  
            }
          
          (RTqPCRBatch,...)
)
              