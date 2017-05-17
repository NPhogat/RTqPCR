#'@name NRQeffsAll
#'@aliases NRQeffsAll
#'@title Compute normalized relative expression of all genes at one time.
#'@description Compute relative expression ration of all genes at one time of an 
#'object of \code{"RTqPCRBatch"}.
#'@param RTqPCRBatch An object of class \code{"RTqPCRBatch"} produced as an output of CombineTechReps function.
#'@param y Name of Reference gene under column Target name.
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
#'## Combine technical replicates - on the basis of efficiency
#'
#'effsreps <- CombineTechReps(Cqeffs.LC480, cRepCq = FALSE)
#'
#'##To compute the relative expression ratio
#'nrqeffsall <- NRQeffsAll(effsreps, y = "hk")
#'nrqeffsall
#'
#'@export 
setMethod("NRQeffsAll", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, y, ...)
            {
              x <- as.data.frame(fData(RTqPCRBatch))
              x1 <- as.data.frame(x[x$"Combined sample and target type" == "Target Unknown",])
              x2 <- as.data.frame(x[x$"Combined sample and target type"== "Target PosCalibrator",])
              x.target <- as.data.frame(rbind(x1,x2))
              Target <- x.target[,"Target name"]
              fac <- unique(Target)
              x.nrq<- c()
              for (i in 1:length(fac)){
                x.nrq <- append(NRQeffs(RTqPCRBatch, Target = fac[i], Ref = y),i)
              }
              x2 <- as.data.frame(x.nrq$ID)
              x3 <- as.data.frame(x.nrq$"Roche Method")
              x4 <- as.data.frame(x.nrq$"Pfaffl Method")
              x5 <- as.data.frame(x.nrq$"delta delta Cq Method")
              nrq <- cbind(x2,x3,x4,x5)
              row.names(nrq) <- 1:length(row.names(nrq))
              names(nrq) <- c("ID","Roche Method", "Pfaffl Method", "delta delta Cq Method")
              nrq[,"ID"] <- factor(nrq[,"ID"],levels = nrq[,"ID"])
              metaData <- data.frame(labelDescription = names(nrq),
                                     row.names = names(nrq), check.names = FALSE,
                                     stringsAsFactors = FALSE)
              fData <- AnnotatedDataFrame(data = nrq, varMetadata = metaData)
              nrq.table <- new("RTqPCRBatch", featureData = fData)
              return(nrq.table)
            }
          (RTqPCRBatch,...)
)

