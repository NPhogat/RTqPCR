#'@name DeltaDeltaCq
#'@aliases DeltaDeltaCq
#'@title Compute delta delta Cq, standard deviation of delta delta Cq and fold concentration.
#'@description Compute delta delta Cq, standard deviation of delta delta Cq and fold concentration of an 
#'object of \code{"RTqPCRBatch"}.
#'@param RTqPCRBatch An object of class \code{"RTqPCRBatch"} produced as an output of DeltaCq function.
#'@param Target Name of target gene under column Target name.
#'@param \dots other parameters to be passed to the downstream methods.
#'@return An object of \code{"data.frame"} with new slots.
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
#'## To compute the delta delta Cq value implement the DeltaDeltaCqAll function. It's an auxilliary function
#'##to DeltaDeltaCqAll function
#'@export 
setMethod("DeltaDeltaCq", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, Target, ...)
            {
              x <- as.data.frame(fData(RTqPCRBatch))
              if (missing(Target)) 
              {
                Target <- as.data.frame(x)
                Cal <- as.matrix(Target[Target$"Combined sample and target type" == "Target PosCalibrator", ])
              }
              
              else 
                
              {
                Target <- as.data.frame(x[x$"Target name" == Target, ])
                Cal <- as.matrix(Target[Target$"Combined sample and target type" == "Target PosCalibrator", ])                   
              }
              
              Cal.Cq <- mean((as.numeric(as.vector(Cal[,"deltaCq"]))), na.rm = TRUE)
              Cal.sd.Cq <- mean((as.numeric(as.vector(Cal[,"sd.deltaCq"]))), na.rm = TRUE)
              target.Cq <- as.numeric(as.matrix(Target[Target$"Combined sample and target type" 
                                                       == "Target Unknown","deltaCq"]))
              target.sd.Cq <- as.numeric(as.matrix(Target[Target$"Combined sample and target type"
                                                      == "Target Unknown","sd.deltaCq"]))
              
              deltadeltaCq <- (target.Cq) - (Cal.Cq)
              is.na(deltadeltaCq) <- NA
              deltadeltaCq <- as.matrix(format(deltadeltaCq, digit = 3))
              
              deltadeltaCq.2 <- as.numeric(as.matrix(deltadeltaCq))
              fold.conc <- log(2^(-deltadeltaCq.2))
              fold.conc <- as.matrix(format(fold.conc, digit = 3))
              
              sd.deltadeltaCq <- sqrt((target.sd.Cq)^2 + (Cal.sd.Cq)^2)
              is.na(sd.deltadeltaCq) <- NA
              sd.deltadeltaCq <- as.matrix(format(sd.deltadeltaCq, digit = 3))
              
              ID <- as.matrix(Target[Target$"Combined sample and target type" == "Target Unknown","ID"])
              
              deltadeltaCq.table <- as.data.frame(cbind(ID, deltadeltaCq, fold.conc, sd.deltadeltaCq))
                                                   
              names(deltadeltaCq.table) <- c("ID", "deltadeltaCq", "Fold Conc.[log]", "sd.deltadeltaCq")
              
              return (deltadeltaCq.table)
            }
          
          (RTqPCRBatch, ...)
)
