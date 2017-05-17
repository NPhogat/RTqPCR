#'@name NRQeffs
#'@aliases NRQeffs
#'@title Compute relative expression ratio
#'@description Perform normalization and compute relative expression ratio of target sample on the 
#'basis of efficiency based methods of Roche, Pfaffl and delta delta Cq from an object of class \code{"RTqPCRBatch"}.
#'@param x An object of class \code{"RTqPCRBatch"} produced as an result of function CombineTechReps.
#'@param Target Name of target sample gene under column Target name.
#'@param Ref Name of reference gene under column Target name.
#'@param \dots Other parameters to be passed to the downstream methods.
#'@return An object of class \code{"data.frame"} with new slots.
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
#'##NRQeffs is an auxiliarry function to NRQeffsAll function. It can be implemented directly through the 
#'##NRQeffsAll function
#'
#'@export 
setMethod("NRQeffs", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, Target, Ref, ...)
            {
              x <- as.data.frame(fData(RTqPCRBatch))
                x.effs.Cq <- as.data.frame(x)  
                
                if (missing(Ref))
                {
                  warning("Please provide Reference gene. Otherwise, all reference genes 
                          will be considered")
                  ref.cal <- as.matrix(x.effs.Cq[x.effs.Cq$
                                      "Combined sample and target type" == "Ref PosCalibrator", ])
                  ref.sample <- as.matrix(x.effs.Cq[x.effs.Cq$
                                      "Combined sample and target type" == "Ref Unknown", ])
                }
                
                else
                {
                x.ref <- as.data.frame(x.effs.Cq[x.effs.Cq$"Target name" == Ref, ])
                ref.cal <- as.matrix(x.ref[x.ref$
                           "Combined sample and target type" == "Ref PosCalibrator", ])
                ref.sample <- as.matrix(x.ref[x.ref$
                              "Combined sample and target type" == "Ref Unknown", ])
                }
                
                if (missing(Target))
                {
                  target.cal <- as.matrix(x.effs.Cq[x.effs.Cq$
                                "Combined sample and target type" == "Target PosCalibrator", ])
                  target.sample <- as.matrix(x.effs.Cq[x.effs.Cq$
                                   "Combined sample and target type" == "Target Unknown", ])
                }
                
                else
                {
                x.target <- as.data.frame(x.effs.Cq[x.effs.Cq$"Target name"== Target, ])
                target.cal <- as.matrix(x.target[x.target$
                              "Combined sample and target type" == "Target PosCalibrator", ])
                target.sample <- as.matrix(x.target[x.target$
                                 "Combined sample and target type" == "Target Unknown", ])
                }
                
                ## Means of Calibrator & sample of Reference 
                ref.cal.effs <- mean(as.numeric(as.vector(ref.cal[,"effs"])), na.rm = TRUE)
                ref.sample.effs <- mean((as.numeric(as.vector(ref.sample[,"effs"]))), na.rm = TRUE)
                ref.cal.Cq <- mean((as.numeric(as.vector(ref.cal[,"Cq"]))), na.rm = TRUE)
                ref.sample.Cq <- mean((as.numeric(as.vector(ref.cal[,"Cq"]))), na.rm = TRUE)
                
                ## Target sample and Mean of calibrator of Target
                target.cal.effs <- mean((as.numeric(as.vector(target.cal[,"effs"]))), 
                                   na.rm = TRUE)
                target.cal.Cq <- mean((as.numeric(as.vector(target.cal[,"Cq"]))), na.rm = TRUE)
                target.sample.effs <- as.numeric(as.vector(target.sample[,"effs"]))
                target.sample.Cq <- as.numeric(as.vector(target.sample[,"Cq"]))

                ## Calculation of effs^Cq
                ref.cal.effs.Cq <- (ref.cal.effs)^(ref.cal.Cq)
                ref.sample.effs.Cq <- (ref.sample.effs)^(ref.sample.Cq)
                target.cal.effs.Cq <- (target.cal.effs)^(target.cal.Cq)
                target.sample.effs.Cq <- (target.sample.effs)^(target.sample.Cq)

                ## Ref Sample/Target Sample
                ref.target.sample <- (ref.sample.effs.Cq)/(target.sample.effs.Cq)

                ## Ref Cal/Target Cal
                ref.target.cal <- (ref.cal.effs.Cq)/(target.cal.effs.Cq)

                ## Final Ratio By Roche Method
                ratio <- log((ref.target.sample)/(ref.target.cal))
                is.na(ratio) <- NA
                ID <- target.sample[,"ID"]
                
                ## Delta Cq of Target and Ref
                deltaCq.target <- (target.cal.Cq) - (target.sample.Cq)
                deltaCq.ref <- (ref.cal.Cq) - (ref.sample.Cq)
                
                ## Final Ratio by M.W.Pfaffl Method
                ratio2 <- log(((target.sample.effs)^(deltaCq.target))/((ref.sample.effs)^(deltaCq.ref)))
                is.na(ratio2) <- NA
                
                ## Final ratio by delta delta Cq relative method
                deltadeltaCq <- (deltaCq.target) - (deltaCq.ref)
                ratio3 <- log(2^(-(deltadeltaCq)))
                is.na(ratio3) <- NA
                
                relative.expression.table <- as.data.frame(cbind(ID, format(ratio, digit = 3),
                                            format(ratio2, digit = 3), format(ratio3, digit = 3)))
                names(relative.expression.table) <- c("ID", "Roche Method",
                                                    "Pfaffl Method", "delta delta Cq Method") 
                
                return(relative.expression.table)
              }
               (x, ...)
            )
                                         
                
                
                
                