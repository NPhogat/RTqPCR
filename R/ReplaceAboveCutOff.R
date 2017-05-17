#'@name ReplaceAboveCutOff
#'@aliases ReplaceAboveCutOff
#'@title Replace the Cq and amplification efficiency above cut off by non detects (\code{NA}).
#'@description Replace the Cq and amplification efficiency above cut off by non detects (\code{NA}) 
#'in object of class \code{"RTqPCRBatch"}.
#'@param Object of class \code{RTqPCRBatch}. It can be output of \code{\link[RTqPCR]{CqValues}} or
#'\code{\link[RTqPCR]{ReplaceValue}} or \code{\link[RTqPCR]{ReplaceNAs}} or \code{\link[RTqPCR]{NonDetects}}.
#'@param NewVal The new value (\code{NA}) to replace the values above cut off.
#'@param Cqcutoff The minimum threshold Cq value above which the values will be replaced by NA.
#'@param effscutoff The minimum efficiency above which the values will be replaced by NA.
#'@param \dots Other parameters to be passed to downstream methods.
#'@return \code{"RTqPCRBatch"} object with new exprs and effs slots.
#'@details This function \code{ReplaceAboveCutOff} is specifically designed to replace Cq values
#'and amplification efficiencies, which are above the cutoffs by \code{NA}. There will be no effect on
#'remaining data. The best point to implement the \code{ReplaceAboveCutOff} function is when some values
#'of Cq and amplification efficiencies are coming too large for their consideration, due to some manual,
#'computational or instrumental error. So, the user can define the cutoffs for Cq values and amplification
#'efficiencies appropriately, based on their experience. For further details, see also vignettes in 
#'\pkg{RTqPCR} package.  
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@examples 
#'
#'## Read in the raw fluorescent data
#'
#'LC480.example <- file.path(path, "LC480_Example.txt") 
#'cycData.LC480 <- read.RTqPCR(LC480.example, PCRtype = "LC480")
#'
#'## Read in the sample information data
#'
#'SampleInfoLC480 <- file.path(path, "LC480_example_SampleInfo.txt")
#'samInfoLC480 <- read.RTqPCRSampleInfo(SampleInfoLC480, PCRtype = "LC480")
#'
#'## Merge the fluorescence and sample information data through merge function
#' 
#'merge.LC480<-merge(cycData.LC480,samInfo.LC480) 
#'
#'## Compute the Cq values and amplification efficiencies
#'
#'Cqeffs.LC480 <- CqValues(merge.LC480, PCRtype = "LC480", Effmethod = "sigfit", baseline = "none")
#'Cqeffs.LC480 #To see the overview of data
#'exprs(Cqeffs.LC480)[1:5]  ##to visualise the first five CqValues
#'effs(Cqeffs.LC480)[1:5]  ##to visualise the first five amplification efficiencies
#'exprs(Cqeffs.LC480)  ##to visualise all Cq values
#'effs(Cqeffs.LC480)   ##to visualise all amplification efficiencies
#' 
#'## Replace the Cq values and amplification efficiencies above the user defined cutoffs by NA
#'
#'Cqeffs.cutoff <- ReplaceAboveCutOff(Cqeffs.LC480, NewVal = NA, Cqcutoff = 38, effscutoff = 1)
#'Cqeffs.cutoff ## To see the overview of data
#'exprs(Cqeffs.cutoff)[1:5]  ## to visualise the first five Cq values 
#'effs(Cqeffs.cutoff)[1:5]  ## to visualise the first five amplification efficiencies
#'exprs(Cqeffs.cutoff)  ## to visualise all Cq values
#'effs(Cqeffs.cutoff)   ## to visualise all amplification efficiencies
#'
#'## This function ReplaceAboveCutOff can also be implemented on the outputs of other functions ReplaceValue, 
#'## ReplaceNAs and NonDetects as well as on the output of itself (ReplaceAboveCutOff), if the user
#'## wants to define a new value of cutoffs smaller than previous one. For further details, see the vignettes
#'## of RTqPCR package.
#'
#'@export
setMethod("ReplaceAboveCutOff", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, NewVal = NA, Cqcutoff = 38, effscutoff = 2, ...)
            {
                if (!(is.na(NewVal)))
                  stop ("NewVal should be NA")
                
                else if ((is.na(Cqcutoff)) && (is.na(effscutoff)))
                  stop ("Both Cqcutoff and effscutoff can not be NA")
                
                else
                  
                  if ((is.null(Cqcutoff)) && (is.null(effscutoff)))
                  {
                    stop ("Please provide the values.")
                  }
                
                    else if ((is.null(Cqcutoff)) && (!(is.null(effscutoff))))
                      {
                        exprs(RTqPCRBatch) <- exprs(RTqPCRBatch)
                        effs(RTqPCRBatch)[effs(RTqPCRBatch) > effscutoff] <- NewVal
                      }               
                        else if ((!(is.null(Cqcutoff))) && (is.null(effscutoff)))
                        {
                          exprs(RTqPCRBatch)[exprs(RTqPCRBatch) > Cqcutoff] <- NewVal
                          effs(RTqPCRBatch) <- effs(RTqPCRBatch)
                        }
                else 
                {
                  exprs(RTqPCRBatch)[exprs(RTqPCRBatch) > Cqcutoff] <- NewVal 
                  effs(RTqPCRBatch)[effs(RTqPCRBatch) > effscutoff] <- NewVal
                }
                
                return(RTqPCRBatch)
              }
                (RTqPCRBatch, ...)      
            )














