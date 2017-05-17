#'@name ReplaceValue
#'@aliases ReplaceValue
#'@title Replace the specific Cq value and amplification efficiency by user defined values or \code{NA}. 
#'@description Replace the specific Cq and amplification efficiency value by user defined values or \code{NA}
#'for Cq and amplification efficiencies in the object of class \code{"RTqPCRBatch"}.
#'@param RTqPCRBatch Object of class \code{RTqPCRBatch}. It can be output of \code{\link[RTqPCR]{CqValues}} or
#'\code{\link[RTqPCR]{ReplaceValue}} or \code{\link[RTqPCR]{ReplaceNAs}} or \code{\link[RTqPCR]{NonDetects}} or
#'\code{\link[RTqPCR]{ReplaceAboveCutOff}}.
#'@param NewCq New User defined value of Cq to replace a specific Cq value. It can be defined as \code{NA} too.
#'@param Neweffs New user defined value of amplification efficiency to replace a specific amplification 
#'efficiency. It can be defined as \code{NA} too.
#'@param Cqrow The specific row of Cq in exprs of object \code{RTqPCRBatch} where the new value of 
#'Cq is to be placed.
#'@param effsrow The specific row of amplification efficiency in effs of object \code{RTqPCRBatch}
#'where the new value of amplification efficiency is to be placed.
#'@param \dots Other parameters to be passed to downstream methods.
#'@return \code{"RTqPCRBatch"} object with a new exprs and effs slots.
#'@details The user can implement the \code{ReplaceValue} function if user wants to change any 
#'specific Cq value or amplification efficiency by another value or \code{NA}, without changing the other data.
#'For further details, see the vignettes in \pkg{RTqPCR} package.
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
#'## Replace the specific Cq value and amplification efficiency by user defined values
#'
#'Cqeffs.value <- ReplaceValue(Cqeffs.LC480, NewCq = 35, Neweffs = 1, Cqrow = 30, effsrow = 30)
#'Cqeffs.value ## To see the overview of data
#'exprs(Cqeffs.value)[1:5]  ## to visualise the first five Cq values 
#'effs(Cqeffs.value)[1:5]  ## to visualise the first five amplification efficiencies
#'exprs(Cqeffs.value)  ## to visualise all Cq values
#'effs(Cqeffs.value)   ## to visualise all amplification efficiencies
#'
#'## This function ReplaceValue can also be implemented on the output of other functions ReplaceNAs, 
#'## ReplaceAboveCutOff and NonDetects as well as on the output of itself (ReplaceValue).
#'## For further details, see the vignettes in RTqPCR package.
#'
#'@export
setMethod("ReplaceValue", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, NewCq = 38, Neweffs = 2, Cqrow = 30, effsrow = 30, ...)
            {
                if ((is.null(Cqrow)) && (is.null(effsrow)) && (is.null(NewCq)) && 
                      (is.null(Neweffs))) 
                  stop ("Please provide the appropriate values")
                
                else if (((!(is.null(NewCq))) && (is.null(Cqrow))) || 
                  ((is.null(NewCq)) && (!(is.null(Cqrow)))))
                  
                {
                  stop ("Please provide the appropriate values")
                }
                
                 else if (((!(is.null(Neweffs))) && (is.null(effsrow))) || 
                  ((is.null(Neweffs)) && (!(is.null(effsrow))))) 
                 {
                   stop ("Please provide the appropriate values.")
                 }
                
                 else
                   
                    if ((!(is.null(NewCq))) && (is.null(Neweffs)))
                     {
                       exprs(RTqPCRBatch)[Cqrow,] <- NewCq
                       effs(RTqPCRBatch) <- effs(RTqPCRBatch)
                     }
                
                    else if ((is.null(NewCq)) && (!(is.null(Neweffs))))
                    {
                      exprs(RTqPCRBatch) <- exprs(RTqPCRBatch)
                      effs(RTqPCRBatch)[effsrow,] <- Neweffs
                    }
                
                    else
                    {
                      exprs(RTqPCRBatch)[Cqrow,] <- NewCq
                      effs(RTqPCRBatch)[effsrow,] <- Neweffs
                    }
                        
                    return(RTqPCRBatch)
              }
                   (RTqPCRBatch, ...)      
            )