#'@name ReplaceNAs
#'@aliases ReplaceNAs
#'@title Replace the all the non detects of Cq and amplification efficiencies by user defined values
#'@description Replace all the non detects of Cq and amplification efficiencies by user defined values 
#'in object of \code{"RTqPCRBatch"}.
#'@param RTqPCRBatch Object of class \code{RTqPCRBatch}. It can be output of \code{\link[RTqPCR]{CqValues}} or
#'\code{\link[RTqPCR]{ReplaceValue}} or \code{\link[RTqPCR]{ReplaceAboveCutOff}} or \code{\link[RTqPCR]{NonDetects}}. 
#'@param NewCqNA New user defined value of Cq for replacing non - detects (\code{NA}) of Cq.
#'@param NeweffsNA New user defined value of amplification efficiency for replacing non detects (\code{NA})
#'of amplification efficiency.
#'@param \dots other parameters to be passed to downstream methods.
#'@return \code{"RTqPCRBatch"} object with new exprs and effs slots.
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@details \code{ReplaceNAs} works on the object of \code{RTqPCRBatch}, which is obtained through implementation
#'of \code{CqValues} or \code{ReplaceAboveCutOff} or \code{ReplaceValue} or \code{NonDetects}. The three 
#'functions \code{ReplaceNAs}, \code{NonDetects} and \code{ReplaceAboveCutOff} can be called among themselves 
#'as well as for implementation on results of \code{CqValues} function. The functions for replacing the values
#'are mainly designed to have a better control of user over the resultant data of Cq values and amplification 
#'efficiencies. For further details, see also \code{\link[RTqPCR]{vignettes}}.   
#'@examples 
#'
#'## i) Read in the raw fluorescent data
#'
#'LC480.example <- file.path(path, "LC480_Example.txt") 
#'cycData.LC480 <- read.RTqPCR(LC480.example, PCRtype = "LC480")
#'
#'## ii) Read in the sample information data
#'
#'SampleInfoLC480 <- file.path(path, "LC480_example_SampleInfo.txt")
#'samInfoLC480 <- read.RTqPCRSampleInfo(SampleInfoLC480, PCRtype = "LC480")
#'
#'## iii) Merge the fluorescence and sample information data through merge function
#' 
#'merge.LC480<-merge(cycData.LC480,samInfo.LC480) 
#'
#'## iv) Compute the Cq values and amplification efficiencies
#'
#'Cqeffs.LC480 <- CqValues(merge.LC480, PCRtype = "LC480", Effmethod = "sigfit", baseline = "none")
#'Cqeffs.LC480 #To see the overview of data
#'exprs(Cqeffs.LC480)[1:5]  ##to visualise the first five CqValues
#'effs(Cqeffs.LC480)[1:5]  ##to visualise the first five amplification efficiencies
#'exprs(Cqeffs.LC480)  ##to visualise all Cq values
#'effs(Cqeffs.LC480)   ##to visualise all amplification efficiencies
#' 
#'## v) Replace the non-detects (NA) by user defined values
#'
#'Cqeffs.ReNA <- ReplaceNAs(Cqeffs.LC480, NewCqNA = 35, NeweffsNA = 1)
#'Cqeffs.ReNA ## To see the overview of data
#'exprs(Cqeffs.ReNA)[1:5]  ## to visualise the first five Cq values 
#'effs(Cqeffs.ReNA)[1:5]  ## to visualise the first five amplification efficiencies
#'exprs(Cqeffs.ReNA)  ## to visualise all Cq values
#'effs(Cqeffs.ReNA)   ## to visualise all amplification efficiencies
#'
#'## This function ReplaceNAs can also be implemented on the outputs of other functions ReplaceValue, 
#'## ReplaceAboveCutOff and NonDetects and for that the user doesn't need to calculate the CqValues 
#'## again and again. For detail, check the vignettes of RTqPCR package.  
#'    
#'@export
setMethod("ReplaceNAs", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, NewCqNA, NeweffsNA, ...)
            {
                if ((is.null(NewCqNA)) && (is.null(NeweffsNA)))
                {
                  stop ("Please provide the appropriate values.")
                }
                
                else if ((is.na(NewCqNA)) && (is.na(NeweffsNA)))
                {
                  stop ("Please provide the appropriate values.")
                }
                
                else
                  
                if ((!(is.null(NewCqNA))) && (is.null(NeweffsNA)))
                {
                  exprs(RTqPCRBatch)[is.na(exprs(RTqPCRBatch))] <- NewCqNA
                  effs(RTqPCRBatch) <- effs(RTqPCRBatch)
                }
                
                else if ((is.null(NewCqNA)) && (!(is.null(NeweffsNA))))
                {
                  exprs(RTqPCRBatch) <- exprs(RTqPCRBatch)
                  effs(RTqPCRBatch)[is.na(effs(RTqPCRBatch))] <- NeweffsNA  
                }    
                 
                else
                {
                 exprs(RTqPCRBatch)[is.na(exprs(RTqPCRBatch))] <- NewCqNA 
                 effs(RTqPCRBatch)[is.na(effs(RTqPCRBatch))] <- NeweffsNA
                }
                 
                return(RTqPCRBatch)
              }

                (RTqPCRBatch, ...)      
            )