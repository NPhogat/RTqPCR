#'@name NonDetects
#'@aliases NonDetects
#'@title Remove the non-detects (\code{NA}) of Cq and amplification efficiency values
#'@description Remove the non detects (\code{NA}) of Cq values and amplification efficiencies within
#'the replicates in the object of class \code{"RTqPCRBatch"} on the basis of statistical concepts
#'of mean and median.
#'@param RTqPCRBatch Object of class \code{RTqPCRBatch}. It can be output of \code{\link[RTqPCR]{CqValues}} or
#'\code{\link[RTqPCR]{ReplaceValue}} or \code{\link[RTqPCR]{NonDetects}} or \code{\link[RTqPCR]{ReplaceAboveCutOff}}
#'functions.
#'@param Calc Which method to be used out of \code{"Mean"} or \code{"Median"}. Default is \code{"Mean"}.
#'@param \dots Other parameters to be passed to downstream methods.
#'@return \code{"RTqPCRBatch"} object with new exprs and effs slots.
#'@details See the vignettes in \pkg{RTqPCR} package. 
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@examples
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
#'## Replace the non-detects (NA) within replicates, based on the concept of mean
#'## default method is mean.
#'
#'Cqeffs.NA.mean <- NonDetects(Cqeffs.LC480)
#'Cqeffs.NA.mean ## To see the overview of data
#'exprs(Cqeffs.NA.mean)[1:5]  ## to visualise the first five Cq values 
#'effs(Cqeffs.NA.mean)[1:5]  ## to visualise the first five amplification efficiencies
#'exprs(Cqeffs.NA.mean)  ## to visualise all Cq values
#'effs(Cqeffs.NA.mean)   ## to visualise all amplification efficiencies
#'
#'## Replace the non-detects (NA) within replicates, based on the concept of median
#'
#'Cqeffs.NA.med <- NonDetects(Cqeffs.LC480, Calc = "Median")
#'Cqeffs.NA.med ## To see the overview of data
#'exprs(Cqeffs.NA.med)[1:5]  ## to visualise the first five Cq values 
#'effs(Cqeffs.NA.med)[1:5]  ## to visualise the first five amplification efficiencies
#'exprs(Cqeffs.NA.med)  ## to visualise all Cq values
#'effs(Cqeffs.NA.med)   ## to visualise all amplification efficiencies
#'
#'## The function NonDetects can also be implemented on the outputs of other functions ReplaceValue, 
#'## ReplaceAboveCutOff and NonDetects. For detail, see the vignettes in RTqPCR package.  

#'@export
setMethod("NonDetects", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, Calc = "Mean", ...)
            {
              expM <- exprs(RTqPCRBatch)
              effsM <- effs(RTqPCRBatch)
              x <- fData(RTqPCRBatch)[,"Replicate of"]
              if ((any(is.na(x))) || (any(is.element("",x))))
              {
                stop ("Please be ensure that there is no NA and vacant row in Replicate of column")
              }
              else
              {
                row.names(expM) <- x
                row.names(effsM) <- x
                x.fac <- factor(x, levels = unique(x))
                
                if (Calc == "Mean")
                {
                  expM.split <- unsplit(lapply(split(expM,x.fac), mean, na.rm= TRUE), x.fac)
                  effsM.split <- unsplit(lapply(split(effsM,x.fac), mean, na.rm= TRUE), x.fac)
                }
                
                else if (Calc == "Median")
                {
                  expM.split <- unsplit(lapply(split(expM,x.fac), median, na.rm= TRUE), x.fac)
                  effsM.split <- unsplit(lapply(split(effsM,x.fac), median, na.rm= TRUE), x.fac)
                }
                
                else 
                {
                  stop ("Be ensure that correct centrality measure Mean, Median
                        has been chosen")
                  }
                  
                  exprs.row <- which(is.na(expM))
                  effs.row <- which(is.na(effsM))
                  expM1 <- as.data.frame(expM.split)
                  effsM1 <- as.data.frame(effsM.split)
                  expM2 <- expM1[exprs.row, ]
                  effsM2 <- effsM1[effs.row, ]
                  expM[exprs.row, ] <- expM2
                  effsM[effs.row, ] <- effsM2
                  exprs1 <- expM
                  effs1 <- effsM
                  exprs1[is.na(exprs1)] <- NA
                  effs1[is.na(effs1)] <- NA
                  row.names(exprs1) <- featureNames(RTqPCRBatch)
                  row.names(effs1) <- featureNames(RTqPCRBatch)
                  RTqPCRBatch <- new("RTqPCRBatch", exprs = exprs1, featureData = featureData(RTqPCRBatch),
                                    effs = effs1)
                }
                  return (RTqPCRBatch)
                
                }
                
                 (RTqPCRBatch, ...)
              )