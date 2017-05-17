#'@name CombineTechReps
#'@aliases CombineTechReps
#'@title Combine technical replicates of Cq and amplification efficiencies
#'@description Combine the technical replicates of Cq values and amplification efficiencies on the basis 
#'of mean, median and geometric mean and compute standard deviation of Cq and efficiency within the replicates
#'in object of class \code{"RTqPCRBatch"}.
#'@param RTqPCRBatch Object of class \code{RTqPCRBatch}. It can be output of \code{\link[RTqPCR]{CqValues}} or
#'\code{\link[RTqPCR]{ReplaceValue}} or \code{\link[RTqPCR]{NonDetects}} or \code{\link[RTqPCR]{ReplaceAboveCutOff}}
#'or \code{\link[RTqPCR]{ReplaceNAs}} functions.
#'@param calc Which methods to be used out of Mean, Median and Geometric mean.
#'@param cRepCq if TRUE, then combine replicates in Cq. If FALSE, then combine replicates in amplification 
#'efficiencies.
#'@param \dots Other parameters to be passed to downstream methods.
#'@return Object of class \code{"RTqPCRBatch"} with new slots. 
#'@details See the vignettes in \pkg{RTqPCR} package.
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
#'##  Implement the functions ReplaceNAs, ReplaceAboveCutOff, ReplaceValue and NonDetects as per requirement.
#'
#'## Combine technical replicates - In this example, we are combining the replicates directly from 
#'## the output of CqValues function
#'
#'##  replicates can be combined through following three methods
#'
#'## 1) based on mean (default method is for Cq values and is based on mean)
#'
#'## i) combine Cq values and compute the standard deviation of Cq values within replicates
#'
#'Cqreps <- CombineTechReps(Cqeffs.LC480)
#'Cqreps   ## to visualise the overview of the rsulting data
#'
#'
#'ii) combine amplification efficiencies and compute the standard deviation of amplification efficiencies
#'within replicates
#'
#'effsreps <- CombineTechReps(Cqeffs.LC480, cRepCq = FALSE)
#'
#'## 2) Based on median
#'
#'## i) Combine Cq values and compute standard deviation of Cq values within replicates
#'
#'Cqreps2 <- CombineTechReps(Cqeffs.LC480, calc = "Median")
#'
#'ii) Combine amplification efficiencies and compute their standard deviation within replicates
#'
#'effsreps2 <- CombineTechReps(Cqeffs.LC480, calc = "Median", cRepCq = FALSE)
#'
#'## 3) Based on Geometric Mean
#'
#'## i) combine Cq values and compute their standard deviation within replicates
#'
#'Cqreps3 <- CombineTechReps(Cqeffs.LC480, calc = "Geom")
#'
#'ii) combine amplification efficiencies and compute their standard deviation within replicates
#'
#'effsreps <- CombineTehcReps(Cqeffs.LC480, calc = "Geom", cRepCq = FALSE)
#'
#'@export
setMethod("CombineTechReps", signature = "RTqPCRBatch", definition = 
            function(RTqPCRBatch, calc = "Mean", cRepCq = TRUE, ...)
            {
                exprs.row <- which(!(is.na(exprs(RTqPCRBatch))))
                effs.row <- which(!(is.na(effs(RTqPCRBatch))))
                fData.exprs <- fData(RTqPCRBatch)[exprs.row, ]
                fData.effs <- fData(RTqPCRBatch)[effs.row, ]
                x.exprs <- fData.exprs[,"Replicate of"]
                x.effs <- fData.effs[,"Replicate of"]
                  
                if ((any(is.na(x.exprs))) || (any(is.na(x.effs))) ||
                      (any(is.element("",x.exprs))) || (any(is.element("",x.effs))))
                {
                  stop ("Please fill up the complete column Replicate of and then
                        try again")
                }
                
                else 
                  
                {  
                 expM <- exprs(RTqPCRBatch)[exprs.row, ] 
                 effsM <- effs(RTqPCRBatch)[effs.row, ]
                 exprs.effs <- exprs(RTqPCRBatch)[effs.row, ]
                 expM <- as.matrix(expM)
                 effsM <- as.matrix(effsM)
                 exprs.effs <- as.matrix(exprs.effs)
                 row.names(expM) <- x.exprs
                 row.names(effsM) <- x.effs
                 row.names(exprs.effs) <- x.effs
                 fac.expM <- factor(x.exprs, levels = unique(x.exprs))
                 fac.effsM <- factor(x.effs, levels = unique(x.effs))
                 sd.expM <- lapply(split(expM,fac.expM), sd, na.rm= TRUE)
                 sd.effsM <- lapply(split(effsM,fac.effsM), sd, na.rm= TRUE)
                 sd.exprs.effs <- lapply(split(exprs.effs,fac.effsM), sd, na.rm= TRUE)
                if (calc == "Mean")
                {
                 expM.split <- lapply(split(expM,fac.expM), mean, na.rm= TRUE)
                 effsM.split <- lapply(split(effsM,fac.effsM), mean, na.rm= TRUE)
                 exprs.effs.split <- lapply(split(exprs.effs,fac.effsM), mean, na.rm= TRUE)
                }
                else if (calc == "Median")
                {
                 expM.split <- lapply(split(expM,fac.expM), median, na.rm= TRUE)
                 effsM.split <- lapply(split(effsM,fac.effsM), median, na.rm= TRUE)
                 exprs.effs.split <- lapply(split(exprs.effs,fac.effsM), median, na.rm= TRUE)
                }
                else if (calc == "Geom")
                {
                 expM.split <- lapply(split(expM,fac.expM), geomMean, na.rm = TRUE)
                 effsM.split <- lapply(split(effsM,fac.effsM), geomMean, na.rm= TRUE)
                 exprs.effs.split <- lapply(split(exprs.effs,fac.effsM), geomMean, na.rm= TRUE)
                }
                else 
                {
                  stop ("ensure you have specified the correct centrality measure,
                        'Mean', Median' or 'Geom'")
                }
                 
                x1 <- as.matrix(x.exprs)
                row.names(x1) <- row.names(fData.exprs)
                x1 <- as.matrix(x1[!duplicated(x1[,1]),])
                ID.Cq <- as.matrix(row.names(x1))
                x2 <- as.matrix(x.effs) 
                row.names(x2) <- row.names(fData.effs)
                x2 <- as.matrix(x2[!duplicated(x2[,1]),])
                ID.effs <- as.matrix(row.names(x2))
                NewCq <- matrix(unlist(expM.split), ncol = 1, byrow = TRUE)
                Neweffs <- matrix(unlist(effsM.split), ncol = 1, byrow = TRUE)
                NewCq.effs <- matrix(unlist(exprs.effs.split), ncol = 1, byrow = TRUE) 
                sd.Cq <- matrix(unlist(sd.expM), ncol = 1, byrow = TRUE)
                sd.effs <- matrix(unlist(sd.effsM), ncol = 1,byrow = TRUE)
                sd.cq.effs <- matrix(unlist(sd.exprs.effs), ncol = 1,byrow = TRUE)
                crepCq <- as.data.frame(cbind(ID.Cq, format(NewCq,digit = 3), 
                          format(sd.Cq, digit = 3),fData.exprs[ID.Cq,"Target name"],
                          fData.exprs[ID.Cq,"Combined sample and target type"]))
                names(crepCq) <- c("ID","Cq","sd.Cq","Target name","Combined sample and target type")
                crepCq[,"ID"] <- factor(crepCq[,"ID"],levels = crepCq[,"ID"])
                metaData <- data.frame(labelDescription = names(crepCq),
                                       row.names = names(crepCq), check.names = FALSE,
                                       stringsAsFactors = FALSE)
                fData <- AnnotatedDataFrame(data = crepCq, varMetadata = metaData)
                rep.Cq <- new("RTqPCRBatch", featureData = fData)
                crepeffs <- as.data.frame(cbind(ID.effs, format(Neweffs,digit =3),
                            format(sd.effs, digit = 3), format(NewCq.effs, digit = 3),
                            format(sd.cq.effs, digit = 3), fData.effs[ID.effs,"Target name"],
                            fData.effs[ID.effs,"Combined sample and target type"]))
                names(crepeffs) <- c("ID","effs","sd.effs", "Cq","sd.Cq","Target name",
                                     "Combined sample and target type")
                crepeffs[,"ID"] <- factor(crepeffs[,"ID"],levels = crepeffs[,"ID"])
                metaData <- data.frame(labelDescription = names(crepeffs),
                                       row.names = names(crepeffs), check.names = FALSE,
                                       stringsAsFactors = FALSE)
                fData <- AnnotatedDataFrame(data = crepeffs, varMetadata = metaData)
                rep.effs <- new("RTqPCRBatch", featureData = fData)
                }
                
                if (cRepCq)
                  return (rep.Cq)
                else
                  return (rep.effs)
                }
              
                 (RTqPCRBatch, ...)
              )




                                    