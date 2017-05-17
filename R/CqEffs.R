#'@name CqEffs
#'@aliases CqEffs
#'@title Compute Cq values and amplification efficiencies
#'@description Compute Cq values from an object of Class \code{"RTqPCRBatch"} and produce an object of Class
#'\code{"RTqPCRBatch"} as a result.
#'@param object An object of Class \code{"RTqPCRBatch"}.
#'@param PCRtype The type of RT-qPCRs ("LC480" or "Mx3005P") for which Cq values and amplification efficiencies
#'are to be computed.
#'@param Effmethod A character vector defining the methods for computing amplification efficiency.
#'@param group A vector containing the grouping for possible replicates.
#'@param model The model to be used for all runs. Default model is \code{l5}.
#'@param check The method for kinetic outlier detection in \code{\link[qpcR]{KOD}}. Method \code{"uni2"}
#' is set as default which is a test on sigmoidal structure.
#'@param checkPAR Parameters to be supplied to the \code{check} method. See \code{\link[qpcR]{parKOD}}.
#'@param remove Indicates which runs to be removed. Either \code{none} of them, those which failed 
#'to \code{fit} or from the outlier methods.
#'@param exclude Indicates samples to be excluded from calculation, either "" for samples with missing column names or a regular expression defining columns (samples); 
#'see Examples in \code{\link[qpcR]{modlist}}
#'@param type The point on the amplification curve which is used for efficiency estimation; 
#'see \code{\link[qpcR]{efficiency}}.
#'@param labels A vector containing labels which define replicate groups. See more details in 
#'\code{\link[qpcR]{pcrbatch}} and \code{\link[qpcR]{ratiobatch}}.
#'@param norm a logical Value which determines whether the raw data should be normalized within [0, 1] 
#'before model fitting or not.
#'@param baselinemethod Type of baseline subtraction. More details in \code{\link[qpcR]{efficiency}}.
#'@param basecyc Cycles to be considered for the baseline subtraction.
#'@param basefac A factor when using averaged baseline cycles, such as \code{0.95}.
#'@param smooth The curve smoothing method. See more details in \code{\link[qpcR]{pcrbatch}}.
#'@param smoothPAR parameters to be supplied to smoothing method in \code{smooth}.
#'@param factor A multiplication factor for the fluorescence response values.
#'@param opt A logical value which determines whether model selection should be applied to each model 
#'or not.
#'@param optPAR Parameters to be supplied for model selection in \code{\link[qpcR]{mselect}}.
#'@param plot A logical value. If \code{TRUE}, the single runs are plotted from the internal 
#'\code{modlist} for diagnostics.
#'@param verbose A logical value. If \code{TRUE}, fitting and tagging results will be displayed in the 
#'console. 
#'@param \dots Other parameters to be passed to the downstream methods.
#'@return Object of Class \code{"RTqPCRBatch"}.
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@examples
#'## 1) To compute the CqValues and amplification efficiencies for LC480 light cycler
#'
#'## Read in the raw fluorescent data of LC480 light cycler
#'
#'LC480.example <- file.path(path, "LC480_Example.txt") 
#'cycData.LC480 <- read.RTqPCR(LC480.example, PCRtype = "LC480")
#'
#'## Read in the sample information data of LC480 light cycler
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
#'## i) sigmoidal model (sigfit - default method)
#'## There are two ways two call this method. One is by defining the Effmethod as "sigfit".
#'## The second is calling CqValues function directly, without defining the Effmethod. As, it
#'## is a default method, so there is no need to define the Effmethod as "sigfit".
#'
#'Cqeffs.LC4801 <- CqValues(merge.LC480, PCRtype = "LC480", Effmethod = "sigfit", baseline = "none")
#'
#'## without mentioning Effmethod (as a default method)
#'Cqeffs.LC4801 <- CqValues(merge.LC480, PCRtype = "LC480", baseline = "none")
#'
#'Cqeffs.LC4801 #To see the overview of data
#'exprs(Cqeffs.LC4801)[1:5]  ##to visualise the first five CqValues
#'effs(Cqeffs.LC4801)[1:5]  ##to visualise the first five amplification efficiencies
#'exprs(Cqeffs.LC4801)  ##to visualise all Cq values
#'effs(Cqeffs.LC4801)   ##to visualise all amplification efficiencies
#'
#'## ii)  fit exponential model (expfit)
#'
#'Cqeffs.LC4802 <- CqValues(merge.LC480, PCRtype = "LC480", Effmethod = "expfit", baseline = "none")
#'Cqeffs.LC4802 
#'exprs(Cqeffs.LC4802)[1:5]  
#'effs(Cqeffs.LC4802)[1:5]  
#'exprs(Cqeffs.LC4802)  
#'effs(Cqeffs.LC4802)   
#'
#'## iii) window of linearity (sliwin)
#'
#'Cqeffs.LC4803 <- CqValues(merge.LC480, PCRtype = "LC480", Effmethod = "sliwin", baseline = "none")
#'Cqeffs.LC4803 
#'exprs(Cqeffs.LC4803)[1:5]  
#'effs(Cqeffs.LC4803)[1:5]  
#'exprs(Cqeffs.LC4803)  
#'effs(Cqeffs.LC4803)
#'
#'## iv) linear regression of efficiency (LRE)
#'
#'Cqeffs.LC4804 <- CqValues(merge.LC480, PCRtype = "LC480", Effmethod = "expfit", baseline = "none")
#'Cqeffs.LC4804 
#'exprs(Cqeffs.LC4804)[1:5]  
#'effs(Cqeffs.LC4804)[1:5]  
#'exprs(Cqeffs.LC4804)  
#'effs(Cqeffs.LC4804)

#'## 2) To compute the Cq values and amplification efficiencies for the Mx3005P RT-qPCR
#'
#'## Read in the raw fluorescent data of Mx3005P RT-qPCR
#'
#'Mx3005P.example <- file.path(path, "Mx3005P_Example.txt") 
#'cycData.Mx <- read.RTqPCR(Mx3005P.example, PCRtype = "Mx3005P")
#'
#'## Read in the sample information data of Mx3005P RT-qPCR
#'
#'SampleInfoMx <- file.path(path, "Mx3005P_example_SampleInfo.txt")
#'samInfoMx <- read.RTqPCRSampleInfo(SampleInfoMx, PCRtype = "Mx3005P")
#'
#'## Merge the fluorescence and sample information data
#'
#'merge.Mx<-merge(cycData.Mx,samInfoMx)
#'
#'## To compute the Cq values and amplification efficiencies
#'
#'## i) SIgmoidal model (sigfit - default method)
#'
#'Cqeffs.Mx<- CqValues(merge.Mx, PCRtype = "Mx3005P", Effmethod = "sigfit", baseline = "none")
#'Cqeffs.Mx #To see the overview of data
#'exprs(Cqeffs.Mx)[1:5]  ##to visualise the first five CqValues
#'effs(Cqeffs.Mx)[1:5]  ##to visualise the first five amplification efficiencies
#'exprs(Cqeffs.Mx)  ##to visualise all Cq values
#'effs(Cqeffs.Mx)   ##to visualise all amplification efficiencies
#'
#'## ii) fit exponential model (expfit)
#'
#'Cqeffs.Mx<- CqValues(merge.Mx, PCRtype = "Mx3005P", Effmethod = "expfit", baseline = "none")
#'Cqeffs.Mx #To see the overview of data
#'exprs(Cqeffs.Mx)[1:5]  
#'effs(Cqeffs.Mx)[1:5]  
#'exprs(Cqeffs.Mx)  
#'effs(Cqeffs.Mx)
#'
#'## iii) Window of linearity (sliwin)
#'
#'Cqeffs.Mx<- CqValues(merge.Mx, PCRtype = "Mx3005P", Effmethod = "sliwin", baseline = "none")
#'Cqeffs.Mx
#'exprs(Cqeffs.Mx)[1:5]
#'effs(Cqeffs.Mx)[1:5]  
#'exprs(Cqeffs.Mx)  
#'effs(Cqeffs.Mx)
#'
#'## iv) linear regression of efficiency
#'
#'Cqeffs.Mx<- CqValues(merge.Mx, PCRtype = "Mx3005P", Effmethod = "LRE", baseline = "none")
#'Cqeffs.Mx 
#'exprs(Cqeffs.Mx)[1:5]  
#'effs(Cqeffs.Mx)[1:5]  
#'exprs(Cqeffs.Mx)  
#'effs(Cqeffs.Mx)
#'
#'@export
setMethod("CqEffs", signature = "RTqPCRBatch", definition =
            function(object, PCRtype, Effmethod = "expfit", group = NULL, model = l5,
                     check = "uni2", checkPAR = parKOD(), remove = "none",
                     exclude = NULL, type = "cpD2", labels = NULL,
                     norm = FALSE, baselinemethod = "mean", basecyc = 1:8, basefac = 1,
                     smooth = NULL, smoothPAR = list(span = 0.1), factor = 1,
                     opt = FALSE, optPAR = list(sig.level = 0.05, crit = "ftest"),
                     plot = FALSE, verbose = FALSE, ...){
              if(missing(Effmethod))
                Effmethod <- "sigfit"
              else{
                if(length(Effmethod) > 1){
                  Effmethod <- Effmethod[1]
                  warning("Only first element of 'Effmethod' is used!")
                }
              }
              if(Effmethod != "sigfit"){
                methods <- c("sigfit", Effmethod)
              }
              else{
                methods <- Effmethod
              }
              if (missing(baselinemethod)){
                baselinemethod <- "none"
              }
              else{
                if(length(baselinemethod) >1){
                  baselinemethod <- baselinemethod[1]
                  warning("only first element of 'baseline' is used!")
                }
              }
                if(baselinemethod != "none"){
                  baseline <- c("none", baselinemethod)
                }
              else {
                baseline <- baselinemethod
              }
              
              fluoData <- exprs(object)
              if (PCRtype == "LC480")
              {
              x <- data.frame(Cycles = fData(object)[,"Cycle number"], fluoData)
              names(x) <- c("Cycles", pData(object)[,"Sample position"])
              }
              else 
              {
                if (PCRtype == "Mx3005P")
                {
                  x <- data.frame(Cycles = fData(object)[,"Cycle#"], fluoData)
                  names(x) <- c("Cycles", pData(object)[,"Well"])
                }
                else
                {
                  stop("Define PCRtype out of LC480 or Mx3005P!")
                }
              }
              if(missing(group))
                res <- pcrbatch(x = x, cyc = 1, fluo = NULL, methods = methods,
                                model = model, check = check, checkPAR = checkPAR,
                                remove = remove, exclude = exclude, type = type,
                                labels = labels,norm = norm, baseline = baseline, 
                                basecyc = basecyc, basefac = basefac,
                                smooth = smooth, smoothPAR = smoothPAR,
                                factor = factor, opt = opt, optPAR = optPAR,
                                group = group,
                                names = "first", plot = plot, verbose = verbose, ...)
              
              ## extract Cq values and efficiencies
              CqVals <- as.matrix(as.numeric(res[res[,"Vars"] == paste("sig", type, sep = "."),-1]))
              
              if(Effmethod == "sigfit"){
                Effs1 <- res[res[,"Vars"] == "sig.eff", -1]
                se.Effs1 <- as.matrix(sqrt(as.numeric(res[res[,"Vars"] == "sig.resVar", -1])))
              }
              if(Effmethod == "sliwin"){
                Effs1 <- res[res[,"Vars"] == "sli.eff", -1]
                se.Effs1 <- NULL
              }
              
              if(Effmethod == "expfit"){
                Effs1 <- res[res[,"Vars"] == "exp.eff", -1]
                se.Effs1 <- as.matrix(sqrt(as.numeric(res[res[,"Vars"] == "exp.resVar", -1])))
              }
              if(Effmethod == "LRE"){
                Effs1 <- res[res[,"Vars"] == "LRE.eff", -1]
                se.Effs1 <- NULL
              }
              Effs1 <- as.matrix(as.numeric(Effs1))
              
              metData <- data.frame(labelDescription = res[,"Vars"],
                                    row.names = res[,"Vars"], check.names = FALSE,
                                    stringsAsFactors = FALSE)
              pD <- t(res[,-1])
              colnames(pD) <- res[,"Vars"]
              metData1 <- varMetadata(phenoData(object))
              metData <- rbind(metData1, metData)
              pData(object) <- cbind(pData(object), pD)
              varMetadata(phenoData(object)) <- metData
              if ((Effmethod == "sigfit") || (Effmethod == "expfit"))
              {
                effs <- assayDataNew("environment", effs = Effs1)
                se.effs <- assayDataNew("environment", se.effs = se.Effs1)
                
                res <- new("RTqPCRBatch", exprs = CqVals, featureData = phenoData(object),
                           effs = Effs1, se.effs = se.Effs1)
                
              }
              
              if ((Effmethod == "sliwin") || (Effmethod == "LRE"))
              {
                effs <- assayDataNew("environment", effs = Effs1)
                
                res <- new("RTqPCRBatch", exprs = CqVals, featureData = phenoData(object),
                           effs = Effs1)
                
              }
              res
            }
          
)
