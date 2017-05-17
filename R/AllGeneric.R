setGeneric("effs", function(object) standardGeneric("effs"))

setGeneric("effs<-", function(object, value) standardGeneric("effs<-"))

setGeneric("se.effs", function(object) standardGeneric("se.effs"))

setGeneric("se.effs<-", function(object, value) standardGeneric("se.effs<-"))

setGeneric("exprs.well.order", function(object) standardGeneric("exprs.well.order"))

setGeneric("exprs.well.order<-", function(object, ..., value) standardGeneric("exprs.well.order<-"))

setGeneric("CqEffs", function(object, PCRtype, ...) standardGeneric("CqEffs"))

setGeneric("ReplaceAboveCutOff", function (RTqPCRBatch,...) standardGeneric("ReplaceAboveCutOff"))

setGeneric("ReplaceNAs", function (RTqPCRBatch,...) standardGeneric("ReplaceNAs"))

setGeneric("ReplaceValue", function (RTqPCRBatch,...) standardGeneric("ReplaceValue"))

setGeneric("RTqPCR.dataframe", function (x, fun, pcrtype, ...) standardGeneric("RTqPCR.dataframe"))

setGeneric("NonDetects", function (RTqPCRBatch, ...) standardGeneric("NonDetects"))

setGeneric("CombineTechReps", function (RTqPCRBatch,...) standardGeneric("CombineTechReps"))

setGeneric("DeltaCq", function (RTqPCRBatch,...) standardGeneric("DeltaCq"))

setGeneric("DeltaDeltaCq", function (RTqPCRBatch,...) standardGeneric("DeltaDeltaCq"))

setGeneric("DeltaDeltaCqAll", function (RTqPCRBatch,...) standardGeneric("DeltaDeltaCqAll"))

setGeneric("NRQeffs", function (RTqPCRBatch,...) standardGeneric("NRQeffs"))

setGeneric("NRQeffsAll", function (RTqPCRBatch, y, ...) standardGeneric("NRQeffsAll"))



