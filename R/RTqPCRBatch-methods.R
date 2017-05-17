setMethod("exprs", signature = "RTqPCRBatch", definition = 
            function (object) assayDataElement(object, "exprs")
)

setReplaceMethod("exprs", signature = "RTqPCRBatch", definition = 
                   function (object, value) assayDataElementReplace(object, "exprs", value)
)

setMethod("se.exprs", signature = "RTqPCRBatch", definition = 
            function (object) assayDataElement(object, "se.exprs")
)

setReplaceMethod("se.exprs", signature = "RTqPCRBatch", definition = 
                   function (object, value) assayDataElementReplace(object, "se.exprs", value)
)

setMethod("exprs.well.order", signature = "RTqPCRBatch", definition = 
            function (object) assayDataElement(object, "exprs.well.order")
)

setReplaceMethod("exprs.well.order", signature = "RTqPCRBatch", definition = 
                   function (object, value) assayDataElementReplace(object, "exprs.well.order", value)
)

setMethod("effs", signature = "RTqPCRBatch", definition = 
            function (object) assayDataElement(object, "effs")
)


setReplaceMethod("effs", signature = "RTqPCRBatch", definition = 
                   function (object, value) assayDataElementReplace(object, "effs", value)
)

setMethod("se.effs", signature = "RTqPCRBatch", definition = 
            function (object) assayDataElement(object, "se.effs")
)

setReplaceMethod("se.effs", signature = "RTqPCRBatch", definition = 
                   function (object, value) assayDataElementReplace(object, "se.effs", value)
)
