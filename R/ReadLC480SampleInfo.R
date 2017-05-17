#'@name Read.LC480SampleInfo
#'@aliases Read.LC480SampleInfo
#'@title Read sample information file of raw data of LC480 RT-qPCR
#'@description Read the .txt sample information file of raw data of experiment run on Light Cycler LC480
#'and use the data to populate an object of \code{"AnnodatedDataFrame"}.
#'@param file Name of the file to read in.
#'@param removeEmptyCols a logical value which indicates whether the empty column(s) should be removed 
#'or not. It should always be considered as \code{"FALSE"} to perform the downstream methods of the work
#'flow.
#'@param header a logical value indicating whether the file contains the names of the variables as its 
#'first line. If missing, the value is determined from the file format: header is set to TRUE if and 
#'only if the first row contains one fewer field than the number of columns.
#'@param sep the field separator character. Values on each line of the file are separated by this 
#'character. If \code{sep = ""} (the default for \code{\link[utils]{read.table}}) the separator is 
#''white space', that is one or more spaces, tabs, newlines or carriage returns.
#'@param quote the set of quoting characters. To disable quoting altogether, use \code{quote = ""}. 
#'See \code{\link[base]{scan}} for the behaviour on quotes embedded in quotes. Quoting is only considered
#'for columns read as character, which is all of them unless colClasses is specified.
#'@param dec logical. If TRUE then in case the rows have unequal length, blank fields are implicitly 
#'added. See \code{\link[utils]{read.table}}.
#'@param fill character: a character vector of length one containing a single character or an empty 
#'string. Use \code{""} to turn off the interpretation of comments altogether.
#'@param comment.char character: a character vector of length one containing a single character or an 
#'empty string. Use \code{""} to turn off the interpretation of comments altogether.
#'@param skip integer: the number of lines of the data file to skip before beginning to read data.
#'@return An object of \code{"AnnodatedDataFrame"}.
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@examples
#'LC480.sampleInfo <- file.path(path, "LC480_example_SampleInfo.txt")
#'LC480.saminfo <- Read.LC480SampleInfo(Mx3005P.sampleInfo)
#'samInfo.Ag [1:5] ##To visualise first five sample information data
#'@export
Read.LC480SampleInfo <- function(file, removeEmptyCols = FALSE,
                                 header = TRUE, sep = "\t", quote = "\"",
                                 dec = ".",fill = TRUE, comment.char = "",
                                 skip = 0){
  sampleInfo <- read.table(file = file, header = header, sep = sep,
                           quote = quote, dec = dec, fill = fill,
                           comment.char = comment.char, skip = skip,
                           stringsAsFactors = FALSE, check.names = FALSE)
  if(!("General:Pos" %in% names(sampleInfo)))
    stop("No 'General:Pos' (i.e. sample position) available in sample information file!")
  ## remove empty columns
  if(removeEmptyCols){
    sampleInfo <- sampleInfo[,colSums(is.na(sampleInfo)) != nrow(sampleInfo)]
  }
  ## possible column names in exported sample information file
  origNames <- c("General:Pos", "General:Sample Name", "General:Repl. Of",
                 "General:Filt. Comb.", "General:Target Name",
                 "General:Color", "General:Subsets", "General:Notes",
                 "General:Sample ID", "General:Sample Prep Notes",
                 "Sample Preferences:Width", "Sample Preferences:Line Style",
                 "Sample Preferences:Point Style", "Sample Preferences:Color",
                 "Color Comp:Dominant Channel",
                 "Endpt. Geno:EndPt Sample Type",
                 "Endpt. Geno:EndPt Genotype",
                 "Abs Quant:Sample Type", "Abs Quant:Concentration",
                 "Abs Quant: Cp Low", "Abs Quant:Cp High",
                 "Melt Geno:Sample Type", "Melt Geno:Genotype",
                 "Rel Quant:Target Type",
                 "Rel Quant:Combined Sample/Target type",
                 "Rel Quant:Efficiency",
                 "Gene Scanning:Scanning Sample Type",
                 "Gene Scanning:Scanning Genotype")
  ## new column names
  newNames <- c("Sample position", "Sample name", "Replicate of",
                "Filter combination", "Target name", "Color",
                "Subsets", "Notes", "Sample ID", "Sample Prep Notes",
                "Sample Pref width", "Sample Pref line style",
                "Sample Pref point style", "Sample Pref color",
                "Dominant channel",
                "Endpoint sample type", "Endpoint genotype",
                "Quantification sample type", "Concentration",
                "Cp Low", "Cp High",
                "Melt Geno sample type", "Melt Geno genotype",
                "Target type", "Combined sample and target type",
                "Efficiency",
                "Scanning sample type", "Scanning genotype")
  ## change column names
  for(i in 1:ncol(sampleInfo)){
    if(names(sampleInfo)[i] %in% origNames)
      names(sampleInfo)[i] <- newNames[which(origNames == names(sampleInfo)[i])]
  }
  sampleInfo[,"Sample position"] <- factor(sampleInfo[,"Sample position"],
                                           levels = sampleInfo[,"Sample position"])
  metaData <- data.frame(labelDescription = names(sampleInfo),
                         row.names = names(sampleInfo), check.names = FALSE,
                         stringsAsFactors = FALSE)
  AnnotatedDataFrame(data = sampleInfo, varMetadata = metaData)
}

