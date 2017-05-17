#'@name Read.Mx3005PSampleInfo
#'@aliases Read.Mx3005PSampleInfo
#'@title Read sample information file of raw data of Mx3005P RT-qPCR
#'@description Read the .txt sample information file of raw data of experiment run on Mx3005P RT-qPCR
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
#'Mx3005P.sampleInfo <- file.path(path, "Mx3005P_example_SampleInfo.txt")
#'Mx3005P.samInfo <- Read.Mx3005PSampleInfo(Mx3005P.SamInfo)
#'Mx3005P.samInfo [1:5]   ## To visualise first five sample information data
#'@keywords \code{"AnnotatedDataFrame"}
#'@export
Read.Mx3005PSampleInfo <- function(file, removeEmptyCols = FALSE,
                                 header = TRUE, sep = "\t", quote = "\"",
                                 dec = ".",fill = TRUE, comment.char = "",
                                 skip = 0){
  sampleInfo <- read.table(file = file, header = header, sep = sep,
                           quote = quote, dec = dec, fill = fill,
                           comment.char = comment.char, skip = skip,
                           stringsAsFactors = FALSE, check.names = FALSE)
  if(!("General:Well" %in% names(sampleInfo)))
    stop("No 'General:Well' (i.e. Well) available in sample information file!")
  ## remove empty columns
  if(removeEmptyCols){
    sampleInfo <- sampleInfo[,colSums(is.na(sampleInfo)) != nrow(sampleInfo)]
  }
  ## possible column names in exported sample information file
  origNames <- c("General:Well", "General:Repl. Of", "General:Target Name",
                 "Rel Quant:Combined Sample/Target type")

  ## new column names
  newNames <- c("Well", "Replicate of", "Target name", "Combined sample and target type")
                
  ## change column names
  for(i in 1:ncol(sampleInfo)){
    if(names(sampleInfo)[i] %in% origNames)
      names(sampleInfo)[i] <- newNames[which(origNames == names(sampleInfo)[i])]
  }
  sampleInfo[,"Well"] <- factor(sampleInfo[,"Well"],
                                           levels = sampleInfo[,"Well"])
  metaData <- data.frame(labelDescription = names(sampleInfo),
                         row.names = names(sampleInfo), check.names = FALSE,
                         stringsAsFactors = FALSE)
  AnnotatedDataFrame(data = sampleInfo, varMetadata = metaData)
}
