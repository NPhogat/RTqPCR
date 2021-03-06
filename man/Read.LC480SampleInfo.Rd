% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadLC480SampleInfo.R
\name{Read.LC480SampleInfo}
\alias{Read.LC480SampleInfo}
\title{Read sample information file of raw data of LC480 RT-qPCR}
\usage{
Read.LC480SampleInfo(file, removeEmptyCols = FALSE, header = TRUE,
  sep = "\\t", quote = "\\"", dec = ".", fill = TRUE,
  comment.char = "", skip = 0)
}
\arguments{
\item{file}{Name of the file to read in.}

\item{removeEmptyCols}{a logical value which indicates whether the empty column(s) should be removed 
or not. It should always be considered as \code{"FALSE"} to perform the downstream methods of the work
flow.}

\item{header}{a logical value indicating whether the file contains the names of the variables as its 
first line. If missing, the value is determined from the file format: header is set to TRUE if and 
only if the first row contains one fewer field than the number of columns.}

\item{sep}{the field separator character. Values on each line of the file are separated by this 
character. If \code{sep = ""} (the default for \code{\link[utils]{read.table}}) the separator is 
'white space', that is one or more spaces, tabs, newlines or carriage returns.}

\item{quote}{the set of quoting characters. To disable quoting altogether, use \code{quote = ""}. 
See \code{\link[base]{scan}} for the behaviour on quotes embedded in quotes. Quoting is only considered
for columns read as character, which is all of them unless colClasses is specified.}

\item{dec}{logical. If TRUE then in case the rows have unequal length, blank fields are implicitly 
added. See \code{\link[utils]{read.table}}.}

\item{fill}{character: a character vector of length one containing a single character or an empty 
string. Use \code{""} to turn off the interpretation of comments altogether.}

\item{comment.char}{character: a character vector of length one containing a single character or an 
empty string. Use \code{""} to turn off the interpretation of comments altogether.}

\item{skip}{integer: the number of lines of the data file to skip before beginning to read data.}
}
\value{
An object of \code{"AnnodatedDataFrame"}.
}
\description{
Read the .txt sample information file of raw data of experiment run on Light Cycler LC480
and use the data to populate an object of \code{"AnnodatedDataFrame"}.
}
\examples{
LC480.sampleInfo <- file.path(path, "LC480_example_SampleInfo.txt")
LC480.saminfo <- Read.LC480SampleInfo(Mx3005P.sampleInfo)
samInfo.Ag [1:5] ##To visualise first five sample information data
}
\author{
Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@stamats.de}
}
