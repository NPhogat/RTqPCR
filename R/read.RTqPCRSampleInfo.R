#'@name read.RTqPCRSampleInfo
#'@aliases read.RTqPCRSampleInfo
#'@title Read sample information file of raw data of RT-qPCRs (LC480 light cycler and Mx3005P)
#'@description Read the .txt sample information file of raw data of experiment run on RT-qPCRs 
#'(LC480 light cycler and Mx3005P) and use the data to populate an object of \code{"AnnodatedDataFrame"}.
#'@param file Name of the file to read in.
#'@param PCRtype The type of RT-qPCRs ("LC480" or "Mx3005P") whose file is to read in.
#'@param \dots Other parameters to be passed to downstream methods.
#'@return An object of \code{"AnnodatedDataFrame"}.
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@examples
#'## To read in the sample information data of LC480 light cycler
#'
#'SampleInfoLC480 <- file.path(path, "LC480_example_SampleInfo.txt")
#'samInfoLC480 <- read.RTqPCRSampleInfo(SampleInfoLC480, PCRtype = "LC480")
#'samInfoLC480[1:5]       ## to visualise the first five sample information data
#'
#'## To read in the sample information data of Mx3005P RT-qPCR
#'
#'SampleInfoMx <- file.path(path, "Mx3005P_example_SampleInfo.txt")
#'samInfoMx <- read.RTqPCRSampleInfo(SampleInfoMx, PCRtype = "Mx3005P")
#'samInfoMx[1:5]     ##to visualise the first five sample information data
#'@export
read.RTqPCRSampleInfo <- function (file, PCRtype, ...)
{
  if (PCRtype == "LC480")
  {
    samInfo <- Read.LC480SampleInfo(file = file)
  }
  else
  {
    if (PCRtype == "Mx3005P")
    {
      samInfo <- Read.Mx3005PSampleInfo(file = file)
    }
    else
    {
      stop("Define appropriate PCRtype out of LC480 or Mx3005P!")
    }
  }
  return(samInfo)
}