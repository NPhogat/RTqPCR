#'@name RTqPCR.gui
#'@aliases RTqPCR.gui
#'@title Run the graphical User Interface
#'@description Run the graphical user interface 
#'@author Navneet Phogat, Matthias Kohl, \email{Matthias.Kohl@@stamats.de}
#'@examples
#'RTqPCR.gui()
#'@export
RTqPCR.gui <- function() {
  runApp(system.file("RTqPCR.gui", package = "RTqPCR"))
}


