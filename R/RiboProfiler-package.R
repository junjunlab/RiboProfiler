#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")
  start_up_mg <- cli::cat_boxx("Welcome to use RiboProfiler package for Ribo-seq analysis.",
                               col = "#8B1874")
  packageStartupMessage(start_up_mg)
  packageStartupMessage(paste("The version of RiboProfiler:",
                              pkgVersion,
                              "\nAny advice or suggestions please contact with me: 3219030654@stu.cpu.edu.cn.",
                              sep = " "))
}
