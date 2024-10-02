#' Convert BAM file to BigWig format
#'
#' This function converts a BAM file to BigWig format, with options for
#' sequencing type, read assignment, normalization, and read length filtering.
#'
#' @param bam_file Character string. Path to the input BAM file.
#' @param output_file Character string. Path for the output BigWig file.
#' @param julia_path The julia program path on the computer.
#' @param XAM_version The XAM version you have installed, 3(0.3.1) or 4(0.4.0).
#' @param seq_type Character string. Sequencing type, either "pairedEnd" or "singleEnd". Default is "pairedEnd".
#' @param assignType Character string. Read assignment type, either "end5" or "end3". Default is "end5".
#' @param normalization Character string. Normalization method, either "rpm" or "counts". Default is "rpm".
#' @param min_length Integer. Minimum read length to include. Default is 28.
#' @param max_length Integer. Maximum read length to include. Default is 31.
#' @param offset Offset to apply to read positions eg "28,29,30,31|-15,-15,-15,-15".
#' Default is "".
#'
#' @return None. The function writes the output to the specified BigWig file.
#'
#' @details
#' This function uses Julia and the XAM and BigWig packages to process the BAM file.
#' It will install Julia and necessary packages if they are not already available.
#' The conversion process considers sequencing type, read assignment strategy,
#' normalization method, and filters reads based on length.
#'
#' @note
#' Ensure that you have sufficient permissions to install Julia packages and
#' write to the output directory.
#'
#' @examples
#' \dontrun{
#' ribo_bam_to_bw(bam_file = "path/to/input.bam",
#'                output_file = "path/to/output.bw",
#'                seq_type = "singleEnd",
#'                assignType = "end5",
#'                normalization = "rpm",
#'                min_length = 25,
#'                max_length = 35)
#' }
#'
#'
#' @export
ribo_bam_to_bw <- function(bam_file = NULL,
                           output_file = NULL,
                           julia_path = NULL,
                           XAM_version = c(3,4),
                           seq_type = c("pairedEnd","singleEnd"),
                           assignType = c("end5","end3"),
                           normalization = c("rpm","counts"),
                           min_length = 28,
                           max_length = 31,
                           offset = ""){
  seq_type <- match.arg(seq_type,c("pairedEnd","singleEnd"))
  assignType <- match.arg(assignType,c("end5","end3"))
  normalization <- match.arg(normalization,c("rpm","counts"))
  XAM_version <- match.arg(XAM_version,c(3,4))

  # check julia
  if(!is.null(julia_path)){
    options(JULIA_HOME = julia_path)
  }

  options(timeout = max(1000, getOption("timeout")))
  JuliaCall::julia_setup(installJulia = TRUE)
  JuliaCall::julia_install_package_if_needed("XAM")
  JuliaCall::julia_library("XAM")
  JuliaCall::julia_install_package_if_needed("BigWig")
  JuliaCall::julia_library("BigWig")

  script_path <- paste0('include("',
                        system.file("extdata", "ribo_bam_to_bw.jl",package = "RiboProfiler"),
                        '")',collapse = "")

  # choose function
  JuliaCall::julia_eval(script_path)

  # check XAM version
  if(XAM_version == 3){
    ribo_bam_to_bw <- JuliaCall::julia_eval("ribo_bam_to_bw")
  }else{
    ribo_bam_to_bw <- JuliaCall::julia_eval("ribo_bam_to_bw2")
  }


  # run
  ribo_bam_to_bw(bam_file = bam_file,
                 output_file = output_file,
                 seq_type = seq_type,
                 assignType = assignType,
                 normalization = normalization,
                 min_length = as.integer(min_length),
                 max_length = as.integer(max_length),
                 offset = offset)

  message(paste0(bam_file,"has been processed!"))
}
