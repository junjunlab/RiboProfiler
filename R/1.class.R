#' ribosomeObj Class
#'
#' A S4 class that represents ribosome profiling data.
#'
#' @slot raw.counts A data.frame containing raw ribosome profiling counts.
#' @slot shifted.data A data.frame containing ribosome profiling data after E/P/A-site shifting.
#' @slot normalized.data A data.frame containing normalized ribosome profiling data.
#' @slot enriched.data A data.frame containing enrichment data for seRP data.
#' @slot longest.annotation A data.frame containing the longest annotated sequences.
#' @slot shifted A character vector indicating whether the data has been E/P/A-site shifted.
#' @slot normalized A character vector indicating whether the data has been normalized.
#' @slot genome.file A character vector indicating the path of genome fasta file.
#' @slot longest.anno.file A character vector indicating the path of longest transcript annotation file.
#' @slot CDS.sequence A character vector containing the coding DNA sequences (CDS).
#' @slot amino.acid.sequence A character vector containing the amino acid sequences derived from the CDS.
#'
#' @details
#' This class is used for storing and manipulating ribosome profiling data.
#' The class contains several slots that store raw, shifted, and normalized data, as well as
#' annotation information related to coding sequences and their corresponding amino acid sequences.
#'
#' @examples
#' \dontrun{
#' # Create an empty ribosomeObj object
#' ribo <- new("ribosomeObj",
#'             raw.counts = data.frame(),
#'             shifted.data = data.frame(),
#'             normalized.data = data.frame(),
#'             longest.annotation = data.frame(),
#'             shifted = "",
#'             normalized = "",
#'             genome.file = "",
#'             longest.anno.file = "",
#'             CDS.sequence = "",
#'             amino.acid.sequence = "")
#' }
#'
#' @export
ribosomeObj <- setClass("ribosomeObj",
                        slots = list("raw.counts" = "data.frame",
                                     "shifted.data" = "data.frame",
                                     "normalized.data" = "data.frame",
                                     "enriched.data" = "data.frame",
                                     "longest.annotation" = "data.frame",
                                     "shifted" = "character",
                                     "normalized" = "character",
                                     "genome.file" = "character",
                                     "longest.anno.file" = "character",
                                     "CDS.sequence" = "character",
                                     "amino.acid.sequence" = "character"))



