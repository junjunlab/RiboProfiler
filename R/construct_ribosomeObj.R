#' Construct a ribosomeObj Object
#'
#' This function constructs an S4 object of class \code{ribosomeObj}, which contains ribosome profiling data
#' and associated gene annotations. It takes as input a gene annotation file, and optionally allows for UTR extension
#' and selection of mapping, sequencing, and assignment types.
#'
#' @param gtf_file A character string specifying the path to the GTF file containing gene annotations.
#' This file is required for preparing gene annotation data.
#' @param genome_file A character string specifying the path of genome sequence.
#' @param extend_utr A logical value indicating whether to extend the UTR regions of genes.
#' Defaults to \code{FALSE}.
#' @param upstream_scale A numeric value specifying the scaling factor for the upstream UTR extension.
#' Used if \code{extend_utr = TRUE}.
#' @param upstream_extend A numeric value specifying how much to extend the upstream UTR.
#' Used if \code{extend_utr = TRUE}.
#' @param downstream_scale A numeric value specifying the scaling factor for the downstream UTR extension.
#' Used if \code{extend_utr = TRUE}.
#' @param downstream_extend A numeric value specifying how much to extend the downstream UTR.
#'  Used if \code{extend_utr = TRUE}.
#' @param has_created_data whether you have created raw data and save time to renalysis.
#'  Defaluts FALSE.
#' @param has_created_annotation whether you have created lognest gene annotation and save time to renalysis.
#'  Defaluts FALSE.
#' @param mapping_type A character string indicating the type of mapping, either "genome" or "transcriptome".
#' Defaults to "genome".
#' @param sam_file A character string specifying the path to the SAM file containing ribosome profiling
#' alignment data. Either \code{sam_file} or \code{bam_file} must be provided.
#' @param bam_file A character string specifying the path to the BAM file containing ribosome profiling
#' alignment data. Either \code{bam_file} or \code{sam_file} must be provided.
#' @param out_file_prefix A character string specifying the prefix for the output files generated during
#' the QC analysis.
#' @param annotation_prefix A character string specifying the prefix for the longest annotation output files.
#' Defaluts NULL.
#' @param rep_name The replicates name. Default NULL.
#' @param group_name The group name. Default NULL.
#' @param assignType A character string specifying the assignment type for mapping data,
#' either "end5" or "end3". Defaults to "end5".
#' @param seq_type A character string specifying the sequencing type, either "pairedEnd" or "singleEnd".
#' Defaults to "singleEnd".
#'
#' @details
#' This function performs several steps to construct a \code{ribosomeObj} object:
#' \itemize{
#'   \item It prepares gene annotation data from the provided GTF file.
#'   \item It optionally extends UTR regions based on the specified parameters.
#'   \item It performs QC analysis on the ribosome profiling data.
#'   \item It loads the annotation and QC data into a \code{ribosomeObj} object.
#' }
#'
#' The function relies on several helper functions, including \code{pre_longest_trans_info}, \code{gene_anno_extend}, \code{pre_qc_data}, and \code{load_qc_data}, to perform various tasks related to annotation processing and QC analysis.
#'
#' @return An S4 object of class \code{ribosomeObj}, which contains the following slots:
#' \item{raw.counts}{A data.frame containing the raw ribosome profiling counts.}
#' \item{shifted.data}{An empty data.frame for future shifted data (not yet used in this function).}
#' \item{normalized.data}{An empty data.frame for future normalized data (not yet used in this function).}
#' \item{longest.annotation}{A data.frame containing the longest annotated transcripts, potentially extended UTRs.}
#'
#' @examples
#' \dontrun{
#' # Example usage of construct_ribosomeObj
#' ribo <- construct_ribosomeObj(
#'   gtf_file = "path/to/gtf_file.gtf",
#'   mapping_type = "genome",
#'   sam_file = "path/to/alignment.sam",
#'   out_file_prefix = "output/ribo"
#' )}
#'
#' @importFrom methods new
#'
#' @export
construct_ribosomeObj <- function(gtf_file = NULL,
                                  genome_file = NULL,
                                  extend_utr = FALSE,
                                  upstream_scale = NULL,
                                  upstream_extend = NULL,
                                  downstream_scale = NULL,
                                  downstream_extend = NULL,
                                  has_created_data = FALSE,
                                  has_created_annotation = FALSE,
                                  mapping_type = c("genome", "transcriptome"),
                                  sam_file = NULL,
                                  bam_file = NULL,
                                  out_file_prefix = NULL,
                                  annotation_prefix = NULL,
                                  rep_name = NULL,
                                  group_name = NULL,
                                  assignType = c("end5", "end3"),
                                  seq_type = c("pairedEnd", "singleEnd")){
  mapping_type <- match.arg(mapping_type,c("genome", "transcriptome"))
  assignType <- match.arg(assignType,c("end5", "end3"))
  seq_type <- match.arg(seq_type,c("pairedEnd", "singleEnd"))

  # ==============================================================================
  # 1_prepare gene annotation file
  # ==============================================================================
  dir.create("annotation_data",showWarnings = F)

  if(is.null(annotation_prefix)){
    out_anno_file <- "annotation_data/longest_info.txt"
  }else{
    out_anno_file <- paste("annotation_data/",annotation_prefix,"_longest_info.txt",sep = "")
  }

  if(has_created_annotation == FALSE){
    pre_longest_trans_info(gtf_file = gtf_file,
                           out_file = out_anno_file)
  }

  # extending utr
  if(extend_utr == TRUE){
    if(is.null(annotation_prefix)){
      longest.annotation <- "annotation_data/longest_info_extend.txt"
    }else{
      longest.annotation <- paste("annotation_data/",annotation_prefix,"_longest_info_extend.txt",sep = "")
    }

    # extend CDS
    df_extend <- gene_anno_extend(longest_trans_file = out_anno_file,
                                  upstream_scale = upstream_scale,
                                  upstream_extend = upstream_extend,
                                  downstream_scale = downstream_scale,
                                  downstream_extend = downstream_extend,
                                  output_file = longest.annotation)


  }else{
    longest.annotation <- out_anno_file
  }

  # ==============================================================================
  # 2_perform QC analysis
  # ==============================================================================
  qc_df <- pre_qc_data(longest_trans_file = longest.annotation,
                       mapping_type = mapping_type,
                       sam_file = sam_file,
                       bam_file = bam_file,
                       out_file_prefix = out_file_prefix,
                       has_created_data = has_created_data,
                       rep_name = rep_name,
                       group_name = group_name,
                       assignType = assignType,
                       seq_type = seq_type)

  # load qc data
  # qc_df <- load_qc_data()

  # ========================================================================================
  # load gene annotation
  # ========================================================================================
  gene_anno <- read.delim(longest.annotation,header = F)

  # add colnames
  colnames(gene_anno) <- c("id","gene_name","gdi","tid","chr","strand","cds_rg","exon_rg","utr5","cds","utr3")

  # ============================================================================
  # create homerResult object
  # ============================================================================
  res <- methods::new("ribosomeObj",
                      "raw.counts" = qc_df,
                      "genome.file" = genome_file,
                      "shifted" = "FALSE",
                      "normalized" = "FALSE",
                      "longest.anno.file" = longest.annotation,
                      "longest.annotation" = gene_anno)

  return(res)
}
