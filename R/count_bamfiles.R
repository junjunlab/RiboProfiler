#' Calculate read counts from BAM files using a GTF annotation
#'
#' This function processes BAM files to calculate read counts for specified features
#' using a GTF file, with options for ribosomal profiling and RNA sequencing.
#'
#' @param gtf_file A character string representing the path to a GTF file for annotation.
#' @param bam_files A character vector of paths to BAM files.
#' @param type A character string specifying the type of data: "ribo" for ribosome profiling or "rna" for RNA sequencing.
#' @param isPairedEnd A logical value indicating whether the BAM files contain paired-end reads. Default is \code{FALSE}.
#' @param exclude_upstream_nt An integer specifying the number of nucleotides to exclude upstream for ribosome profiling. Default is 60.
#' @param exclude_downstream_nt An integer specifying the number of nucleotides to exclude downstream for ribosome profiling. Default is 60.
#'
#' @return A matrix of read counts with genes as rows and samples as columns.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import.gff export.gff
#'
#' @details For ribosome profiling, the function adjusts the CDS regions by excluding specified nucleotides upstream and downstream before calculating counts.
#'
#' @export
#' @examples
#' \dontrun{
#' gtf_file <- "path/to/annotation.gtf"
#' bam_files <- c("sample1.bam", "sample2.bam")
#' counts <- get_counts_from_bam(gtf_file, bam_files, type = "rna", isPairedEnd = TRUE)
#' }
get_counts_from_bam <- function(gtf_file = NULL,
                                bam_files = NULL,
                                type = c("ribo","rna"),
                                isPairedEnd = FALSE,
                                exclude_upstream_nt = 60,
                                exclude_downstream_nt = 60){
  # check type
  type <- match.arg(type,c("ribo","rna"))

  # ==============================================================================================
  # calculate counts
  # ==============================================================================================

  if(type == "ribo"){
    gtf <- rtracklayer::import.gff(gtf_file,format = "gtf") %>%
      data.frame()

    gtf_ft <- gtf %>%
      dplyr::filter(type == "CDS") %>%
      dplyr::arrange(transcript_id,start,end)

    extend_gtf <- gtf_ft %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::mutate(start = dplyr::if_else(dplyr::row_number() == 1, start + exclude_upstream_nt, start),
                    end = dplyr::if_else(dplyr::row_number() == dplyr::n(), end - exclude_downstream_nt, end)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(width = end - start + 1) %>%
      # filter regions out of range
      dplyr::filter(end >= start)

    extend_gtf_gr <- GenomicRanges::GRanges(extend_gtf)

    # output
    rtracklayer::export.gff(con = "ribo_anno.gtf",object = extend_gtf_gr,format = "gtf")

    gtf_file <- "ribo_anno.gtf"
    feature <- "CDS"
  }else{
    feature <- "exon"
  }

  # count
  count_mtx <- Rsubread::featureCounts(files = bam_files,
                                       annot.ext = gtf_file,
                                       isGTFAnnotationFile = T,
                                       GTF.featureType = feature,
                                       GTF.attrType = "gene_id",
                                       isPairedEnd = isPairedEnd,
                                       GTF.attrType.extra = c("gene_name","gene_biotype"),
                                       isPairedEnd = FALSE,
                                       nthreads = 12)

  return(count_mtx)
}
