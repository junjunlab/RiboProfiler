#' Get Normalized Counts for Transcript Positions
#'
#' This function normalizes the counts at individual positions to the average counts of the transcript,
#' with options to exclude a specified number of codons from the start and end of the coding sequence.
#'
#' @param longest_trans_file A string specifying the path to the file containing the longest transcript information.
#'        This file should contain columns for transcript ID, gene name, gene ID, chromosome, strand,
#'        CDS range, exon range, 5' UTR length, CDS length, and 3' UTR length.
#' @param qc_df A data frame containing quality control data, including columns for sample, transcript ID,
#'        transcript position, and counts.
#' @param exclude_upstreamCodon An integer specifying the number of codons to exclude from the start of the coding sequence.
#'        Default is 30.
#' @param exclude_downstreamCodon An integer specifying the number of codons to exclude from the end of the coding sequence.
#'        Default is 0.
#'
#' @return A data frame with normalized counts for each position in each transcript across all samples.
#'         The returned data frame includes all original columns from qc_df plus a new column 'norm_exp'
#'         representing the normalized expression.
#'
#' @details The function performs the following steps:
#'   1. Reads and processes the gene annotation file.
#'   2. Calculates total reads and average expression for each transcript, excluding specified regions.
#'   3. Normalizes the counts at each position by dividing by the average expression of the transcript.
#'
#' @export
get_nomalized_counts <- function(longest_trans_file = NULL,
                                 qc_df = NULL,
                                 exclude_upstreamCodon = 30,
                                 exclude_downstreamCodon = 0){
  # ============================================================================
  # at individual position was nomalized to the average counts of the transcript
  # ============================================================================
  # loop read file
  # gene annotation
  ganao <- read.delim(longest_trans_file,header = F)
  colnames(ganao) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                       "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")

  gene_lenth <- ganao %>%
    # dplyr::mutate(trancript_len = `5UTR_length` + CDS_length + `3UTR_length`) %>%
    # dplyr::select(trans_id,trancript_len)
    dplyr::select(trans_id,`5UTR_length`,CDS_length,`3UTR_length`)

  # loop
  sp <- unique(qc_df$sample)

  # x = 1
  purrr::map_df(seq_along(sp),function(x){
    tmp <- subset(qc_df,sample == sp[x])
    total_exp <- tmp %>%
      dplyr::left_join(y = gene_lenth,by = "trans_id") %>%
      dplyr::filter(trans_pos >= `5UTR_length` + 1 + exclude_upstreamCodon*3 & trans_pos <= CDS_length + `5UTR_length` - exclude_downstreamCodon*3) %>%
      dplyr::group_by(trans_id) %>%
      dplyr::summarise(total_reads = sum(counts))

    gene_exp <- gene_lenth %>%
      dplyr::filter(trans_id %in% unique(total_exp$trans_id)) %>%
      dplyr::left_join(y = total_exp,by = "trans_id")

    gene_exp$average_exp <- gene_exp$total_reads/gene_exp$CDS_length

    qc_df_new <- tmp %>%
      dplyr::left_join(y = gene_exp[,c("trans_id","average_exp")],by = "trans_id") %>%
      dplyr::mutate(norm_exp = counts/average_exp) %>%
      dplyr::select(-average_exp)

    return(qc_df_new)
  }) -> qc_df_nomalized
}
