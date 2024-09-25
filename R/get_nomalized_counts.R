#' Get Normalized Counts for Transcript Positions
#'
#' This function normalizes the counts at individual positions to the average counts of the transcript,
#' with options to exclude a specified number of codons from the start and end of the coding sequence.
#'
#' @param object ribosomeObj object.
#' @param exclude_upstreamCodon An integer specifying the number of codons to exclude from the start of the coding sequence.
#'        Default is 5
#' @param exclude_downstreamCodon An integer specifying the number of codons to exclude from the end of the coding sequence.
#'        Default is 5.
#' @param norm_type The nomalization methods for ribosome density. "average" is calculated by
#' the count at each position divided by mean density across cds region. "rpm"
#' is calculated by the count at each position divided by the total counts and multiplied with 10^6.
#' Default is "average".
#' @param ... Useless args.
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
#'
#' @export
setGeneric("get_nomalized_counts",
           function(object,
                    exclude_upstreamCodon = 5,
                    exclude_downstreamCodon = 5,
                    norm_type = c("average","rpm"),...) standardGeneric("get_nomalized_counts"))





#' method for get_nomalized_counts
#'
#' @rdname get_nomalized_counts
#' @exportMethod get_nomalized_counts
setMethod("get_nomalized_counts",
          signature(object = "ribosomeObj"),
          function(object,
                   exclude_upstreamCodon = 5,
                   exclude_downstreamCodon = 5,
                   norm_type = c("average","rpm"),...){
            norm_type <- match.arg(norm_type,c("average","rpm"))
            # ============================================================================
            # at individual position was nomalized to the average counts of the transcript
            # ============================================================================
            # loop read file
            # gene annotation
            ganao <- object@longest.annotation
            colnames(ganao) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                                 "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")

            gene_lenth <- ganao %>%
              # dplyr::mutate(trancript_len = `5UTR_length` + CDS_length + `3UTR_length`) %>%
              # dplyr::select(trans_id,trancript_len)
              dplyr::select(trans_id,`5UTR_length`,CDS_length,`3UTR_length`)

            # loop
            if(length(object@shifted.data) != 0){
              qc_df <- object@shifted.data
            }else{
              qc_df <- object@raw.counts
            }

            sp <- unique(qc_df$sample)

            # x = 1
            purrr::map_df(seq_along(sp),function(x){
              tmp <- subset(qc_df,sample == sp[x])

              # check normalization methods
              if(norm_type == "average"){
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
                  dplyr::filter(trans_id %in% unique(gene_exp$trans_id)) %>%
                  dplyr::mutate(norm_exp = counts/average_exp) %>%
                  dplyr::select(-average_exp)
              }else{
                total_counts <- sum(tmp$counts)
                tmp$norm_exp <- (tmp$counts/total_counts)*10^6
                qc_df_new <- tmp
              }

              return(qc_df_new)
            }) -> qc_df_nomalized

            object@normalized.data <- qc_df_nomalized
            object@normalized <- "TRUE"

            return(object)
          }
)
